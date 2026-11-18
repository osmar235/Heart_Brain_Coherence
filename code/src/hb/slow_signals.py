"""Step 06: slow (4 Hz) heart and EEG envelope series for the pilot, and their trial locked averages.

Purpose: give the null model benchmark real response shapes instead of assumed ones. This step does not
compute coherence and does not compare span with span+1. Averages are pooled per sequence length only.

PROVISIONAL EEG handling (pilot only, not the final pipeline): decimation to 200 Hz, common MEDIAN
reference across the EEG channels, zero phase band pass, Hilbert amplitude, mean over the region of
interest, decimation to 4 Hz, and interpolation of envelope samples above median + 8 robust SD.
Bands and regions follow Kotyusov et al. 2026 (preprint): frontal midline theta, posterior alpha,
left centroparietal beta.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gc
import platform
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.signal import resample_poly

from hb import bidsio, hbcoh as H
from hb.download_ds008302 import PROJECT_ROOT, append_execution_log, load_json
from hb.pilot_qc import write_tsv

FS = H.FS_SLOW
ROIS = {
    "theta": ((4.0, 8.0), ["AFz", "Fz", "F1", "F2"]),
    "alpha": ((9.0, 14.0), ["P1", "P2", "P3", "P4", "P5", "P6", "PO3", "PO4"]),
    "beta": ((16.0, 22.0), ["C3", "C1", "CP3", "CP5", "C5", "CP1"]),
}
SIGNALS = ["nn", "theta", "alpha", "beta"]


def nn_series(beats_tsv: Path, n: int) -> tuple[np.ndarray, np.ndarray]:
    """Corrected NN intervals on the 4 Hz grid; marked intervals are bridged by linear interpolation.
    Returns (series, bridged mask)."""
    t, nn, mark = [], [], []
    with open(beats_tsv, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["nn_s"] == "n/a":
                continue
            t.append(float(r["beat_time_s"])); nn.append(float(r["nn_s"])); mark.append(r["interval_marked"] == "1")
    t, nn, mark = np.array(t), np.array(nn), np.array(mark)
    grid = np.arange(n) / FS
    y = np.interp(grid, t[~mark], nn[~mark])
    bridged = np.zeros(n, bool)
    for tb, d in zip(t[mark], nn[mark]):
        bridged |= (grid >= tb - d) & (grid <= tb)
    return y, bridged


def roi_envelopes(rec: dict, log, clean=None) -> tuple[dict, dict]:
    """clean: optional function(D, names) -> dict, applied to the 200 Hz median referenced channels
    before envelopes (used only by robustness step 10; the frozen pipeline passes None)."""
    labels, fs = rec["all_labels"], rec["fs"]
    q = int(round(fs / 200.0))
    eeg_idx = [i for i, lab in enumerate(labels) if lab not in ("ECG", "PPG")]
    n200 = int(np.ceil(rec["n_samples"] / q))
    D = np.empty((len(eeg_idx), n200), np.float32)
    for k, i in enumerate(eeg_idx):
        x = rec["data"][i].astype(np.float64)
        D[k] = resample_poly(x - np.median(x), 1, q)[:n200]
    step = 200000
    for a in range(0, n200, step):
        D[:, a:a + step] -= np.median(D[:, a:a + step], axis=0)
    names = [labels[i] for i in eeg_idx]
    out, info = {}, {}
    if clean is not None:
        info.update(clean(D, names))
    for band, (rng, chans) in ROIS.items():
        use = [c for c in chans if c in names]
        if not use:
            continue
        env = None
        for c in use:
            _, e = H.band_envelope(D[names.index(c)].astype(np.float64), 200.0, rng)
            env = e if env is None else env + e
        env /= len(use)
        med = np.median(env)
        rsd = 1.4826 * np.median(np.abs(env - med))
        bad = env > med + 8 * rsd
        bad = np.convolve(bad.astype(float), np.ones(int(2 * FS) + 1), mode="same") > 0      # widen by 1 s each side
        y = env.copy()
        if bad.any() and (~bad).any():
            idx = np.arange(env.size)
            y[bad] = np.interp(idx[bad], idx[~bad], env[~bad])
        out[band] = (y, bad)
        info[f"{band}_channels"] = ",".join(use)
        info[f"{band}_pct_interpolated"] = round(100 * float(bad.mean()), 2)
        log(f"    {band}: {len(use)} channels, {info[f'{band}_pct_interpolated']}% of samples interpolated")
    return out, info


def trial_profiles(sig: dict, parsed: dict, min_trials: int) -> dict:
    """Mean and SEM per sequence length, from 2 s before trial start to 12 s after the recall prompt.
    Signals are z scored over all in block samples first."""
    n = len(next(iter(sig.values())))
    t = np.arange(n) / FS
    inblock = np.zeros(n, bool)
    for b in parsed["blocks"]:
        if b["trials"]:
            inblock |= (t >= b["t0"]) & (t < b["t1"])
    by_len = {}
    for tr in parsed["trials"]:
        by_len.setdefault(tr["length"], []).append(tr)
    prof = {}
    for length, trials in sorted(by_len.items()):
        if len(trials) < min_trials:
            continue
        enc = float(np.median([tr["t_probe"] - tr["t_start"] for tr in trials]))
        m = int(round((2.0 + enc + 12.0) * FS))
        entry = {"n_trials": len(trials), "t": np.arange(m) / FS - 2.0, "probe_at": enc}
        for name, x in sig.items():
            z = (x - x[inblock].mean()) / x[inblock].std()
            rows = []
            for tr in trials:
                i0 = int(round((tr["t_start"] - 2.0) * FS))
                if i0 >= 0 and i0 + m <= n:
                    rows.append(z[i0:i0 + m])
            rows = np.array(rows)
            mean = rows.mean(axis=0)
            entry[name + "_mean"] = mean
            entry[name + "_sem"] = rows.std(axis=0, ddof=1) / np.sqrt(len(rows))
            entry[name + "_evoked_var_fraction"] = float(mean.var() / rows.var())
            a, b = int(round((2.0 + 2.5) * FS)), int(round((2.0 + enc) * FS))             # digit presentation part
            segm = mean[a:b] - np.polyval(np.polyfit(np.arange(b - a), mean[a:b], 2), np.arange(b - a))
            spec = np.abs(np.fft.rfft(segm * np.hanning(segm.size))) ** 2
            f = np.fft.rfftfreq(segm.size, 1 / FS)
            entry[name + "_digit_rate_power_fraction"] = float(spec[(f > 0.42) & (f < 0.58)].sum() / spec[f > 0.05].sum())
        prof[length] = entry
    return prof


def process(root: Path, sub: str, task: str, outdir: Path, sdir: Path, bdir: Path, log, min_trials=20) -> list[dict]:
    p = bidsio.rec_paths(root, sub, task)
    t0 = time.time()
    outdir.mkdir(parents=True, exist_ok=True)
    sdir.mkdir(parents=True, exist_ok=True)
    rec = bidsio.load_set(p["set"])
    dur = rec["n_samples"] / rec["fs"]
    env, info = roi_envelopes(rec, log)
    del rec
    gc.collect()
    n = min(len(v[0]) for v in env.values())
    nn, bridged = nn_series(bdir / f"{sub}_task-{task}_beats.tsv", n)
    sig = {"nn": nn, **{k: v[0][:n] for k, v in env.items()}}
    sdir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(sdir / f"{sub}_task-{task}_slow.npz", fs=FS, nn=nn, nn_bridged=bridged,
                        **{k: v[0][:n] for k, v in env.items()}, **{k + "_interpolated": v[1][:n] for k, v in env.items()})
    rows = []
    base = {"sub": sub, "task": task, "duration_s": round(dur, 1), "nn_pct_bridged": round(100 * float(bridged.mean()), 2), **info}
    if task != "memory":
        rows.append(base)
        log(f"    saved slow series ({time.time() - t0:.0f} s)")
        return rows
    parsed = bidsio.parse_memory_events(bidsio.read_tsv(p["events"]), dur)
    prof = trial_profiles(sig, parsed, min_trials)
    np.savez_compressed(outdir / f"{sub}_trial_profiles.npz",
                        **{f"L{L}_{k}": v for L, e in prof.items() for k, v in e.items()})
    names = [s for s in SIGNALS if s in sig]
    fig, axes = plt.subplots(len(names), max(len(prof), 1), figsize=(4.2 * max(len(prof), 1), 2.4 * len(names)), squeeze=False, sharey="row")
    for j, (L, e) in enumerate(prof.items()):
        r = {**base, "length": L, "n_trials": e["n_trials"], "trial_start_to_prompt_s": round(e["probe_at"], 2)}
        for i, name in enumerate(names):
            ax = axes[i, j]
            ax.fill_between(e["t"], e[name + "_mean"] - e[name + "_sem"], e[name + "_mean"] + e[name + "_sem"], alpha=0.3)
            ax.plot(e["t"], e[name + "_mean"], lw=1)
            ax.axvline(0, color="k", lw=0.6)
            ax.axvline(e["probe_at"], color="C3", lw=0.6)
            if i == 0:
                ax.set_title(f"{L} digits, {e['n_trials']} trials", fontsize=9)
            if j == 0:
                ax.set_ylabel(f"{name} (z)")
            if i == len(names) - 1:
                ax.set_xlabel("s from trial start; red = recall prompt")
            r[f"{name}_evoked_var_fraction"] = round(e[name + "_evoked_var_fraction"], 4)
            r[f"{name}_digit_rate_power_fraction"] = round(e[name + "_digit_rate_power_fraction"], 4)
            r[f"{name}_profile_range_z"] = round(float(np.ptp(e[name + "_mean"])), 3)
        rows.append(r)
    fig.suptitle(f"{sub}: trial locked averages per sequence length (pooled, no condition contrast)", fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / f"{sub}_trial_profiles.png", dpi=110)
    plt.close(fig)
    log(f"    saved slow series and profiles for lengths {list(prof)} ({time.time() - t0:.0f} s)")
    return rows


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--subjects", nargs="*")
    ap.add_argument("--tasks", nargs="*", default=["memory", "rest", "video"])
    a = ap.parse_args(argv)
    cfg = load_json(PROJECT_ROOT / "config" / "download_config.json")
    root = bidsio.raw_root()
    outdir = PROJECT_ROOT / "results" / "slow_signals"
    sdir = PROJECT_ROOT / "derivatives" / "slow"
    bdir = PROJECT_ROOT / "derivatives" / "beats"
    outdir.mkdir(parents=True, exist_ok=True)
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    (PROJECT_ROOT / "logs").mkdir(exist_ok=True)
    logf = open(PROJECT_ROOT / "logs" / f"slow_signals_{stamp}.log", "w", encoding="utf-8")

    def log(m):
        print(m, flush=True)
        logf.write(m + "\n")
        logf.flush()

    log(f"slow signals {stamp} | python {platform.python_version()} numpy {np.__version__} scipy {scipy.__version__}")
    rows, failed = [], 0
    for sub in a.subjects or cfg["pilot"]["subjects"]:
        for task in a.tasks:
            log(f"{sub} task-{task}")
            try:
                rows += process(root, sub, task, outdir, sdir, bdir, log)
            except Exception as e:
                import traceback
                failed += 1
                log("  ERROR: " + "".join(traceback.format_exception(e)).strip())
    write_tsv(outdir / "slow_signals_summary.tsv", rows)
    append_execution_log([
        f"## {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} slow signals and trial profiles (step 06)",
        f"- Command: python -m hb.slow_signals {' '.join(sys.argv[1:])}",
        f"- Python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}, {platform.platform()}",
        f"- Errors: {failed}. EEG handling is PROVISIONAL (median reference, envelope outlier interpolation).",
        f"- Outputs: derivatives/slow/*.npz; results/slow_signals/*_trial_profiles.png and .npz, slow_signals_summary.tsv; log logs/slow_signals_{stamp}.log",
    ])
    log("DONE" if failed == 0 else f"FINISHED WITH {failed} ERROR(S)")
    return 0 if failed == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
