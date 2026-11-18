"""Step 10 (robustness checks 1 and 3, docs/ROBUSTNESS_PLAN_v1.md): second pass over the raw memory recordings.

For every participant:
  check 1  cardiac field artifact: in each ROI channel (200 Hz, after the frozen median reference) the
           R locked artifact is estimated as a moving average over 80 neighbouring beats of the 1 Hz high
           passed signal, window -0.10 to +0.50 s around each detected R peak, tapered, and subtracted.
           Envelopes are then built exactly as in the frozen step 06.
  check 3  ECG derived respiration (EDR): R wave amplitude at each detected beat, interpolated to 4 Hz,
           band pass 0.08 to 0.6 Hz.
Beat times come from the frozen beat tables (inserted beats are not used). The frozen heart series is
copied unchanged. Output: derivatives/slow_cfa/<sub>_task-memory_slow.npz and a QC table. Resumable.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import gc
import platform
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import scipy
from scipy.signal import butter, sosfiltfilt
from scipy.signal.windows import tukey

from hb import bidsio, cardiac, slow_signals
from hb.download_ds008302 import PROJECT_ROOT, append_execution_log
from hb.pilot_qc import write_tsv

FS200 = 200.0
PRE, POST = 0.10, 0.50
N_AVG = 80
TAPER = 0.25
ROI_CHANNELS = {band: chans for band, (_, chans) in slow_signals.ROIS.items()}


def load_detected_beats(path: Path) -> np.ndarray:
    t = []
    with open(path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["beat_inserted"] != "1":
                t.append(float(r["beat_time_s"]))
    return np.array(t)


def subtract_cfa(x, beats, fs=FS200, pre=PRE, post=POST, n_avg=N_AVG):
    """Return (cleaned copy of x, R locked peak to peak before, after), both measured on the 1 Hz high passed signal."""
    sos = butter(2, 1.0, btype="highpass", fs=fs, output="sos")
    hp = sosfiltfilt(sos, x)
    a, b = int(round(pre * fs)), int(round(post * fs))
    ok = beats[(beats - a >= 0) & (beats + b <= x.size)]
    if ok.size < n_avg:
        return x.copy(), float("nan"), float("nan")
    idx = ok[:, None] + np.arange(-a, b)[None, :]
    E = hp[idx]
    c = np.vstack([np.zeros((1, E.shape[1])), np.cumsum(E, axis=0)])
    lo = np.clip(np.arange(ok.size) - n_avg // 2, 0, ok.size - n_avg)
    T = (c[lo + n_avg] - c[lo]) / n_avg
    y = x.copy()
    np.add.at(y, idx, -(T * tukey(a + b, TAPER)[None, :]))
    after = sosfiltfilt(sos, y)[idx].mean(axis=0)
    return y, float(np.ptp(E.mean(axis=0))), float(np.ptp(after))


def make_cleaner(beat_times_s):
    beats = np.round(np.asarray(beat_times_s) * FS200).astype(int)
    wanted = sorted({c for chans in ROI_CHANNELS.values() for c in chans})

    def clean(D, names):
        info = {}
        for c in wanted:
            if c in names:
                k = names.index(c)
                y, p0, p1 = subtract_cfa(D[k].astype(np.float64), beats)
                D[k] = y.astype(np.float32)
                info[f"_p2p_{c}"] = (p0, p1)
        for band, chans in ROI_CHANNELS.items():
            v = [info[f"_p2p_{c}"] for c in chans if f"_p2p_{c}" in info]
            info[f"cfa_p2p_before_{band}_uV"] = round(float(np.nanmean([p for p, _ in v])), 3) if v else ""
            info[f"cfa_p2p_after_{band}_uV"] = round(float(np.nanmean([q for _, q in v])), 3) if v else ""
        return {k: v for k, v in info.items() if not k.startswith("_")}
    return clean


def edr_series(ecg, fs, beat_times_s, n4, fs_out=4.0):
    """ECG derived respiration from R wave amplitude. Returns (series at fs_out, median R amplitude)."""
    x = cardiac.detect_rpeaks(ecg, fs)["filtered"]
    w = int(round(0.02 * fs))
    s = np.round(np.asarray(beat_times_s) * fs).astype(int)
    s = s[(s - w >= 0) & (s + w < x.size)]
    amp = x[s[:, None] + np.arange(-w, w + 1)[None, :]].max(axis=1)
    grid = np.arange(n4) / fs_out
    y = np.interp(grid, s / fs, amp)
    y = sosfiltfilt(butter(2, [0.08, 0.6], btype="band", fs=fs_out, output="sos"), y - y.mean())
    return y, float(np.median(amp))


def process(root: Path, sub: str, log) -> dict:
    sdir = PROJECT_ROOT / "derivatives" / "slow"
    odir = PROJECT_ROOT / "derivatives" / "slow_cfa"
    odir.mkdir(parents=True, exist_ok=True)
    out = odir / f"{sub}_task-memory_slow.npz"
    if out.exists():
        return {"sub": sub, "status": "skipped (already done)"}
    beats = load_detected_beats(PROJECT_ROOT / "derivatives" / "beats" / f"{sub}_task-memory_beats.tsv")
    frozen = np.load(sdir / f"{sub}_task-memory_slow.npz")
    n = frozen["nn"].size
    rec = bidsio.load_set(bidsio.rec_paths(root, sub, "memory")["set"])
    fs = rec["fs"]
    ecg = rec["data"][rec["all_labels"].index("ECG")].astype(float)
    env, info = slow_signals.roi_envelopes(rec, log, clean=make_cleaner(beats))
    del rec
    gc.collect()
    edr, ramp = edr_series(ecg, fs, beats, n)
    del ecg
    arrays = {"fs": 4.0, "nn": frozen["nn"], "nn_bridged": frozen["nn_bridged"], "edr": edr}
    for k, (y, bad) in env.items():
        if y.size < n:
            raise ValueError(f"{k} envelope shorter than the frozen series ({y.size} < {n})")
        arrays[k], arrays[k + "_interpolated"] = y[:n], bad[:n]
    tmp = out.with_name(out.stem + ".part.npz")
    np.savez_compressed(tmp, **arrays)
    tmp.replace(out)
    row = {"sub": sub, "status": "ok", "n_detected_beats": int(beats.size), "r_amplitude_median": round(ramp, 1), **info}
    for band in ROI_CHANNELS:
        a, b = info.get(f"cfa_p2p_before_{band}_uV"), info.get(f"cfa_p2p_after_{band}_uV")
        row[f"cfa_reduction_{band}_pct"] = round(100 * (1 - b / a), 1) if a not in ("", None) and a > 0 else ""
    return row


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--subjects", nargs="*")
    a = ap.parse_args(argv)
    root = bidsio.raw_root()
    subs = a.subjects or sorted(f.name.split("_")[0] for f in (PROJECT_ROOT / "derivatives" / "slow").glob("sub-*_task-memory_slow.npz"))
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    (PROJECT_ROOT / "logs").mkdir(exist_ok=True)
    logf = open(PROJECT_ROOT / "logs" / f"robust_reprocess_{stamp}.log", "w", encoding="utf-8")

    def log(m):
        print(m, flush=True)
        logf.write(m + "\n")
        logf.flush()

    log(f"robustness reprocess {stamp} | {len(subs)} participants | python {platform.python_version()} numpy {np.__version__} scipy {scipy.__version__}")
    t0, rows, failed = time.time(), [], 0
    qdir = PROJECT_ROOT / "results" / "robustness"
    qdir.mkdir(parents=True, exist_ok=True)
    for k, sub in enumerate(subs, 1):
        t1 = time.time()
        try:
            r = process(root, sub, lambda m: None)
        except Exception as e:
            failed += 1
            r = {"sub": sub, "status": "FAILED: " + "".join(traceback.format_exception(e)).strip().replace("\n", " | ")}
        if r["status"] == "ok":
            (qdir / "qc").mkdir(exist_ok=True)
            write_tsv(qdir / "qc" / f"{sub}_reprocess_qc.tsv", [r])
        rows.append(r)
        extra = "" if r["status"] != "ok" else " | " + ", ".join(f"{b} CFA {r.get(f'cfa_p2p_before_{b}_uV')} -> {r.get(f'cfa_p2p_after_{b}_uV')} uV" for b in ROI_CHANNELS)
        log(f"[{k}/{len(subs)}] {sub}: {r['status'][:300]}{extra} ({time.time() - t1:.0f} s)")
        gc.collect()
    log(f"finished: {len(subs) - failed} ok or skipped, {failed} failed, {(time.time() - t0) / 60:.1f} min")
    append_execution_log([
        f"## {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} robustness reprocess (step 10, checks 1 and 3)",
        f"- Command: python -m hb.robust_reprocess {' '.join(sys.argv[1:])}; Python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}",
        f"- {len(subs)} participants, {failed} failed, {(time.time() - t0) / 60:.1f} min. Outputs: derivatives/slow_cfa, results/robustness/qc; log logs/robust_reprocess_{stamp}.log",
    ])
    return 0 if failed == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
