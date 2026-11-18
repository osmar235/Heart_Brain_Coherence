"""Step 05b: which control separates task rhythm from ongoing coupling when trial timing is irregular?

SIMULATION ONLY. The real trial timing of the pilot participants (events.tsv, public metadata) is used
as a scaffold; every signal is synthetic and no recording is read. For each scenario the truth is known:

  task_only     heart and brain both respond to every trial, with independent trial to trial gains;
                no ongoing coupling. Any "excess coherence" here is a false positive.
  ongoing       task_only plus a shared ongoing 0.07 to 0.13 Hz process. Excess should be detected.
  shared_gain   task_only, but the trial to trial gains of heart and brain are identical. No ongoing
                coupling; shows how each control treats trial level covariation of evoked responses.

Controls compared, each giving excess = same block coherence minus control coherence:
  clock         EEG from another block of the same sequence length, paired by clock time
  aligned       EEG from another block of the same length, time warped so its trials line up
  resid_clock   trial locked averages removed from both signals first, then the clock control
  residfit_clock  same, with the average scaled to each trial by least squares before removal
"""
from __future__ import annotations

import argparse
import datetime as dt
import platform
import sys
import time
from pathlib import Path

import numpy as np
import scipy
from scipy import stats

from hb import bidsio, hbcoh as H
from hb import synth as S
from hb.download_ds008302 import PROJECT_ROOT, append_execution_log, load_json
from hb.vendor.emgsim_v3_3_0 import functions as emg

FS = H.FS_SLOW
BANDS = {"LF": H.BANDS["LF"], "LFsub": (0.07, 0.13), "HF": H.BANDS["HF"]}
CONTROLS = ["clock", "aligned", "resid_clock", "residfit_clock"]
SCENARIOS = ["task_only", "ongoing", "shared_gain"]


def load_timing(root: Path, sub: str) -> dict:
    p = bidsio.rec_paths(root, sub, "memory")
    dur = float(bidsio.sidecar(p["json"])["RecordingDuration"])
    parsed = bidsio.parse_memory_events(bidsio.read_tsv(p["events"]), dur)
    blocks = []
    for b in parsed["blocks"]:
        tr = b["trials"]
        if len(tr) != 10:
            continue
        ends = [t["t_start"] for t in tr[1:]] + [b["t1"]]
        blocks.append(dict(t0=b["t0"], t1=b["t1"], length=b["length"],
                           start=np.array([t["t_start"] for t in tr]), probe=np.array([t["t_probe"] for t in tr]),
                           end=np.array(ends)))
    return dict(sub=sub, dur=dur, blocks=blocks)


def gamma_kernel(tau):
    t = np.arange(0, 8 * tau, 1 / FS)
    k = t * np.exp(-t / tau)
    return k / k.sum()


def simulate(timing: dict, scenario: str, rng) -> tuple[np.ndarray, np.ndarray]:
    n = int(timing["dur"] * FS)
    t = np.arange(n) / FS
    task = {"heart": np.zeros(n), "brain": np.zeros(n)}
    inblock = np.zeros(n, bool)
    for b in timing["blocks"]:
        inblock |= (t >= b["t0"]) & (t < b["t1"])
        for k in range(10):
            z_h, z_b = rng.standard_normal(2)
            if scenario == "shared_gain":
                z_b = z_h
            enc = (t >= b["start"][k] + 2.5) & (t < b["probe"][k])
            rec = (t >= b["probe"][k]) & (t < b["end"][k])
            task["heart"] += (1 + 0.3 * z_h) * enc - 0.6 * (1 + 0.3 * rng.standard_normal()) * rec
            task["brain"] += -(1 + 0.3 * z_b) * enc + 0.8 * (1 + 0.3 * rng.standard_normal()) * rec
    out = []
    shared = S.band_noise(n, FS, (0.07, 0.13), rng) if scenario == "ongoing" else np.zeros(n)
    for name, tau in (("heart", 3.0), ("brain", 1.5)):
        x = np.convolve(task[name], gamma_kernel(tau))[:n]
        x = (x - x[inblock].mean()) / x[inblock].std()
        out.append(x + S.band_noise(n, FS, (0.01, 1.0), rng) + 0.7 * shared)
    return out[0], out[1]


def seg(x, b, n=None):
    i0 = int(round(b["t0"] * FS))
    i1 = int(round(b["t1"] * FS)) if n is None else i0 + n
    return x[i0:i1]


def warp_to(x, src, dst):
    """Sample block src of signal x on the timeline of block dst, lining up trial starts and recall prompts."""
    def anchors(b):
        a = [b["t0"]]
        for k in range(10):
            a += [b["start"][k], b["probe"][k]]
        return np.maximum.accumulate(np.r_[a, b["t1"]])
    td = np.arange(int(round(dst["t0"] * FS)), int(round(dst["t1"] * FS))) / FS
    ts = np.interp(td, anchors(dst), anchors(src))
    return np.interp(ts, np.arange(x.size) / FS, x)


def remove_evoked(x, blocks, fit_gain):
    """Subtract trial locked averages (per sequence length): trial start to recall prompt, then the first 12 s of recall."""
    y = x.copy()
    by_len = {}
    for b in blocks:
        by_len.setdefault(b["length"], []).append(b)
    rec_n = int(12 * FS)
    for length, bl in by_len.items():
        trials = [(b["start"][k], b["probe"][k], b["end"][k]) for b in bl for k in range(10)]
        enc_n = int(round(np.median([p - s for s, p, _ in trials]) * FS))
        enc = np.full((len(trials), enc_n), np.nan)
        rec = np.full((len(trials), rec_n), np.nan)
        for i, (s, p, e) in enumerate(trials):
            i0 = int(round(s * FS))
            m = min(enc_n, x.size - i0)
            enc[i, :m] = x[i0:i0 + m]
            j0 = int(round(p * FS))
            m = max(0, min(rec_n, int(round(e * FS)) - j0, x.size - j0))
            rec[i, :m] = x[j0:j0 + m]
        enc_t = np.nanmean(enc, axis=0)
        cnt = np.isfinite(rec).sum(axis=0)
        rec_t = np.where(cnt >= 5, np.nanmean(np.where(cnt >= 5, rec, 0.0), axis=0), 0.0)
        for i, (s, p, e) in enumerate(trials):
            for row, tpl, start in ((enc[i], enc_t, s), (rec[i], rec_t, p)):
                ok = np.isfinite(row)
                if not ok.any():
                    continue
                g = 1.0
                if fit_gain:
                    tc = tpl[ok] - tpl[ok].mean()
                    g = float(np.dot(row[ok] - row[ok].mean(), tc) / (np.dot(tc, tc) + 1e-12))
                i0 = int(round(start * FS))
                y[i0:i0 + ok.sum()] -= g * tpl[:ok.sum()]
    return y


def band_values(cx, cy):
    coh = H.coherence(cx, cy)
    return {k: H.band_summary(coh, v)[0] for k, v in BANDS.items()}


def one_session(timing, scenario, rng):
    heart, brain = simulate(timing, scenario, rng)
    blocks = timing["blocks"]
    groups = {}
    for i, b in enumerate(blocks):
        groups.setdefault(b["length"], []).append(i)
    partner = {}
    for idx in groups.values():
        if len(idx) >= 2:
            p = emg.random_derangement(len(idx), rng)
            for a, q in zip(idx, p):
                partner[a] = idx[q]
    res = {"plain": (heart, brain),
           "resid": (remove_evoked(heart, blocks, False), remove_evoked(brain, blocks, False)),
           "residfit": (remove_evoked(heart, blocks, True), remove_evoked(brain, blocks, True))}
    rows = []
    for i, b in enumerate(blocks):
        if i not in partner:
            continue
        j = partner[i]
        n = min(int(round((b["t1"] - b["t0"]) * FS)), int(round((blocks[j]["t1"] - blocks[j]["t0"]) * FS)))
        r = {"length": b["length"]}
        same = band_values(H.cwt_slow(seg(heart, b)), H.cwt_slow(seg(brain, b)))
        aligned = band_values(H.cwt_slow(seg(heart, b)), H.cwt_slow(warp_to(brain, blocks[j], b)))
        for band in BANDS:
            r[f"same_{band}"] = same[band]
            r[f"excess_aligned_{band}"] = same[band] - aligned[band]
        for key, name in (("plain", "clock"), ("resid", "resid_clock"), ("residfit", "residfit_clock")):
            h, e = res[key]
            s_n = band_values(H.cwt_slow(seg(h, b, n)), H.cwt_slow(seg(e, b, n)))
            c_n = band_values(H.cwt_slow(seg(h, b, n)), H.cwt_slow(seg(e, blocks[j], n)))
            for band in BANDS:
                r[f"excess_{name}_{band}"] = s_n[band] - c_n[band]
                if key != "plain":
                    r[f"same_{key}_{band}"] = s_n[band]
        rows.append(r)
    out = {"n_blocks": len(rows)}
    for k in rows[0]:
        if k == "length":
            continue
        v = np.array([r[k] for r in rows], float)
        v = v[np.isfinite(v)]
        out[k] = float(v.mean())
        if k.startswith("excess"):
            out["p_" + k] = float(stats.ttest_1samp(v, 0.0, alternative="greater").pvalue)
    return out


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-sims", type=int, default=30)
    ap.add_argument("--seed", type=int, default=20261109)
    a = ap.parse_args(argv)
    cfg = load_json(PROJECT_ROOT / "config" / "download_config.json")
    root = bidsio.raw_root()
    outdir = PROJECT_ROOT / "results" / "null_benchmark"
    outdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    lines = [f"null benchmark | n_sims {a.n_sims} per participant and scenario | seed {a.seed} | python {platform.python_version()} "
             f"numpy {np.__version__} scipy {scipy.__version__}", "SIMULATED SIGNALS ON REAL TRIAL TIMING. No recording was read."]
    print("\n".join(lines), flush=True)
    sessions = []
    for si, sub in enumerate(cfg["pilot"]["subjects"]):
        timing = load_timing(root, sub)
        for ci, scenario in enumerate(SCENARIOS):
            for k in range(a.n_sims):
                rng = np.random.default_rng([a.seed, si, ci, k])
                r = one_session(timing, scenario, rng)
                r.update(sub=sub, scenario=scenario, sim=k)
                sessions.append(r)
            print(f"  {sub} {scenario}: {a.n_sims} sessions done ({time.time() - t0:.0f} s)", flush=True)
    cols = ["sub", "scenario", "sim", "n_blocks"] + sorted(k for k in sessions[0] if k not in ("sub", "scenario", "sim", "n_blocks"))
    with open(outdir / "null_benchmark_sessions.tsv", "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sessions:
            fh.write("\t".join(f"{r[c]:.5f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
    table = ["scenario\tband\tcontrol\tmean_same\tmean_excess\tsd_excess\tshare_sessions_p_below_0.05"]
    for scenario in SCENARIOS:
        ss = [r for r in sessions if r["scenario"] == scenario]
        for band in BANDS:
            for c in CONTROLS:
                e = np.array([r[f"excess_{c}_{band}"] for r in ss])
                p = np.array([r[f"p_excess_{c}_{band}"] for r in ss])
                same = np.mean([r[f"same_{band}"] for r in ss])
                table.append(f"{scenario}\t{band}\t{c}\t{same:.3f}\t{e.mean():+.4f}\t{e.std():.4f}\t{(p < 0.05).mean():.3f}")
    (outdir / "null_benchmark_summary.tsv").write_text("\n".join(table) + "\n", encoding="utf-8")
    print("\n".join(table), flush=True)
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    (PROJECT_ROOT / "logs").mkdir(exist_ok=True)
    (PROJECT_ROOT / "logs" / f"null_benchmark_{stamp}.log").write_text("\n".join(lines + table) + "\n", encoding="utf-8")
    append_execution_log([
        f"## {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} null benchmark (step 05b, simulation only)",
        f"- Command: python -m hb.null_benchmark {' '.join(sys.argv[1:])}",
        f"- Python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}, {platform.platform()}",
        f"- {len(sessions)} simulated sessions, seed {a.seed}, {time.time() - t0:.0f} s",
        f"- Outputs: results/null_benchmark/null_benchmark_summary.tsv, null_benchmark_sessions.tsv; log logs/null_benchmark_{stamp}.log",
    ])
    print(f"DONE in {time.time() - t0:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
