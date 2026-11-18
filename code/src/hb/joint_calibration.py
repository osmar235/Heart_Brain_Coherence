"""Step 13 (review C4 and C5, docs/REVIEW_RESPONSE_PLAN_v1.md): joint pipeline with two null models.

Per participant, the observed analysis is run ONCE with every adjustment applied together:
cardiac field artifact subtracted (step 10 envelopes) -> trial locked averages removed -> ECG derived
respiration regressed out within each block -> frozen block estimand. The same analysis, with the same
exclusions, covariate and partners, is run on two kinds of no coupling surrogates:
  phase  the construction of step 12 (heart and EEG residuals phase randomized within each block,
         brain trial gains shuffled), 16 copies
  shift  the REAL heart series unchanged; the EEG residual circularly shifted within each block by a
         random lag of at least 20 s from either end, brain trial gains shuffled, 16 copies. Each
         signal keeps its own non stationarity and amplitude distribution; only their alignment is broken.
Summary: calibrated H1 (observed minus the participant's own surrogate mean, t test across participants,
Holm over 6 cells) for each null; sensitivity to the number of surrogates (8 against 16), also for the
step 12 surrogates of the frozen analysis.
"""
from __future__ import annotations

import argparse
import datetime as dt
import os
import platform
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import scipy
from scipy import stats

from hb import bidsio
from hb.cohort import COLS, coherence_tables
from hb.cohort_stats import fmt, holm, participant_table
from hb.download_ds008302 import PROJECT_ROOT, append_execution_log
from hb.null_benchmark_real import FS, block_slices, prepare_parts, surrogate_session, synthesize
from hb.null_calibration import PCOLS, load_rows
from hb.pilot_coherence import BANDS, EEG, all_blocks, marked_share
from hb.pilot_qc import write_tsv

SEED = 20261113
MIN_SHIFT_S = 20.0


def shift_session(parts, blocks10, rng):
    """Real heart; EEG = trial locked model with shuffled gains + residual circularly shifted within each block."""
    n = parts["n"]
    heart = synthesize(n, blocks10, parts["heart_model"]) + parts["heart_resid"]
    brain = synthesize(n, blocks10, parts["brain_model"], rng)
    res = parts["brain_resid"]
    for s in block_slices(blocks10, n):
        seg = res[s]
        m = int(MIN_SHIFT_S * FS)
        lag = int(rng.integers(m, seg.size - m)) if seg.size > 2 * m + 1 else seg.size // 2
        brain[s] = brain[s] + np.roll(seg, lag)
    return heart, brain


def table(sub, heart, brain_by_eeg, interp, blocks, marked, covariate):
    rows = []
    for e, b in brain_by_eeg.items():
        n = min(heart.size, b.size)
        rows += coherence_tables(sub, heart[:n], {e: b[:n]}, {e: interp[e][:n]}, blocks, marked,
                                 covariate=None if covariate is None else covariate[:n])
    return participant_table({sub: [{c: r.get(c, "") for c in COLS} for r in rows]})


def run_subject(args):
    root, sub, k = args
    out = PROJECT_ROOT / "derivatives" / "joint_calibration" / f"{sub}.tsv"
    if out.exists():
        return sub, "skipped (already done)"
    t0 = time.time()
    blocks, _ = all_blocks(Path(root), sub)
    blocks10 = [b for b in blocks if b["n_trials"] == 10]
    fz = np.load(PROJECT_ROOT / "derivatives" / "slow" / f"{sub}_task-memory_slow.npz")
    cz = np.load(PROJECT_ROOT / "derivatives" / "slow_cfa" / f"{sub}_task-memory_slow.npz")
    marked = marked_share(PROJECT_ROOT / "derivatives" / "beats" / f"{sub}_task-memory_beats.tsv", blocks)
    nn, edr = fz["nn"], cz["edr"]
    env = {e: cz[e] for e in EEG}
    interp = {e: cz[e + "_interpolated"] for e in EEG}
    rows = []
    for r in table(sub, nn, env, interp, blocks, marked, edr):
        r.update(scenario="observed", k=0)
        rows.append(r)
    parts = {e: prepare_parts(nn, env[e], blocks10) for e in EEG}
    sid = int(sub.split("-")[1])
    for null_id, null in enumerate(("phase", "shift")):
        for kk in range(k):
            for ei, e in enumerate(EEG):
                rng = np.random.default_rng([SEED, sid, null_id, kk, ei])
                h, b = surrogate_session(parts[e], blocks10, "task_only", rng) if null == "phase" else shift_session(parts[e], blocks10, rng)
                for r in table(sub, h, {e: b}, interp, blocks, marked, edr):
                    r.update(scenario=null, k=kk)
                    rows.append(r)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".part")
    write_tsv(tmp, fmt([{c: r.get(c, "") for c in PCOLS} for r in rows]))
    tmp.replace(out)
    return sub, f"ok ({len(rows)} rows, {time.time() - t0:.0f} s)"


def calibrated(obs, null_rows, k_max=None):
    """Per variant and cell: observed participant excess minus the mean of its first k_max null rows."""
    out = []
    for v in ("full", "encret"):
        fam = []
        for e in EEG:
            for band in BANDS:
                null = {}
                for r in null_rows:
                    if r["variant"] == v and r["eeg"] == e and r["band"] == band and np.isfinite(r["mean_excess"]) and (k_max is None or r["k"] < k_max):
                        null.setdefault(r["sub"], []).append(r["mean_excess"])
                x = np.array([(r["mean_excess"], np.mean(null[r["sub"]])) for r in obs
                              if r["variant"] == v and r["eeg"] == e and r["band"] == band and r["in_H1"] and r["sub"] in null and np.isfinite(r["mean_excess"])])
                d = x[:, 0] - x[:, 1]
                t = stats.ttest_1samp(d, 0.0)
                h = stats.t.ppf(0.975, d.size - 1) * d.std(ddof=1) / np.sqrt(d.size)
                fam.append({"variant": v, "eeg": e, "band": band, "n": d.size, "observed_mean": x[:, 0].mean(), "null_mean": x[:, 1].mean(),
                            "calibrated_excess": d.mean(), "ci95_low": d.mean() - h, "ci95_high": d.mean() + h, "cohen_d": d.mean() / d.std(ddof=1),
                            "p": float(t.pvalue)})
        for r, a in zip(fam, holm([r["p"] for r in fam])):
            r["p_holm"] = float(a)
        out += fam
    return out


def summarize(log):
    outdir = PROJECT_ROOT / "results" / "joint_calibration"
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for f in sorted((PROJECT_ROOT / "derivatives" / "joint_calibration").glob("sub-*.tsv")):
        rows += load_rows(f)
    obs = [r for r in rows if r["scenario"] == "observed"]
    lines = [f"joint calibration summary | {len({r['sub'] for r in obs})} participants"]
    res = []
    for null in ("phase", "shift"):
        for kmax in (8, None):
            for r in calibrated(obs, [r for r in rows if r["scenario"] == null], kmax):
                r.update(analysis="joint", null=null, k_used=kmax or 16)
                res.append(r)
    step12 = []
    for f in sorted((PROJECT_ROOT / "derivatives" / "null_calibration").glob("sub-*.tsv")):
        step12 += [r for r in load_rows(f) if r["scenario"] == "none"]
    if step12:
        frozen_obs = load_rows(PROJECT_ROOT / "results" / "cohort_stats" / "participant_level.tsv")
        for kmax in (8, None):
            for r in calibrated(frozen_obs, step12, kmax):
                r.update(analysis="frozen (step 12 surrogates)", null="phase", k_used=kmax or 16)
                res.append(r)
    write_tsv(outdir / "joint_calibrated_H1.tsv", fmt(res))
    for r in res:
        if r["variant"] == "full":
            lines.append(f"{r['analysis']:28s} {r['null']:5s} K={r['k_used']:2d} {r['eeg']:5s} {r['band']}: observed {r['observed_mean']:+.4f} null {r['null_mean']:+.4f} "
                         f"calibrated {r['calibrated_excess']:+.4f} [{r['ci95_low']:+.4f}, {r['ci95_high']:+.4f}] Holm p {r['p_holm']:.3g}")
    for ln in lines:
        log(ln)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--subjects", nargs="*")
    ap.add_argument("--k", type=int, default=16)
    ap.add_argument("--workers", type=int, default=max(1, min(6, (os.cpu_count() or 2) - 1)))
    ap.add_argument("--summary-only", action="store_true")
    a = ap.parse_args(argv)
    root = bidsio.raw_root()
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    (PROJECT_ROOT / "logs").mkdir(exist_ok=True)
    logf = open(PROJECT_ROOT / "logs" / f"joint_calibration_{stamp}.log", "w", encoding="utf-8")

    def log(m):
        print(m, flush=True)
        logf.write(m + "\n")
        logf.flush()

    log(f"joint calibration {stamp} | k {a.k} workers {a.workers} seed {SEED} | python {platform.python_version()} numpy {np.__version__} scipy {scipy.__version__}")
    t0, failed = time.time(), 0
    if not a.summary_only:
        subs = a.subjects or sorted(f.name.split("_")[0] for f in (PROJECT_ROOT / "derivatives" / "slow_cfa").glob("sub-*_task-memory_slow.npz"))
        done = 0
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            futs = {ex.submit(run_subject, (str(root), s, a.k)): s for s in subs}
            for f in as_completed(futs):
                done += 1
                try:
                    sub, status = f.result()
                except Exception as e:
                    failed += 1
                    sub, status = futs[f], "FAILED: " + "".join(traceback.format_exception(e)).strip().replace("\n", " | ")
                log(f"[{done}/{len(subs)}] {sub}: {status[:300]} ({(time.time() - t0) / 60:.1f} min elapsed)")
    try:
        summarize(log)
    except Exception as e:
        failed += 1
        log("SUMMARY FAILED: " + "".join(traceback.format_exception(e)).strip())
    append_execution_log([
        f"## {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} joint calibration (step 13, review C4 and C5)",
        f"- Command: python -m hb.joint_calibration {' '.join(sys.argv[1:])}; seed {SEED}; Python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}",
        f"- Failures: {failed}; {(time.time() - t0) / 60:.1f} min. Outputs: derivatives/joint_calibration, results/joint_calibration; log logs/joint_calibration_{stamp}.log",
    ])
    return 0 if failed == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
