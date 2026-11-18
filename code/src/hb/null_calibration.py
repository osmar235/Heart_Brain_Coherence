"""Step 12 (robustness checks 2 and 4, docs/ROBUSTNESS_PLAN_v1.md): cohort scale surrogates.

For every participant, surrogate copies of the frozen 4 Hz series are built with the null_benchmark_real
construction (trial locked part kept with shuffled brain gains, residuals phase randomized independently
within each block), so heart and EEG share NO ongoing coupling. Each surrogate is analysed with the frozen
estimand, the participant's real block exclusions and real sequence lengths.
  scenario "none"   (default 16 copies)  null distribution of H1 and of H2
  scenario "equal"  (default 4 copies)   a shared 0.3 SD process in LF (0.07 to 0.13 Hz) and in HF
                                         (0.20 to 0.30 Hz), identical in span and span+1 blocks:
                                         H2 must stay at zero here
Statistics:
  check 2  per cell, observed participant excess minus the mean of that participant's "none" surrogates,
           paired t test across participants, Holm over the 6 cells of each variant
  check 4  the frozen H1 and H2 tests repeated on every surrogate cohort: share of cohorts significant
"""
from __future__ import annotations

import argparse
import csv
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
from hb import synth as S
from hb.cohort import COLS, coherence_tables
from hb.cohort_stats import fmt, group_tests, holm, participant_table
from hb.download_ds008302 import PROJECT_ROOT, append_execution_log
from hb.null_benchmark_real import prepare_parts, surrogate_session
from hb.pilot_coherence import BANDS, EEG, all_blocks, marked_share
from hb.pilot_qc import write_tsv

FS = 4.0
SEED = 20261112
PCOLS = ["scenario", "k", "sub", "eeg", "variant", "band", "span", "n_blocks", "mean_excess", "in_H1",
         "n_span", "n_span1", "excess_span", "excess_span1", "in_H2", "diff_span1_minus_span"]


def surrogate_rows(sub, nn, env, interp, blocks, marked, scenario, k):
    blocks10 = [b for b in blocks if b["n_trials"] == 10]
    rows = []
    for ei, e in enumerate(EEG):
        parts = prepare_parts(nn, env[e], blocks10)
        rng = np.random.default_rng([SEED, int(sub.split("-")[1]), 0 if scenario == "none" else 1, k, ei])
        h, b = surrogate_session(parts, blocks10, "task_only", rng)
        if scenario == "equal":
            n = parts["n"]
            s = 0.3 * S.band_noise(n, FS, (0.07, 0.13), rng) + 0.3 * S.band_noise(n, FS, (0.20, 0.30), rng)
            h, b = h + s, b + s
        rows += coherence_tables(sub, h, {e: b}, {e: interp[e][:parts["n"]]}, blocks, marked)
    ptab = participant_table({sub: [{c: r.get(c, "") for c in COLS} for r in rows]})
    for r in ptab:
        r.update(scenario=scenario, k=k)
    return ptab


def run_subject(args):
    root, sub, k_none, k_equal = args
    out = PROJECT_ROOT / "derivatives" / "null_calibration" / f"{sub}.tsv"
    if out.exists():
        return sub, "skipped (already done)"
    t0 = time.time()
    blocks, _ = all_blocks(Path(root), sub)
    z = np.load(PROJECT_ROOT / "derivatives" / "slow" / f"{sub}_task-memory_slow.npz")
    marked = marked_share(PROJECT_ROOT / "derivatives" / "beats" / f"{sub}_task-memory_beats.tsv", blocks)
    env = {e: z[e] for e in EEG}
    interp = {e: z[e + "_interpolated"] for e in EEG}
    rows = []
    for scenario, K in (("none", k_none), ("equal", k_equal)):
        for k in range(K):
            rows += surrogate_rows(sub, z["nn"], env, interp, blocks, marked, scenario, k)
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".part")
    write_tsv(tmp, fmt([{c: r.get(c, "") for c in PCOLS} for r in rows]))
    tmp.replace(out)
    return sub, f"ok ({len(rows)} rows, {time.time() - t0:.0f} s)"


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def load_rows(path):
    with open(path, encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        for c in ("mean_excess", "diff_span1_minus_span", "excess_span", "excess_span1"):
            r[c] = num(r.get(c))
        for c in ("in_H1", "in_H2", "n_blocks", "k"):          # the observed table has no "k" column
            r[c] = int(num(r.get(c))) if np.isfinite(num(r.get(c))) else 0
    return rows


def summarize(log):
    outdir = PROJECT_ROOT / "results" / "null_calibration"
    outdir.mkdir(parents=True, exist_ok=True)
    sur = []
    for f in sorted((PROJECT_ROOT / "derivatives" / "null_calibration").glob("sub-*.tsv")):
        sur += load_rows(f)
    obs = load_rows(PROJECT_ROOT / "results" / "cohort_stats" / "participant_level.tsv")
    lines = [f"null calibration summary | {len({r['sub'] for r in sur})} participants with surrogates"]

    # check 2: observed minus participant's own null mean
    calib = []
    for v in ("full", "encret"):
        fam = []
        for e in EEG:
            for band in BANDS:
                null = {}
                for r in sur:
                    if r["scenario"] == "none" and r["variant"] == v and r["eeg"] == e and r["band"] == band and np.isfinite(r["mean_excess"]):
                        null.setdefault(r["sub"], []).append(r["mean_excess"])
                pairs = [(r["mean_excess"], np.mean(null[r["sub"]])) for r in obs
                         if r["variant"] == v and r["eeg"] == e and r["band"] == band and r["in_H1"] and r["sub"] in null and np.isfinite(r["mean_excess"])]
                x = np.array(pairs)
                d = x[:, 0] - x[:, 1]
                t = stats.ttest_1samp(d, 0.0)
                h = stats.t.ppf(0.975, d.size - 1) * d.std(ddof=1) / np.sqrt(d.size)
                fam.append({"variant": v, "eeg": e, "band": band, "n": d.size, "observed_mean": x[:, 0].mean(), "null_mean": x[:, 1].mean(),
                            "null_mean_ci95_halfwidth": stats.t.ppf(0.975, d.size - 1) * x[:, 1].std(ddof=1) / np.sqrt(d.size),
                            "calibrated_excess": d.mean(), "ci95_low": d.mean() - h, "ci95_high": d.mean() + h,
                            "cohen_d": d.mean() / d.std(ddof=1), "p": float(t.pvalue)})
        for r, a in zip(fam, holm([r["p"] for r in fam])):
            r["p_holm"] = float(a)
        calib += fam
    write_tsv(outdir / "check2_H1_calibrated.tsv", fmt(calib))
    for r in calib:
        lines.append(f"check 2 {r['variant']:6s} {r['eeg']:5s} {r['band']}: observed {r['observed_mean']:+.4f}, null {r['null_mean']:+.4f} "
                     f"(+/- {r['null_mean_ci95_halfwidth']:.4f}), calibrated {r['calibrated_excess']:+.4f} [{r['ci95_low']:+.4f}, {r['ci95_high']:+.4f}], Holm p {r['p_holm']:.3g}")

    # check 4: frozen tests on each surrogate cohort
    rate = []
    for scenario in ("none", "equal"):
        ks = sorted({r["k"] for r in sur if r["scenario"] == scenario})
        res = []
        for k in ks:
            tab = [r for r in sur if r["scenario"] == scenario and r["k"] == k]
            for t in group_tests(tab):
                t["k"] = k
                res.append(t)
        for hyp in ("H1_existence", "H2_span1_minus_span"):
            for v in ("full", "encret"):
                fam_hits = []
                for k in ks:
                    fam = [t for t in res if t["k"] == k and t["hypothesis"] == hyp and t["variant"] == v and "p_holm" in t]
                    fam_hits.append(any(t["p_holm"] < 0.05 for t in fam))
                for e in EEG:
                    for band in BANDS:
                        c = [t for t in res if t["hypothesis"] == hyp and t["variant"] == v and t["eeg"] == e and t["band"] == band and "mean" in t]
                        if not c:
                            continue
                        rate.append({"scenario": scenario, "hypothesis": hyp, "variant": v, "eeg": e, "band": band, "n_cohorts": len(c),
                                     "mean_of_group_means": float(np.mean([t["mean"] for t in c])),
                                     "sd_of_group_means": float(np.std([t["mean"] for t in c], ddof=1)) if len(c) > 1 else float("nan"),
                                     "share_p_below_0.05": float(np.mean([t["p"] < 0.05 for t in c])),
                                     "share_holm_below_0.05": float(np.mean([t["p_holm"] < 0.05 for t in c])),
                                     "share_cohorts_any_cell_holm_significant_in_family": float(np.mean(fam_hits)) if fam_hits else float("nan")})
    write_tsv(outdir / "check4_surrogate_cohort_tests.tsv", fmt(rate))
    for r in rate:
        if r["variant"] == "full":
            lines.append(f"check 4 {r['scenario']:5s} {r['hypothesis']:20s} {r['eeg']:5s} {r['band']}: group mean {r['mean_of_group_means']:+.4f} "
                         f"(sd {r['sd_of_group_means']:.4f}), p<0.05 in {100 * r['share_p_below_0.05']:.0f}% of {r['n_cohorts']} cohorts, "
                         f"any Holm hit in family {100 * r['share_cohorts_any_cell_holm_significant_in_family']:.0f}%")
    for ln in lines:
        log(ln)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--subjects", nargs="*")
    ap.add_argument("--k-none", type=int, default=16)
    ap.add_argument("--k-equal", type=int, default=4)
    ap.add_argument("--workers", type=int, default=max(1, min(6, (os.cpu_count() or 2) - 1)))
    ap.add_argument("--summary-only", action="store_true")
    a = ap.parse_args(argv)
    root = bidsio.raw_root()
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    (PROJECT_ROOT / "logs").mkdir(exist_ok=True)
    logf = open(PROJECT_ROOT / "logs" / f"null_calibration_{stamp}.log", "w", encoding="utf-8")

    def log(m):
        print(m, flush=True)
        logf.write(m + "\n")
        logf.flush()

    log(f"null calibration {stamp} | k_none {a.k_none} k_equal {a.k_equal} workers {a.workers} seed {SEED} | "
        f"python {platform.python_version()} numpy {np.__version__} scipy {scipy.__version__}")
    t0, failed = time.time(), 0
    if not a.summary_only:
        subs = a.subjects or sorted(f.name.split("_")[0] for f in (PROJECT_ROOT / "derivatives" / "cohort_coherence").glob("sub-*_blocks.tsv"))
        jobs = [(str(root), s, a.k_none, a.k_equal) for s in subs]
        done = 0
        with ProcessPoolExecutor(max_workers=a.workers) as ex:
            futs = {ex.submit(run_subject, j): j[1] for j in jobs}
            for f in as_completed(futs):
                done += 1
                try:
                    sub, status = f.result()
                except Exception as e:
                    failed += 1
                    sub, status = futs[f], "FAILED: " + "".join(traceback.format_exception(e)).strip().replace("\n", " | ")
                log(f"[{done}/{len(jobs)}] {sub}: {status[:300]} ({(time.time() - t0) / 60:.1f} min elapsed)")
    summarize(log)
    append_execution_log([
        f"## {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} null calibration (step 12, checks 2 and 4)",
        f"- Command: python -m hb.null_calibration {' '.join(sys.argv[1:])}; seed {SEED}; Python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}",
        f"- Failed participants: {failed}; {(time.time() - t0) / 60:.1f} min. Outputs: derivatives/null_calibration, results/null_calibration; log logs/null_calibration_{stamp}.log",
    ])
    return 0 if failed == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
