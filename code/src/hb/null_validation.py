"""Step 14 (review C1, C2 and C5, docs/REVIEW_RESPONSE_PLAN_v1.md): validation of the calibrated test on
synthetic data from a generative model that is INDEPENDENT of the surrogate construction.

For every real participant, a synthetic session is built on that participant's real block and trial timing
(public events files only; no recording is read). Heart and three EEG envelope like signals have NO ongoing
coupling. Each responds to every trial with:
  - an encoding response (trial start + 2.5 s to recall prompt) and a recall response that lasts the real,
    self paced recall period and grows with its duration;
  - trial to trial amplitude variability (gain sd 0.4) and latency jitter (sd 1.0 s, clipped at 2.5 s);
  - non stationary background noise (0.01 to 1 Hz, variance modulated over minutes).
Scenarios:
  independent    heart and EEG trial gains and latencies independent: the calibrated test must hold 5%
  shared_trial   heart and EEG trial gains correlated (0.7) and latencies shared: task driven covariation
                 with no ongoing coupling; shows how much of it the estimand counts as coupling
Each synthetic participant is analysed with the full frozen estimand and calibrated with 8 phase surrogates
built from its own synthetic series, exactly as step 12. Running hundreds of full cohorts is not feasible
(about 9 CPU hours each), so error rates are derived from one cohort per scenario: the mean and covariance
of the calibrated statistic across the 6 cells give, by Monte Carlo, the per cell error rate and the Holm
family wise error rate at n = 122, at the estimated bias and at the edge of its 95% interval.
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
from hb.cohort_stats import fmt, holm, participant_table
from hb.download_ds008302 import PROJECT_ROOT, append_execution_log
from hb.null_benchmark import gamma_kernel
from hb.null_benchmark_real import FS, prepare_parts, surrogate_session
from hb.pilot_coherence import BANDS, EEG, all_blocks
from hb.pilot_qc import write_tsv

SEED = 20261114
SCENARIOS = ["independent", "shared_trial"]
TAU = {"heart": 3.0, "theta": 1.0, "alpha": 1.5, "beta": 2.0}
EVOKED_SHARE = {"heart": 0.30, "theta": 0.20, "alpha": 0.20, "beta": 0.20}
GAIN_SD, LAT_SD, LAT_MAX, RHO = 0.4, 1.0, 2.5, 0.7


def slow_modulation(n, rng):
    t = np.arange(n) / FS
    f = rng.uniform(0.002, 0.008, 3)
    m = sum(np.sin(2 * np.pi * fi * t + rng.uniform(0, 2 * np.pi)) for fi in f) / np.sqrt(1.5)
    return np.exp(0.35 * m)


def simulate(blocks10, n, scenario, rng):
    t = np.arange(n) / FS
    inblock = np.zeros(n, bool)
    trials = []
    for b in blocks10:
        inblock |= (t >= b["t0"]) & (t < b["t1"])
        for k in range(10):
            trials.append((b["start"][k], b["probe"][k], b["end"][k]))
    dur = np.array([e - p for _, p, e in trials])
    dmed = np.median(dur)
    z_h = rng.standard_normal((len(trials), 2))
    lat_h = np.clip(LAT_SD * rng.standard_normal(len(trials)), -LAT_MAX, LAT_MAX)
    out = {}
    for name in ["heart"] + EEG:
        if name == "heart":
            z, lat = z_h, lat_h
        elif scenario == "shared_trial":
            z = RHO * z_h + np.sqrt(1 - RHO ** 2) * rng.standard_normal((len(trials), 2))
            lat = lat_h
        else:
            z = rng.standard_normal((len(trials), 2))
            lat = np.clip(LAT_SD * rng.standard_normal(len(trials)), -LAT_MAX, LAT_MAX)
        sign = 1.0 if name == "heart" else -1.0
        box = np.zeros(n)
        for i, (s, p, e) in enumerate(trials):
            enc = (t >= s + 2.5 + lat[i]) & (t < p + lat[i])
            rec = (t >= p + lat[i]) & (t < e + lat[i])
            box += sign * (1 + GAIN_SD * z[i, 0]) * enc
            box -= sign * 0.7 * (1 + GAIN_SD * z[i, 1]) * np.sqrt(dur[i] / dmed) * rec
        ev = np.convolve(box, gamma_kernel(TAU[name]))[:n]
        ev = (ev - ev[inblock].mean()) / ev[inblock].std()
        noise = S.band_noise(n, FS, (0.01, 1.0), rng) * slow_modulation(n, rng)
        noise = (noise - noise[inblock].mean()) / noise[inblock].std()
        a = np.sqrt(EVOKED_SHARE[name])
        out[name] = a * ev + np.sqrt(1 - a ** 2) * noise
    return out


def analyse(sub, heart, brains, blocks, n):
    zeros_i = {e: np.zeros(n, bool) for e in brains}
    rows = []
    for e, b in brains.items():
        rows += coherence_tables(sub, heart, {e: b}, {e: zeros_i[e]}, blocks, [0.0] * len(blocks))
    return participant_table({sub: [{c: r.get(c, "") for c in COLS} for r in rows]})


def run_subject(args):
    root, sub, k = args
    out = PROJECT_ROOT / "derivatives" / "null_validation" / f"{sub}.tsv"
    if out.exists():
        return sub, "skipped (already done)"
    t0 = time.time()
    p = bidsio.rec_paths(Path(root), sub, "memory")
    dur = float(bidsio.sidecar(p["json"])["RecordingDuration"])
    blocks, _ = all_blocks(Path(root), sub)
    blocks10 = [b for b in blocks if b["n_trials"] == 10]
    n = int(dur * FS)
    sid = int(sub.split("-")[1])
    rows = []
    for si, scenario in enumerate(SCENARIOS):
        rng = np.random.default_rng([SEED, sid, si])
        sig = simulate(blocks10, n, scenario, rng)
        for r in analyse(sub, sig["heart"], {e: sig[e] for e in EEG}, blocks, n):
            r.update(scenario=scenario, kind="observed", k=0)
            rows.append(r)
        for ei, e in enumerate(EEG):
            parts = prepare_parts(sig["heart"], sig[e], blocks10)
            for kk in range(k):
                h, b = surrogate_session(parts, blocks10, "task_only", np.random.default_rng([SEED, sid, si, ei, kk, 1]))
                for r in analyse(sub, h, {e: b}, blocks, parts["n"]):
                    r.update(scenario=scenario, kind="surrogate", k=kk)
                    rows.append(r)
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = ["scenario", "kind", "k", "sub", "eeg", "variant", "band", "n_blocks", "mean_excess", "in_H1"]
    tmp = out.with_suffix(".part")
    write_tsv(tmp, fmt([{c: r.get(c, "") for c in cols} for r in rows]))
    tmp.replace(out)
    return sub, f"ok ({len(rows)} rows, {time.time() - t0:.0f} s)"


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def error_rates(D, n, rng, reps=20000, alpha=0.05):
    """D: participants x 6 statistic matrix (complete rows). Monte Carlo of cohorts of size n from MVN(mean, cov)."""
    mu, cov = D.mean(axis=0), np.cov(D, rowvar=False)
    sd = np.sqrt(np.diag(cov))
    half = stats.t.ppf(0.975, D.shape[0] - 1) * sd / np.sqrt(D.shape[0])
    worst = np.where(mu >= 0, mu + half, mu - half)
    res = {}
    for label, m in (("at_zero_bias", np.zeros_like(mu)), ("at_estimated_bias", mu), ("at_95pct_edge_of_bias", worst)):
        means = rng.multivariate_normal(m, cov / n, size=reps)
        tt = means / (sd / np.sqrt(n))
        p = 2 * stats.t.sf(np.abs(tt), n - 1)
        per_cell = (p < alpha).mean(axis=0)
        fam = np.mean([(holm(row) < alpha).any() for row in p])
        res[label] = (per_cell, fam)
    return mu, half, res


def summarize(log):
    outdir = PROJECT_ROOT / "results" / "null_validation"
    outdir.mkdir(parents=True, exist_ok=True)
    rows = []
    for f in sorted((PROJECT_ROOT / "derivatives" / "null_validation").glob("sub-*.tsv")):
        with open(f, encoding="utf-8", newline="") as fh:
            rows += list(csv.DictReader(fh, delimiter="\t"))
    cells = [(e, b) for e in EEG for b in BANDS]
    lines, table = [f"null validation summary | {len({r['sub'] for r in rows})} synthetic participants"], []
    rng = np.random.default_rng(SEED)
    for scenario in SCENARIOS:
        for v in ("full", "encret"):
            obs, sur = {}, {}
            for r in rows:
                if r["scenario"] != scenario or r["variant"] != v or not np.isfinite(num(r["mean_excess"])):
                    continue
                key = (r["sub"], r["eeg"], r["band"])
                if r["kind"] == "observed" and r["in_H1"] == "1":
                    obs[key] = num(r["mean_excess"])
                elif r["kind"] == "surrogate":
                    sur.setdefault(key, []).append(num(r["mean_excess"]))
            subs = sorted({k[0] for k in obs})
            Dcal, Draw = [], []
            for s in subs:
                if all((s, e, b) in obs and (s, e, b) in sur for e, b in cells):
                    Draw.append([obs[(s, e, b)] for e, b in cells])
                    Dcal.append([obs[(s, e, b)] - np.mean(sur[(s, e, b)]) for e, b in cells])
            for label, D in (("calibrated", np.array(Dcal)), ("uncalibrated", np.array(Draw))):
                n = D.shape[0]
                mu, half, res = error_rates(D, 122, rng)
                direct = [stats.ttest_1samp(D[:, j], 0.0).pvalue for j in range(len(cells))]
                for j, (e, b) in enumerate(cells):
                    table.append({"scenario": scenario, "variant": v, "statistic": label, "eeg": e, "band": b, "n_synthetic": n,
                                  "mean": mu[j], "ci95_halfwidth": half[j], "p_this_cohort": direct[j],
                                  "error_rate_at_zero_bias": res["at_zero_bias"][0][j], "error_rate_at_estimated_bias": res["at_estimated_bias"][0][j],
                                  "error_rate_at_95pct_edge": res["at_95pct_edge_of_bias"][0][j]})
                lines.append(f"{scenario:12s} {v:6s} {label:12s}: Holm family error {100 * res['at_estimated_bias'][1]:.1f}% at estimated bias, "
                             f"{100 * res['at_95pct_edge_of_bias'][1]:.1f}% at the 95% edge, {100 * res['at_zero_bias'][1]:.1f}% at zero bias; "
                             f"cell means " + ", ".join(f"{e} {b} {mu[j]:+.4f}" for j, (e, b) in enumerate(cells)))
    write_tsv(outdir / "null_validation_cells.tsv", fmt(table))
    for ln in lines:
        log(ln)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--subjects", nargs="*")
    ap.add_argument("--k", type=int, default=8)
    ap.add_argument("--workers", type=int, default=max(1, min(6, (os.cpu_count() or 2) - 1)))
    ap.add_argument("--summary-only", action="store_true")
    a = ap.parse_args(argv)
    root = bidsio.raw_root()
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    (PROJECT_ROOT / "logs").mkdir(exist_ok=True)
    logf = open(PROJECT_ROOT / "logs" / f"null_validation_{stamp}.log", "w", encoding="utf-8")

    def log(m):
        print(m, flush=True)
        logf.write(m + "\n")
        logf.flush()

    log(f"null validation {stamp} | k {a.k} workers {a.workers} seed {SEED} | SYNTHETIC SIGNALS ON REAL TIMING | python {platform.python_version()} numpy {np.__version__} scipy {scipy.__version__}")
    t0, failed = time.time(), 0
    if not a.summary_only:
        subs = a.subjects or sorted(f.name.split("_")[0] for f in (PROJECT_ROOT / "derivatives" / "cohort_coherence").glob("sub-*_blocks.tsv"))
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
        f"## {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} null validation (step 14, review C1, C2, C5; synthetic only)",
        f"- Command: python -m hb.null_validation {' '.join(sys.argv[1:])}; seed {SEED}; Python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}",
        f"- Failures: {failed}; {(time.time() - t0) / 60:.1f} min. Outputs: derivatives/null_validation, results/null_validation; log logs/null_validation_{stamp}.log",
    ])
    return 0 if failed == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
