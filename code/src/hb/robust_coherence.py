"""Step 11 (robustness checks 1 and 3, docs/ROBUSTNESS_PLAN_v1.md): frozen estimand and statistics on
  cfa  envelopes with the cardiac field artifact subtracted (step 10)
  edr  frozen envelopes, with ECG derived respiration (lags -2 to +2 s) regressed out of both residual
       series inside every block before coherence
Everything else is identical to the frozen analysis: blocks, exclusions, partners, variants, tests.
Also reports, per participant, the HF coherence between the heart series and the EDR (validity of EDR),
and across participants the rank correlation between the frozen excess and the size of the cardiac
field artifact (descriptive).
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import platform
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import scipy
from scipy import stats

from hb import bidsio, hbcoh as H
from hb.cohort import COLS, coherence_tables
from hb.cohort_stats import fmt, group_tests, participant_table
from hb.download_ds008302 import PROJECT_ROOT, append_execution_log
from hb.pilot_coherence import BANDS, EEG, MAX_MARKED, all_blocks, marked_share
from hb.pilot_qc import write_tsv

ROBUST = ["cfa", "edr"]
FS = H.FS_SLOW


def rsa_check(nn, edr, blocks, marked):
    """Mean HF squared coherence between heart series and EDR over usable 10 trial blocks."""
    v = []
    for b, ms in zip(blocks, marked):
        if b["n_trials"] != 10 or b["boundary"] or ms > MAX_MARKED:
            continue
        s = slice(int(round(b["t0"] * FS)), min(int(round(b["t1"] * FS)), nn.size, edr.size))
        try:
            coh = H.coherence(H.cwt_slow(nn[s]), H.cwt_slow(edr[s]))
            v.append(H.band_summary(coh, BANDS["HF"])[0])
        except ValueError:
            continue
    v = [x for x in v if np.isfinite(x)]
    return float(np.mean(v)) if v else float("nan")


def run_subject(root, sub):
    odir = PROJECT_ROOT / "derivatives" / "robust"
    odir.mkdir(parents=True, exist_ok=True)
    outs = {v: odir / f"{sub}_{v}_blocks.tsv" for v in ROBUST}
    if all(p.exists() for p in outs.values()):
        return "skipped (already done)", None
    blocks, _ = all_blocks(root, sub)
    fz = np.load(PROJECT_ROOT / "derivatives" / "slow" / f"{sub}_task-memory_slow.npz")
    cz = np.load(PROJECT_ROOT / "derivatives" / "slow_cfa" / f"{sub}_task-memory_slow.npz")
    marked = marked_share(PROJECT_ROOT / "derivatives" / "beats" / f"{sub}_task-memory_beats.tsv", blocks)
    rows = {
        "cfa": coherence_tables(sub, fz["nn"], {e: cz[e] for e in EEG}, {e: cz[e + "_interpolated"] for e in EEG}, blocks, marked),
        "edr": coherence_tables(sub, fz["nn"], {e: fz[e] for e in EEG}, {e: fz[e + "_interpolated"] for e in EEG}, blocks, marked,
                                covariate=cz["edr"]),
    }
    for v, r in rows.items():
        tmp = outs[v].with_suffix(".part")
        write_tsv(tmp, [{c: x.get(c, "") for c in COLS} for x in r])
        tmp.replace(outs[v])
    return "ok", {"sub": sub, "rsa_hf_coherence_nn_edr": round(rsa_check(fz["nn"], cz["edr"], blocks, marked), 4)}


def load_variant(v):
    tables = {}
    for f in sorted((PROJECT_ROOT / "derivatives" / "robust").glob(f"sub-*_{v}_blocks.tsv")):
        with open(f, encoding="utf-8", newline="") as fh:
            tables[f.name.split("_")[0]] = list(csv.DictReader(fh, delimiter="\t"))
    return tables


def read_tsv(p):
    with open(p, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def paired(frozen, robust, col, flag):
    """Paired comparison frozen minus robust per cell, over participants included in both."""
    out = []
    key = lambda r: (r["variant"], r["eeg"], r["band"])
    fz = {(r["sub"],) + key(r): r for r in frozen}
    for r in robust:
        f = fz.get((r["sub"],) + key(r))
        if f is None:
            continue
        a, b = num(f[col]), num(r[col])
        if int(num(f[flag]) or 0) and int(num(r[flag]) or 0) and np.isfinite(a) and np.isfinite(b):
            out.append((key(r), a, b))
    res = []
    for k in sorted({k for k, _, _ in out}):
        x = np.array([(a, b) for kk, a, b in out if kk == k])
        d = x[:, 0] - x[:, 1]
        t = stats.ttest_1samp(d, 0.0)
        h = stats.t.ppf(0.975, d.size - 1) * d.std(ddof=1) / np.sqrt(d.size)
        res.append({"variant": k[0], "eeg": k[1], "band": k[2], "n": d.size, "mean_frozen": x[:, 0].mean(), "mean_robust": x[:, 1].mean(),
                    "frozen_minus_robust": d.mean(), "ci95_low": d.mean() - h, "ci95_high": d.mean() + h, "p": float(t.pvalue),
                    "share_of_frozen_effect_remaining": x[:, 1].mean() / x[:, 0].mean() if x[:, 0].mean() != 0 else float("nan")})
    return res


def stats_all(log):
    outdir = PROJECT_ROOT / "results" / "robustness"
    outdir.mkdir(parents=True, exist_ok=True)
    frozen_p = read_tsv(PROJECT_ROOT / "results" / "cohort_stats" / "participant_level.tsv")
    lines = []
    for v in ROBUST:
        tables = load_variant(v)
        if len(tables) < 10:
            log(f"{v}: only {len(tables)} participant tables, statistics skipped")
            continue
        ptab = participant_table(tables)
        write_tsv(outdir / f"{v}_participant_level.tsv", fmt(ptab))
        tests = group_tests(ptab)
        write_tsv(outdir / f"{v}_group_tests.tsv", fmt(tests))
        comp = paired(frozen_p, [{k: str(x) for k, x in r.items()} for r in ptab], "mean_excess", "in_H1")
        write_tsv(outdir / f"{v}_vs_frozen_H1.tsv", fmt(comp))
        lines.append(f"== {v}: {len(tables)} participants")
        for r in tests:
            if "mean" in r:
                lines.append(f"{v:4s} {r['hypothesis']:20s} {r['variant']:6s} {r['eeg']:5s} {r['band']}  n={r['n']:3d}  mean {r['mean']:+.4f} "
                             f"[{r['ci95_low']:+.4f}, {r['ci95_high']:+.4f}]  Holm p={r['p_holm']:.4g}")
        for r in comp:
            lines.append(f"{v:4s} frozen minus {v} {r['variant']:6s} {r['eeg']:5s} {r['band']}  {r['frozen_minus_robust']:+.4f} "
                         f"[{r['ci95_low']:+.4f}, {r['ci95_high']:+.4f}]  share remaining {r['share_of_frozen_effect_remaining']:.2f}")
    qc = []
    for f in sorted((outdir / "qc").glob("sub-*_reprocess_qc.tsv")):
        qc += read_tsv(f)
    if qc:
        amp = {r["sub"]: r for r in qc}
        corr = []
        for r in frozen_p:
            if r["variant"] != "full" or not int(num(r["in_H1"]) or 0):
                continue
            a = num(amp.get(r["sub"], {}).get(f"cfa_p2p_before_{r['eeg']}_uV"))
            if np.isfinite(a) and np.isfinite(num(r["mean_excess"])):
                corr.append((r["eeg"], r["band"], a, num(r["mean_excess"])))
        rows = []
        for e in EEG:
            for band in BANDS:
                x = np.array([(a, m) for ee, bb, a, m in corr if ee == e and bb == band])
                if len(x) > 10:
                    rho = stats.spearmanr(x[:, 0], x[:, 1])
                    rows.append({"eeg": e, "band": band, "n": len(x), "spearman_rho": float(rho.statistic), "p": float(rho.pvalue)})
                    lines.append(f"descriptive: frozen excess vs CFA size, {e} {band}: rho {rho.statistic:+.3f} (p {rho.pvalue:.3g}, n {len(x)})")
        write_tsv(outdir / "excess_vs_cfa_size.tsv", fmt(rows))
        red = {b: [num(r.get(f"cfa_reduction_{b}_pct")) for r in qc] for b in ("theta", "alpha", "beta")}
        for b, v in red.items():
            v = [x for x in v if np.isfinite(x)]
            if v:
                lines.append(f"CFA reduction {b}: median {np.median(v):.0f}% (range {min(v):.0f} to {max(v):.0f}%)")
    rsa = read_tsv(outdir / "rsa_check.tsv") if (outdir / "rsa_check.tsv").exists() else []
    if rsa:
        v = [num(r["rsa_hf_coherence_nn_edr"]) for r in rsa]
        v = [x for x in v if np.isfinite(x)]
        lines.append(f"EDR validity: HF coherence between heart series and EDR, median {np.median(v):.3f} (IQR {np.percentile(v, 25):.3f} to {np.percentile(v, 75):.3f}, n {len(v)})")
    for ln in lines:
        log(ln)
    return lines


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--subjects", nargs="*")
    ap.add_argument("--stats-only", action="store_true")
    a = ap.parse_args(argv)
    root = bidsio.raw_root()
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    (PROJECT_ROOT / "logs").mkdir(exist_ok=True)
    logf = open(PROJECT_ROOT / "logs" / f"robust_coherence_{stamp}.log", "w", encoding="utf-8")

    def log(m):
        print(m, flush=True)
        logf.write(m + "\n")
        logf.flush()

    log(f"robustness coherence {stamp} | python {platform.python_version()} numpy {np.__version__} scipy {scipy.__version__}")
    t0, failed = time.time(), 0
    if not a.stats_only:
        subs = a.subjects or sorted(f.name.split("_")[0] for f in (PROJECT_ROOT / "derivatives" / "slow_cfa").glob("sub-*_task-memory_slow.npz"))
        rsa_path = PROJECT_ROOT / "results" / "robustness" / "rsa_check.tsv"
        rsa = {r["sub"]: r for r in read_tsv(rsa_path)} if rsa_path.exists() else {}
        for k, sub in enumerate(subs, 1):
            t1 = time.time()
            try:
                status, r = run_subject(root, sub)
                if r:
                    rsa[sub] = r
            except Exception as e:
                failed += 1
                status = "FAILED: " + "".join(traceback.format_exception(e)).strip().replace("\n", " | ")
            log(f"[{k}/{len(subs)}] {sub}: {status[:300]} ({time.time() - t1:.0f} s)")
        rsa_path.parent.mkdir(parents=True, exist_ok=True)
        write_tsv(rsa_path, [rsa[s] for s in sorted(rsa)])
    stats_all(log)
    append_execution_log([
        f"## {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} robustness coherence and statistics (step 11, checks 1 and 3)",
        f"- Command: python -m hb.robust_coherence {' '.join(sys.argv[1:])}; Python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}",
        f"- Failed participants: {failed}; {(time.time() - t0) / 60:.1f} min. Outputs: derivatives/robust, results/robustness; log logs/robust_coherence_{stamp}.log",
    ])
    return 0 if failed == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
