"""Step 09: group statistics exactly as frozen in docs/ANALYSIS_FREEZE_v1.md. Reads only the block tables.

H1 (existence): participant mean excess coherence differs from zero.
H2 (condition): participant mean excess in span+1 blocks minus span blocks differs from zero.
Each hypothesis is tested in six cells (theta, alpha, beta x LF, HF) per variant, two sided one sample
t test with Holm correction over the six cells; Wilcoxon signed rank reported as sensitivity.
The "full" variant is the primary family, "encret" the secondary family.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import platform
import sys
from pathlib import Path

import numpy as np
import scipy
from scipy import stats

from hb.cohort import VARIANTS
from hb.download_ds008302 import PROJECT_ROOT, append_execution_log, load_json
from hb.pilot_coherence import BANDS, EEG
from hb.pilot_qc import write_tsv

MIN_BLOCKS_H1 = 6
MIN_BLOCKS_H2 = 3


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def span_of(block_rows):
    """Individual span from the blocks visited: the lower of two adjacent lengths that each have at least
    six blocks. With several such pairs, the one with most blocks; ties go to the pair reached later."""
    seen = {}
    for r in block_rows:
        seen.setdefault(int(r["block"]), (int(r["length"]), fnum(r["t0_s"])))
    count, last = {}, {}
    for _, (L, t) in seen.items():
        count[L] = count.get(L, 0) + 1
        last[L] = max(last.get(L, 0.0), t)
    pairs = [(count[L] + count[L + 1], max(last[L], last[L + 1]), L) for L in count if count[L] >= 6 and count.get(L + 1, 0) >= 6]
    return max(pairs)[2] if pairs else None


def participant_table(tables: dict) -> list[dict]:
    out = []
    for sub, rows in sorted(tables.items()):
        span = span_of(rows)
        for e in EEG:
            er = [r for r in rows if r["eeg"] == e and not r["excluded"]]
            for v in VARIANTS:
                for band in BANDS:
                    vals = [(int(r["length"]), fnum(r[f"excess_{v}_{band}"])) for r in er]
                    vals = [(L, x) for L, x in vals if np.isfinite(x)]
                    x = np.array([q for _, q in vals])
                    row = {"sub": sub, "eeg": e, "variant": v, "band": band, "span": span if span is not None else "",
                           "n_blocks": x.size, "mean_excess": float(x.mean()) if x.size else float("nan"),
                           "in_H1": int(x.size >= MIN_BLOCKS_H1)}
                    a = np.array([q for L, q in vals if span is not None and L == span])
                    b = np.array([q for L, q in vals if span is not None and L == span + 1])
                    row.update(n_span=a.size, n_span1=b.size, excess_span=float(a.mean()) if a.size else float("nan"),
                               excess_span1=float(b.mean()) if b.size else float("nan"))
                    row["in_H2"] = int(a.size >= MIN_BLOCKS_H2 and b.size >= MIN_BLOCKS_H2)
                    row["diff_span1_minus_span"] = row["excess_span1"] - row["excess_span"] if row["in_H2"] else float("nan")
                    out.append(row)
    return out


def holm(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    adj = np.empty(p.size)
    run = 0.0
    for rank, i in enumerate(order):
        run = max(run, (p.size - rank) * p[i])
        adj[i] = min(1.0, run)
    return adj


def group_tests(ptab: list[dict], exclude=()) -> list[dict]:
    out = []
    for hyp, col, flag in (("H1_existence", "mean_excess", "in_H1"), ("H2_span1_minus_span", "diff_span1_minus_span", "in_H2")):
        for v in VARIANTS:
            fam = []
            for e in EEG:
                for band in BANDS:
                    x = np.array([r[col] for r in ptab if r["variant"] == v and r["eeg"] == e and r["band"] == band
                                  and r[flag] and r["sub"] not in exclude and np.isfinite(r[col])])
                    row = {"hypothesis": hyp, "variant": v, "family": "primary" if v == "full" else "secondary", "eeg": e, "band": band, "n": x.size}
                    if x.size >= 5:
                        t = stats.ttest_1samp(x, 0.0)
                        sem = x.std(ddof=1) / np.sqrt(x.size)
                        h = stats.t.ppf(0.975, x.size - 1) * sem
                        try:
                            w = stats.wilcoxon(x).pvalue
                        except ValueError:
                            w = float("nan")
                        row.update(mean=float(x.mean()), ci95_low=float(x.mean() - h), ci95_high=float(x.mean() + h),
                                   sd=float(x.std(ddof=1)), cohen_d=float(x.mean() / x.std(ddof=1)), t=float(t.statistic),
                                   p=float(t.pvalue), p_wilcoxon=float(w), share_positive=float((x > 0).mean()))
                    else:
                        row.update(p=float("nan"))
                    fam.append(row)
            ok = [r for r in fam if np.isfinite(r["p"])]
            for r, a in zip(ok, holm([r["p"] for r in ok])):
                r["p_holm"] = float(a)
            out += fam
    return out


def load_tables(cdir: Path) -> dict:
    tables = {}
    for f in sorted(cdir.glob("sub-*_blocks.tsv")):
        with open(f, encoding="utf-8", newline="") as fh:
            tables[f.name.split("_")[0]] = list(csv.DictReader(fh, delimiter="\t"))
    return tables


def fmt(rows):
    return [{k: (round(v, 5) if isinstance(v, float) else v) for k, v in r.items()} for r in rows]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.parse_args(argv)
    cfg = load_json(PROJECT_ROOT / "config" / "download_config.json")
    tables = load_tables(PROJECT_ROOT / "derivatives" / "cohort_coherence")
    if len(tables) < 10:
        print(f"Only {len(tables)} participant tables found. Run scripts\\08_run_cohort.cmd first.")
        return 3
    outdir = PROJECT_ROOT / "results" / "cohort_stats"
    outdir.mkdir(parents=True, exist_ok=True)
    ptab = participant_table(tables)
    write_tsv(outdir / "participant_level.tsv", fmt(ptab))
    main_tests = group_tests(ptab)
    write_tsv(outdir / "group_tests.tsv", fmt(main_tests))
    write_tsv(outdir / "group_tests_without_pilot_participants.tsv", fmt(group_tests(ptab, exclude=set(cfg["pilot"]["subjects"]))))
    excl = {}
    for rows in tables.values():
        for r in rows:
            key = (r["eeg"], r["excluded"] or "retained")
            excl[key] = excl.get(key, 0) + 1
    write_tsv(outdir / "block_exclusions.tsv", [{"eeg": k[0], "status": k[1], "n_block_rows": v} for k, v in sorted(excl.items())])
    spans = [span_of(r) for r in tables.values()]
    lines = [f"cohort statistics | {len(tables)} participants with tables | span determined for {sum(s is not None for s in spans)} | "
             f"python {platform.python_version()} numpy {np.__version__} scipy {scipy.__version__}"]
    for r in main_tests:
        if "mean" in r:
            lines.append(f"{r['hypothesis']:22s} {r['family']:9s} {r['eeg']:5s} {r['band']}  n={r['n']:3d}  mean {r['mean']:+.4f} "
                         f"[{r['ci95_low']:+.4f}, {r['ci95_high']:+.4f}]  d={r['cohen_d']:+.2f}  p={r['p']:.4g}  Holm p={r['p_holm']:.4g}  Wilcoxon p={r['p_wilcoxon']:.4g}")
    print("\n".join(lines))
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    (PROJECT_ROOT / "logs").mkdir(exist_ok=True)
    (PROJECT_ROOT / "logs" / f"cohort_stats_{stamp}.log").write_text("\n".join(lines) + "\n", encoding="utf-8")
    append_execution_log([
        f"## {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} cohort statistics (step 09)",
        f"- Command: python -m hb.cohort_stats; Python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}",
        f"- Participants with block tables: {len(tables)}",
        f"- Outputs: results/cohort_stats/participant_level.tsv, group_tests.tsv, group_tests_without_pilot_participants.tsv, block_exclusions.tsv; log logs/cohort_stats_{stamp}.log",
    ])
    return 0


if __name__ == "__main__":
    sys.exit(main())
