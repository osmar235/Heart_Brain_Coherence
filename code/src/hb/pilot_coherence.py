"""Step 07: first coherence run on the pilot recordings. FEASIBILITY ONLY (docs/PILOT_COHERENCE_PLAN.md).

Approved by O. Pinto Neto on 2026-11-10 (D-009 and the plan, thresholds of 10%, GO criteria as written).

Per block of 10 trials: remove trial locked averages from the heart and EEG series, compute squared
wavelet coherence, and subtract the mean coherence obtained with the EEG residuals of the other blocks
of the same sequence length (same truncated duration for both terms). Reports block counts, exclusions,
cone of influence coverage and the pooled distribution of the excess.

Not done here, by design: any span against span+1 comparison, any p value, rest or video recordings.
Sequence length is used only to choose pairing partners and is not shown in the figures.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import platform
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy

from hb import bidsio, hbcoh as H
from hb.download_ds008302 import PROJECT_ROOT, append_execution_log, load_json
from hb.null_benchmark import remove_evoked
from hb.pilot_qc import write_tsv

FS = H.FS_SLOW
BANDS = {"LF": (0.04, 0.15), "HF": (0.15, 0.40)}
EEG = ["theta", "alpha", "beta"]
MAX_MARKED = 0.10          # share of NN intervals marked in a block
MAX_INTERP = 0.10          # share of EEG envelope samples interpolated in a block
MIN_COI = 0.50             # share of the band inside the cone of influence
SEED = 20261111
N_BOOT = 2000
CODE_NOTE = ("Code added 2026-11-10 by Claude: src/hb/pilot_coherence.py, tests/test_pilot_coherence.py, "
             "scripts/07_pilot_coherence.cmd. Tested on synthetic inputs in Claude's Linux workspace before delivery.")


def all_blocks(root: Path, sub: str) -> tuple[list[dict], float]:
    p = bidsio.rec_paths(root, sub, "memory")
    dur = float(bidsio.sidecar(p["json"])["RecordingDuration"])
    parsed = bidsio.parse_memory_events(bidsio.read_tsv(p["events"]), dur)
    out = []
    for bi, b in enumerate(parsed["blocks"]):
        tr = b["trials"]
        if not tr:
            continue
        d = dict(index=bi + 1, t0=b["t0"], t1=b["t1"], length=b["length"], n_trials=len(tr),
                 boundary=any(b["t0"] < q < b["t1"] for q in parsed["boundaries"]),
                 start=np.array([t["t_start"] for t in tr]), probe=np.array([t["t_probe"] for t in tr]),
                 end=np.array([t["t_start"] for t in tr[1:]] + [b["t1"]]))
        out.append(d)
    return out, dur


def marked_share(beats_tsv: Path, blocks: list[dict]) -> list[float]:
    t, m = [], []
    with open(beats_tsv, encoding="utf-8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["nn_s"] != "n/a":
                t.append(float(r["beat_time_s"])); m.append(r["interval_marked"] == "1")
    t, m = np.array(t), np.array(m)
    out = []
    for b in blocks:
        k = (t >= b["t0"]) & (t < b["t1"])
        out.append(float(m[k].mean()) if k.any() else 1.0)
    return out


def zscore_over(x, blocks):
    k = np.zeros(x.size, bool)
    for b in blocks:
        k[int(round(b["t0"] * FS)):int(round(b["t1"] * FS))] = True
    return (x - x[k].mean()) / x[k].std()


def analyse_pair(heart, brain, blocks, reasons):
    """blocks: all task blocks; reasons: list of exclusion strings ('' = keep) for this signal pair.
    Returns one row per block with excess per band."""
    keep = [i for i, r in enumerate(reasons) if not r]
    kb = [blocks[i] for i in keep]
    rows = [{"block": b["index"], "length": b["length"], "excluded": reasons[i]} for i, b in enumerate(blocks)]
    if len(kb) < 2:
        return rows
    hr = remove_evoked(zscore_over(heart, kb), kb, False)
    er = remove_evoked(zscore_over(brain, kb), kb, False)
    cache = {}

    def cwt(name, x, b, n):
        key = (name, b["index"], n)
        if key not in cache:
            i0 = int(round(b["t0"] * FS))
            cache[key] = H.cwt_slow(x[i0:i0 + n])
        return cache[key]

    def nlen(b):
        return int(round((b["t1"] - b["t0"]) * FS))

    for i in keep:
        b = blocks[i]
        partners = [blocks[j] for j in keep if j != i and blocks[j]["length"] == b["length"]]
        row = rows[i]
        if not partners:
            row["excluded"] = "no other block of the same length"
            continue
        full = H.coherence(cwt("h", hr, b, nlen(b)), cwt("e", er, b, nlen(b)))
        acc = {band: [] for band in BANDS}
        for band, rng in BANDS.items():
            v, frac = H.band_summary(full, rng, min_valid_fraction=MIN_COI)
            row[f"same_{band}"] = v
            row[f"coi_share_{band}"] = round(frac, 3)
        for pb in partners:
            n = min(nlen(b), nlen(pb))
            same = H.coherence(cwt("h", hr, b, n), cwt("e", er, b, n))
            cross = H.coherence(cwt("h", hr, b, n), cwt("e", er, pb, n))
            for band, rng in BANDS.items():
                s, _ = H.band_summary(same, rng, min_valid_fraction=MIN_COI)
                c, _ = H.band_summary(cross, rng, min_valid_fraction=MIN_COI)
                if np.isfinite(s) and np.isfinite(c):
                    acc[band].append(s - c)
        row["n_partners"] = len(partners)
        for band in BANDS:
            row[f"excess_{band}"] = float(np.mean(acc[band])) if acc[band] else float("nan")
    return rows


def boot_ci(v, rng):
    v = np.asarray(v, float)
    if v.size < 3:
        return float("nan"), float("nan")
    m = rng.choice(v, (N_BOOT, v.size), replace=True).mean(axis=1)
    return float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def run_participant(sub, heart, eeg, interp, blocks, marked, log):
    """eeg and interp: dicts by EEG signal name. Returns (block rows, summary rows)."""
    base = []
    for b, ms in zip(blocks, marked):
        r = ""
        if b["n_trials"] != 10:
            r = "fewer than 10 trials"
        elif b["boundary"]:
            r = "recording discontinuity inside block"
        elif ms > MAX_MARKED:
            r = f"more than {int(100 * MAX_MARKED)}% of NN intervals marked"
        base.append(r)
    brows, srows = [], []
    rng = np.random.default_rng(SEED)
    for name in eeg:
        reasons = list(base)
        for i, b in enumerate(blocks):
            if not reasons[i]:
                s = slice(int(round(b["t0"] * FS)), int(round(b["t1"] * FS)))
                if interp[name][s].mean() > MAX_INTERP:
                    reasons[i] = f"more than {int(100 * MAX_INTERP)}% of envelope samples interpolated"
        rows = analyse_pair(heart, eeg[name], blocks, reasons)
        for r, ms in zip(rows, marked):
            r.update(sub=sub, eeg=name, nn_marked_share=round(ms, 4))
        brows += rows
        for band in BANDS:
            ok = [r for r in rows if not r["excluded"] and np.isfinite(r.get(f"excess_{band}", np.nan))]
            v = np.array([r[f"excess_{band}"] for r in ok])
            lo, hi = boot_ci(v, rng)
            srows.append({
                "sub": sub, "eeg": name, "band": band, "blocks_total": len(rows), "blocks_retained": len(ok),
                "pct_retained": round(100 * len(ok) / len(rows), 1),
                "mean_coi_share": round(float(np.mean([r[f"coi_share_{band}"] for r in rows if f"coi_share_{band}" in r])), 3) if ok else "",
                "mean_same_block_coherence": round(float(np.mean([r[f"same_{band}"] for r in ok])), 4) if ok else "",
                "mean_excess": round(float(v.mean()), 4) if v.size else "", "median_excess": round(float(np.median(v)), 4) if v.size else "",
                "boot95_low": round(lo, 4), "boot95_high": round(hi, 4),
                "share_blocks_positive": round(float((v > 0).mean()), 3) if v.size else "",
            })
        log(f"  {sub} {name}: retained " + ", ".join(f"{s['band']} {s['blocks_retained']}/{s['blocks_total']}" for s in srows[-2:]))
    return brows, srows


def go_check(srows):
    lines, passing = [], 0
    for sub in sorted({s["sub"] for s in srows}):
        good = [e for e in EEG if all(s["pct_retained"] >= 80 for s in srows if s["sub"] == sub and s["eeg"] == e)]
        coi = [s["mean_coi_share"] for s in srows if s["sub"] == sub and s["band"] == "LF" and s["mean_coi_share"] != ""]
        lf_ok = bool(coi) and float(np.mean(coi)) >= 0.5
        ok = len(good) >= 2 and lf_ok
        passing += ok
        lines.append(f"{sub}: EEG signals with at least 80% of blocks retained: {len(good)} ({', '.join(good) or 'none'}); "
                     f"mean LF share inside cone of influence {float(np.mean(coi)) if coi else float('nan'):.2f}; passes: {ok}")
    n = len(lines)
    verdict = "GO criteria met" if passing == n else ("NO-GO criteria met" if passing < 2 else "criteria met by some participants only: owners decide between GO and MODIFY")
    return lines, f"{passing} of {n} participants pass. {verdict}. The excess values are not part of these criteria."


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--subjects", nargs="*")
    a = ap.parse_args(argv)
    cfg = load_json(PROJECT_ROOT / "config" / "download_config.json")
    root = bidsio.raw_root()
    outdir = PROJECT_ROOT / "results" / "pilot_coherence"
    outdir.mkdir(parents=True, exist_ok=True)
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    (PROJECT_ROOT / "logs").mkdir(exist_ok=True)
    logf = open(PROJECT_ROOT / "logs" / f"pilot_coherence_{stamp}.log", "w", encoding="utf-8")

    def log(m):
        print(m, flush=True)
        logf.write(m + "\n")
        logf.flush()

    t0 = time.time()
    log(f"pilot coherence {stamp} | FEASIBILITY ONLY | python {platform.python_version()} numpy {np.__version__} scipy {scipy.__version__}")
    brows, srows = [], []
    for sub in a.subjects or cfg["pilot"]["subjects"]:
        blocks, _ = all_blocks(root, sub)
        z = np.load(PROJECT_ROOT / "derivatives" / "slow" / f"{sub}_task-memory_slow.npz")
        marked = marked_share(PROJECT_ROOT / "derivatives" / "beats" / f"{sub}_task-memory_beats.tsv", blocks)
        b, s = run_participant(sub, z["nn"], {e: z[e] for e in EEG if e in z}, {e: z[e + "_interpolated"] for e in EEG if e in z}, blocks, marked, log)
        brows += b
        srows += s
    cols = ["sub", "eeg", "block", "length", "excluded", "nn_marked_share", "n_partners"] + [f"{k}_{b}" for b in BANDS for k in ("coi_share", "same", "excess")]
    write_tsv(outdir / "pilot_coherence_blocks.tsv", [{c: r.get(c, "") for c in cols} for r in brows])
    write_tsv(outdir / "pilot_coherence_summary.tsv", srows)
    lines, verdict = go_check(srows)
    (outdir / "go_criteria_check.txt").write_text("\n".join(lines + [verdict]) + "\n", encoding="utf-8")

    subs = sorted({s["sub"] for s in srows})
    fig, axes = plt.subplots(1, len(subs), figsize=(5.2 * len(subs), 4.6), sharey=True, squeeze=False)
    rng = np.random.default_rng(1)
    for ax, sub in zip(axes[0], subs):
        k = 0
        ticks = []
        for e in EEG:
            for band in BANDS:
                v = np.array([r[f"excess_{band}"] for r in brows if r["sub"] == sub and r["eeg"] == e and not r["excluded"]
                              and np.isfinite(r.get(f"excess_{band}", np.nan))])
                s = next(x for x in srows if x["sub"] == sub and x["eeg"] == e and x["band"] == band)
                ax.plot(k + rng.uniform(-0.15, 0.15, v.size), v, ".", color="0.5", ms=5)
                if v.size:
                    ax.errorbar(k, s["mean_excess"], yerr=[[s["mean_excess"] - s["boot95_low"]], [s["boot95_high"] - s["mean_excess"]]],
                                fmt="o", color="C3", capsize=4)
                ticks.append(f"{e}\n{band}")
                k += 1
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks(range(k))
        ax.set_xticklabels(ticks, fontsize=8)
        ax.set_title(sub)
    axes[0, 0].set_ylabel("excess squared coherence per block\n(own block minus other blocks of same length)")
    fig.suptitle("Pilot feasibility run: grey = blocks, red = mean with bootstrap 95% interval. Not a hypothesis test.", fontsize=10)
    fig.tight_layout()
    fig.savefig(outdir / "pilot_coherence_excess.png", dpi=120)
    plt.close(fig)

    for s in srows:
        log(f"  {s['sub']} {s['eeg']} {s['band']}: retained {s['blocks_retained']}/{s['blocks_total']}, COI share {s['mean_coi_share']}, "
            f"same block coherence {s['mean_same_block_coherence']}, mean excess {s['mean_excess']} [{s['boot95_low']}, {s['boot95_high']}]")
    for ln in lines + [verdict]:
        log(ln)
    append_execution_log([
        f"## {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} pilot coherence, feasibility only (step 07)",
        f"- {CODE_NOTE}",
        f"- Command: python -m hb.pilot_coherence {' '.join(sys.argv[1:])}",
        f"- Python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}, {platform.platform()}; bootstrap seed {SEED}",
        f"- {verdict}",
        f"- Outputs: results/pilot_coherence/pilot_coherence_summary.tsv, pilot_coherence_blocks.tsv, go_criteria_check.txt, pilot_coherence_excess.png; log logs/pilot_coherence_{stamp}.log ({time.time() - t0:.0f} s)",
    ])
    log("DONE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
