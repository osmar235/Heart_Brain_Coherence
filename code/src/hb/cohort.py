"""Step 08: whole cohort run, as frozen in docs/ANALYSIS_FREEZE_v1.md (GO decision D-010, 2026-11-10).

Per participant, resumable (finished participants are skipped):
  stage 1  beat correction (as step 03b) and slow 4 Hz series (as step 06) for the memory task
  stage 2  per block excess coherence in two prespecified variants:
             full    band mean over the whole block (the pilot estimand, D-009)
             encret  band mean restricted to the fixed part of each trial, from trial start to the
                     recall prompt (encoding and retention); recall periods are left out
Writes one block table per participant to derivatives/cohort_coherence. No group statistics here:
those are in hb.cohort_stats, which implements the frozen tests.
"""
from __future__ import annotations

import argparse
import datetime as dt
import gc
import platform
import sys
import time
import traceback
from pathlib import Path

import numpy as np
import scipy

from hb import beat_correction, bidsio, hbcoh as H, slow_signals
from hb.download_ds008302 import PROJECT_ROOT, append_execution_log
from hb.null_benchmark import remove_evoked
from hb.pilot_coherence import BANDS, EEG, MAX_INTERP, MAX_MARKED, MIN_COI, all_blocks, marked_share, zscore_over
from hb.pilot_qc import write_tsv

FS = H.FS_SLOW
VARIANTS = ["full", "encret"]


def band_mean(coh, band, cols):
    """Mean squared coherence in band over the time samples selected by cols, inside the cone of influence.
    Returns (value, share of selected samples that are valid); NaN when the share is below MIN_COI."""
    rows = (coh["freqs"] >= band[0]) & (coh["freqs"] < band[1])
    if not cols.any():
        return float("nan"), 0.0
    R = coh["Rsq"][np.ix_(rows, cols)]
    ok = coh["inside"][np.ix_(rows, cols)] & np.isfinite(R)
    frac = float(ok.mean())
    return (float(R[ok].mean()) if frac >= MIN_COI else float("nan")), frac


def masks(b, n):
    """Time masks on the block's own timeline: everything, and trial start to recall prompt only."""
    t = b["t0"] + np.arange(n) / FS
    enc = np.zeros(n, bool)
    for s, p in zip(b["start"], b["probe"]):
        enc |= (t >= s) & (t < p)
    return {"full": np.ones(n, bool), "encret": enc}


LAGS = tuple(range(-8, 9, 2))      # covariate lags in 4 Hz samples: -2 s to +2 s in 0.5 s steps


def regress_lagged(x, c, lags=LAGS):
    """Remove from x the part linearly predictable from c at the given lags (edge padded shifts)."""
    cols = []
    for lag in lags:
        if lag > 0:
            cols.append(np.r_[np.full(lag, c[0]), c[:-lag]])
        elif lag < 0:
            cols.append(np.r_[c[-lag:], np.full(-lag, c[-1])])
        else:
            cols.append(c)
    X = np.column_stack(cols + [np.ones(x.size)])
    beta, *_ = np.linalg.lstsq(X, x, rcond=None)
    return x - X @ beta


def analyse_pair(heart, brain, blocks, reasons, covariate=None):
    keep = [i for i, r in enumerate(reasons) if not r]
    kb = [blocks[i] for i in keep]
    rows = [{"block": b["index"], "length": b["length"], "n_trials": b["n_trials"], "t0_s": round(b["t0"], 1), "excluded": reasons[i]}
            for i, b in enumerate(blocks)]
    if len(kb) < 2:
        for i in keep:
            rows[i]["excluded"] = "fewer than two usable blocks"
        return rows
    hr = remove_evoked(zscore_over(heart, kb), kb, False)
    er = remove_evoked(zscore_over(brain, kb), kb, False)
    if covariate is not None:            # robustness check 3 only; the frozen analysis passes None
        for b in kb:
            s = slice(int(round(b["t0"] * FS)), min(int(round(b["t1"] * FS)), hr.size, er.size, covariate.size))
            c = covariate[s] - covariate[s].mean()
            c = c / (c.std() + 1e-12)
            hr[s] = regress_lagged(hr[s], c)
            er[s] = regress_lagged(er[s], c)
    cache = {}

    def cwt(name, x, b, n):
        key = (name, b["index"], n)
        if key not in cache:
            i0 = int(round(b["t0"] * FS))
            cache[key] = H.cwt_slow(x[i0:i0 + n])
        return cache[key]

    def nlen(b):
        return min(int(round((b["t1"] - b["t0"]) * FS)), heart.size - int(round(b["t0"] * FS)), brain.size - int(round(b["t0"] * FS)))

    for i in keep:
        b, row = blocks[i], rows[i]
        partners = [blocks[j] for j in keep if j != i and blocks[j]["length"] == b["length"]]
        if not partners:
            row["excluded"] = "no other block of the same length"
            continue
        n0 = nlen(b)
        full = H.coherence(cwt("h", hr, b, n0), cwt("e", er, b, n0))
        m0 = masks(b, n0)
        for v in VARIANTS:
            for band, rng in BANDS.items():
                val, frac = band_mean(full, rng, m0[v])
                row[f"same_{v}_{band}"] = val
                row[f"valid_share_{v}_{band}"] = round(frac, 3)
        acc = {(v, band): [] for v in VARIANTS for band in BANDS}
        for pb in partners:
            n = min(n0, nlen(pb))
            same = H.coherence(cwt("h", hr, b, n), cwt("e", er, b, n))
            cross = H.coherence(cwt("h", hr, b, n), cwt("e", er, pb, n))
            m = masks(b, n)
            for v in VARIANTS:
                for band, rng in BANDS.items():
                    s, _ = band_mean(same, rng, m[v])
                    c, _ = band_mean(cross, rng, m[v])
                    if np.isfinite(s) and np.isfinite(c):
                        acc[(v, band)].append(s - c)
        row["n_partners"] = len(partners)
        for (v, band), vals in acc.items():
            row[f"excess_{v}_{band}"] = float(np.mean(vals)) if vals else float("nan")
    return rows


def coherence_tables(sub, heart, eeg, interp, blocks, marked, covariate=None):
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
    out = []
    for name in eeg:
        reasons = list(base)
        for i, b in enumerate(blocks):
            if not reasons[i]:
                s = slice(int(round(b["t0"] * FS)), int(round(b["t1"] * FS)))
                if interp[name][s].size == 0 or interp[name][s].mean() > MAX_INTERP:
                    reasons[i] = f"more than {int(100 * MAX_INTERP)}% of envelope samples interpolated"
        rows = analyse_pair(heart, eeg[name], blocks, reasons, covariate)
        for r, ms in zip(rows, marked):
            r.update(sub=sub, eeg=name, nn_marked_share=round(ms, 4))
        out += rows
    return out


COLS = (["sub", "eeg", "block", "length", "n_trials", "t0_s", "excluded", "nn_marked_share", "n_partners"]
        + [f"{k}_{v}_{b}" for v in VARIANTS for b in BANDS for k in ("valid_share", "same", "excess")])


def run_subject(root: Path, sub: str, log) -> str:
    bdir = PROJECT_ROOT / "derivatives" / "beats"
    sdir = PROJECT_ROOT / "derivatives" / "slow"
    cdir = PROJECT_ROOT / "derivatives" / "cohort_coherence"
    qdir = PROJECT_ROOT / "results" / "cohort_qc"
    out = cdir / f"{sub}_blocks.tsv"
    if out.exists():
        return "skipped (already done)"
    p = bidsio.rec_paths(root, sub, "memory")
    if not p["set"].exists():
        return "MISSING recording"
    beats = bdir / f"{sub}_task-memory_beats.tsv"
    slow = sdir / f"{sub}_task-memory_slow.npz"
    if not beats.exists():
        row, _ = beat_correction.process(root, sub, "memory", qdir / "beats", bdir)
        write_tsv(qdir / "beats" / f"{sub}_beat_summary.tsv", [row])
        gc.collect()
    if not slow.exists():
        rows = slow_signals.process(root, sub, "memory", qdir / "slow", sdir, bdir, lambda m: None)
        write_tsv(qdir / "slow" / f"{sub}_slow_summary.tsv", rows)
        gc.collect()
    blocks, _ = all_blocks(root, sub)
    z = np.load(slow)
    marked = marked_share(beats, blocks)
    rows = coherence_tables(sub, z["nn"], {e: z[e] for e in EEG if e in z}, {e: z[e + "_interpolated"] for e in EEG if e in z}, blocks, marked)
    cdir.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".part")
    write_tsv(tmp, [{c: r.get(c, "") for c in COLS} for r in rows])
    tmp.replace(out)
    kept = sum(1 for r in rows if not r["excluded"])
    return f"ok ({len(blocks)} blocks, {kept} of {len(rows)} block and signal rows retained)"


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--subjects", nargs="*", help="default: every participant folder with a memory recording on disk")
    a = ap.parse_args(argv)
    root = bidsio.raw_root()
    subs = a.subjects or sorted(d.name for d in root.glob("sub-*") if bidsio.rec_paths(root, d.name, "memory")["set"].exists())
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    (PROJECT_ROOT / "logs").mkdir(exist_ok=True)
    logf = open(PROJECT_ROOT / "logs" / f"cohort_{stamp}.log", "w", encoding="utf-8")

    def log(m):
        print(m, flush=True)
        logf.write(m + "\n")
        logf.flush()

    log(f"cohort run {stamp} | {len(subs)} participants | python {platform.python_version()} numpy {np.__version__} scipy {scipy.__version__}")
    t0 = time.time()
    counts = {"ok": 0, "skipped": 0, "failed": 0}
    for k, sub in enumerate(subs, 1):
        t1 = time.time()
        try:
            status = run_subject(root, sub, log)
        except MemoryError:
            status = "FAILED: not enough memory for this recording"
        except Exception as e:
            status = "FAILED: " + "".join(traceback.format_exception(e)).strip().replace("\n", " | ")
        key = "ok" if status.startswith("ok") else ("skipped" if status.startswith("skipped") else "failed")
        counts[key] += 1
        log(f"[{k}/{len(subs)}] {sub}: {status} ({time.time() - t1:.0f} s)")
        gc.collect()
    log(f"finished: {counts} in {(time.time() - t0) / 60:.1f} min")
    append_execution_log([
        f"## {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} cohort run (step 08)",
        f"- Command: python -m hb.cohort {' '.join(sys.argv[1:])}",
        f"- Python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}, {platform.platform()}",
        f"- Participants on disk: {len(subs)}; result: {counts}; {(time.time() - t0) / 60:.1f} min",
        f"- Outputs: derivatives/beats, derivatives/slow, derivatives/cohort_coherence/*_blocks.tsv, results/cohort_qc; log logs/cohort_{stamp}.log",
    ])
    return 0 if counts["failed"] == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
