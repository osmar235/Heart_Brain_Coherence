"""Step 03b: apply the proposed beat correction policy (D-007) to the pilot recordings.

For each recording: detect R peaks, drop detections that do not match the participant's QRS shape,
classify the remaining intervals with adaptive thresholds, remove extra beats, insert single missed
beats, and mark abnormal intervals for interpolation. Writes derivatives/beats, a summary table, a
per block table and one before/after figure per recording. Feasibility only; no condition contrast.
"""
from __future__ import annotations

import argparse
import datetime as dt
import gc
import platform
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy

from hb import bidsio, cardiac
from hb.download_ds008302 import PROJECT_ROOT, append_execution_log, load_json
from hb.pilot_qc import interp_series, write_tsv


def process(root: Path, sub: str, task: str, outdir: Path, bdir: Path) -> tuple[dict, list[dict]]:
    p = bidsio.rec_paths(root, sub, task)
    rec = bidsio.load_set(p["set"], channels=["ECG", "PPG"])
    fs, dur = rec["fs"], rec["n_samples"] / rec["fs"]
    ecg, ppg = rec["data"][0].astype(float), rec["data"][1].astype(float)
    del rec
    gc.collect()
    det = cardiac.detect_rpeaks(ecg, fs)
    pk = det["peaks"]
    keep, corr = cardiac.template_screen(det["filtered"], pk, fs)
    res = cardiac.correct_beats(pk[keep] / fs)
    t, nn, art = res["times"], res["nn"], res["artifact"]
    row = {"sub": sub, "task": task, "detected": int(pk.size), "template_rejected": int((~keep).sum()),
           **res["counts"], "final_beats": int(t.size), "intervals_marked": int(art.sum()),
           "pct_intervals_marked": round(100 * float(art.mean()), 3),
           "pct_beats_changed": round(100 * ((~keep).sum() + res["counts"]["extra_removed"] + res["counts"]["missed_inserted"]) / pk.size, 3)}
    good = ~art
    row["mean_hr_bpm"] = round(60.0 / float(np.median(nn[good])), 1)
    row["rmssd_ms"] = round(float(1000 * np.sqrt(np.mean(np.diff(nn[good]) ** 2))), 1)
    row["sdnn_ms"] = round(float(1000 * nn[good].std()), 1)
    # longest stretch without any marked interval
    edges = np.r_[0.0, t[1:][art], dur]
    row["longest_unmarked_run_s"] = round(float(np.diff(edges).max()), 1)
    # agreement with PPG, before and after
    pp = cardiac.ppg_peaks(ppg, fs)
    grid = np.arange(5, dur - 5, 0.25)
    if pp.size > 10:
        b = interp_series(pp[1:] / fs, np.diff(pp) / fs, grid)
        raw = interp_series(pk[1:] / fs, np.diff(pk) / fs, grid)
        cor = interp_series(t[1:][good], nn[good], grid)
        for name, a in (("raw", raw), ("corrected", cor)):
            ok = np.isfinite(a) & np.isfinite(b)
            row[f"ppg_corr_{name}"] = round(float(np.corrcoef(a[ok], b[ok])[0, 1]), 3)
            row[f"ppg_median_abs_diff_ms_{name}"] = round(float(1000 * np.median(np.abs(a[ok] - b[ok]))), 1)

    blocks = []
    parsed = None
    if task == "memory":
        parsed = bidsio.parse_memory_events(bidsio.read_tsv(p["events"]), dur)
        for bi, bl in enumerate(parsed["blocks"]):
            if not bl["trials"]:
                continue
            inb = (t[1:] >= bl["t0"]) & (t[1:] < bl["t1"])
            blocks.append({"sub": sub, "block": bi + 1, "length": bl["length"], "n_intervals": int(inb.sum()),
                           "n_marked": int(art[inb].sum()), "pct_marked": round(100 * float(art[inb].mean()), 2) if inb.any() else "",
                           "template_rejected": int(((~keep) & (pk / fs >= bl["t0"]) & (pk / fs < bl["t1"])).sum())})

    bdir.mkdir(parents=True, exist_ok=True)
    with open(bdir / f"{sub}_task-{task}_beats.tsv", "w", encoding="utf-8", newline="\n") as fh:
        fh.write("beat_time_s\tnn_s\tinterval_marked\tbeat_inserted\n")
        fh.write(f"{t[0]:.3f}\tn/a\tn/a\t{int(res['inserted'][0])}\n")
        for i in range(nn.size):
            fh.write(f"{t[i + 1]:.3f}\t{nn[i]:.4f}\t{int(art[i])}\t{int(res['inserted'][i + 1])}\n")

    fig, axes = plt.subplots(3, 1, figsize=(15, 9), sharex=True)
    rr_raw = np.diff(pk) / fs
    lo, hi = np.percentile(nn[good], [0.5, 99.5])
    axes[0].plot(pk[1:] / fs, rr_raw, lw=0.5, color="0.3")
    axes[0].plot(pk[~keep] / fs, np.full((~keep).sum(), lo - 0.1), "|", color="C3", ms=10, label=f"rejected by QRS shape ({(~keep).sum()})")
    axes[0].set_title("before: all detections")
    axes[0].legend(fontsize=8, loc="upper right")
    axes[1].plot(t[1:][good], nn[good], lw=0.5, color="C0")
    axes[1].plot(t[1:][art], np.clip(nn[art], lo - 0.12, hi + 0.12), "x", color="C3", ms=4, label=f"marked intervals ({art.sum()})")
    axes[1].plot(t[res["inserted"]], np.full(res["inserted"].sum(), hi + 0.1), "v", color="C2", ms=4, label=f"inserted beats ({res['inserted'].sum()})")
    axes[1].set_title("after: corrected series (marked intervals are interpolated, never used as data)")
    axes[1].legend(fontsize=8, loc="upper right", ncol=2)
    if pp.size > 10:
        axes[2].plot(grid, 1000 * (cor - b), lw=0.5, color="C4")
        axes[2].set_ylim(-150, 150)
    axes[2].set_title("corrected ECG interval minus PPG pulse interval (ms)")
    axes[2].set_xlabel("time (s)")
    for ax in axes[:2]:
        ax.set_ylim(lo - 0.15, hi + 0.15)
        ax.set_ylabel("interval (s)")
        if parsed:
            for bl in parsed["blocks"]:
                if bl["trials"]:
                    ax.axvspan(bl["t0"], bl["t1"], color="0.9", lw=0)
    axes[0].set_xlim(0, dur)
    fig.suptitle(f"{sub} task-{task} | {row['pct_beats_changed']}% beats changed, {row['pct_intervals_marked']}% intervals marked, "
                 f"RMSSD {row['rmssd_ms']} ms", fontsize=11)
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{sub}_task-{task}_beats.png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    return row, blocks


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--subjects", nargs="*")
    ap.add_argument("--tasks", nargs="*", default=["video", "rest", "memory"])
    a = ap.parse_args(argv)
    cfg = load_json(PROJECT_ROOT / "config" / "download_config.json")
    subs = a.subjects or cfg["pilot"]["subjects"]
    root = bidsio.raw_root()
    outdir = PROJECT_ROOT / "results" / "beat_correction"
    bdir = PROJECT_ROOT / "derivatives" / "beats"
    outdir.mkdir(parents=True, exist_ok=True)
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    (PROJECT_ROOT / "logs").mkdir(exist_ok=True)
    lines, rows, blocks, failed = [], [], [], 0
    print(f"beat correction {stamp} | python {platform.python_version()} numpy {np.__version__} scipy {scipy.__version__}")
    for sub in subs:
        for task in a.tasks:
            try:
                r, b = process(root, sub, task, outdir, bdir)
                msg = (f"{sub} {task}: detected {r['detected']}, shape rejected {r['template_rejected']}, extra {r['extra_removed']}, "
                       f"missed {r['missed_inserted']}, marked {r['pct_intervals_marked']}%, RMSSD {r['rmssd_ms']} ms, "
                       f"PPG corr {r.get('ppg_corr_raw', 'n/a')} -> {r.get('ppg_corr_corrected', 'n/a')}")
            except Exception as e:
                import traceback
                failed += 1
                r, b = {"sub": sub, "task": task, "error": f"{type(e).__name__}: {e}"}, []
                msg = f"{sub} {task}: ERROR " + "".join(traceback.format_exception(e)).strip()
            print(msg, flush=True)
            lines.append(msg)
            rows.append(r)
            blocks += b
    write_tsv(outdir / "beat_correction_summary.tsv", rows)
    if blocks:
        write_tsv(outdir / "beat_correction_blocks.tsv", blocks)
    (PROJECT_ROOT / "logs" / f"beat_correction_{stamp}.log").write_text("\n".join(lines) + "\n", encoding="utf-8")
    append_execution_log([
        f"## {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} beat correction (step 03b)",
        f"- Command: python -m hb.beat_correction {' '.join(sys.argv[1:])}",
        f"- Python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}, {platform.platform()}",
        f"- Recordings: {len(rows)}; errors: {failed}",
        f"- Outputs: results/beat_correction/*.png, beat_correction_summary.tsv, beat_correction_blocks.tsv; derivatives/beats/*.tsv; log logs/beat_correction_{stamp}.log",
    ])
    print("DONE" if failed == 0 else f"FINISHED WITH {failed} ERROR(S)")
    return 0 if failed == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
