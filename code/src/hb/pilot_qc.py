"""Step 03: inventory and ECG quality check of the pilot recordings.

Reads each recording once, checks file content against the BIDS text files, detects R peaks,
flags suspicious beats, compares with PPG, measures heartbeat artifact in the EEG, and for the
memory task tabulates blocks. Writes figures and tables to results/pilot_qc and beat times to
derivatives/rpeaks. Feasibility information only: no coherence and no condition contrast here.
"""
from __future__ import annotations

import argparse
import datetime as dt
import gc
import platform
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy

from hb import bidsio, cardiac
from hb.download_ds008302 import PROJECT_ROOT, append_execution_log, load_json

EEG_PROBE = ["Fz", "Cz", "Pz", "Oz", "T7", "T8"]


def interp_series(t_beats, values, t):
    return np.interp(t, t_beats, values, left=np.nan, right=np.nan)


def qc_recording(root: Path, sub: str, task: str, outdir: Path, rdir: Path, log) -> tuple[dict, list[dict]]:
    p = bidsio.rec_paths(root, sub, task)
    row = {"sub": sub, "task": task}
    if not p["set"].exists():
        row["status"] = "MISSING FILE"
        return row, []
    t0 = time.time()
    side = bidsio.sidecar(p["json"])
    chan_rows = bidsio.read_tsv(p["channels"])
    ev_rows = bidsio.read_tsv(p["events"])
    rec = bidsio.load_set(p["set"])
    fs, n = rec["fs"], rec["n_samples"]
    dur = n / fs
    labels = rec["all_labels"]
    log(f"  loaded {p['set'].name}: {len(labels)} channels, {n} samples, {fs:g} Hz, {dur:.1f} s ({time.time() - t0:.0f} s to read)")

    # ---- file content against BIDS text files
    row.update(fs=fs, duration_s=round(dur, 3), n_channels=len(labels), ref=rec["ref"])
    row["labels_match_channels_tsv"] = labels == [r["name"] for r in chan_rows]
    row["fs_matches_sidecar"] = abs(fs - float(side["SamplingFrequency"])) < 1e-9
    row["duration_minus_sidecar_s"] = round(dur - float(side["RecordingDuration"]), 3)
    set_ev = [(t, lat) for t, lat, _ in rec["events"] if t != "boundary"]
    tsv_ev = [(r["value"], float(r["onset"])) for r in ev_rows if r["value"] != "boundary"]
    row["n_events_set"], row["n_events_tsv"] = len(set_ev), len(tsv_ev)
    if len(set_ev) == len(tsv_ev) and set_ev:
        row["event_types_match"] = all(a[0].strip() == b[0].strip() for a, b in zip(set_ev, tsv_ev))
        row["max_event_time_diff_ms"] = round(1000 * max(abs(a[1] - b[1]) for a, b in zip(set_ev, tsv_ev)), 2)
    bnd = [(lat, d) for t, lat, d in rec["events"] if t == "boundary" and 1.0 < lat < dur - 1.0]
    row["interior_boundaries"] = len({round(b[0], 2) for b in bnd})
    row["boundary_gap_s"] = ";".join(f"{lat:.1f}:{(d / fs if np.isfinite(d) else float('nan')):.1f}" for lat, d in bnd)

    # ---- signals
    ecg = rec["data"][labels.index("ECG")].astype(float)
    ppg = rec["data"][labels.index("PPG")].astype(float) if "PPG" in labels else None
    probe = [c for c in EEG_PROBE if c in labels]
    eeg = {c: rec["data"][labels.index(c)].astype(float) for c in probe}
    eeg_sd = np.array([rec["data"][i, ::10].std() for i, lab in enumerate(labels) if lab not in ("ECG", "PPG")])
    row["eeg_flat_channels"] = int((eeg_sd < 0.5).sum())
    row["eeg_median_sd_uV"] = round(float(np.median(eeg_sd)), 1)
    del rec
    gc.collect()

    row["ecg_dc_offset"] = round(float(np.median(ecg)), 0)
    det = cardiac.detect_rpeaks(ecg, fs)
    pk, x = det["peaks"], det["filtered"]
    rr, flag = cardiac.rr_flags(pk, fs)
    tb = pk[1:] / fs
    row["ecg_polarity"] = det["polarity"]
    row["n_beats"] = int(pk.size)
    row["mean_hr_bpm"] = round(60.0 / np.median(rr), 1)
    row["r_amp_median"] = round(float(np.median(x[pk])), 0)
    noise = np.median(np.abs(x - np.median(x))) * 1.4826
    row["r_amp_over_noise"] = round(float(np.median(x[pk]) / noise), 1)
    row["pct_rr_out_of_range"] = round(100 * float((flag == 1).mean()), 3)
    row["pct_rr_deviant"] = round(100 * float((flag == 2).mean()), 3)
    row["pct_rr_ok"] = round(100 * float((flag == 0).mean()), 3)
    bad_t = tb[flag > 0]
    edges = np.r_[0.0, bad_t, dur]
    row["longest_clean_run_s"] = round(float(np.diff(edges).max()), 1)
    row["rmssd_ms_ok_beats"] = round(float(1000 * np.sqrt(np.mean(np.diff(rr[flag == 0]) ** 2))), 1)

    # ---- PPG cross check
    tp = ibi = None
    if ppg is not None:
        pp = cardiac.ppg_peaks(ppg, fs)
        if pp.size > 10:
            ibi, tp = np.diff(pp) / fs, pp[1:] / fs
            row["ppg_beats_over_ecg_beats"] = round(pp.size / pk.size, 4)
            row["ppg_minus_ecg_median_ibi_ms"] = round(1000 * float(np.median(ibi) - np.median(rr)), 1)
            grid = np.arange(5, dur - 5, 0.25)
            a, b = interp_series(tb, rr, grid), interp_series(tp, ibi, grid)
            ok = np.isfinite(a) & np.isfinite(b)
            row["rr_ppg_interval_correlation"] = round(float(np.corrcoef(a[ok], b[ok])[0, 1]), 3)

    # ---- heartbeat artifact in EEG (R locked average, -200 to +400 ms)
    pre, post = int(0.2 * fs), int(0.4 * fs)
    good = pk[(pk > pre) & (pk < n - post)]
    erp = {}
    for c in probe:
        seg = np.stack([eeg[c][g - pre:g + post] for g in good[:: max(1, good.size // 2000)]])
        m = seg.mean(axis=0)
        erp[c] = m - m[:int(0.1 * fs)].mean()
    row["eeg_rlocked_p2p_uV"] = ";".join(f"{c}:{np.ptp(erp[c]):.1f}" for c in probe)
    ecg_avg = np.stack([x[g - pre:g + post] for g in good[:: max(1, good.size // 2000)]]).mean(axis=0)

    # ---- memory task blocks
    block_rows = []
    parsed = None
    if task == "memory":
        parsed = bidsio.parse_memory_events(ev_rows, dur)
        tr = parsed["trials"]
        row["n_trials"], row["n_blocks"] = len(tr), len(parsed["blocks"])
        for bi, b in enumerate(parsed["blocks"]):
            if not b["trials"]:
                continue
            starts = np.array([t["t_start"] for t in b["trials"]])
            per = np.diff(starts)
            inb = (tb >= b["t0"]) & (tb < b["t1"])
            known = [t["correct"] for t in b["trials"] if t["correct"] is not None]
            block_rows.append({
                "sub": sub, "block": bi + 1, "length": b["length"], "t0_s": round(b["t0"], 2), "t1_s": round(b["t1"], 2),
                "duration_s": round(b["t1"] - b["t0"], 1), "n_trials": len(b["trials"]),
                "n_correct": sum(known) if known else "", "trial_period_mean_s": round(float(per.mean()), 2) if per.size else "",
                "trial_period_sd_s": round(float(per.std()), 2) if per.size else "",
                "trial_freq_hz": round(float(1 / per.mean()), 4) if per.size else "",
                "boundary_inside": int(any(b["t0"] < q < b["t1"] for q in parsed["boundaries"])),
                "n_beats": int(inb.sum()), "pct_beats_flagged": round(100 * float((flag[inb] > 0).mean()), 2) if inb.any() else "",
                "mean_hr_bpm": round(60.0 / float(np.median(rr[inb])), 1) if inb.any() else "",
            })

    # ---- beat table
    rdir.mkdir(parents=True, exist_ok=True)
    with open(rdir / f"{sub}_task-{task}_rpeaks.tsv", "w", encoding="utf-8", newline="\n") as fh:
        fh.write("r_time_s\trr_s\tflag\n")
        fh.write(f"{pk[0] / fs:.3f}\tn/a\tn/a\n")
        for t, r, f in zip(tb, rr, flag):
            fh.write(f"{t:.3f}\t{r:.4f}\t{f}\n")

    # ---- figure
    fig = plt.figure(figsize=(15, 11))
    gs = fig.add_gridspec(4, 3, hspace=0.45, wspace=0.25)
    for j, (name, c) in enumerate((("start", 30.0), ("middle", dur / 2), ("end", dur - 30.0))):
        ax = fig.add_subplot(gs[0, j])
        a, b = int(max(c - 4, 0) * fs), int(min(c + 4, dur) * fs)
        tt = np.arange(a, b) / fs
        ax.plot(tt, x[a:b], lw=0.7, color="0.2")
        sel = pk[(pk >= a) & (pk < b)]
        ax.plot(sel / fs, x[sel], "v", ms=5, color="C3")
        ax.set_title(f"ECG 5 to 25 Hz, polarity corrected ({name})", fontsize=9)
        ax.set_xlabel("time (s)")
    ax = fig.add_subplot(gs[1, :])
    ax.plot(tb, rr, lw=0.5, color="C0", label="RR from ECG")
    if tp is not None:
        ax.plot(tp, ibi, lw=0.4, color="C2", alpha=0.6, label="pulse interval from PPG")
    ax.plot(tb[flag > 0], rr[flag > 0], "x", color="C3", ms=5, label=f"flagged beats ({(flag > 0).sum()})")
    if parsed:
        for b in parsed["blocks"]:
            if b["trials"]:
                ax.axvspan(b["t0"], b["t1"], color="0.85", lw=0)
                ax.text((b["t0"] + b["t1"]) / 2, 0.02, str(b["length"]), transform=ax.get_xaxis_transform(), ha="center", fontsize=7)
        for q in parsed["boundaries"]:
            ax.axvline(q, color="k", ls=":", lw=1)
    lo_, hi_ = np.percentile(rr, [0.5, 99.5])
    ax.set_ylim(lo_ - 0.15, hi_ + 0.15)
    ax.set_xlim(0, dur)
    ax.set_ylabel("interval (s)")
    ax.set_xlabel("time (s); grey = task blocks with sequence length, dotted = recording discontinuity")
    ax.legend(fontsize=8, loc="upper right", ncol=3)
    ax = fig.add_subplot(gs[2, 0])
    tl = (np.arange(-pre, post)) / fs * 1000
    ax.plot(tl, ecg_avg, color="0.2")
    ax.set_title("R locked average ECG (filtered)", fontsize=9)
    ax.set_xlabel("ms from R")
    ax = fig.add_subplot(gs[2, 1])
    for c in probe:
        ax.plot(tl, erp[c], lw=0.9, label=c)
    ax.set_title("R locked average EEG, online reference (uV)", fontsize=9)
    ax.set_xlabel("ms from R")
    ax.legend(fontsize=7, ncol=3)
    ax = fig.add_subplot(gs[2, 2])
    ax.hist(rr[flag == 0], bins=60, color="C0")
    ax.set_title("RR distribution, unflagged beats", fontsize=9)
    ax.set_xlabel("s")
    ax = fig.add_subplot(gs[3, :])
    if parsed and block_rows:
        per_all = [(t["t_start"], t["length"]) for t in parsed["trials"]]
        ts = np.array([a for a, _ in per_all])
        ax.plot(ts[1:], np.diff(ts), ".", ms=3, color="C1")
        ax.set_ylim(0, 80)
        ax.axhspan(1 / 0.15, 1 / 0.04, color="C0", alpha=0.12)
        ax.set_ylabel("trial period (s)")
        ax.set_xlabel("time (s); blue band = periods inside the LF range (6.7 to 25 s)")
        ax.set_xlim(0, dur)
    else:
        grid = np.arange(0, dur, 0.25)
        ax.plot(grid, 60.0 / interp_series(tb, rr, grid), lw=0.7)
        ax.set_ylabel("heart rate (bpm)")
        ax.set_xlabel("time (s)")
        ax.set_xlim(0, dur)
    fig.suptitle(f"{sub} task-{task} | {row['n_beats']} beats, {row['mean_hr_bpm']} bpm, {row['pct_rr_ok']}% intervals unflagged, "
                 f"R/noise {row['r_amp_over_noise']}, polarity {det['polarity']:+d}", fontsize=11)
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / f"{sub}_task-{task}_qc.png", dpi=110, bbox_inches="tight")
    plt.close(fig)
    row["status"] = "ok"
    row["seconds"] = round(time.time() - t0)
    return row, block_rows


def write_tsv(path: Path, rows: list[dict]) -> None:
    cols = []
    for r in rows:
        cols += [k for k in r if k not in cols]
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in cols) + "\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--subjects", nargs="*")
    ap.add_argument("--tasks", nargs="*", default=["video", "rest", "memory"])
    a = ap.parse_args(argv)
    cfg = load_json(PROJECT_ROOT / "config" / "download_config.json")
    subs = a.subjects or cfg["pilot"]["subjects"]
    root = bidsio.raw_root()
    outdir = PROJECT_ROOT / "results" / "pilot_qc"
    rdir = PROJECT_ROOT / "derivatives" / "rpeaks"
    outdir.mkdir(parents=True, exist_ok=True)
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    (PROJECT_ROOT / "logs").mkdir(exist_ok=True)
    logf = open(PROJECT_ROOT / "logs" / f"pilot_qc_{stamp}.log", "w", encoding="utf-8")

    def log(msg):
        print(msg, flush=True)
        logf.write(msg + "\n")
        logf.flush()

    log(f"pilot QC {stamp} | python {platform.python_version()} numpy {np.__version__} scipy {scipy.__version__} matplotlib {matplotlib.__version__}")
    log(f"raw data root: {root}")
    rows, blocks, failed = [], [], 0
    for sub in subs:
        for task in a.tasks:
            log(f"{sub} task-{task}")
            try:
                r, b = qc_recording(root, sub, task, outdir, rdir, log)
            except Exception as e:  # keep going, report at the end
                import traceback
                log("  ERROR: " + "".join(traceback.format_exception(e)).strip())
                r, b = {"sub": sub, "task": task, "status": f"ERROR {type(e).__name__}: {e}"}, []
            if r.get("status") != "ok":
                failed += 1
            else:
                log(f"  beats {r['n_beats']}, HR {r['mean_hr_bpm']} bpm, unflagged {r['pct_rr_ok']}%, polarity {r['ecg_polarity']:+d}, "
                    f"R/noise {r['r_amp_over_noise']}, PPG r {r.get('rr_ppg_interval_correlation', 'n/a')}, "
                    f"boundaries {r['interior_boundaries']}, events match {r.get('event_types_match', 'n/a')}")
            rows.append(r)
            blocks += b
    write_tsv(outdir / "pilot_qc_summary.tsv", rows)
    if blocks:
        write_tsv(outdir / "pilot_blocks.tsv", blocks)
    log(f"wrote {outdir / 'pilot_qc_summary.tsv'}" + (f" and pilot_blocks.tsv ({len(blocks)} blocks)" if blocks else ""))
    append_execution_log([
        f"## {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} pilot QC (step 03)",
        f"- Command: python -m hb.pilot_qc {' '.join(sys.argv[1:])}",
        f"- Python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}, matplotlib {matplotlib.__version__}, {platform.platform()}",
        f"- Recordings: {len(rows)}; not ok: {failed}",
        f"- Outputs: results/pilot_qc/pilot_qc_summary.tsv, pilot_blocks.tsv, one PNG per recording; derivatives/rpeaks/*.tsv; log logs/pilot_qc_{stamp}.log",
    ])
    log("DONE" if failed == 0 else f"FINISHED WITH {failed} RECORDING(S) NOT OK")
    return 0 if failed == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
