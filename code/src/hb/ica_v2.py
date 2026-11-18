"""Step 16a (docs/ANALYSIS_FREEZE_v2_EEG.md): EEG cleaning by ICA, fit stage only.

Per participant, from the raw memory task recording:
  1. EEG channels (ECG and PPG dropped), median removed per channel, resampled 1000 -> 250 Hz (polyphase)
  2. notch at 50 and 100 Hz; band pass 1 to 100 Hz
  3. bad channels by fixed automatic criteria (flat, deviation, high frequency noise, low correlation),
     spherical spline interpolation, then common average reference
  4. ICA fitted on the task blocks only, 2 s windows above 1 mV peak to peak skipped, every 2nd sample:
       picard  extended infomax (Picard solver, ortho=False, extended=True), n_components = good channels - 1
  5. ICLabel probabilities for every component (7 classes)
Nothing is removed from the data here. Removal (rule fixed in the freeze) and envelopes are step 16b.
Outputs: derivatives/ica_v2/<method>/<sub>-ica.fif, <sub>_iclabel.tsv, <sub>_qc.json; review figures for the
pilot participants in results/ica_v2/review/.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import platform
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import scipy
from scipy.signal import resample_poly

from hb import bidsio
from hb.download_ds008302 import PROJECT_ROOT, append_execution_log
from hb.pilot_coherence import all_blocks

SEED = 20261118
FS_ICA = 250.0
ICLABEL_CLASSES = ["brain", "muscle", "eye", "heart", "line_noise", "channel_noise", "other"]
PILOT = ["sub-124", "sub-125", "sub-128"]
CRIT = {"flat_uV": 0.5, "deviation_z": 5.0, "hf_noise_z": 5.0, "corr_min": 0.4, "corr_bad_window_share": 0.01,
        "max_bad_channels": 12, "reject_p2p_uV": 1000.0, "decim": 2}
OUT = lambda *p: PROJECT_ROOT.joinpath(*p)


def robust_z(x):
    med = np.median(x)
    mad = 1.4826 * np.median(np.abs(x - med))
    return (x - med) / (mad if mad > 0 else 1.0)


def find_bad_channels(data_uV, fs, names, inblock):
    """PREP style criteria (Bigdely-Shamlo et al. 2015 defaults) on 1 to 100 Hz data, task blocks only.
    Returns dict channel -> reason, at most CRIT['max_bad_channels'] channels, ranked flat, deviation,
    high frequency noise, low correlation (the participant is flagged when more were found, see prepare)."""
    from scipy.signal import butter, sosfiltfilt
    X = data_uV[:, inblock]
    bad = {}
    sd = 1.4826 * np.median(np.abs(X - np.median(X, axis=1, keepdims=True)), axis=1)
    for i in np.where(sd < CRIT["flat_uV"])[0]:
        bad[names[i]] = "flat"
    iqr_sd = (np.percentile(X, 75, axis=1) - np.percentile(X, 25, axis=1)) * 0.7413
    zd = robust_z(iqr_sd)
    for i in np.argsort(-zd):
        if zd[i] > CRIT["deviation_z"]:
            bad.setdefault(names[i], "deviation")
    lo = sosfiltfilt(butter(4, 50.0, btype="low", fs=fs, output="sos"), X, axis=1)
    hf = np.median(np.abs(X - lo), axis=1) / np.maximum(np.median(np.abs(lo), axis=1), 1e-12)
    zh = robust_z(hf)
    for i in np.argsort(-zh):
        if zh[i] > CRIT["hf_noise_z"]:
            bad.setdefault(names[i], "high frequency noise")
    w = int(fs)
    starts = np.arange(0, X.shape[1] - w, w)[::10]
    low = np.zeros(X.shape[0])
    for s in starts:
        seg = lo[:, s:s + w]
        seg = seg - seg.mean(axis=1, keepdims=True)
        nrm = np.linalg.norm(seg, axis=1)
        ok = nrm > 0
        C = np.zeros((X.shape[0], X.shape[0]))
        C[np.ix_(ok, ok)] = np.abs((seg[ok] / nrm[ok, None]) @ (seg[ok] / nrm[ok, None]).T)
        np.fill_diagonal(C, 0.0)
        low += C.max(axis=1) < CRIT["corr_min"]
    share = low / max(len(starts), 1)
    for i in np.argsort(-share):
        if share[i] > CRIT["corr_bad_window_share"]:
            bad.setdefault(names[i], f"low correlation in {100 * share[i]:.1f}% of windows")
    find_bad_channels.n_found = len(bad)
    return dict(list(bad.items())[:CRIT["max_bad_channels"]])


def prepare(root, sub):
    """Returns (raw 1 to 100 Hz average referenced, blocks, qc dict)."""
    import mne
    p = bidsio.rec_paths(Path(root), sub, "memory")
    rec = bidsio.load_set(p["set"])
    fs = rec["fs"]
    q = int(round(fs / FS_ICA))
    idx = [i for i, lab in enumerate(rec["all_labels"]) if lab not in ("ECG", "PPG")]
    names = [rec["all_labels"][i] for i in idx]
    n = int(np.ceil(rec["n_samples"] / q))
    D = np.empty((len(idx), n), np.float64)
    for k, i in enumerate(idx):
        x = rec["data"][i].astype(np.float64)
        D[k] = resample_poly(x - np.median(x), 1, q)[:n]
    del rec
    info = mne.create_info(names, FS_ICA, "eeg")
    raw = mne.io.RawArray(D * 1e-6, info, verbose="error")
    del D
    name = "colin27_1005" if "colin27_1005" in mne.channels.get_builtin_montages() else "standard_1005"   # renamed in MNE 1.13 (deprecation warning)
    raw.set_montage(mne.channels.make_standard_montage(name), on_missing="warn", verbose="error")
    raw.notch_filter([50.0, 100.0], verbose="error")
    raw.filter(1.0, 100.0, fir_design="firwin", verbose="error")
    blocks, _ = all_blocks(Path(root), sub)
    t = raw.times
    inblock = np.zeros(t.size, bool)
    for b in blocks:
        inblock |= (t >= b["t0"]) & (t < b["t1"])
    bad = find_bad_channels(raw.get_data() * 1e6, FS_ICA, names, inblock)
    n_found = find_bad_channels.n_found
    qc = {"sub": sub, "n_channels": len(names), "bad_channels": bad, "n_bad": len(bad), "n_bad_found": n_found,
          "too_many_bad": n_found > CRIT["max_bad_channels"], "task_seconds": float(inblock.sum() / FS_ICA)}
    raw.info["bads"] = sorted(bad)
    if bad:
        raw.interpolate_bads(reset_bads=True, verbose="error")
    raw.set_eeg_reference("average", verbose="error")
    onsets, durs = [], []
    edges = np.flatnonzero(np.diff(np.r_[0, (~inblock).astype(int), 0]))
    for a, b in zip(edges[::2], edges[1::2]):
        onsets.append(a / FS_ICA)
        durs.append((b - a) / FS_ICA)
    raw.set_annotations(mne.Annotations(onsets, durs, ["BAD_outside_task"] * len(onsets)))
    return raw, blocks, qc


def fit_subject(args):
    root, sub, method = args
    odir = OUT("derivatives", "ica_v2", method)
    out = odir / f"{sub}-ica.fif"
    if out.exists():
        return sub, "skipped (already done)"
    import mne
    from mne_icalabel.iclabel import iclabel_label_components
    t0 = time.time()
    raw, blocks, qc = prepare(root, sub)
    n_comp = len(raw.ch_names) - len(qc["bad_channels"]) - 1
    import inspect
    seed_kw = {"rng": SEED} if "rng" in inspect.signature(mne.preprocessing.ICA.__init__).parameters else {"random_state": SEED}
    if method == "picard":
        ica = mne.preprocessing.ICA(n_components=n_comp, method="picard", fit_params=dict(ortho=False, extended=True),
                                    max_iter=500, **seed_kw)
    else:
        raise ValueError(f"unknown method {method}")
    ica.fit(raw, reject=dict(eeg=CRIT["reject_p2p_uV"] * 1e-6), reject_by_annotation=True, decim=CRIT["decim"], verbose="error")
    proba = iclabel_label_components(raw, ica, inplace=False)
    odir.mkdir(parents=True, exist_ok=True)
    rows = ["component\t" + "\t".join(ICLABEL_CLASSES) + "\tlabel\tprobability"]
    for k, pr in enumerate(proba):
        j = int(np.argmax(pr))
        rows.append(f"{k}\t" + "\t".join(f"{v:.4f}" for v in pr) + f"\t{ICLABEL_CLASSES[j]}\t{pr[j]:.4f}")
    (odir / f"{sub}_iclabel.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")
    qc.update(method=method, n_components=int(ica.n_components_), n_iter=int(getattr(ica, "n_iter_", -1)),
              fit_samples=int(ica.n_samples_), fit_seconds=round(time.time() - t0, 1), seed=SEED, criteria=CRIT,
              mne=mne.__version__, labels_ge_0p8={c: int(((proba.argmax(1) == i) & (proba.max(1) >= 0.8)).sum()) for i, c in enumerate(ICLABEL_CLASSES)})
    (odir / f"{sub}_qc.json").write_text(json.dumps(qc, indent=1), encoding="utf-8")
    tmp = odir / f"{sub}-part-ica.fif"
    ica.save(tmp, overwrite=True, verbose="error")
    tmp.replace(out)
    if sub in PILOT:
        try:
            review_figures(raw, ica, proba, sub, method)
        except Exception as e:  # figures are a convenience; never fail the fit for them
            qc["review_error"] = str(e)
            (odir / f"{sub}_qc.json").write_text(json.dumps(qc, indent=1), encoding="utf-8")
    lab = ", ".join(f"{c} {v}" for c, v in qc["labels_ge_0p8"].items() if v)
    return sub, f"ok ({qc['n_bad']} bad channels, {qc['n_components']} components; >=0.8: {lab}; {qc['fit_seconds']:.0f} s)"


def review_figures(raw, ica, proba, sub, method):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    rdir = OUT("results", "ica_v2", "review")
    rdir.mkdir(parents=True, exist_ok=True)
    figs = ica.plot_components(show=False, verbose="error")
    figs = figs if isinstance(figs, list) else [figs]
    for fi, fig in enumerate(figs):
        for ax in fig.axes:
            title = ax.get_title()
            if title.startswith("ICA"):
                k = int(title[3:])
                j = int(np.argmax(proba[k]))
                ax.set_title(f"{k}: {ICLABEL_CLASSES[j]} {proba[k, j]:.2f}", fontsize=8)
        fig.savefig(rdir / f"{sub}_{method}_components_{fi + 1}.png", dpi=110)
        plt.close(fig)


def run(args_list, workers, log, t0):
    failed = 0
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(fit_subject, a): a[1] for a in args_list}
        for i, f in enumerate(as_completed(futs), 1):
            try:
                sub, status = f.result()
            except Exception as e:
                failed += 1
                sub, status = futs[f], "FAILED: " + "".join(traceback.format_exception(e)).strip().replace("\n", " | ")
            log(f"[{i}/{len(args_list)}] {sub}: {status[:400]} ({(time.time() - t0) / 60:.1f} min elapsed)")
    return failed


def summarize(method, log):
    rows = []
    for f in sorted(OUT("derivatives", "ica_v2", method).glob("sub-*_qc.json")):
        rows.append(json.loads(f.read_text(encoding="utf-8")))
    if not rows:
        return
    nb = np.array([r["n_bad"] for r in rows])
    log(f"summary {method}: {len(rows)} participants; bad channels median {np.median(nb):.0f} (max {nb.max()}); "
        f"more than {CRIT['max_bad_channels']}: {sum(r['too_many_bad'] for r in rows)}")
    for c in ICLABEL_CLASSES:
        v = np.array([r["labels_ge_0p8"][c] for r in rows])
        log(f"  components labelled {c} with probability >= 0.8: median {np.median(v):.0f} per participant (range {v.min()} to {v.max()})")
    OUT("results", "ica_v2").mkdir(parents=True, exist_ok=True)
    with open(OUT("results", "ica_v2", f"fit_summary_{method}.tsv"), "w", encoding="utf-8") as fh:
        fh.write("sub\tn_bad\tbad_channels\tn_components\tfit_seconds\t" + "\t".join(f"n_{c}_ge_0p8" for c in ICLABEL_CLASSES) + "\n")
        for r in rows:
            fh.write(f"{r['sub']}\t{r['n_bad']}\t{';'.join(f'{k}:{v}' for k, v in r['bad_channels'].items())}\t{r['n_components']}\t{r['fit_seconds']}\t"
                     + "\t".join(str(r["labels_ge_0p8"][c]) for c in ICLABEL_CLASSES) + "\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["fit"])
    ap.add_argument("--method", default="picard", choices=["picard"])
    ap.add_argument("--subjects", nargs="*")
    ap.add_argument("--workers", type=int, default=3, help="each worker holds about 5 GB of memory")
    ap.add_argument("--summary-only", action="store_true")
    a = ap.parse_args(argv)
    root = bidsio.raw_root()
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    OUT("logs").mkdir(exist_ok=True)
    logf = open(OUT("logs", f"ica_v2_{a.command}_{a.method}_{stamp}.log"), "w", encoding="utf-8")

    def log(m):
        print(m, flush=True)
        logf.write(m + "\n")
        logf.flush()

    import mne
    import mne_icalabel
    import picard
    log(f"ica v2 {a.command} {a.method} {stamp} | workers {a.workers} seed {SEED} | python {platform.python_version()} numpy {np.__version__} "
        f"scipy {scipy.__version__} mne {mne.__version__} mne-icalabel {mne_icalabel.__version__} picard {picard.__version__}")
    log(f"criteria {json.dumps(CRIT)}")
    t0, failed = time.time(), 0
    subs = a.subjects or sorted(f.name.split("_")[0] for f in OUT("derivatives", "slow").glob("sub-*_task-memory_slow.npz"))
    subs = [s for s in PILOT if s in subs] + [s for s in subs if s not in PILOT]          # pilot first: review figures early
    if not a.summary_only:
        failed = run([(str(root), s, a.method) for s in subs], a.workers, log, t0)
    try:
        summarize(a.method, log)
    except Exception as e:
        failed += 1
        log("SUMMARY FAILED: " + "".join(traceback.format_exception(e)).strip())
    append_execution_log([
        f"## {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} ICA v2 {a.command} {a.method} (step 16a, docs/ANALYSIS_FREEZE_v2_EEG.md)",
        f"- Command: python -m hb.ica_v2 {' '.join(sys.argv[1:])}; seed {SEED}; Python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}, mne {mne.__version__}, mne-icalabel {mne_icalabel.__version__}, picard {picard.__version__}",
        f"- Failures: {failed}; {(time.time() - t0) / 60:.1f} min. Outputs: derivatives/ica_v2/{a.method}, results/ica_v2; log logs/ica_v2_{a.command}_{a.method}_{stamp}.log",
    ])
    return 0 if failed == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
