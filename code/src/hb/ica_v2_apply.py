"""Steps 16b and 16c (docs/ANALYSIS_FREEZE_v2_EEG.md).

16b  apply     per participant: the 16a preprocessing again (deterministic), ICA components removed by the fixed
               ICLabel rule, region envelopes as v1 (from 200 Hz data), the same envelopes without removal
               ("noica" control), and a 55 to 95 Hz muscle proxy per region (from the 250 Hz cleaned data).
16c  analyse   frozen estimand on three envelope variants, trial shuffle null of step 15 (K = 16):
                 v2          ICA cleaned (primary)
                 v2_noica    same pipeline, nothing removed
                 v2_emgadj   v2 with the muscle proxy regressed out of the envelope within each block
               claim rule = rule 2 of docs/STEP15_PLAN_v1.md with the step 15b synthetic bias edges.
"""
from __future__ import annotations

import argparse
import csv
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
from scipy import stats
from scipy.signal import butter, hilbert, resample_poly, sosfiltfilt

from hb import bidsio
from hb import hbcoh as H
from hb import ica_v2 as IV
from hb import trial_shuffle as TS
from hb.cohort import COLS, coherence_tables, regress_lagged
from hb.cohort_stats import fmt, holm, participant_table
from hb.download_ds008302 import PROJECT_ROOT, append_execution_log
from hb.pilot_coherence import BANDS, EEG, all_blocks, marked_share
from hb.pilot_qc import write_tsv
from hb.slow_signals import ROIS

SEED = 20261119
REMOVE = ("muscle", "eye", "heart", "line_noise", "channel_noise")
P_MIN = 0.8
EMG_BAND = (55.0, 95.0)
VARIANTS_V2 = ["v2", "v2_noica", "v2_emgadj"]
OUT = lambda *p: PROJECT_ROOT.joinpath(*p)


# ---------------------------------------------------------------- 16b
def components_to_remove(tsv_path):
    with open(tsv_path, encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    out = []
    for r in rows:
        p = np.array([float(r[c]) for c in IV.ICLABEL_CLASSES])
        j = int(np.argmax(p))
        if IV.ICLABEL_CLASSES[j] in REMOVE and p[j] >= P_MIN:
            out.append((int(r["component"]), IV.ICLABEL_CLASSES[j]))
    return out


def clean_envelope(env):
    """Outlier handling identical to v1 (slow_signals.roi_envelopes): above median + 8 robust SD, 1 s margins."""
    med = np.median(env)
    rsd = 1.4826 * np.median(np.abs(env - med))
    bad = env > med + 8 * rsd
    bad = np.convolve(bad.astype(float), np.ones(int(2 * H.FS_SLOW) + 1), mode="same") > 0
    y = env.copy()
    if bad.any() and (~bad).any():
        idx = np.arange(env.size)
        y[bad] = np.interp(idx[bad], idx[~bad], env[~bad])
    return y, bad


def region_envelopes(data_uV, names):
    """data_uV at 250 Hz -> 200 Hz -> band envelopes at 4 Hz exactly as in v1."""
    out = {}
    D = resample_poly(data_uV, 4, 5, axis=1)
    for e, (rng, chans) in ROIS.items():
        use = [c for c in chans if c in names]
        if not use:
            raise ValueError(f"no channel of region {e} in the recording")
        env = None
        for c in use:
            _, x = H.band_envelope(D[names.index(c)], 200.0, rng)
            env = x if env is None else env + x
        out[e] = clean_envelope(env / len(use))
    return out


def muscle_proxy(data_uV, names):
    sos = butter(4, EMG_BAND, btype="band", fs=IV.FS_ICA, output="sos")
    out = {}
    for e, (_, chans) in ROIS.items():
        use = [names.index(c) for c in chans if c in names]
        env = np.abs(hilbert(sosfiltfilt(sos, data_uV[use], axis=1), axis=1)).mean(axis=0)
        out[e] = resample_poly(env, 2, 125)                           # 250 Hz -> 4 Hz
    return out


def apply_subject(args):
    root, sub, method = args
    odir = OUT("derivatives", "slow_v2", method)
    out = odir / f"{sub}_task-memory_slow.npz"
    if out.exists():
        return sub, "skipped (already done)"
    import mne
    t0 = time.time()
    idir = OUT("derivatives", "ica_v2", method)
    qc_fit = json.loads((idir / f"{sub}_qc.json").read_text(encoding="utf-8"))
    raw, blocks, qc = IV.prepare(root, sub)
    if qc["bad_channels"] != qc_fit["bad_channels"]:
        raise RuntimeError(f"bad channels differ from the 16a fit: {qc['bad_channels']} vs {qc_fit['bad_channels']}")
    ica = mne.preprocessing.read_ica(idir / f"{sub}-ica.fif", verbose="error")
    remove = components_to_remove(idir / f"{sub}_iclabel.tsv")
    names = raw.ch_names
    noica = raw.get_data() * 1e6
    clean = ica.apply(raw.copy(), exclude=[k for k, _ in remove], verbose="error").get_data() * 1e6
    del raw
    env_ica, env_no, emg = region_envelopes(clean, names), region_envelopes(noica, names), muscle_proxy(clean, names)
    frozen = np.load(OUT("derivatives", "slow", f"{sub}_task-memory_slow.npz"))
    n = frozen["nn"].size
    arrays = {"fs": 4.0, "nn": frozen["nn"], "nn_bridged": frozen["nn_bridged"]}
    for e in EEG:
        for tag, src in (("ica", env_ica), ("noica", env_no)):
            y, bad = src[e]
            if y.size < n:
                raise ValueError(f"{e} envelope shorter than the frozen series ({y.size} < {n})")
            arrays[f"{e}_{tag}"], arrays[f"{e}_{tag}_interpolated"] = y[:n], bad[:n]
        m = emg[e]
        arrays[f"emg_{e}"] = np.r_[m, np.full(max(0, n - m.size), m[-1])][:n]
        chans = ROIS[e][1]
        arrays[f"region_excluded_{e}"] = np.array(sum(c in qc["bad_channels"] for c in chans) > len(chans) / 2)
    odir.mkdir(parents=True, exist_ok=True)
    tmp = odir / f"{sub}_task-memory_slow.part.npz"
    np.savez_compressed(tmp, **arrays)
    tmp.replace(out)
    counts = {c: sum(1 for _, k in remove if k == c) for c in REMOVE}
    info = {"sub": sub, "method": method, "n_components": int(ica.n_components_), "removed": [k for k, _ in remove], "removed_by_class": counts,
            "bad_channels": qc["bad_channels"], "region_excluded": {e: bool(arrays[f"region_excluded_{e}"]) for e in EEG},
            "pct_interpolated": {f"{e}_{t}": round(100 * float(arrays[f"{e}_{t}_interpolated"].mean()), 2) for e in EEG for t in ("ica", "noica")},
            "seconds": round(time.time() - t0, 1)}
    (odir / f"{sub}_apply.json").write_text(json.dumps(info, indent=1), encoding="utf-8")
    return sub, f"ok (removed {len(remove)}: " + ", ".join(f"{c} {v}" for c, v in counts.items() if v) + f"; {info['seconds']:.0f} s)"


# ---------------------------------------------------------------- 16c
def emg_adjust(env, emg, blocks):
    y = env.astype(float).copy()
    for b in blocks:
        i0, i1 = int(round(b["t0"] * H.FS_SLOW)), min(int(round(b["t1"] * H.FS_SLOW)), y.size, emg.size)
        if i1 - i0 < 20:
            continue
        c = emg[i0:i1] - emg[i0:i1].mean()
        c = c / (c.std() + 1e-12)
        seg = y[i0:i1]
        y[i0:i1] = regress_lagged(seg, c) + seg.mean()
    return y


def variant_signal(z, e, variant, blocks):
    if variant == "v2":
        env, interp = z[f"{e}_ica"], z[f"{e}_ica_interpolated"]
    elif variant == "v2_noica":
        env, interp = z[f"{e}_noica"], z[f"{e}_noica_interpolated"]
    else:
        env, interp = emg_adjust(z[f"{e}_ica"], z[f"emg_{e}"], blocks), z[f"{e}_ica_interpolated"]
    interp = np.asarray(interp, bool)
    if bool(z[f"region_excluded_{e}"]):
        interp = np.ones_like(interp)
    return np.asarray(env, float), interp


def analyse_subject(args):
    root, sub, k, method = args
    out = OUT("derivatives", "trial_shuffle_v2", method, f"{sub}.tsv")
    if out.exists():
        return sub, "skipped (already done)"
    t0 = time.time()
    blocks, b10, _ = TS.timing(root, sub)
    z = np.load(OUT("derivatives", "slow_v2", method, f"{sub}_task-memory_slow.npz"))
    marked = marked_share(OUT("derivatives", "beats", f"{sub}_task-memory_beats.tsv"), blocks)
    sid = int(sub.split("-")[1])
    nn = z["nn"]
    rows, realism = [], []
    for vi, variant in enumerate(VARIANTS_V2):
        for ei, e in enumerate(EEG):
            env, interp = variant_signal(z, e, variant, blocks)
            r2, vmi = TS.realism_metrics(env, b10)
            realism.append({"sub": sub, "variant": variant, "signal": e, "R2": r2, "VMI": vmi})

            def tab(b):
                n = min(nn.size, b.size)
                r = coherence_tables(sub, nn[:n], {e: b[:n]}, {e: interp[:n]}, blocks, marked)
                return [x for x in participant_table({sub: [{c: q.get(c, "") for c in COLS} for q in r]}) if x["eeg"] == e]

            for r in tab(env):
                r.update(cohort=variant, kind="observed", null="", k=0)
                rows.append(r)
            for kk in range(k):
                b = TS.trial_shuffle(env, b10, np.random.default_rng([SEED, sid, vi, ei, kk]))
                for r in tab(b):
                    r.update(cohort=variant, kind="surrogate", null="shuffle", k=kk)
                    rows.append(r)
    TS.save(out, rows)
    write_tsv(out.with_name(f"{sub}_realism.tsv"), fmt(realism))
    return sub, f"ok ({len(rows)} rows, {time.time() - t0:.0f} s)"


def summarize_analysis(method, log):
    folder = OUT("derivatives", "trial_shuffle_v2", method)
    rows, realism = [], []
    for f in sorted(folder.glob("sub-*.tsv")):
        target = realism if f.name.endswith("_realism.tsv") else rows
        with open(f, encoding="utf-8", newline="") as fh:
            target += list(csv.DictReader(fh, delimiter="\t"))
    edges = json.loads(OUT("results", "trial_shuffle", "gate.json").read_text(encoding="utf-8"))["bias_upper_edges"]
    res = []
    for variant in VARIANTS_V2:
        for v in ("full", "encret"):
            per_k = {}
            for k_max in (8, None):
                _, per = TS.calibrated_matrix(rows, variant, "shuffle", v, k_max)
                fam = []
                for e, b in TS.CELLS:
                    d = per[(e, b)]
                    if d.size < 3:
                        fam.append({"variant_eeg": variant, "variant": v, "k_used": k_max or 16, "eeg": e, "band": b, "n": d.size,
                                    "calibrated_excess": np.nan, "ci95_low": np.nan, "ci95_high": np.nan, "cohen_d": np.nan, "p": np.nan})
                        continue
                    t = stats.ttest_1samp(d, 0.0)
                    h = stats.t.ppf(0.975, d.size - 1) * d.std(ddof=1) / np.sqrt(d.size)
                    fam.append({"variant_eeg": variant, "variant": v, "k_used": k_max or 16, "eeg": e, "band": b, "n": d.size,
                                "calibrated_excess": d.mean(), "ci95_low": d.mean() - h, "ci95_high": d.mean() + h,
                                "cohen_d": d.mean() / d.std(ddof=1), "p": float(t.pvalue)})
                ps = np.array([r["p"] for r in fam])
                adj = np.full(ps.size, np.nan)
                ok = np.isfinite(ps)
                adj[ok] = holm(ps[ok])
                for r, a in zip(fam, adj):
                    r["p_holm"] = float(a)
                per_k[k_max or 16] = fam
            for r16, r8 in zip(per_k[16], per_k[8]):
                edge = edges.get(f"{r16['eeg']} {r16['band']}", np.nan)
                claim = (v == "full" and r16["p_holm"] < 0.05 and r8["p_holm"] < 0.05
                         and r16["ci95_low"] > edge and r8["ci95_low"] > edge)
                for r in (r16, r8):
                    r["synthetic_bias_upper_edge"] = edge
                    r["meets_claim_rule"] = int(bool(claim))
                log(f"{r16['variant_eeg']:9s} {v:6s} {r16['eeg']:5s} {r16['band']}: n {r16['n']} K16 {r16['calibrated_excess']:+.4f} "
                    f"[{r16['ci95_low']:+.4f}, {r16['ci95_high']:+.4f}] d {r16['cohen_d']:.2f} Holm p {r16['p_holm']:.3g} | "
                    f"K8 Holm p {r8['p_holm']:.3g} low {r8['ci95_low']:+.4f} | edge {edge:+.4f}" + (" | CLAIM" if claim else ""))
                res += [r16, r8]
    outdir = OUT("results", "ica_v2")
    outdir.mkdir(parents=True, exist_ok=True)
    write_tsv(outdir / f"shuffle_calibrated_H1_{method}.tsv", fmt(res))
    summ = []
    for variant in VARIANTS_V2:
        for e in EEG:
            r = [x for x in realism if x["variant"] == variant and x["signal"] == e]
            if r:
                r2 = np.median([float(x["R2"]) for x in r])
                vmi = np.median([float(x["VMI"]) for x in r])
                summ.append({"variant": variant, "signal": e, "R2_median": r2, "VMI_median": vmi})
                log(f"realism {variant:9s} {e:5s}: R2 median {r2:.3f}, VMI median {vmi:.3f}")
    write_tsv(outdir / f"realism_{method}.tsv", fmt(summ))


# ---------------------------------------------------------------- driver
def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["apply", "analyse"])
    ap.add_argument("--method", default="picard", choices=["picard"])
    ap.add_argument("--subjects", nargs="*")
    ap.add_argument("--k", type=int, default=16)
    ap.add_argument("--workers", type=int)
    ap.add_argument("--summary-only", action="store_true")
    a = ap.parse_args(argv)
    workers = a.workers or (3 if a.command == "apply" else max(1, min(6, (os.cpu_count() or 2) - 1)))
    root = bidsio.raw_root()
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    OUT("logs").mkdir(exist_ok=True)
    logf = open(OUT("logs", f"ica_v2_{a.command}_{a.method}_{stamp}.log"), "w", encoding="utf-8")

    def log(m):
        print(m, flush=True)
        logf.write(m + "\n")
        logf.flush()

    import mne
    log(f"ica v2 {a.command} {a.method} {stamp} | workers {workers} k {a.k} seed {SEED} | remove {REMOVE} at p >= {P_MIN} | "
        f"python {platform.python_version()} numpy {np.__version__} scipy {scipy.__version__} mne {mne.__version__}")
    t0, failed = time.time(), 0
    subs = a.subjects or sorted(f.name.split("-ica")[0] for f in OUT("derivatives", "ica_v2", a.method).glob("sub-*-ica.fif"))
    fn, arg = ((apply_subject, lambda s: (str(root), s, a.method)) if a.command == "apply"
               else (analyse_subject, lambda s: (str(root), s, a.k, a.method)))
    if not a.summary_only:
        with ProcessPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(fn, arg(s)): s for s in subs}
            for i, f in enumerate(as_completed(futs), 1):
                try:
                    sub, status = f.result()
                except Exception as e:
                    failed += 1
                    sub, status = futs[f], "FAILED: " + "".join(traceback.format_exception(e)).strip().replace("\n", " | ")
                log(f"[{i}/{len(subs)}] {sub}: {status[:400]} ({(time.time() - t0) / 60:.1f} min elapsed)")
    try:
        if a.command == "apply":
            infos = [json.loads(f.read_text(encoding="utf-8")) for f in sorted(OUT("derivatives", "slow_v2", a.method).glob("sub-*_apply.json"))]
            for c in REMOVE:
                v = np.array([i["removed_by_class"][c] for i in infos])
                log(f"removed {c}: median {np.median(v):.0f} per participant (range {v.min()} to {v.max()})")
            for e in EEG:
                log(f"region {e} excluded (more than half of its channels interpolated): {sum(i['region_excluded'][e] for i in infos)} participants")
        else:
            summarize_analysis(a.method, log)
    except Exception as e:
        failed += 1
        log("SUMMARY FAILED: " + "".join(traceback.format_exception(e)).strip())
    append_execution_log([
        f"## {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} ICA v2 {a.command} {a.method} (step 16{'b' if a.command == 'apply' else 'c'}, docs/ANALYSIS_FREEZE_v2_EEG.md)",
        f"- Command: python -m hb.ica_v2_apply {' '.join(sys.argv[1:])}; seed {SEED}; Python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}, mne {mne.__version__}",
        f"- Failures: {failed}; {(time.time() - t0) / 60:.1f} min. Log logs/ica_v2_{a.command}_{a.method}_{stamp}.log",
    ])
    return 0 if failed == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
