"""Step 15 (docs/STEP15_PLAN_v1.md): realism, decomposition and a trial shuffle null.

Subcommands
  realism     15a part 1: trial locked variance share (R2) and residual variance modulation (VMI),
              real data and step 14 synthetic data; writes the matched evoked shares
  decompose   15a part 2: step 14 bias with model components switched off (61 participants, K = 4)
  validate    15b: trial shuffle calibrated test on synthetic data (independent and matched), K = 8
  real        15c: trial shuffle calibrated test on the real data, K = 16 (only if 15b passed)
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

from hb import bidsio
from hb import synth as S
from hb.cohort import COLS, coherence_tables
from hb.cohort_stats import fmt, holm, participant_table
from hb.download_ds008302 import PROJECT_ROOT, append_execution_log
from hb.null_benchmark import gamma_kernel, remove_evoked
from hb.null_benchmark_real import FS, prepare_parts, surrogate_session
from hb.null_validation import EVOKED_SHARE, GAIN_SD, LAT_MAX, LAT_SD, TAU, analyse, error_rates, num, slow_modulation
from hb.null_validation import SEED as SEED14
from hb.pilot_coherence import BANDS, EEG, all_blocks, marked_share, zscore_over
from hb.pilot_qc import write_tsv

SEED_DECOMP, SEED_VALID, SEED_REAL = 20261117, 20261115, 20261116
SIGNALS = ["heart"] + EEG
CELLS = [(e, b) for e in EEG for b in BANDS]
VARIANTS = {
    "base": {},
    "no_gain": {"gain_sd": 0.0},
    "no_latency": {"lat_sd": 0.0},
    "stationary": {"modulate": False},
    "no_recall_scaling": {"recall_scale": False},
    "mean_only": {"gain_sd": 0.0, "lat_sd": 0.0, "modulate": False, "recall_scale": False},
}
OUT = lambda *p: PROJECT_ROOT.joinpath(*p)


# ---------------------------------------------------------------- synthetic model (step 14, parameterised)
def simulate(blocks10, n, rng, share=None, gain_sd=GAIN_SD, lat_sd=LAT_SD, modulate=True, recall_scale=True):
    """Independent scenario of hb.null_validation.simulate with switchable parts.
    With the defaults it reproduces that function draw for draw; every variant consumes the same draws."""
    share = share or EVOKED_SHARE
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
    for name in SIGNALS:
        if name == "heart":
            z, lat = z_h, lat_h
        else:
            z = rng.standard_normal((len(trials), 2))
            lat = np.clip(LAT_SD * rng.standard_normal(len(trials)), -LAT_MAX, LAT_MAX)
        lat = lat * (lat_sd / LAT_SD)
        sign = 1.0 if name == "heart" else -1.0
        box = np.zeros(n)
        for i, (s, p, e) in enumerate(trials):
            enc = (t >= s + 2.5 + lat[i]) & (t < p + lat[i])
            rec = (t >= p + lat[i]) & (t < e + lat[i])
            box += sign * (1 + gain_sd * z[i, 0]) * enc
            box -= sign * 0.7 * (1 + gain_sd * z[i, 1]) * (np.sqrt(dur[i] / dmed) if recall_scale else 1.0) * rec
        ev = np.convolve(box, gamma_kernel(TAU[name]))[:n]
        ev = (ev - ev[inblock].mean()) / ev[inblock].std()
        noise = S.band_noise(n, FS, (0.01, 1.0), rng)
        mod = slow_modulation(n, rng)
        if modulate:
            noise = noise * mod
        noise = (noise - noise[inblock].mean()) / noise[inblock].std()
        a = np.sqrt(share[name])
        out[name] = a * ev + np.sqrt(1 - a ** 2) * noise
    return out


def timing(root, sub):
    p = bidsio.rec_paths(Path(root), sub, "memory")
    dur = float(bidsio.sidecar(p["json"])["RecordingDuration"])
    blocks, _ = all_blocks(Path(root), sub)
    return blocks, [b for b in blocks if b["n_trials"] == 10], int(dur * FS)


# ---------------------------------------------------------------- trial shuffle surrogate
def derangement(m, rng):
    while True:
        p = rng.permutation(m)
        if not np.any(p == np.arange(m)):
            return p


def trial_shuffle(x, blocks10, rng):
    """Within each block, the signal of each trial (start to next start) is replaced by that of another trial of
    the same block, warped so that start, recall prompt and end line up with the target trial."""
    y = x.copy()
    tt = np.arange(x.size) / FS
    for b in blocks10:
        perm = derangement(10, rng)
        for k in range(10):
            j = perm[k]
            s, p, e = b["start"][k], b["probe"][k], b["end"][k]
            s2, p2, e2 = b["start"][j], b["probe"][j], b["end"][j]
            i0, i1 = int(round(s * FS)), min(int(round(e * FS)), x.size)
            if i1 <= i0:
                continue
            ti = tt[i0:i1]
            src = np.where(ti < p, s2 + (ti - s) * (p2 - s2) / max(p - s, 1e-9), p2 + (ti - p) * (e2 - p2) / max(e - p, 1e-9))
            y[i0:i1] = np.interp(src, tt, x)
    return y


# ---------------------------------------------------------------- realism metrics
def realism_metrics(x, blocks10):
    """R2 of the trial locked averages and variance modulation index of the residual."""
    z = zscore_over(x, blocks10)
    r = remove_evoked(z, blocks10, False)
    k = np.zeros(x.size, bool)
    for b in blocks10:
        k[int(round(b["t0"] * FS)):int(round(b["t1"] * FS))] = True
    r2 = 1.0 - r[k].var() / z[k].var()
    bins = {}
    by_len = {}
    for b in blocks10:
        by_len.setdefault(b["length"], []).append(b)
    for length, bl in by_len.items():
        enc_s = np.median([p - s for b in bl for s, p in zip(b["start"], b["probe"])])
        for b in bl:
            for s, p, e in zip(b["start"], b["probe"], b["end"]):
                for i in range(int(enc_s)):
                    a = int(round((s + i) * FS))
                    bins.setdefault(("enc", length, i), []).append(np.mean(r[a:a + int(FS)] ** 2))
                for i in range(12):
                    if p + i + 1 > e:
                        break
                    a = int(round((p + i) * FS))
                    bins.setdefault(("rec", length, i), []).append(np.mean(r[a:a + int(FS)] ** 2))
    m = np.array([np.mean(v) for v in bins.values() if len(v) >= 5])
    return float(r2), float(m.std() / m.mean())


def realism_subject(args):
    root, sub = args
    blocks, b10, n = timing(root, sub)
    z = np.load(OUT("derivatives", "slow", f"{sub}_task-memory_slow.npz"))
    real = {"heart": z["nn"]} | {e: z[e] for e in EEG}
    syn = simulate(b10, n, np.random.default_rng([SEED14, int(sub.split("-")[1]), 0]))
    rows = []
    for name in SIGNALS:
        for source, x in (("real", real[name]), ("synthetic", syn[name])):
            r2, vmi = realism_metrics(np.asarray(x, float), b10)
            rows.append({"sub": sub, "signal": name, "source": source, "R2": r2, "VMI": vmi})
    return sub, rows


def summarize_realism(rows, log):
    outdir = OUT("results", "trial_shuffle")
    outdir.mkdir(parents=True, exist_ok=True)
    table, matched = [], {}
    for name in SIGNALS:
        rec = {"signal": name}
        for metric in ("R2", "VMI"):
            for source in ("real", "synthetic"):
                v = np.array([r[metric] for r in rows if r["signal"] == name and r["source"] == source and np.isfinite(r[metric])])
                rec[f"{metric}_{source}_median"] = np.median(v)
                rec[f"{metric}_{source}_q1"], rec[f"{metric}_{source}_q3"] = np.percentile(v, [25, 75])
            rec[f"{metric}_ratio_real_over_synthetic"] = rec[f"{metric}_real_median"] / rec[f"{metric}_synthetic_median"]
        matched[name] = float(min(0.9, EVOKED_SHARE[name] * rec["R2_ratio_real_over_synthetic"]))
        rec["step14_share"], rec["matched_share"] = EVOKED_SHARE[name], matched[name]
        table.append(rec)
        log(f"realism {name:6s}: R2 real {rec['R2_real_median']:.3f} synthetic {rec['R2_synthetic_median']:.3f} ratio {rec['R2_ratio_real_over_synthetic']:.2f} | "
            f"VMI real {rec['VMI_real_median']:.3f} synthetic {rec['VMI_synthetic_median']:.3f} ratio {rec['VMI_ratio_real_over_synthetic']:.2f} | "
            f"share {EVOKED_SHARE[name]:.2f} -> matched {matched[name]:.3f}")
    write_tsv(outdir / "realism_summary.tsv", fmt(table))
    write_tsv(outdir / "realism_participants.tsv", fmt(rows))
    (outdir / "matched_shares.json").write_text(json.dumps(matched, indent=1), encoding="utf-8")


# ---------------------------------------------------------------- per participant runs
def cal_rows(sub, sig, b10, blocks, n, k, seed, tag, null):
    """Observed and k surrogate participant tables for one synthetic or real session."""
    rows = []
    for r in analyse(sub, sig["heart"], {e: sig[e] for e in EEG}, blocks, n):
        r.update(cohort=tag, kind="observed", null="", k=0)
        rows.append(r)
    for ei, e in enumerate(EEG):
        parts = prepare_parts(sig["heart"], sig[e], b10) if null == "phase" else None
        for kk in range(k):
            rng = np.random.default_rng(seed + [ei, kk])
            if null == "phase":
                h, b = surrogate_session(parts, b10, "task_only", rng)
                nn = parts["n"]
            else:
                h, b, nn = sig["heart"], trial_shuffle(sig[e], b10, rng), n
            for r in analyse(sub, h, {e: b}, blocks, nn):
                r.update(cohort=tag, kind="surrogate", null=null, k=kk)
                rows.append(r)
    return rows


RCOLS = ["cohort", "kind", "null", "k", "sub", "eeg", "variant", "band", "n_blocks", "mean_excess", "in_H1"]


def save(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".part")
    write_tsv(tmp, fmt([{c: r.get(c, "") for c in RCOLS} for r in rows]))
    tmp.replace(path)


def decompose_subject(args):
    root, sub, k = args
    out = OUT("derivatives", "trial_shuffle", "decompose", f"{sub}.tsv")
    if out.exists():
        return sub, "skipped (already done)"
    t0 = time.time()
    blocks, b10, n = timing(root, sub)
    sid = int(sub.split("-")[1])
    rows = []
    for vi, (v, kw) in enumerate(VARIANTS.items()):
        sig = simulate(b10, n, np.random.default_rng([SEED14, sid, 0]), **kw)
        rows += cal_rows(sub, sig, b10, blocks, n, k, [SEED_DECOMP, sid, vi], v, "phase")
    save(out, rows)
    return sub, f"ok ({len(rows)} rows, {time.time() - t0:.0f} s)"


def validate_subject(args):
    root, sub, k, matched = args
    out = OUT("derivatives", "trial_shuffle", "validate", f"{sub}.tsv")
    if out.exists():
        return sub, "skipped (already done)"
    t0 = time.time()
    blocks, b10, n = timing(root, sub)
    sid = int(sub.split("-")[1])
    rows = []
    sig = simulate(b10, n, np.random.default_rng([SEED14, sid, 0]))
    rows += cal_rows(sub, sig, b10, blocks, n, k, [SEED_VALID, sid, 0], "independent", "shuffle")
    sig = simulate(b10, n, np.random.default_rng([SEED_VALID, sid, 1]), share=matched)
    rm = cal_rows(sub, sig, b10, blocks, n, k, [SEED_VALID, sid, 1], "matched", "shuffle")
    rp = cal_rows(sub, sig, b10, blocks, n, k, [SEED_VALID, sid, 2], "matched", "phase")
    rows += rm + [r for r in rp if r["kind"] == "surrogate"]
    save(out, rows)
    return sub, f"ok ({len(rows)} rows, {time.time() - t0:.0f} s)"


def real_subject(args):
    root, sub, k = args
    out = OUT("derivatives", "trial_shuffle", "real", f"{sub}.tsv")
    if out.exists():
        return sub, "skipped (already done)"
    t0 = time.time()
    blocks, b10, _ = timing(root, sub)
    z = np.load(OUT("derivatives", "slow", f"{sub}_task-memory_slow.npz"))
    marked = marked_share(OUT("derivatives", "beats", f"{sub}_task-memory_beats.tsv"), blocks)
    sid = int(sub.split("-")[1])
    nn = z["nn"]
    rows = []

    def tab(heart, e, b):
        n = min(heart.size, b.size)
        r = coherence_tables(sub, heart[:n], {e: b[:n]}, {e: z[e + "_interpolated"][:n]}, blocks, marked)
        return participant_table({sub: [{c: x.get(c, "") for c in COLS} for x in r]})

    for ei, e in enumerate(EEG):
        for r in tab(nn, e, z[e]):
            if r["eeg"] == e:
                r.update(cohort="real", kind="observed", null="", k=0)
                rows.append(r)
        for kk in range(k):
            b = trial_shuffle(np.asarray(z[e], float), b10, np.random.default_rng([SEED_REAL, sid, ei, kk]))
            for r in tab(nn, e, b):
                if r["eeg"] == e:
                    r.update(cohort="real", kind="surrogate", null="shuffle", k=kk)
                    rows.append(r)
    save(out, rows)
    return sub, f"ok ({len(rows)} rows, {time.time() - t0:.0f} s)"


# ---------------------------------------------------------------- summaries
def load(folder):
    rows = []
    for f in sorted(OUT("derivatives", "trial_shuffle", folder).glob("sub-*.tsv")):
        with open(f, encoding="utf-8", newline="") as fh:
            rows += list(csv.DictReader(fh, delimiter="\t"))
    return rows


def calibrated_matrix(rows, cohort, null, variant, k_max=None):
    obs, sur = {}, {}
    for r in rows:
        if r["cohort"] != cohort or r["variant"] != variant or not np.isfinite(num(r["mean_excess"])):
            continue
        key = (r["sub"], r["eeg"], r["band"])
        if r["kind"] == "observed" and r["in_H1"] in ("1", "True"):
            obs[key] = num(r["mean_excess"])
        elif r["kind"] == "surrogate" and r["null"] == null and (k_max is None or int(r["k"]) < k_max):
            sur.setdefault(key, []).append(num(r["mean_excess"]))
    D = [[obs[(s, e, b)] - np.mean(sur[(s, e, b)]) for e, b in CELLS]
         for s in sorted({k[0] for k in obs}) if all((s, e, b) in obs and (s, e, b) in sur for e, b in CELLS)]
    per_cell = {(e, b): np.array([obs[k] - np.mean(sur[k]) for k in obs if k[1] == e and k[2] == b and k in sur]) for e, b in CELLS}
    return np.array(D), per_cell


def synthetic_summary(rows, folder, specs, log):
    """specs: list of (cohort, null, k_max). Writes per cell bias and Monte Carlo error rates."""
    rng = np.random.default_rng(SEED_VALID)
    table = []
    for cohort, null, k_max in specs:
        for v in ("full", "encret"):
            D, _ = calibrated_matrix(rows, cohort, null, v, k_max)
            if D.size == 0:
                continue
            mu, half, res = error_rates(D, 122, rng)
            p_bias = [float(stats.ttest_1samp(D[:, j], 0.0).pvalue) for j in range(len(CELLS))]
            p_bias_holm = holm(p_bias)
            for j, (e, b) in enumerate(CELLS):
                table.append({"cohort": cohort, "null": null, "variant": v, "eeg": e, "band": b, "n_synthetic": D.shape[0],
                              "bias": mu[j], "ci95_low": mu[j] - half[j], "ci95_high": mu[j] + half[j],
                              "p_bias": p_bias[j], "p_bias_holm": float(p_bias_holm[j]),
                              "error_rate_at_zero_bias": res["at_zero_bias"][0][j],
                              "error_rate_at_estimated_bias": res["at_estimated_bias"][0][j],
                              "error_rate_at_95pct_edge": res["at_95pct_edge_of_bias"][0][j],
                              "holm_family_at_zero_bias": res["at_zero_bias"][1],
                              "holm_family_at_estimated_bias": res["at_estimated_bias"][1],
                              "holm_family_at_95pct_edge": res["at_95pct_edge_of_bias"][1]})
            log(f"{cohort:18s} {null:7s} {v:6s} n {D.shape[0]}: bias Holm p min {min(p_bias_holm):.3g} | Holm family error {res['at_estimated_bias'][1]:.3f} at estimated bias, "
                f"{res['at_95pct_edge_of_bias'][1]:.3f} at edge, {res['at_zero_bias'][1]:.3f} at zero | bias "
                + ", ".join(f"{e} {b} {mu[j]:+.4f} [{mu[j] - half[j]:+.4f}, {mu[j] + half[j]:+.4f}]" for j, (e, b) in enumerate(CELLS)))
    outdir = OUT("results", "trial_shuffle")
    outdir.mkdir(parents=True, exist_ok=True)
    write_tsv(outdir / f"{folder}_cells.tsv", fmt(table))
    return table


def gate(table):
    """Rule 1 of docs/STEP15_PLAN_v1.md: in both synthetic cohorts, no cell of the shuffle calibrated bias (primary
    variant) is Holm significant against zero."""
    for cohort in ("independent", "matched"):
        r = [t for t in table if t["cohort"] == cohort and t["null"] == "shuffle" and t["variant"] == "full"]
        if len(r) != len(CELLS) or min(t["p_bias_holm"] for t in r) < 0.05:
            return False
    return True


def bias_edges(table):
    """Upper edge of the shuffle calibrated synthetic bias per cell, the larger of the two cohorts (rule 2)."""
    out = {}
    for t in table:
        if t["null"] == "shuffle" and t["variant"] == "full" and t["cohort"] in ("independent", "matched"):
            out[(t["eeg"], t["band"])] = max(out.get((t["eeg"], t["band"]), -np.inf), t["ci95_high"])
    return out


def real_summary(rows, log, edges=None):
    out = []
    for k_max in (8, None):
        for v in ("full", "encret"):
            fam = []
            for e, b in CELLS:
                _, per = calibrated_matrix(rows, "real", "shuffle", v, k_max)
                d = per[(e, b)]
                t = stats.ttest_1samp(d, 0.0)
                h = stats.t.ppf(0.975, d.size - 1) * d.std(ddof=1) / np.sqrt(d.size)
                fam.append({"variant": v, "k_used": k_max or 16, "eeg": e, "band": b, "n": d.size, "calibrated_excess": d.mean(),
                            "ci95_low": d.mean() - h, "ci95_high": d.mean() + h, "cohen_d": d.mean() / d.std(ddof=1), "p": float(t.pvalue)})
            for r, a in zip(fam, holm([r["p"] for r in fam])):
                r["p_holm"] = float(a)
                edge = (edges or {}).get(f"{r['eeg']} {r['band']}", float("nan"))
                r["synthetic_bias_upper_edge"] = edge
                r["lower_bound_above_edge"] = int(np.isfinite(edge) and r["ci95_low"] > edge)
                log(f"real shuffle K={r['k_used']:2d} {v:6s} {r['eeg']:5s} {r['band']}: calibrated {r['calibrated_excess']:+.4f} "
                    f"[{r['ci95_low']:+.4f}, {r['ci95_high']:+.4f}] d {r['cohen_d']:.2f} Holm p {r['p_holm']:.3g} | synthetic bias edge {edge:+.4f}")
            out += fam
    frozen = {}
    with open(OUT("results", "cohort_stats", "participant_level.tsv"), encoding="utf-8", newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            frozen[(r["sub"], r["eeg"], r["variant"], r["band"])] = num(r["mean_excess"])
    diffs = [abs(num(r["mean_excess"]) - frozen[(r["sub"], r["eeg"], r["variant"], r["band"])]) for r in rows
             if r["kind"] == "observed" and (r["sub"], r["eeg"], r["variant"], r["band"]) in frozen
             and np.isfinite(num(r["mean_excess"])) and np.isfinite(frozen[(r["sub"], r["eeg"], r["variant"], r["band"])])]
    log(f"check: recomputed observed against frozen participant_level.tsv, max abs difference {max(diffs) if diffs else float('nan'):.3g} over {len(diffs)} values")
    OUT("results", "trial_shuffle").mkdir(parents=True, exist_ok=True)
    write_tsv(OUT("results", "trial_shuffle", "real_shuffle_calibrated_H1.tsv"), fmt(out))


# ---------------------------------------------------------------- driver
def run_pool(fn, args_list, workers, log, t0):
    failed, done, results = 0, 0, []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(fn, a): a[1] for a in args_list}
        for f in as_completed(futs):
            done += 1
            try:
                sub, status = f.result()
                results.append(status)
                if isinstance(status, list):
                    status = f"ok ({len(status)} rows)"
            except Exception as e:
                failed += 1
                sub, status = futs[f], "FAILED: " + "".join(traceback.format_exception(e)).strip().replace("\n", " | ")
            log(f"[{done}/{len(args_list)}] {sub}: {str(status)[:300]} ({(time.time() - t0) / 60:.1f} min elapsed)")
    return failed, results


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["realism", "decompose", "validate", "real"])
    ap.add_argument("--subjects", nargs="*")
    ap.add_argument("--k", type=int)
    ap.add_argument("--workers", type=int, default=max(1, min(6, (os.cpu_count() or 2) - 1)))
    ap.add_argument("--summary-only", action="store_true")
    ap.add_argument("--force", action="store_true", help="real: run even if the 15b gate failed (not allowed by the plan)")
    a = ap.parse_args(argv)
    root = bidsio.raw_root()
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    OUT("logs").mkdir(exist_ok=True)
    logf = open(OUT("logs", f"trial_shuffle_{a.command}_{stamp}.log"), "w", encoding="utf-8")

    def log(m):
        print(m, flush=True)
        logf.write(m + "\n")
        logf.flush()

    k = a.k or {"realism": 0, "decompose": 4, "validate": 8, "real": 16}[a.command]
    log(f"trial shuffle {a.command} {stamp} | k {k} workers {a.workers} | python {platform.python_version()} numpy {np.__version__} scipy {scipy.__version__}")
    t0, failed = time.time(), 0
    subs = a.subjects or sorted(f.name.split("_")[0] for f in OUT("derivatives", "slow").glob("sub-*_task-memory_slow.npz"))
    try:
        if a.command == "realism":
            failed, res = run_pool(realism_subject, [(str(root), s) for s in subs], a.workers, log, t0)
            summarize_realism([r for rows in res if isinstance(rows, list) for r in rows], log)
        elif a.command == "decompose":
            if not a.summary_only:
                failed, _ = run_pool(decompose_subject, [(str(root), s, k) for s in subs[::2]], a.workers, log, t0)
            synthetic_summary(load("decompose"), "decompose", [(v, "phase", None) for v in VARIANTS], log)
        elif a.command == "validate":
            matched = json.loads(OUT("results", "trial_shuffle", "matched_shares.json").read_text(encoding="utf-8"))
            log(f"matched shares (from 15a): {matched}")
            if not a.summary_only:
                failed, _ = run_pool(validate_subject, [(str(root), s, k, matched) for s in subs], a.workers, log, t0)
            table = synthetic_summary(load("validate"), "validate",
                                      [("independent", "shuffle", None), ("matched", "shuffle", None), ("matched", "phase", None)], log)
            ok = gate(table)
            edges = {f"{e} {b}": v for (e, b), v in bias_edges(table).items()}
            OUT("results", "trial_shuffle", "gate.json").write_text(json.dumps({"passed": ok, "written": stamp, "bias_upper_edges": edges}, indent=1), encoding="utf-8")
            log(f"GATE (rule 1 of docs/STEP15_PLAN_v1.md): {'PASSED, step 15c may be run' if ok else 'FAILED, step 15c is not to be run'}")
        elif a.command == "real":
            g = OUT("results", "trial_shuffle", "gate.json")
            passed = g.exists() and json.loads(g.read_text(encoding="utf-8"))["passed"]
            if not passed and not a.force:
                log("Step 15b gate has not passed (results/trial_shuffle/gate.json). Not running, as fixed in docs/STEP15_PLAN_v1.md.")
                return 3
            if not a.summary_only:
                failed, _ = run_pool(real_subject, [(str(root), s, k) for s in subs], a.workers, log, t0)
            real_summary(load("real"), log, json.loads(g.read_text(encoding="utf-8")).get("bias_upper_edges", {}) if g.exists() else {})
    except Exception as e:
        failed += 1
        log("SUMMARY FAILED: " + "".join(traceback.format_exception(e)).strip())
    append_execution_log([
        f"## {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} trial shuffle {a.command} (step 15, docs/STEP15_PLAN_v1.md)",
        f"- Command: python -m hb.trial_shuffle {' '.join(sys.argv[1:])}; Python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}",
        f"- Failures: {failed}; {(time.time() - t0) / 60:.1f} min. Outputs: derivatives/trial_shuffle, results/trial_shuffle; log logs/trial_shuffle_{a.command}_{stamp}.log",
    ])
    return 0 if failed == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
