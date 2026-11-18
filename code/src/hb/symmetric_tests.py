"""Step 18 (docs/STEP18_PLAN_v1.md): symmetric tests.

Every test so far asked whether an apparent coupling could be explained without coupling. Step 18 asks the
same questions in the other direction, with the same weight:

  pilot      synthetic only: planted coupling level against the resulting observed excess, used once to fix
             the planted levels of the positive control (written into the plan before any other part runs)
  synthetic  (a) cross participant null validated on synthetic data without coupling (independent and matched
                 cohorts, K = 8), same gate as the trial shuffle null
             (b) positive control: matched synthetic cohort with a planted ongoing coupling at three levels,
                 calibrated with the trial shuffle null and with the cross participant null (K = 8 each)
  real       cross participant null on the real data (first pipeline v1, ICA v2, ICA v2 with muscle proxy
             regressed), K = 16, only if the gate of (a) passed
  inject     positive control on real data: the middle planted level added to the real ICA v2 heart period and
             envelopes; trial shuffle (same seeds as step 16c) and cross participant (same partners and seeds as
             'real') surrogates, K = 8

Cross participant surrogate ("heart of one person, brain of another"): heart unchanged. For surrogate k a
partner participant is drawn; every complete block of the participant receives the EEG envelope of a block of
the partner with the same sequence length (nearest length if the partner has none), warped piecewise linearly
so that block start, every trial start, every recall prompt and block end line up with the participant's own.
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
from hb import ica_v2_apply as IA
from hb import synth as S
from hb import trial_shuffle as TS
from hb.cohort import COLS, coherence_tables
from hb.cohort_stats import fmt, holm, participant_table
from hb.download_ds008302 import PROJECT_ROOT, append_execution_log
from hb.null_benchmark_real import FS
from hb.null_validation import SEED as SEED14
from hb.pilot_coherence import EEG, marked_share, zscore_over
from hb.pilot_qc import write_tsv

SEED18 = 20261120
PLANT_BANDS = ((0.04, 0.15), (0.15, 0.40))
PILOT_LEVELS = [0.0, 0.005, 0.01, 0.02, 0.04, 0.08]
# Fixed after the synthetic pilot, before any other part of step 18 was run (docs/STEP18_PLAN_v1.md).
LEVELS = {"low": None, "mid": None, "high": None}
PIPELINES = ["v1", "v2", "v2_emgadj"]
CELLS = TS.CELLS
OUT = lambda *p: PROJECT_ROOT.joinpath(*p)


def levels():
    f = OUT("results", "symmetric", "levels.json")
    if f.exists():
        return json.loads(f.read_text(encoding="utf-8"))["levels"]
    if any(v is None for v in LEVELS.values()):
        raise RuntimeError("planted levels not fixed yet: run the pilot and write results/symmetric/levels.json")
    return LEVELS


# ---------------------------------------------------------------- planted coupling
def shared_process(n, b10, rng):
    """Ongoing slow process with equal LF and HF power, unit variance inside the complete blocks."""
    s = sum(S.band_noise(n, FS, band, rng) for band in PLANT_BANDS)
    return zscore_over(s, b10)


def plant(x, s, c, b10):
    """x z scored inside the complete blocks, mixed with the shared process so that it carries a share c of
    the variance. c = 0 returns the z scored signal."""
    z = zscore_over(np.asarray(x, float), b10)
    return np.sqrt(1.0 - c) * z + np.sqrt(c) * s


def plant_all(sig, c, b10, rng):
    s = shared_process(len(sig["heart"]), b10, rng)
    return {k: plant(v, s, c, b10) for k, v in sig.items()}


# ---------------------------------------------------------------- cross participant surrogate
def anchors(b):
    """Block start, then trial start and recall prompt of every trial, then block end."""
    a = [b["t0"]]
    for s, p in zip(b["start"], b["probe"]):
        a += [s, p]
    a.append(b["t1"])
    return np.array(a, float)


def warp_block(y, b_self, x_src, b_src):
    """Replace the samples of block b_self in y by x_src warped from block b_src (in place)."""
    a_self, a_src = anchors(b_self), anchors(b_src)
    keep = np.r_[True, np.diff(a_self) > 0] & np.r_[True, np.diff(a_src) > 0]
    a_self, a_src = a_self[keep], a_src[keep]
    i0, i1 = int(round(b_self["t0"] * FS)), min(int(round(b_self["t1"] * FS)), y.size)
    if i1 <= i0:
        return
    t = np.arange(i0, i1) / FS
    src = np.interp(t, a_self, a_src)
    y[i0:i1] = np.interp(src, np.arange(x_src.size) / FS, x_src)


def choose_partner_blocks(b10_self, partner_b10, rng):
    """For each complete block of the participant, one complete block of the partner, same length if possible."""
    out, exact = [], 0
    lengths = np.array([b["length"] for b in partner_b10])
    for b in b10_self:
        same = np.flatnonzero(lengths == b["length"])
        if same.size:
            exact += 1
            out.append(int(rng.choice(same)))
        else:
            d = np.abs(lengths - b["length"])
            out.append(int(rng.choice(np.flatnonzero(d == d.min()))))
    return out, exact


def cross_surrogate(x_self, b10_self, x_partner, b10_partner, rng):
    y = np.asarray(x_self, float).copy()
    pick, exact = choose_partner_blocks(b10_self, b10_partner, rng)
    for b, j in zip(b10_self, pick):
        warp_block(y, b, np.asarray(x_partner, float), b10_partner[j])
    return y, exact


def draw_partner(sub, subjects, rng):
    others = [s for s in subjects if s != sub]
    return others[int(rng.integers(len(others)))]


# ---------------------------------------------------------------- shared helpers
def index_timing(root, subs):
    idx = {}
    for s in subs:
        blocks, b10, n = TS.timing(root, s)
        idx[s] = {"blocks": blocks, "b10": b10, "n": n}
    return idx


def table_rows(sub, heart, e, b, interp, blocks, marked):
    n = min(heart.size, b.size, interp.size)
    r = coherence_tables(sub, heart[:n], {e: b[:n]}, {e: np.asarray(interp[:n], bool)}, blocks, marked)
    return [x for x in participant_table({sub: [{c: q.get(c, "") for c in COLS} for q in r]}) if x["eeg"] == e]


def tag(rows, **kw):
    for r in rows:
        r.update(kw)
    return rows


# ---------------------------------------------------------------- synthetic signals
def synth_signals(cohort, sub, idx, matched, lv):
    """Signals of one synthetic participant, deterministic from the participant id and the cohort."""
    t = idx[sub]
    sid = int(sub.split("-")[1])
    if cohort == "independent":
        return TS.simulate(t["b10"], t["n"], np.random.default_rng([SEED14, sid, 0]))
    sig = TS.simulate(t["b10"], t["n"], np.random.default_rng([TS.SEED_VALID, sid, 1]), share=matched)
    if cohort == "matched":
        return sig
    li = list(lv).index(cohort.replace("plant_", ""))
    return plant_all(sig, lv[cohort.replace("plant_", "")], t["b10"], np.random.default_rng([SEED18, sid, 100 + li]))


def pilot_subject(args):
    root, sub, idx, matched = args
    t = idx[sub]
    sid = int(sub.split("-")[1])
    base = TS.simulate(t["b10"], t["n"], np.random.default_rng([TS.SEED_VALID, sid, 1]), share=matched)
    rows = []
    for li, c in enumerate(PILOT_LEVELS):
        sig = plant_all(base, c, t["b10"], np.random.default_rng([SEED18, sid, 900]))
        zeros = np.zeros(t["n"], bool)
        for e in EEG:
            rows += tag(table_rows(sub, sig["heart"], e, sig[e], zeros, t["blocks"], [0.0] * len(t["blocks"])),
                        cohort=f"pilot_{c}", kind="observed", null="", k=0)
    return sub, rows


def synthetic_subject(args):
    root, sub, idx, matched, lv, k, cohorts = args
    out = OUT("derivatives", "symmetric", "synthetic", f"{sub}.tsv")
    if out.exists():
        return sub, "skipped (already done)"
    t0 = time.time()
    t = idx[sub]
    sid = int(sub.split("-")[1])
    subjects = sorted(idx)
    zeros = np.zeros(t["n"], bool)
    marked = [0.0] * len(t["blocks"])
    rows, exact, total = [], 0, 0
    cache = {}
    for ci, cohort in enumerate(cohorts):
        sig = synth_signals(cohort, sub, idx, matched, lv)
        nulls = ["cross"] if cohort in ("independent", "matched") else ["shuffle", "cross"]
        for ei, e in enumerate(EEG):
            rows += tag(table_rows(sub, sig["heart"], e, sig[e], zeros, t["blocks"], marked), cohort=cohort, kind="observed", null="", k=0)
            for kk in range(k):
                if "shuffle" in nulls:
                    b = TS.trial_shuffle(sig[e], t["b10"], np.random.default_rng([SEED18, sid, ci, ei, kk, 1]))
                    rows += tag(table_rows(sub, sig["heart"], e, b, zeros, t["blocks"], marked), cohort=cohort, kind="surrogate", null="shuffle", k=kk)
                rng = np.random.default_rng([SEED18, sid, ci, ei, kk, 2])
                partner = draw_partner(sub, subjects, rng)
                if (cohort, partner) not in cache:
                    cache = {key: v for key, v in cache.items() if key[0] == cohort}
                    cache[(cohort, partner)] = synth_signals(cohort, partner, idx, matched, lv)
                b, ex = cross_surrogate(sig[e], t["b10"], cache[(cohort, partner)][e], idx[partner]["b10"], rng)
                exact, total = exact + ex, total + len(t["b10"])
                rows += tag(table_rows(sub, sig["heart"], e, b, zeros, t["blocks"], marked), cohort=cohort, kind="surrogate", null="cross", k=kk)
    TS.save(out, rows)
    return sub, f"ok ({len(rows)} rows, partner blocks of the same length {exact}/{total}, {time.time() - t0:.0f} s)"


# ---------------------------------------------------------------- real signals
_NPZ = {}


def npz(pipeline, sub):
    key = ("v1" if pipeline == "v1" else "v2", sub)
    if key not in _NPZ:
        if len(_NPZ) > 24:
            _NPZ.clear()
        f = OUT("derivatives", "slow", f"{sub}_task-memory_slow.npz") if key[0] == "v1" else \
            OUT("derivatives", "slow_v2", "picard", f"{sub}_task-memory_slow.npz")
        z = np.load(f)
        _NPZ[key] = {k: z[k] for k in z.files}
    return _NPZ[key]


def real_series(pipeline, sub, e, blocks):
    """Heart period, envelope and interpolation mask of one real participant for one pipeline."""
    z = npz(pipeline, sub)
    if pipeline == "v1":
        return np.asarray(z["nn"], float), np.asarray(z[e], float), np.asarray(z[e + "_interpolated"], bool)
    env, interp = IA.variant_signal(z, e, pipeline, blocks)
    return np.asarray(z["nn"], float), env, interp


def real_subject(args):
    root, sub, idx, k, pipelines = args
    out = OUT("derivatives", "symmetric", "real", f"{sub}.tsv")
    if out.exists():
        return sub, "skipped (already done)"
    t0 = time.time()
    t = idx[sub]
    sid = int(sub.split("-")[1])
    subjects = sorted(idx)
    marked = marked_share(OUT("derivatives", "beats", f"{sub}_task-memory_beats.tsv"), t["blocks"])
    rows, exact, total = [], 0, 0
    for pi, pipe in enumerate(pipelines):
        for ei, e in enumerate(EEG):
            nn, env, interp = real_series(pipe, sub, e, t["blocks"])
            rows += tag(table_rows(sub, nn, e, env, interp, t["blocks"], marked), cohort=pipe, kind="observed", null="", k=0)
            for kk in range(k):
                rng = np.random.default_rng([SEED18, sid, 50 + pi, ei, kk])
                partner = draw_partner(sub, subjects, rng)
                _, penv, _ = real_series(pipe, partner, e, idx[partner]["blocks"])
                b, ex = cross_surrogate(env, t["b10"], penv, idx[partner]["b10"], rng)
                exact, total = exact + ex, total + len(t["b10"])
                rows += tag(table_rows(sub, nn, e, b, interp, t["blocks"], marked), cohort=pipe, kind="surrogate", null="cross", k=kk)
    TS.save(out, rows)
    return sub, f"ok ({len(rows)} rows, partner blocks of the same length {exact}/{total}, {time.time() - t0:.0f} s)"


def inject_subject(args):
    root, sub, idx, k, c = args
    out = OUT("derivatives", "symmetric", "inject", f"{sub}.tsv")
    if out.exists():
        return sub, "skipped (already done)"
    t0 = time.time()
    t = idx[sub]
    sid = int(sub.split("-")[1])
    subjects = sorted(idx)
    marked = marked_share(OUT("derivatives", "beats", f"{sub}_task-memory_beats.tsv"), t["blocks"])
    pipe, pi, vi = "v2", PIPELINES.index("v2"), IA.VARIANTS_V2.index("v2")
    s_self = None
    rows = []
    for ei, e in enumerate(EEG):
        nn, env, interp = real_series(pipe, sub, e, t["blocks"])
        n = min(nn.size, env.size)
        if s_self is None:
            s_self = shared_process(n, t["b10"], np.random.default_rng([SEED18, sid, 700]))
        h_inj, e_inj = plant(nn[:n], s_self, c, t["b10"]), plant(env[:n], s_self, c, t["b10"])
        rows += tag(table_rows(sub, h_inj, e, e_inj, interp, t["blocks"], marked), cohort="v2_injected", kind="observed", null="", k=0)
        for kk in range(k):
            b = TS.trial_shuffle(e_inj, t["b10"], np.random.default_rng([IA.SEED, sid, vi, ei, kk]))
            rows += tag(table_rows(sub, h_inj, e, b, interp, t["blocks"], marked), cohort="v2_injected", kind="surrogate", null="shuffle", k=kk)
            rng = np.random.default_rng([SEED18, sid, 50 + pi, ei, kk])
            partner = draw_partner(sub, subjects, rng)
            pn, penv, _ = real_series(pipe, partner, e, idx[partner]["blocks"])
            m = min(pn.size, penv.size)
            sp = shared_process(m, idx[partner]["b10"], np.random.default_rng([SEED18, int(partner.split("-")[1]), 700]))
            b, _ = cross_surrogate(e_inj, t["b10"], plant(penv[:m], sp, c, idx[partner]["b10"]), idx[partner]["b10"], rng)
            rows += tag(table_rows(sub, h_inj, e, b, interp, t["blocks"], marked), cohort="v2_injected", kind="surrogate", null="cross", k=kk)
    TS.save(out, rows)
    return sub, f"ok ({len(rows)} rows, {time.time() - t0:.0f} s)"


# ---------------------------------------------------------------- summaries
def load(*folder):
    rows = []
    for f in sorted(OUT("derivatives", *folder).glob("sub-*.tsv")):
        if f.name.endswith("_realism.tsv"):
            continue
        with open(f, encoding="utf-8", newline="") as fh:
            rows += list(csv.DictReader(fh, delimiter="\t"))
    return rows


def cell_family(rows, cohort, null, variant="full", k_max=None):
    _, per = TS.calibrated_matrix(rows, cohort, null, variant, k_max)
    fam = []
    for e, b in CELLS:
        d = per[(e, b)]
        if d.size < 3:
            fam.append({"cohort": cohort, "null": null, "variant": variant, "k_used": k_max or "all", "eeg": e, "band": b, "n": d.size,
                        "mean": np.nan, "ci95_low": np.nan, "ci95_high": np.nan, "cohen_d": np.nan, "p": np.nan})
            continue
        h = stats.t.ppf(0.975, d.size - 1) * d.std(ddof=1) / np.sqrt(d.size)
        fam.append({"cohort": cohort, "null": null, "variant": variant, "k_used": k_max or "all", "eeg": e, "band": b, "n": d.size,
                    "mean": d.mean(), "ci95_low": d.mean() - h, "ci95_high": d.mean() + h,
                    "cohen_d": d.mean() / d.std(ddof=1), "p": float(stats.ttest_1samp(d, 0.0).pvalue)})
    ps = np.array([r["p"] for r in fam])
    adj = np.full(ps.size, np.nan)
    ok = np.isfinite(ps)
    adj[ok] = holm(ps[ok])
    for r, a in zip(fam, adj):
        r["p_holm"] = float(a)
    return fam


def shuffle_limits():
    return json.loads(OUT("results", "trial_shuffle", "gate.json").read_text(encoding="utf-8"))["bias_upper_edges"]


def summarize_pilot(rows, log):
    out = []
    obs = {}
    for r in rows:
        if r["variant"] == "full" and r["in_H1"] in ("1", "True"):
            obs.setdefault((r["cohort"], r["sub"], r["eeg"], r["band"]), float(r["mean_excess"]))
    subs = sorted({k[1] for k in obs})
    for c in PILOT_LEVELS:
        for e, b in CELLS:
            d = [obs[(f"pilot_{c}", s, e, b)] - obs[(f"pilot_{PILOT_LEVELS[0]}", s, e, b)] for s in subs
                 if (f"pilot_{c}", s, e, b) in obs and (f"pilot_{PILOT_LEVELS[0]}", s, e, b) in obs]
            out.append({"level": c, "eeg": e, "band": b, "n": len(d), "added_excess_mean": float(np.mean(d)), "added_excess_sd": float(np.std(d, ddof=1)) if len(d) > 1 else np.nan})
            log(f"pilot level {c:.3f} {e:5s} {b}: added observed excess {np.mean(d):+.4f} (n {len(d)})")
    OUT("results", "symmetric").mkdir(parents=True, exist_ok=True)
    write_tsv(OUT("results", "symmetric", "pilot_levels.tsv"), fmt(out))


def summarize_synthetic(rows, lv, log):
    table = []
    for cohort in ("independent", "matched"):
        table += cell_family(rows, cohort, "cross")
    limits = {}
    for r in table:
        key = f"{r['eeg']} {r['band']}"
        limits[key] = max(limits.get(key, -np.inf), r["ci95_high"])
    passed = all(min(r["p_holm"] for r in table if r["cohort"] == c) >= 0.05 for c in ("independent", "matched"))
    sh = shuffle_limits()
    for name in lv:
        for null, lim in (("shuffle", sh), ("cross", limits)):
            for r in cell_family(rows, f"plant_{name}", null):
                edge = lim.get(f"{r['eeg']} {r['band']}", np.nan)
                r["upper_bias_limit"] = edge
                r["detected_by_rule"] = int(r["p_holm"] < 0.05 and r["ci95_low"] > edge)
                table.append(r)
    for r in table:
        log(f"{r['cohort']:12s} {r['null']:7s} {r['eeg']:5s} {r['band']}: {r['mean']:+.4f} [{r['ci95_low']:+.4f}, {r['ci95_high']:+.4f}] "
            f"Holm p {r['p_holm']:.3g}" + (f" | limit {r['upper_bias_limit']:+.4f} detected {r['detected_by_rule']}" if "upper_bias_limit" in r else ""))
    OUT("results", "symmetric").mkdir(parents=True, exist_ok=True)
    write_tsv(OUT("results", "symmetric", "synthetic_cells.tsv"), fmt(table))
    OUT("results", "symmetric", "gate_cross.json").write_text(json.dumps(
        {"passed": bool(passed), "written": dt.datetime.now().strftime("%Y%m%d_%H%M%S"), "cross_upper_limits": limits}, indent=1), encoding="utf-8")
    log(f"GATE cross participant null (rule 1 of docs/STEP18_PLAN_v1.md): {'PASSED' if passed else 'FAILED'}")
    return passed


def summarize_real(rows, log):
    limits = json.loads(OUT("results", "symmetric", "gate_cross.json").read_text(encoding="utf-8"))["cross_upper_limits"]
    res = []
    for pipe in PIPELINES:
        f16, f8 = cell_family(rows, pipe, "cross"), cell_family(rows, pipe, "cross", k_max=8)
        for r16, r8 in zip(f16, f8):
            edge = limits.get(f"{r16['eeg']} {r16['band']}", np.nan)
            claim = r16["p_holm"] < 0.05 and r8["p_holm"] < 0.05 and r16["ci95_low"] > edge and r8["ci95_low"] > edge
            for r in (r16, r8):
                r["upper_bias_limit"], r["meets_claim_rule"] = edge, int(bool(claim))
            log(f"real {pipe:9s} {r16['eeg']:5s} {r16['band']}: cross K16 {r16['mean']:+.4f} [{r16['ci95_low']:+.4f}, {r16['ci95_high']:+.4f}] "
                f"Holm p {r16['p_holm']:.3g} | K8 Holm p {r8['p_holm']:.3g} low {r8['ci95_low']:+.4f} | limit {edge:+.4f}" + (" | CLAIM" if claim else ""))
            res += [r16, r8]
    write_tsv(OUT("results", "symmetric", "real_cross_calibrated_H1.tsv"), fmt(res))


def summarize_inject(log):
    inj = load("symmetric", "inject")
    base = {"shuffle": load("trial_shuffle_v2", "picard"), "cross": load("symmetric", "real")}
    lim = {"shuffle": shuffle_limits(), "cross": json.loads(OUT("results", "symmetric", "gate_cross.json").read_text(encoding="utf-8"))["cross_upper_limits"]}
    res = []
    for null in ("shuffle", "cross"):
        fam = cell_family(inj, "v2_injected", null, k_max=8)
        oi, _ = obs_cal(inj, "v2_injected", null)
        ob, _ = obs_cal(base[null], "v2", null)
        for r in fam:
            key = (r["eeg"], r["band"])
            common = sorted(set(oi[key]) & set(ob[key]))
            d = np.array([oi[key][s] - ob[key][s] for s in common])
            h = stats.t.ppf(0.975, d.size - 1) * d.std(ddof=1) / np.sqrt(d.size) if d.size > 2 else np.nan
            edge = lim[null].get(f"{r['eeg']} {r['band']}", np.nan)
            r.update(upper_bias_limit=edge, detected_by_rule=int(r["p_holm"] < 0.05 and r["ci95_low"] > edge),
                     delta_vs_real=d.mean() if d.size else np.nan, delta_ci95_low=d.mean() - h if d.size else np.nan,
                     delta_ci95_high=d.mean() + h if d.size else np.nan, n_paired=d.size)
            log(f"inject {null:7s} {r['eeg']:5s} {r['band']}: injected {r['mean']:+.4f} [{r['ci95_low']:+.4f}, {r['ci95_high']:+.4f}] Holm p {r['p_holm']:.3g} "
                f"limit {edge:+.4f} detected {r['detected_by_rule']} | added over real {r['delta_vs_real']:+.4f} [{r['delta_ci95_low']:+.4f}, {r['delta_ci95_high']:+.4f}] n {d.size}")
            res.append(r)
    write_tsv(OUT("results", "symmetric", "inject_cells.tsv"), fmt(res))


def obs_cal(rows, cohort, null, k_max=8):
    """Per participant calibrated excess {cell: {sub: value}} using the first k_max surrogates."""
    obs, sur = {}, {}
    for r in rows:
        if r["cohort"] != cohort or r["variant"] != "full":
            continue
        try:
            v = float(r["mean_excess"])
        except ValueError:
            continue
        if not np.isfinite(v):
            continue
        key = (r["eeg"], r["band"])
        if r["kind"] == "observed" and r["in_H1"] in ("1", "True"):
            obs.setdefault(key, {})[r["sub"]] = v
        elif r["kind"] == "surrogate" and r["null"] == null and int(r["k"]) < k_max:
            sur.setdefault(key, {}).setdefault(r["sub"], []).append(v)
    return {key: {s: obs[key][s] - np.mean(sur[key][s]) for s in obs.get(key, {}) if s in sur.get(key, {})} for key in CELLS}, None


# ---------------------------------------------------------------- driver
def run_pool(fn, args_list, workers, log, t0):
    failed, results = 0, []
    with ProcessPoolExecutor(max_workers=workers) as ex:
        futs = {ex.submit(fn, a): a[1] for a in args_list}
        for i, f in enumerate(as_completed(futs), 1):
            try:
                sub, status = f.result()
                results.append(status)
                if isinstance(status, list):
                    status = f"ok ({len(status)} rows)"
            except Exception as e:
                failed += 1
                sub, status = futs[f], "FAILED: " + "".join(traceback.format_exception(e)).strip().replace("\n", " | ")
            log(f"[{i}/{len(args_list)}] {sub}: {str(status)[:300]} ({(time.time() - t0) / 60:.1f} min elapsed)")
    return failed, results


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("command", choices=["pilot", "synthetic", "real", "inject"])
    ap.add_argument("--subjects", nargs="*")
    ap.add_argument("--k", type=int)
    ap.add_argument("--workers", type=int, default=max(1, min(6, (os.cpu_count() or 2) - 1)))
    ap.add_argument("--summary-only", action="store_true")
    a = ap.parse_args(argv)
    root = bidsio.raw_root()
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    OUT("logs").mkdir(exist_ok=True)
    logf = open(OUT("logs", f"symmetric_{a.command}_{stamp}.log"), "w", encoding="utf-8")

    def log(m):
        print(m, flush=True)
        logf.write(m + "\n")
        logf.flush()

    k = a.k or {"pilot": 0, "synthetic": 8, "real": 16, "inject": 8}[a.command]
    log(f"symmetric tests {a.command} {stamp} | k {k} workers {a.workers} seed {SEED18} | python {platform.python_version()} numpy {np.__version__} scipy {scipy.__version__}")
    t0, failed = time.time(), 0
    everyone = sorted(f.name.split("_")[0] for f in OUT("derivatives", "slow").glob("sub-*_task-memory_slow.npz"))
    subs = a.subjects or everyone
    try:
        idx = index_timing(root, everyone)
        matched = json.loads(OUT("results", "trial_shuffle", "matched_shares.json").read_text(encoding="utf-8"))
        if a.command == "pilot":
            subs = a.subjects or everyone[::10]
            failed, res = run_pool(pilot_subject, [(str(root), s, idx, matched) for s in subs], a.workers, log, t0)
            summarize_pilot([r for rows in res if isinstance(rows, list) for r in rows], log)
        elif a.command == "synthetic":
            lv = levels()
            log(f"planted levels (share of variance): {lv}")
            cohorts = ["independent", "matched"] + [f"plant_{n}" for n in lv]
            if not a.summary_only:
                failed, _ = run_pool(synthetic_subject, [(str(root), s, idx, matched, lv, k, cohorts) for s in subs], a.workers, log, t0)
            summarize_synthetic(load("symmetric", "synthetic"), lv, log)
        elif a.command in ("real", "inject"):
            g = OUT("results", "symmetric", "gate_cross.json")
            if not (g.exists() and json.loads(g.read_text(encoding="utf-8"))["passed"]) and a.command == "real":
                log("The cross participant gate has not passed (results/symmetric/gate_cross.json). Not running, as fixed in docs/STEP18_PLAN_v1.md.")
                return 3
            if a.command == "real":
                if not a.summary_only:
                    failed, _ = run_pool(real_subject, [(str(root), s, idx, k, PIPELINES) for s in subs], a.workers, log, t0)
                summarize_real(load("symmetric", "real"), log)
            else:
                c = levels()["mid"]
                log(f"injected level (share of variance): {c}")
                if not a.summary_only:
                    failed, _ = run_pool(inject_subject, [(str(root), s, idx, k, c) for s in subs], a.workers, log, t0)
                summarize_inject(log)
    except Exception as e:
        failed += 1
        log("FAILED: " + "".join(traceback.format_exception(e)).strip())
    append_execution_log([
        f"## {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} symmetric tests {a.command} (step 18, docs/STEP18_PLAN_v1.md)",
        f"- Command: python -m hb.symmetric_tests {' '.join(sys.argv[1:])}; seed {SEED18}; Python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}",
        f"- Failures: {failed}; {(time.time() - t0) / 60:.1f} min. Outputs: derivatives/symmetric, results/symmetric; log logs/symmetric_{a.command}_{stamp}.log",
    ])
    return 0 if failed == 0 else 2


if __name__ == "__main__":
    sys.exit(main())
