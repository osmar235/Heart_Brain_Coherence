"""Step 05c: null model benchmark rebuilt on the pilot's own response shapes and noise.

SURROGATE DATA WITH KNOWN TRUTH. For each pilot participant the real 4 Hz series (step 06) are split
into a trial locked part (average per sequence length, scaled to each trial) and a residual. Surrogate
sessions are then built so that heart and brain share NO ongoing coupling:
  - the trial locked parts are put back on the real trial timing, the brain trial gains shuffled across
    trials of the same length (breaks any shared trial to trial gain);
  - each residual is phase randomized within its block, independently for heart and brain (keeps its
    spectrum, destroys any real coupling).
Scenarios: task_only (as above), ongoing_LF (adds a shared 0.07 to 0.13 Hz process, 0.3 SD),
ongoing_HF (adds a shared 0.20 to 0.30 Hz process, 0.3 SD).
Responses locked to single digits are part of the real trial locked averages, so they are included.

The real, unshuffled heart and brain coherence is NOT computed here. No result about real coupling
comes out of this step, only error rates and sensitivity of the candidate controls.
"""
from __future__ import annotations

import argparse
import datetime as dt
import platform
import sys
import time

import numpy as np
import scipy
from scipy import stats

from hb import bidsio, hbcoh as H
from hb import synth as S
from hb.download_ds008302 import PROJECT_ROOT, append_execution_log, load_json
from hb.null_benchmark import load_timing, remove_evoked, seg, warp_to
from hb.vendor.emgsim_v3_3_0 import functions as emg

FS = H.FS_SLOW
BANDS = {"LF": (0.04, 0.15), "HF": (0.15, 0.40), "HFlow": (0.15, 0.33)}
CONTROLS = ["clock", "aligned", "resid_clock"]
SCENARIOS = {"task_only": None, "ongoing_LF": (0.07, 0.13), "ongoing_HF": (0.20, 0.30)}
REC_N = int(12 * FS)


def decompose(x, blocks):
    """Trial locked model of x: per sequence length an encoding and a recall average, one gain per trial each."""
    by_len = {}
    for bi, b in enumerate(blocks):
        by_len.setdefault(b["length"], []).append(bi)
    model = {}
    for length, idx in by_len.items():
        trials = [(bi, k) for bi in idx for k in range(10)]
        enc_n = int(round(np.median([blocks[bi]["probe"][k] - blocks[bi]["start"][k] for bi, k in trials]) * FS))
        enc = np.full((len(trials), enc_n), np.nan)
        rec = np.full((len(trials), REC_N), np.nan)
        for i, (bi, k) in enumerate(trials):
            b = blocks[bi]
            i0 = int(round(b["start"][k] * FS))
            m = min(enc_n, x.size - i0)
            enc[i, :m] = x[i0:i0 + m]
            j0 = int(round(b["probe"][k] * FS))
            m = max(0, min(REC_N, int(round(b["end"][k] * FS)) - j0, x.size - j0))
            rec[i, :m] = x[j0:j0 + m]
        enc_t = np.nanmean(enc, axis=0)
        cnt = np.isfinite(rec).sum(axis=0)
        rec_t = np.where(cnt >= 5, np.nanmean(np.where(cnt >= 5, rec, 0.0), axis=0), 0.0)
        gains = np.ones((len(trials), 2))
        for i in range(len(trials)):
            for j, (row, tpl) in enumerate(((enc[i], enc_t), (rec[i], rec_t))):
                ok = np.isfinite(row)
                tc = tpl[ok] - tpl[ok].mean()
                if ok.sum() > 8 and np.dot(tc, tc) > 1e-9:
                    gains[i, j] = np.clip(np.dot(row[ok] - row[ok].mean(), tc) / np.dot(tc, tc), -2.0, 4.0)
        model[length] = dict(trials=trials, enc_t=enc_t, rec_t=rec_t, gains=gains)
    return model


def synthesize(n, blocks, model, rng=None):
    """Put the trial locked model back on the timeline. With rng, gains are shuffled across trials of a length."""
    y = np.zeros(n)
    for length, m in model.items():
        g = m["gains"] if rng is None else m["gains"][rng.permutation(len(m["trials"]))]
        for i, (bi, k) in enumerate(m["trials"]):
            b = blocks[bi]
            i0 = int(round(b["start"][k] * FS))
            e = min(m["enc_t"].size, n - i0)
            y[i0:i0 + e] += g[i, 0] * m["enc_t"][:e]
            j0 = int(round(b["probe"][k] * FS))
            r = max(0, min(REC_N, int(round(b["end"][k] * FS)) - j0, n - j0))
            y[j0:j0 + r] += g[i, 1] * m["rec_t"][:r]
    return y


def block_slices(blocks, n):
    return [slice(int(round(b["t0"] * FS)), min(n, int(round(b["t1"] * FS)))) for b in blocks]


def zscore_inblock(x, blocks):
    m = np.zeros(x.size, bool)
    for s in block_slices(blocks, x.size):
        m[s] = True
    return (x - x[m].mean()) / x[m].std()


def band_values(cx, cy):
    coh = H.coherence(cx, cy)
    return {k: H.band_summary(coh, v)[0] for k, v in BANDS.items()}


def surrogate_session(parts, blocks, scenario, rng):
    n = parts["n"]
    sl = block_slices(blocks, n)
    shared = np.zeros(n)
    if SCENARIOS[scenario]:
        shared = 0.3 * S.band_noise(n, FS, SCENARIOS[scenario], rng)
    out = {}
    for name in ("heart", "brain"):
        ev = synthesize(n, blocks, parts[name + "_model"], rng if name == "brain" else None)
        res = np.zeros(n)
        for s in sl:
            res[s] = emg.phase_randomize(parts[name + "_resid"][s], rng)
        out[name] = ev + res + shared
    return out["heart"], out["brain"]


def evaluate(heart, brain, blocks, rng):
    groups = {}
    for i, b in enumerate(blocks):
        groups.setdefault(b["length"], []).append(i)
    partner = {}
    for idx in groups.values():
        if len(idx) >= 2:
            p = emg.random_derangement(len(idx), rng)
            for a, q in zip(idx, p):
                partner[a] = idx[q]
    hr, br = remove_evoked(heart, blocks, False), remove_evoked(brain, blocks, False)
    rows = []
    for i, b in enumerate(blocks):
        if i not in partner:
            continue
        j = partner[i]
        n = min(int(round((b["t1"] - b["t0"]) * FS)), int(round((blocks[j]["t1"] - blocks[j]["t0"]) * FS)))
        r = {}
        ch = H.cwt_slow(seg(heart, b))
        same = band_values(ch, H.cwt_slow(seg(brain, b)))
        aligned = band_values(ch, H.cwt_slow(warp_to(brain, blocks[j], b)))
        chn = H.cwt_slow(seg(heart, b, n))
        s_n = band_values(chn, H.cwt_slow(seg(brain, b, n)))
        c_n = band_values(chn, H.cwt_slow(seg(brain, blocks[j], n)))
        crn = H.cwt_slow(seg(hr, b, n))
        rs = band_values(crn, H.cwt_slow(seg(br, b, n)))
        rc = band_values(crn, H.cwt_slow(seg(br, blocks[j], n)))
        for band in BANDS:
            r[f"same_{band}"] = same[band]
            r[f"excess_aligned_{band}"] = same[band] - aligned[band]
            r[f"excess_clock_{band}"] = s_n[band] - c_n[band]
            r[f"excess_resid_clock_{band}"] = rs[band] - rc[band]
        rows.append(r)
    out = {"n_blocks": len(rows)}
    for k in rows[0]:
        v = np.array([r[k] for r in rows], float)
        v = v[np.isfinite(v)]
        out[k] = float(v.mean())
        if k.startswith("excess"):
            out["p_" + k] = float(stats.ttest_1samp(v, 0.0, alternative="greater").pvalue)
    return out


def prepare_parts(heart_real, brain_real, blocks):
    n = min(heart_real.size, brain_real.size)
    parts = {"n": n}
    for name, x in (("heart", heart_real[:n]), ("brain", brain_real[:n])):
        z = zscore_inblock(x, blocks)
        model = decompose(z, blocks)
        fitted = synthesize(n, blocks, model)
        parts[name + "_model"] = model
        parts[name + "_resid"] = z - fitted
        m = np.zeros(n, bool)
        for s in block_slices(blocks, n):
            m[s] = True
        parts[name + "_evoked_var_fraction"] = float(fitted[m].var() / z[m].var())
    return parts


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-sims", type=int, default=8)
    ap.add_argument("--seed", type=int, default=20261110)
    ap.add_argument("--eeg", nargs="*", default=["theta", "alpha", "beta"])
    a = ap.parse_args(argv)
    cfg = load_json(PROJECT_ROOT / "config" / "download_config.json")
    root = bidsio.raw_root()
    outdir = PROJECT_ROOT / "results" / "null_benchmark_real"
    outdir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    head = [f"null benchmark on real response shapes | n_sims {a.n_sims} | seed {a.seed} | python {platform.python_version()} "
            f"numpy {np.__version__} scipy {scipy.__version__}",
            "SURROGATE SESSIONS WITH NO ONGOING COUPLING BY CONSTRUCTION. Real heart and brain coherence is not computed."]
    print("\n".join(head), flush=True)
    sessions, evf = [], []
    for si, sub in enumerate(cfg["pilot"]["subjects"]):
        timing = load_timing(root, sub)
        z = np.load(PROJECT_ROOT / "derivatives" / "slow" / f"{sub}_task-memory_slow.npz")
        for ei, eeg in enumerate(a.eeg):
            parts = prepare_parts(z["nn"], z[eeg], timing["blocks"])
            evf.append(f"{sub}\t{eeg}\theart {parts['heart_evoked_var_fraction']:.3f}\tbrain {parts['brain_evoked_var_fraction']:.3f}")
            for ci, scenario in enumerate(SCENARIOS):
                for k in range(a.n_sims):
                    rng = np.random.default_rng([a.seed, si, ei, ci, k])
                    h, b = surrogate_session(parts, timing["blocks"], scenario, rng)
                    r = evaluate(h, b, timing["blocks"], rng)
                    r.update(sub=sub, eeg=eeg, scenario=scenario, sim=k)
                    sessions.append(r)
            print(f"  {sub} {eeg}: done ({time.time() - t0:.0f} s) | trial locked share of variance: heart "
                  f"{parts['heart_evoked_var_fraction']:.3f}, {eeg} {parts['brain_evoked_var_fraction']:.3f}", flush=True)
    idc = ["sub", "eeg", "scenario", "sim", "n_blocks"]
    cols = idc + sorted(k for k in sessions[0] if k not in idc)
    with open(outdir / "null_benchmark_real_sessions.tsv", "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in sessions:
            fh.write("\t".join(f"{r[c]:.5f}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
    table = ["scenario\teeg\tband\tcontrol\tn_sessions\tmean_same\tmean_excess\tsd_excess\tshare_sessions_p_below_0.05"]
    for scenario in SCENARIOS:
        for eeg in a.eeg + ["all"]:
            ss = [r for r in sessions if r["scenario"] == scenario and (eeg == "all" or r["eeg"] == eeg)]
            for band in BANDS:
                for c in CONTROLS:
                    e = np.array([r[f"excess_{c}_{band}"] for r in ss])
                    p = np.array([r[f"p_excess_{c}_{band}"] for r in ss])
                    same = np.mean([r[f"same_{band}"] for r in ss])
                    table.append(f"{scenario}\t{eeg}\t{band}\t{c}\t{len(ss)}\t{same:.3f}\t{e.mean():+.4f}\t{e.std():.4f}\t{(p < 0.05).mean():.3f}")
    (outdir / "null_benchmark_real_summary.tsv").write_text("\n".join(table) + "\n", encoding="utf-8")
    (outdir / "trial_locked_variance_share.tsv").write_text("sub\teeg\theart\tbrain\n" + "\n".join(evf) + "\n", encoding="utf-8")
    print("\n".join(t for t in table if "\tall\t" in t or t.startswith("scenario")), flush=True)
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    (PROJECT_ROOT / "logs").mkdir(exist_ok=True)
    (PROJECT_ROOT / "logs" / f"null_benchmark_real_{stamp}.log").write_text("\n".join(head + evf + table) + "\n", encoding="utf-8")
    append_execution_log([
        f"## {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} null benchmark on real response shapes (step 05c, surrogates only)",
        f"- Command: python -m hb.null_benchmark_real {' '.join(sys.argv[1:])}",
        f"- Python {platform.python_version()}, numpy {np.__version__}, scipy {scipy.__version__}, {platform.platform()}",
        f"- {len(sessions)} surrogate sessions, seed {a.seed}, {time.time() - t0:.0f} s. Real heart and brain coherence not computed.",
        f"- Outputs: results/null_benchmark_real/*.tsv; log logs/null_benchmark_real_{stamp}.log",
    ])
    print(f"DONE in {time.time() - t0:.0f} s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
