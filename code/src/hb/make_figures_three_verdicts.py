"""Step 19: figures of the Three Verdicts manuscript, drawn only from result tables (no recording is read).

Figure 1  the three verdicts on one axis: initial analysis; trial shuffle after ICA and after the muscle proxy;
          cross participant null after ICA and after the muscle proxy (upper bias limits as dashes)
Figure 2  synthetic bias without coupling (phase against trial shuffle calibration) and its decomposition
Figure 3  positive control: planted coupling at three levels under both nulls, limits, and the observed values
Outputs: results/figures_three_verdicts/figure_<n>.png (300 dpi) and .pdf, figure_<n>_source_data.tsv (every
plotted value) and figure_sources.tsv (source tables with SHA256 prefixes).
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import json
import sys

import numpy as np

from hb import make_figures as MF
from hb.download_ds008302 import PROJECT_ROOT, append_execution_log
from hb.make_figures import ACCENT, CELLS, DARK, GREY, LABEL, f, pick



def OUT(*p):
    """Project path. HB_RESULTS_ROOT and HB_LOGS_ROOT redirect the 'results' and 'logs' trees (used by the release
    package, where the frozen result tables live in results/tables)."""
    import os
    from pathlib import Path
    for head, var in (("results", "HB_RESULTS_ROOT"), ("logs", "HB_LOGS_ROOT")):
        if p and p[0] == head and os.environ.get(var):
            return Path(os.environ[var]).joinpath(*p[1:])
    return PROJECT_ROOT.joinpath(*p)

ORANGE = "#c05621"
SRC = dict(MF.SRC)
SRC.update({
    "sym_synthetic": ("results", "symmetric", "synthetic_cells.tsv"),
    "sym_real": ("results", "symmetric", "real_cross_calibrated_H1.tsv"),
    "gate_cross": ("results", "symmetric", "gate_cross.json"),
})


def read(key):
    with open(OUT(*SRC[key]), encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def val(r):
    for k in ("mean", "calibrated_excess", "bias"):
        if k in r and r[k] not in ("", None):
            return f(r[k])
    raise KeyError("no value column")


def limits(which):
    if which == "shuffle":
        return json.loads(OUT(*SRC["gate"]).read_text(encoding="utf-8"))["bias_upper_edges"]
    return json.loads(OUT(*SRC["gate_cross"]).read_text(encoding="utf-8"))["cross_upper_limits"]


def bar(ax, x, r, data, panel, series, **kw):
    m, lo, hi = val(r), f(r["ci95_low"]), f(r["ci95_high"])
    ax.errorbar(x, m, yerr=[[m - lo], [hi - m]], capsize=0, lw=1.2, ms=4.5, **kw)
    data.append({"panel": panel, "series": series, "eeg": r["eeg"], "band": r["band"], "value": m, "ci95_low": lo, "ci95_high": hi})


def limit_marks(ax, lim, data, panel, name, half=0.36):
    for j, (e, b) in enumerate(CELLS):
        v = lim[f"{e} {b}"]
        ax.plot([j - half, j + half], [v, v], color=ORANGE, lw=1.0, ls="--")
        data.append({"panel": panel, "series": name, "eeg": e, "band": b, "value": v, "ci95_low": "", "ci95_high": ""})


def figure_1(plt):
    gt, ica, cr = read("group_tests"), read("ica"), read("sym_real")
    data = []
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.3), sharey=True)
    for j, (e, b) in enumerate(CELLS):
        bar(axes[0], j, pick(gt, hypothesis="H1_existence", variant="full", family="primary", eeg=e, band=b), data,
            "A", "initial analysis, between block reference", fmt="o", color=DARK)
        bar(axes[1], j - 0.15, pick(ica, variant_eeg="v2", variant="full", k_used="16", eeg=e, band=b), data,
            "B", "trial shuffle, ICA", fmt="o", color=ACCENT)
        bar(axes[1], j + 0.15, pick(ica, variant_eeg="v2_emgadj", variant="full", k_used="16", eeg=e, band=b), data,
            "B", "trial shuffle, ICA, muscle proxy", fmt="s", color=ACCENT, alpha=0.5)
        bar(axes[2], j - 0.15, pick(cr, cohort="v2", k_used="all", eeg=e, band=b), data,
            "C", "cross participant, ICA", fmt="o", color=ACCENT)
        bar(axes[2], j + 0.15, pick(cr, cohort="v2_emgadj", k_used="all", eeg=e, band=b), data,
            "C", "cross participant, ICA, muscle proxy", fmt="s", color=ACCENT, alpha=0.5)
    limit_marks(axes[1], limits("shuffle"), data, "B", "upper bias limit, trial shuffle")
    limit_marks(axes[2], limits("cross"), data, "C", "upper bias limit, cross participant")
    for ax, title in zip(axes, ("A  Verdict 1: initial analysis", "B  Verdict 2: trial shuffle null", "C  Verdict 3: cross participant null")):
        MF.cell_ticks(ax)
        ax.tick_params(axis="x", labelsize=7)
        ax.set_title(title, loc="left", fontsize=9)
    axes[0].set_ylabel("Excess squared coherence")
    h = [plt.Line2D([], [], marker="o", color=DARK, lw=0, ms=4.5, label="First pipeline, uncalibrated"),
         plt.Line2D([], [], marker="o", color=ACCENT, lw=0, ms=4.5, label="After ICA"),
         plt.Line2D([], [], marker="s", color=ACCENT, alpha=0.5, lw=0, ms=4.5, label="After ICA and muscle proxy"),
         plt.Line2D([], [], color=ORANGE, lw=1.0, ls="--", label="Upper bias limit (synthetic, no coupling)")]
    fig.legend(handles=h, loc="lower center", ncol=4, fontsize=7.5, bbox_to_anchor=(0.5, -0.06))
    fig.tight_layout()
    return fig, ["group_tests", "ica", "sym_real", "gate", "gate_cross"], data


def figure_2(plt):
    nv, va, ob, dc = read("null_validation"), read("validate"), read("phase_observed"), read("decompose")
    data = []
    fig = plt.figure(figsize=(7.2, 6.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.15, 1.0], hspace=0.45, wspace=0.12)
    ax = fig.add_subplot(gs[0, :])
    series = [("Phase randomized, independent model", "o", ACCENT, 1.0, lambda e, b: pick(nv, scenario="independent", variant="full", statistic="calibrated", eeg=e, band=b)),
              ("Phase randomized, matched model", "s", ACCENT, 0.55, lambda e, b: pick(va, cohort="matched", null="phase", variant="full", eeg=e, band=b)),
              ("Trial shuffle, independent model", "o", GREY, 1.0, lambda e, b: pick(va, cohort="independent", null="shuffle", variant="full", eeg=e, band=b)),
              ("Trial shuffle, matched model", "s", GREY, 0.6, lambda e, b: pick(va, cohort="matched", null="shuffle", variant="full", eeg=e, band=b)),
              ("Observed, phase calibrated", "D", DARK, 1.0, lambda e, b: pick(ob, variant="full", eeg=e, band=b))]
    off = np.linspace(-0.27, 0.27, len(series))
    for k, (name, mk, col, alpha, get) in enumerate(series):
        for j, (e, b) in enumerate(CELLS):
            r = dict(get(e, b))
            if "ci95_halfwidth" in r:
                r["ci95_low"], r["ci95_high"] = f(r["mean"]) - f(r["ci95_halfwidth"]), f(r["mean"]) + f(r["ci95_halfwidth"])
            r.setdefault("eeg", e), r.setdefault("band", b)
            kw = dict(fmt=mk, color=col, alpha=alpha)
            if mk == "D":
                kw.update(mfc="white", mec=DARK)
            bar(ax, j + off[k], r, data, "A", name, label=name if j == 0 else None, **kw)
    MF.cell_ticks(ax)
    ax.set_ylabel("Calibrated excess coherence")
    ax.set_title("A  Synthetic participants without coupling, real trial timing", loc="left", fontsize=9)
    ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=7.5)
    variants = [("base", "Full model"), ("no_gain", "No gain variability"), ("no_latency", "No latency jitter"),
                ("stationary", "Stationary noise"), ("no_recall_scaling", "No recall scaling"), ("mean_only", "Identical trials,\nstationary noise")]
    first = None
    for c, e in enumerate(("theta", "alpha", "beta")):
        a2 = fig.add_subplot(gs[1, c], sharex=first, sharey=first)
        first = first or a2
        for i, (v, name) in enumerate(variants):
            r = pick(dc, cohort=v, null="phase", variant="full", eeg=e, band="LF")
            m = f(r["bias"])
            a2.errorbar(m, i, xerr=[[m - f(r["ci95_low"])], [f(r["ci95_high"]) - m]], fmt="o", color=ACCENT if v in ("base", "mean_only") else GREY, ms=4, lw=1.1, capsize=0)
            data.append({"panel": "B", "series": name.replace("\n", " "), "eeg": e, "band": "LF", "value": m, "ci95_low": f(r["ci95_low"]), "ci95_high": f(r["ci95_high"])})
        a2.axvline(0, color="#555555", lw=0.8)
        a2.set_title(("B  " if c == 0 else "") + f"{LABEL[e]} LF", loc="left", fontsize=9)
        a2.set_xlabel("Phase calibrated bias")
        a2.set_yticks(range(len(variants)))
        if c == 0:
            a2.set_yticklabels([n for _, n in variants], fontsize=8)
            a2.invert_yaxis()
        else:
            plt.setp(a2.get_yticklabels(), visible=False)
    return fig, ["null_validation", "validate", "phase_observed", "decompose"], data


def figure_3(plt):
    sy, ica, cr = read("sym_synthetic"), read("ica"), read("sym_real")
    data = []
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.5), sharey=True)
    lv = [("plant_low", "Planted, low", 0.35), ("plant_mid", "Planted, middle", 0.65), ("plant_high", "Planted, high", 1.0)]
    for ax, null, panel, title in ((axes[0], "shuffle", "A", "A  Trial shuffle null"), (axes[1], "cross", "B", "B  Cross participant null")):
        for k, (cohort, name, alpha) in enumerate(lv):
            for j, (e, b) in enumerate(CELLS):
                bar(ax, j - 0.3 + 0.15 * k, pick(sy, cohort=cohort, null=null, variant="full", eeg=e, band=b), data, panel, name, fmt="o", color=GREY if alpha < 1 else DARK, alpha=max(alpha, 0.45))
        for j, (e, b) in enumerate(CELLS):
            if null == "shuffle":
                r1 = pick(ica, variant_eeg="v2", variant="full", k_used="16", eeg=e, band=b)
                r2 = pick(ica, variant_eeg="v2_emgadj", variant="full", k_used="16", eeg=e, band=b)
            else:
                r1, r2 = pick(cr, cohort="v2", k_used="all", eeg=e, band=b), pick(cr, cohort="v2_emgadj", k_used="all", eeg=e, band=b)
            bar(ax, j + 0.18, r1, data, panel, "observed, ICA", fmt="D", color=ACCENT, mfc="white", mec=ACCENT)
            bar(ax, j + 0.33, r2, data, panel, "observed, ICA, muscle proxy", fmt="s", color=ACCENT, alpha=0.5)
        limit_marks(ax, limits(null), data, panel, f"upper bias limit, {null}", half=0.42)
        MF.cell_ticks(ax)
        ax.tick_params(axis="x", labelsize=7.5)
        ax.set_title(title, loc="left", fontsize=9)
    axes[0].set_ylabel("Calibrated excess coherence")
    h = [plt.Line2D([], [], marker="o", color=GREY, alpha=0.45, lw=0, ms=4.5, label="Planted, low"),
         plt.Line2D([], [], marker="o", color=GREY, alpha=0.65, lw=0, ms=4.5, label="Planted, middle"),
         plt.Line2D([], [], marker="o", color=DARK, lw=0, ms=4.5, label="Planted, high"),
         plt.Line2D([], [], marker="D", mfc="white", mec=ACCENT, color=ACCENT, lw=0, ms=4.5, label="Observed, ICA"),
         plt.Line2D([], [], marker="s", color=ACCENT, alpha=0.5, lw=0, ms=4.5, label="Observed, ICA, muscle proxy"),
         plt.Line2D([], [], color=ORANGE, lw=1.0, ls="--", label="Upper bias limit")]
    fig.legend(handles=h, loc="lower center", ncol=3, fontsize=7.5, bbox_to_anchor=(0.5, -0.12))
    fig.tight_layout()
    return fig, ["sym_synthetic", "ica", "sym_real", "gate", "gate_cross"], data


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.parse_args(argv)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    MF.style(plt)
    import os
    from pathlib import Path
    odir = Path(os.environ["HB_OUT_DIR"]) if os.environ.get("HB_OUT_DIR") else OUT("results", "figures_three_verdicts")
    odir.mkdir(parents=True, exist_ok=True)
    rows = []
    for n, fn in ((1, figure_1), (2, figure_2), (3, figure_3)):
        fig, keys, data = fn(plt)
        for ext in ("png", "pdf"):
            fig.savefig(odir / f"figure_{n}.{ext}", dpi=300)
        plt.close(fig)
        with open(odir / f"figure_{n}_source_data.tsv", "w", encoding="utf-8", newline="") as fh:
            w = csv.DictWriter(fh, ["panel", "series", "eeg", "band", "value", "ci95_low", "ci95_high"], delimiter="\t")
            w.writeheader()
            w.writerows(data)
        for k in keys:
            rows.append({"figure": n, "source": "/".join(SRC[k]), "sha256_16": MF.sha(OUT(*SRC[k]))})
        print(f"figure {n}: " + ", ".join("/".join(SRC[k]) for k in keys))
    with open(odir / "figure_sources.tsv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, ["figure", "source", "sha256_16"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    append_execution_log([f"## {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} figures of the Three Verdicts manuscript (step 19)",
                          "- Command: python -m hb.make_figures_three_verdicts; outputs results/figures_three_verdicts (png, pdf, source data, figure_sources.tsv)"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
