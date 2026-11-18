"""Step 17: manuscript figures, drawn only from result tables (no recording is read).

Figure 1  synthetic bias without coupling, per cell: phase calibration (step 14 model, matched model) against
          trial shuffle calibration (both models), with the observed phase calibrated values for comparison
Figure 2  decomposition of the phase calibrated synthetic bias, LF cells, model components switched off
Figure 3  real data per cell across analysis versions, with the trial shuffle bias edge
Outputs: results/figures/figure_<n>.png (300 dpi) and .pdf, results/figures/figure_sources.tsv
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

from hb.download_ds008302 import PROJECT_ROOT, append_execution_log

OUT = lambda *p: PROJECT_ROOT.joinpath(*p)
CELLS = [("theta", "LF"), ("theta", "HF"), ("alpha", "LF"), ("alpha", "HF"), ("beta", "LF"), ("beta", "HF")]
LABEL = {"theta": "Frontal theta", "alpha": "Posterior alpha", "beta": "Centroparietal beta"}
ACCENT, GREY, DARK = "#2b6cb0", "#8a8f98", "#222222"
SRC = {
    "group_tests": ("results", "cohort_stats", "group_tests.tsv"),
    "phase_observed": ("results", "null_calibration", "check2_H1_calibrated.tsv"),
    "null_validation": ("results", "null_validation", "null_validation_cells.tsv"),
    "validate": ("results", "trial_shuffle", "validate_cells.tsv"),
    "decompose": ("results", "trial_shuffle", "decompose_cells.tsv"),
    "real_shuffle": ("results", "trial_shuffle", "real_shuffle_calibrated_H1.tsv"),
    "ica": ("results", "ica_v2", "shuffle_calibrated_H1_picard.tsv"),
    "gate": ("results", "trial_shuffle", "gate.json"),
}


def read(key):
    with open(OUT(*SRC[key]), encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def f(x):
    return float(x)


def pick(rows, **kw):
    out = [r for r in rows if all(r.get(k) == v for k, v in kw.items())]
    if len(out) != 1:
        raise ValueError(f"expected one row for {kw}, found {len(out)}")
    return out[0]


def style(plt):
    plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans", "axes.spines.top": False, "axes.spines.right": False,
                         "axes.edgecolor": "#555555", "axes.linewidth": 0.8, "xtick.color": "#333333", "ytick.color": "#333333",
                         "legend.frameon": False, "savefig.bbox": "tight"})


def cell_ticks(ax):
    ax.set_xticks(range(len(CELLS)))
    ax.set_xticklabels([f"{LABEL[e].split()[-1]}\n{b}" for e, b in CELLS])
    ax.axhline(0, color="#555555", lw=0.8)


def figure_1(plt):
    nv, va, ob = read("null_validation"), read("validate"), read("phase_observed")
    series = []
    m, lo, hi = [], [], []
    for e, b in CELLS:
        r = pick(nv, scenario="independent", variant="full", statistic="calibrated", eeg=e, band=b)
        m.append(f(r["mean"])); lo.append(f(r["mean"]) - f(r["ci95_halfwidth"])); hi.append(f(r["mean"]) + f(r["ci95_halfwidth"]))
    series.append(("Phase calibration, step 14 model", m, lo, hi, ACCENT, "o", 1.0))
    for cohort, null, name, col, mk, alpha in (("matched", "phase", "Phase calibration, matched model", ACCENT, "s", 0.55),
                                                ("independent", "shuffle", "Trial shuffle, step 14 model", GREY, "o", 1.0),
                                                ("matched", "shuffle", "Trial shuffle, matched model", GREY, "s", 0.6)):
        rows = [pick(va, cohort=cohort, null=null, variant="full", eeg=e, band=b) for e, b in CELLS]
        series.append((name, [f(r["bias"]) for r in rows], [f(r["ci95_low"]) for r in rows], [f(r["ci95_high"]) for r in rows], col, mk, alpha))
    obs = [pick(ob, variant="full", eeg=e, band=b) for e, b in CELLS]
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    off = np.linspace(-0.27, 0.27, len(series) + 1)
    for k, (name, mm, ll, hh, col, mk, alpha) in enumerate(series):
        x = np.arange(len(CELLS)) + off[k]
        ax.errorbar(x, mm, yerr=[np.array(mm) - ll, np.array(hh) - mm], fmt=mk, color=col, alpha=alpha, ms=4.5, lw=1.1, capsize=0, label=name)
    x = np.arange(len(CELLS)) + off[-1]
    ax.errorbar(x, [f(r["calibrated_excess"]) for r in obs],
                yerr=[[f(r["calibrated_excess"]) - f(r["ci95_low"]) for r in obs], [f(r["ci95_high"]) - f(r["calibrated_excess"]) for r in obs]],
                fmt="D", mfc="white", mec=DARK, color=DARK, ms=4.5, lw=1.1, capsize=0, label="Observed, phase calibrated")
    cell_ticks(ax)
    ax.set_ylabel("Calibrated excess coherence")
    ax.set_title("Without coupling, phase calibration leaves a timing bias in LF; trial shuffling removes it", loc="left", fontsize=9.5)
    ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0), fontsize=8)
    return fig, ["null_validation", "validate", "phase_observed"]


def figure_2(plt):
    dc = read("decompose")
    variants = [("base", "Full model"), ("no_gain", "No gain variability"), ("no_latency", "No latency jitter"),
                ("stationary", "Stationary noise"), ("no_recall_scaling", "No recall scaling"), ("mean_only", "Identical trials, stationary noise")]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.9), sharey=True, sharex=True)
    for ax, e in zip(axes, ("theta", "alpha", "beta")):
        for i, (v, name) in enumerate(variants):
            r = pick(dc, cohort=v, null="phase", variant="full", eeg=e, band="LF")
            m = f(r["bias"])
            col = ACCENT if v in ("base", "mean_only") else GREY
            ax.errorbar(m, i, xerr=[[m - f(r["ci95_low"])], [f(r["ci95_high"]) - m]], fmt="o", color=col, ms=4, lw=1.1, capsize=0)
        ax.axvline(0, color="#555555", lw=0.8)
        ax.set_title(f"{LABEL[e]} LF", fontsize=9)
        ax.set_yticks(range(len(variants)))
        ax.set_yticklabels([n for _, n in variants])
        ax.invert_yaxis()
        ax.set_xlabel("Calibrated bias")
    fig.suptitle("The bias remains with identical trial responses: it follows from self paced trial timing", x=0.0, ha="left", fontsize=9.5)
    fig.tight_layout()
    return fig, ["decompose"]


def figure_3(plt):
    gt, ob, rs, ica = read("group_tests"), read("phase_observed"), read("real_shuffle"), read("ica")
    edges = json.loads(OUT(*SRC["gate"]).read_text(encoding="utf-8"))["bias_upper_edges"]
    versions = [("Standard, uncalibrated", GREY, 0.5), ("Phase calibrated", GREY, 0.8), ("Trial shuffle", DARK, 1.0),
                ("Trial shuffle, ICA", ACCENT, 1.0), ("Trial shuffle, ICA, muscle proxy", ACCENT, 0.5)]
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.6), sharey=True)
    for j, (e, b) in enumerate(CELLS):
        ax = axes[j % 2, j // 2]
        rows = [pick(gt, hypothesis="H1_existence", variant="full", family="primary", eeg=e, band=b),
                pick(ob, variant="full", eeg=e, band=b),
                pick(rs, variant="full", k_used="16", eeg=e, band=b),
                pick(ica, variant_eeg="v2", variant="full", k_used="16", eeg=e, band=b),
                pick(ica, variant_eeg="v2_emgadj", variant="full", k_used="16", eeg=e, band=b)]
        for i, ((name, col, alpha), r) in enumerate(zip(versions, rows)):
            m = f(r["mean"]) if "mean" in r else f(r["calibrated_excess"])
            ax.errorbar(i, m, yerr=[[m - f(r["ci95_low"])], [f(r["ci95_high"]) - m]], fmt="o", color=col, alpha=alpha, ms=4.5, lw=1.2, capsize=0)
        edge = edges[f"{e} {b}"]
        ax.plot([1.6, 4.4], [edge, edge], color="#c05621", lw=1.0, ls="--")
        ax.axhline(0, color="#555555", lw=0.8)
        ax.set_title(f"{LABEL[e]} {b}", fontsize=9)
        ax.set_xticks(range(len(versions)))
        ax.set_xticklabels(["Std", "Phase", "Shuffle", "ICA", "ICA+\nmuscle"], fontsize=7.5)
    axes[0, 0].set_ylabel("Excess coherence")
    axes[1, 0].set_ylabel("Excess coherence")
    fig.suptitle("Effects shrink as timing bias and artifacts are controlled; dashed line: trial shuffle bias edge", x=0.0, ha="left", fontsize=9.5)
    fig.tight_layout()
    return fig, ["group_tests", "phase_observed", "real_shuffle", "ica", "gate"]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.parse_args(argv)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    style(plt)
    odir = OUT("results", "figures")
    odir.mkdir(parents=True, exist_ok=True)
    rows = []
    for n, fn in ((1, figure_1), (2, figure_2), (3, figure_3)):
        fig, keys = fn(plt)
        for ext in ("png", "pdf"):
            fig.savefig(odir / f"figure_{n}.{ext}", dpi=300)
        plt.close(fig)
        for k in keys:
            rows.append({"figure": n, "source": "/".join(SRC[k]), "sha256_16": sha(OUT(*SRC[k]))})
        print(f"figure {n}: " + ", ".join("/".join(SRC[k]) for k in keys))
    with open(odir / "figure_sources.tsv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, ["figure", "source", "sha256_16"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    append_execution_log([f"## {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} figures (step 17)",
                          "- Command: python -m hb.make_figures; outputs results/figures/figure_1 to figure_3 (.png, .pdf) and figure_sources.tsv (source tables with SHA256 prefixes)"])
    return 0


if __name__ == "__main__":
    sys.exit(main())
