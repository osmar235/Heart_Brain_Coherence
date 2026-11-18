"""Step 19: Table 3 and numerical claim checks of the Three Verdicts manuscript, from result tables and logs only.

Outputs in results/tables_three_verdicts/: table_3.tsv, table_3.md, claims_check.tsv
"""
from __future__ import annotations

import csv
import json
import re
import sys

from hb.download_ds008302 import PROJECT_ROOT
from hb.make_figures import CELLS, LABEL, f, pick



def OUT(*p):
    """Project path. HB_RESULTS_ROOT and HB_LOGS_ROOT redirect the 'results' and 'logs' trees (used by the release
    package, where the frozen result tables live in results/tables)."""
    import os
    from pathlib import Path
    for head, var in (("results", "HB_RESULTS_ROOT"), ("logs", "HB_LOGS_ROOT")):
        if p and p[0] == head and os.environ.get(var):
            return Path(os.environ[var]).joinpath(*p[1:])
    return PROJECT_ROOT.joinpath(*p)

MINUS = "−"


def read(*p):
    with open(OUT(*p), encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def num(x, d=4):
    return f"{x:+.{d}f}".replace("-", MINUS)


def ci(r, key):
    return f"{num(f(r[key]))} ({num(f(r['ci95_low']))} to {num(f(r['ci95_high']))})"


def main(argv=None) -> int:
    gt = read("results", "cohort_stats", "group_tests.tsv")
    ica = read("results", "ica_v2", "shuffle_calibrated_H1_picard.tsv")
    cr = read("results", "symmetric", "real_cross_calibrated_H1.tsv")
    sy = read("results", "symmetric", "synthetic_cells.tsv")
    lim_s = json.loads(OUT("results", "trial_shuffle", "gate.json").read_text(encoding="utf-8"))["bias_upper_edges"]
    lim_c = json.loads(OUT("results", "symmetric", "gate_cross.json").read_text(encoding="utf-8"))["cross_upper_limits"]
    rows, claims = [], []
    for e, b in CELLS:
        key = f"{e} {b}"
        g = pick(gt, hypothesis="H1_existence", variant="full", family="primary", eeg=e, band=b)
        s1 = pick(ica, variant_eeg="v2", variant="full", k_used="16", eeg=e, band=b)
        s2 = pick(ica, variant_eeg="v2_emgadj", variant="full", k_used="16", eeg=e, band=b)
        c1 = pick(cr, cohort="v2", k_used="all", eeg=e, band=b)
        c2 = pick(cr, cohort="v2_emgadj", k_used="all", eeg=e, band=b)
        c2k8 = pick(cr, cohort="v2_emgadj", k_used="8", eeg=e, band=b)
        smallest = "none"
        for lv in ("low", "mid", "high"):
            p = pick(sy, cohort=f"plant_{lv}", null="shuffle", variant="full", eeg=e, band=b)
            if p["detected_by_rule"] in ("1", "1.0"):
                smallest = f"{f(p['mean']):.4f} ({lv})"
                break
        met = {"shuffle ICA": s1["meets_claim_rule"], "shuffle ICA+proxy": s2["meets_claim_rule"],
               "cross ICA": c1["meets_claim_rule"], "cross ICA+proxy": c2["meets_claim_rule"]}
        any_met = any(v in ("1", "1.0") for v in met.values())
        if any_met:
            status = "Unresolved candidate: " + ", ".join(k for k, v in met.items() if v in ("1", "1.0")) + " met the rule"
        elif "low" in smallest:
            status = "No supported effect; absence informative at the low planted size"
        else:
            status = "Neither supported nor excluded (below the detected size)"
        rows.append({"Cell": f"{LABEL[e]} {b}", "Initial analysis": ci(g, "mean"),
                     "Trial shuffle, ICA": ci(s1, "calibrated_excess"), "Trial shuffle, ICA and muscle proxy": ci(s2, "calibrated_excess"),
                     "Limit, trial shuffle": num(lim_s[key]),
                     "Cross participant, ICA": ci(c1, "mean"), "Cross participant, ICA and muscle proxy": ci(c2, "mean"),
                     "Limit, cross participant": num(lim_c[key]), "Smallest planted effect detected": smallest, "Status": status})
        claims.append({"claim": f"{key}: cross ICA+proxy lower limit minus limit, K16", "value": f(c2["ci95_low"]) - lim_c[key]})
        claims.append({"claim": f"{key}: cross ICA+proxy lower limit minus limit, K8", "value": f(c2k8["ci95_low"]) - lim_c[key]})
        claims.append({"claim": f"{key}: trial shuffle ICA+proxy lower limit minus limit, K16", "value": f(s2["ci95_low"]) - lim_s[key]})
        pl = pick(sy, cohort="plant_low", null="cross", variant="full", eeg=e, band=b)
        claims.append({"claim": f"{key}: planted low, cross, lower limit minus limit", "value": f(pl["ci95_low"]) - lim_c[key]})
        for ve in ("v2", "v2_emgadj"):
            for v in ("full", "encret"):
                r = pick(ica, variant_eeg=ve, variant=v, k_used="16", eeg=e, band=b)
                claims.append({"claim": f"{key}: trial shuffle {ve} {v} calibrated excess", "value": f(r["calibrated_excess"])})
    ex = tot = 0
    logs = sorted(OUT("logs").glob("symmetric_real_*.log"))
    if logs:
        for m in re.finditer(r"same length (\d+)/(\d+)", logs[-1].read_text(encoding="utf-8")):
            ex, tot = ex + int(m.group(1)), tot + int(m.group(2))
        claims.append({"claim": f"cross participant, real data: share of blocks with a partner block of the same length ({logs[-1].name})", "value": ex / tot if tot else float("nan")})
    import os
    from pathlib import Path
    odir = Path(os.environ["HB_OUT_DIR"]) if os.environ.get("HB_OUT_DIR") else OUT("results", "tables_three_verdicts")
    odir.mkdir(parents=True, exist_ok=True)
    cols = list(rows[0])
    with open(odir / "table_3.tsv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    md = ["| " + " | ".join(cols) + " |", "| " + " | ".join("---" for _ in cols) + " |"] + ["| " + " | ".join(r[c] for c in cols) + " |" for r in rows]
    (odir / "table_3.md").write_bytes(("\n".join(md) + "\n").encode("utf-8"))  # LF on every system, so the file is byte identical across platforms
    with open(odir / "claims_check.tsv", "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, ["claim", "value"], delimiter="\t")
        w.writeheader()
        w.writerows([{"claim": c["claim"], "value": f"{c['value']:.6f}"} for c in claims])
    print(f"table 3 and {len(claims)} claim checks written to {odir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
