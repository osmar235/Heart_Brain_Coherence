"""Assemble the Three Verdicts manuscript: numbered citations in order of first use, reference list, Table 3.

Input : docs/manuscript/manuscript_three_verdicts.src.md, docs/manuscript/references.tsv,
        results/tables_three_verdicts/table_3.md
Output: docs/manuscript/manuscript_three_verdicts.md
Usage : python -m hb.build_manuscript [--tests N]   (N = number of automated tests, written into the Methods)
"""
from __future__ import annotations

import argparse
import csv
import re
import sys

from hb.download_ds008302 import PROJECT_ROOT

OUT = lambda *p: PROJECT_ROOT.joinpath(*p)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tests", default="85")
    ap.add_argument("--src", help="source markdown (default docs/manuscript/manuscript_three_verdicts.src.md)")
    ap.add_argument("--refs", help="references table (default docs/manuscript/references.tsv)")
    ap.add_argument("--table3", help="Table 3 markdown (default results/tables_three_verdicts/table_3.md)")
    ap.add_argument("--out", help="output markdown (default docs/manuscript/manuscript_three_verdicts.md)")
    a = ap.parse_args(argv)
    from pathlib import Path
    f_src = Path(a.src) if a.src else OUT("docs", "manuscript", "manuscript_three_verdicts.src.md")
    f_refs = Path(a.refs) if a.refs else OUT("docs", "manuscript", "references.tsv")
    f_t3 = Path(a.table3) if a.table3 else OUT("results", "tables_three_verdicts", "table_3.md")
    f_out = Path(a.out) if a.out else OUT("docs", "manuscript", "manuscript_three_verdicts.md")
    src = f_src.read_text(encoding="utf-8")
    with open(f_refs, encoding="utf-8", newline="") as fh:
        refs = {r["key"]: r["reference"] for r in csv.DictReader(fh, delimiter="\t")}
    src = src.replace("{{TABLE3}}", f_t3.read_text(encoding="utf-8").strip())
    src = src.replace("{{NTESTS}}", str(a.tests))
    order = []

    def cite(m):
        nums = []
        for k in re.findall(r"@([A-Za-z0-9]+)", m.group(0)):
            if k not in refs:
                raise KeyError(f"citation key without reference: {k}")
            if k not in order:
                order.append(k)
            nums.append(order.index(k) + 1)
        return "\\[" + ", ".join(str(n) for n in nums) + "\\]"

    body = re.sub(r"\[@[^\]]+\]", cite, src)
    unused = [k for k in refs if k not in order]
    lines = [f"{i}. " + refs[k].replace("<", "\\<").replace(">", "\\>") for i, k in enumerate(order, 1)]
    body = body.replace("{{REFERENCES}}", "\n".join(lines))
    if "{{" in body:
        raise ValueError("unresolved placeholder")
    f_out.parent.mkdir(parents=True, exist_ok=True)
    f_out.write_bytes(body.encode("utf-8"))
    print(f"manuscript built: {len(order)} references cited" + (f"; not cited and left out: {unused}" if unused else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
