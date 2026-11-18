"""Verify a release folder built by hb.build_release.

Usage:  python -m hb.verify_release --root <release folder>      (Windows: scripts\\verify_release.cmd)

Checks
  1. required files and folders are present
  2. every file listed in SHA256SUMS.txt exists and matches; no unlisted files (except volatile and runtime ones)
  3. privacy: no absolute paths of the authors' computers, no private e-mail addresses, in any text file
  4. no raw recordings or large files (> 25 MB)
  5. launchers and code contain no fixed computer specific paths
  6. regeneration: Table 3, claim checks, figure source data and the manuscript are rebuilt from the stored
     result tables in a temporary folder and compared byte for byte with the stored files
Writes results/validation/verify_report_latest.md. Exit code 0 only if every check passes.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import os
import platform
import re
import subprocess
import sys
import tempfile
from pathlib import Path

REQUIRED = [
    "README.md", "pyproject.toml", "examples/coherence_example.py", "LICENSE", "LICENSE-CC-BY-4.0.md", ".gitattributes", "CITATION.cff", "requirements.txt", "requirements-lock.txt", ".gitignore", "SHA256SUMS.txt",
    "config/example_config.json", "config/download_config.json",
    "code/src/hb/__init__.py", "code/src/hb/vendor/emgsim_v3_3_0/functions.py", "code/src/hb/vendor/emgsim_v3_3_0/PROVENANCE.md",
    "code/tests",
    "scripts/setup_windows.cmd", "scripts/run_tests.cmd", "scripts/run_analysis.cmd", "scripts/run_figures.cmd",
    "scripts/verify_release.cmd", "scripts/pipeline/_env.cmd",
    "data_manifest/README.md", "data_manifest/dataset_manifest.tsv",
    "results/tables", "results/figure_source_data/figure_sources.tsv", "results/validation/claims_check.tsv",
    "figures/png/figure_1.png", "figures/png/figure_2.png", "figures/png/figure_3.png",
    "figures/pdf/figure_1.pdf", "figures/pdf/figure_2.pdf", "figures/pdf/figure_3.pdf",
    "manuscript/manuscript_three_verdicts.md", "manuscript/manuscript_three_verdicts.src.md", "manuscript/references.tsv",
    "docs/METHODS.md", "docs/REPRODUCIBILITY.md", "docs/RESULTS_INDEX.md",
]
VOLATILE = {"SHA256SUMS.txt", "docs/RELEASE_STATUS.md", "results/validation/verify_report_latest.md"}
RUNTIME_DIRS = {".git", ".venv", "regenerated", "derivatives", "logs", "data", "tmp", "__pycache__", ".pytest_cache", "build"}
TEXT_EXT = {".md", ".txt", ".tsv", ".json", ".py", ".cmd", ".cff", ".log", ".csv", ".yml", ".yaml", ""}
BINARY_BAN = {".set", ".fdt", ".edf", ".bdf", ".fif", ".npz", ".zip", ".mat"}
PRIVATE = [
    ("Windows user folder", re.compile(r"[A-Za-z]:[\\/]+Users[\\/]+[A-Za-z0-9]", re.I)),
    ("Linux workspace folder", re.compile(r"/sessions/[A-Za-z0-9]")),
    ("home folder", re.compile(r"/home/[a-z0-9]")),
    ("private e-mail", re.compile("arena" + "235|@" + "gmail" + r"\.com", re.I)),
]


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def files(root: Path):
    for p in sorted(root.rglob("*")):
        rel = p.relative_to(root)
        if p.is_dir() or (rel.parts and rel.parts[0] in RUNTIME_DIRS) or "__pycache__" in rel.parts or any(part.endswith(".egg-info") for part in rel.parts):
            continue
        if rel.parts[0] == "results" and len(rel.parts) > 1 and rel.parts[1] not in {"tables", "figure_source_data", "validation"}:
            continue  # outputs of new analysis runs in the release folder
        yield p, rel.as_posix()


def check_required(root):
    missing = [r for r in REQUIRED if not (root / r).exists()]
    return not missing, "all present" if not missing else "missing: " + ", ".join(missing)


def check_sums(root):
    f = root / "SHA256SUMS.txt"
    if not f.exists():
        return False, "SHA256SUMS.txt missing"
    listed, bad = {}, []
    for line in f.read_text(encoding="utf-8").splitlines():
        if line.strip():
            h, rel = line.split("  ", 1)
            listed[rel] = h
    for rel, h in listed.items():
        p = root / rel
        if not p.exists():
            bad.append(f"missing {rel}")
        elif sha256(p) != h:
            bad.append(f"changed {rel}")
    extra = [rel for _, rel in files(root) if rel not in listed and rel not in VOLATILE]
    detail = f"{len(listed)} files listed"
    if bad:
        detail += "; " + "; ".join(bad[:10]) + (f" (and {len(bad) - 10} more)" if len(bad) > 10 else "")
    if extra:
        detail += "; not listed: " + ", ".join(extra[:10])
    return not bad and not extra, detail


def check_privacy(root):
    hits = []
    for p, rel in files(root):
        if p.suffix.lower() not in TEXT_EXT or rel == "SHA256SUMS.txt":
            continue
        try:
            text = p.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        for name, rx in PRIVATE:
            m = rx.search(text)
            if m:
                line = text.count("\n", 0, m.start()) + 1
                hits.append(f"{rel}:{line} ({name})")
    return not hits, "none found" if not hits else "; ".join(hits[:15])


def check_size(root):
    bad = []
    for p, rel in files(root):
        if p.suffix.lower() in BINARY_BAN:
            bad.append(f"{rel} (banned type)")
        elif p.stat().st_size > 25 * 1024 * 1024:
            bad.append(f"{rel} ({p.stat().st_size / 1e6:.0f} MB)")
    total = sum(p.stat().st_size for p, _ in files(root))
    return not bad, f"total {total / 1e6:.1f} MB" + ("; " + "; ".join(bad) if bad else "")


def check_portable(root):
    bad = []
    rx = re.compile(r"(?<![%\w])[A-Za-z]:[\\/](?!\\)")
    for p, rel in files(root):
        if not (rel.startswith("scripts/") or rel.startswith("code/") or rel.startswith("config/")):
            continue
        if rel == "config/example_config.json":  # an example path is the purpose of this file
            continue
        if p.suffix.lower() not in {".cmd", ".py", ".json"}:
            continue
        for k, line in enumerate(p.read_text(encoding="utf-8").splitlines(), 1):
            s = line.strip()
            if s.lower().startswith(("rem ", "#", "echo ")) or '"_comment"' in s or "example" in s.lower():
                continue
            if rx.search(line):
                bad.append(f"{rel}:{k}")
    return not bad, "no fixed drive paths in code, launchers or configuration" if not bad else "; ".join(bad[:15])


def run(cmd, env, cwd):
    r = subprocess.run(cmd, env=env, cwd=cwd, capture_output=True, text=True)
    return r.returncode, (r.stdout + r.stderr)[-2000:]


def check_regeneration(root):
    tmp = Path(tempfile.mkdtemp(prefix="hb_verify_"))
    env = dict(os.environ, PYTHONPATH=str(root / "code" / "src"), HB_PROJECT_ROOT=str(root), PYTHONUTF8="1",
               HB_RESULTS_ROOT=str(root / "results" / "tables"), HB_LOGS_ROOT=str(root / "results" / "validation" / "logs"),
               HB_NO_EXECUTION_LOG="1", MPLBACKEND="Agg")
    py = sys.executable
    pairs, notes = [], []
    rc, out = run([py, "-W", "ignore", "-m", "hb.make_tables_three_verdicts"], dict(env, HB_OUT_DIR=str(tmp / "tables")), root)
    if rc:
        return False, "table script failed: " + out
    for name in ("table_3.tsv", "table_3.md", "claims_check.tsv"):
        pairs.append((tmp / "tables" / name, root / "results" / "tables" / "tables_three_verdicts" / name))
    pairs.append((tmp / "tables" / "claims_check.tsv", root / "results" / "validation" / "claims_check.tsv"))
    rc, out = run([py, "-W", "ignore", "-m", "hb.make_figures_three_verdicts"], dict(env, HB_OUT_DIR=str(tmp / "figures")), root)
    if rc:
        return False, "figure script failed: " + out
    for n in (1, 2, 3):
        pairs.append((tmp / "figures" / f"figure_{n}_source_data.tsv", root / "results" / "figure_source_data" / f"figure_{n}_source_data.tsv"))
    pairs.append((tmp / "figures" / "figure_sources.tsv", root / "results" / "figure_source_data" / "figure_sources.tsv"))
    rc, out = run([py, "-W", "ignore", "-m", "hb.build_manuscript", "--src", str(root / "manuscript" / "manuscript_three_verdicts.src.md"),
                   "--refs", str(root / "manuscript" / "references.tsv"), "--table3", str(tmp / "tables" / "table_3.md"),
                   "--out", str(tmp / "manuscript.md")], env, root)
    if rc:
        return False, "manuscript build failed: " + out
    pairs.append((tmp / "manuscript.md", root / "manuscript" / "manuscript_three_verdicts.md"))
    diff = [f"{b.relative_to(root).as_posix()}" for a, b in pairs if not a.exists() or not b.exists() or a.read_bytes() != b.read_bytes()]
    notes.append(f"{len(pairs) - len(diff)} of {len(pairs)} regenerated files identical")
    notes.append("figure images not compared (depend on Matplotlib version and fonts)")
    if diff:
        notes.append("different: " + ", ".join(diff))
    return not diff, "; ".join(notes)


def verify(root: Path, log=print, regenerate=True) -> dict:
    root = Path(root).resolve()
    checks = []
    for name, fn in (("Required files", check_required), ("Checksums", check_sums), ("Privacy scan", check_privacy),
                     ("No raw data or large files", check_size), ("Portable paths", check_portable)):
        ok, detail = fn(root)
        checks.append({"name": name, "ok": ok, "detail": detail})
        log(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}")
    if regenerate:
        ok, detail = check_regeneration(root)
        checks.append({"name": "Regeneration from stored tables", "ok": ok, "detail": detail})
        log(f"{'PASS' if ok else 'FAIL'}  Regeneration from stored tables: {detail}")
    passed = all(c["ok"] for c in checks)
    rep = root / "results" / "validation" / "verify_report_latest.md"
    rep.parent.mkdir(parents=True, exist_ok=True)
    rows = "\n".join(f"| {c['name']} | {'pass' if c['ok'] else 'FAIL'} | {c['detail']} |" for c in checks)
    rep.write_bytes((f"# Verification report\n\n{dt.datetime.now():%Y-%m-%d %H:%M}, Python {platform.python_version()} on "
                     f"{platform.system()} {platform.release()}.\n\n**{'PASSED' if passed else 'FAILED'}**\n\n"
                     f"| Check | Result | Detail |\n| --- | --- | --- |\n{rows}\n").encode("utf-8"))
    return {"passed": passed, "checks": checks}


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=os.environ.get("HB_PROJECT_ROOT", "."))
    ap.add_argument("--no-regenerate", action="store_true", help="skip the regeneration check")
    a = ap.parse_args(argv)
    r = verify(Path(a.root), regenerate=not a.no_regenerate)
    print("VERIFICATION " + ("PASSED" if r["passed"] else "FAILED"))
    return 0 if r["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
