"""Build the portable release folder submission_ready/latest/ from the active project.

Usage:  python -m hb.build_release [--no-freeze] [--skip-verify] [--snapshot NAME]
Windows launcher: scripts\\19_build_release.cmd

What it does
  - deletes and recreates ONLY the contents of <project>/submission_ready/latest (a local .git there is kept);
    nothing outside that folder is changed, moved or deleted
  - copies (never moves) code, tests, result tables, figures, manuscript, documents and manifests
  - writes the release launchers from release_template/ and adapts the numbered pipeline launchers
  - replaces absolute paths of the authors' computers by placeholders in copied text files
  - writes requirements-lock.txt (pip freeze of the running interpreter; with --no-freeze, the versions recorded
    in the logs, marked as incomplete)
  - writes SHA256SUMS.txt, runs hb.verify_release on the result, and writes docs/RELEASE_STATUS.md
  - with --snapshot NAME, copies the verified folder to submission_ready/NAME, refusing to overwrite
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import os
import platform
import re
import shutil
import subprocess
import sys
from pathlib import Path

from hb.download_ds008302 import PROJECT_ROOT, append_execution_log

TEXT_EXT = {".md", ".txt", ".tsv", ".json", ".py", ".cmd", ".cff", ".log", ".csv", ".yml", ".yaml"}
VOLATILE = {"SHA256SUMS.txt", "docs/RELEASE_STATUS.md", "results/validation/verify_report_latest.md"}
SKIP_DIRS = {"__pycache__", ".pytest_cache"}
LOGGED_VERSIONS = [
    ("python", "3.14.3"), ("numpy", "2.4.3"), ("scipy", "1.17.1"), ("mne", "1.13.2"),
    ("mne-icalabel", "0.10.0"), ("python-picard", "0.8.2"), ("onnxruntime", "1.31.1"),
]
REDACTIONS = [
    (re.compile(r"[A-Za-z]:[\\/]+Users[\\/]+[^\\/\s\"']+", re.I), "<USER_HOME>"),
    (re.compile(r"/sessions/[^/\s\"']+"), "<WORKSPACE_HOME>"),
    (re.compile(r"/home/[^/\s\"']+"), "<WORKSPACE_HOME>"),
    (re.compile(r"/mnt/user-data/[^\s\"']*"), "<WORKSPACE_FILE>"),
]
CONFIRM = (
    'if /i "%HB_CONFIRM_DOWNLOAD%"=="YES" goto :confirmed\r\n'
    "echo This step downloads large recordings from OpenNeuro (see the comment at the top of this file).\r\n"
    'set /p "ANS=Type YES to start the download: "\r\n'
    'if /i not "%ANS%"=="YES" (echo Download not started. & popd & exit /b 4)\r\n'
    ":confirmed\r\n"
)
PIPELINE_ENV = (
    "@echo off\r\n"
    "rem Shared setup of the pipeline launchers inside the release folder. Not meant to be run alone.\r\n"
    'set "HB_ROOT=%~dp0..\\.."\r\n'
    'for %%I in ("%HB_ROOT%") do set "HB_ROOT=%%~fI"\r\n'
    'set "PY=%HB_ROOT%\\.venv\\Scripts\\python.exe"\r\n'
    'if not exist "%PY%" if defined HB_PYTHON set "PY=%HB_PYTHON%"\r\n'
    'if not exist "%PY%" set "PY=python"\r\n'
    'set "PYTHONPATH=%HB_ROOT%\\code\\src"\r\n'
    'set "PYTHONUTF8=1"\r\n'
    'set "HB_PROJECT_ROOT=%HB_ROOT%"\r\n'
    'if not exist "%HB_ROOT%\\logs" mkdir "%HB_ROOT%\\logs"\r\n'
    "exit /b 0\r\n"
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def redact(text: str) -> str:
    for rx, rep in REDACTIONS:
        text = rx.sub(rep, text)
    return text


class Builder:
    def __init__(self, src: Path, dest: Path, log):
        self.src, self.dest, self.log = src, dest, log
        self.copied = 0

    def put(self, src_file: Path, rel: str, text_transform=None, crlf=None):
        out = self.dest / rel
        out.parent.mkdir(parents=True, exist_ok=True)
        if src_file.suffix.lower() in TEXT_EXT or text_transform or crlf is not None:
            data = src_file.read_bytes().decode("utf-8")
            if src_file.suffix.lower() != ".py":  # code is copied unchanged; the privacy scan checks it
                data = redact(data)
            if text_transform:
                data = text_transform(data)
            if crlf is True:
                data = data.replace("\r\n", "\n").replace("\n", "\r\n")
            out.write_bytes(data.encode("utf-8"))
        else:
            shutil.copyfile(src_file, out)
        self.copied += 1

    def write(self, rel: str, text: str, crlf=False):
        out = self.dest / rel
        out.parent.mkdir(parents=True, exist_ok=True)
        if crlf:
            text = text.replace("\r\n", "\n").replace("\n", "\r\n")
        out.write_bytes(text.encode("utf-8"))
        self.copied += 1

    def tree(self, src_dir: Path, rel_dir: str, keep=lambda p: True):
        for p in sorted(src_dir.rglob("*")):
            if p.is_dir() or any(part in SKIP_DIRS for part in p.parts) or not keep(p):
                continue
            self.put(p, f"{rel_dir}/{p.relative_to(src_dir).as_posix()}")


def clean_dest(dest: Path, project: Path):
    dest = dest.resolve()
    expected = (project / "submission_ready" / "latest").resolve()
    if dest != expected:
        raise RuntimeError(f"refusing to clean {dest}: only {expected} may be rebuilt")
    if dest.exists():
        for child in dest.iterdir():
            if child.name == ".git":
                continue
            shutil.rmtree(child) if child.is_dir() else child.unlink()
    dest.mkdir(parents=True, exist_ok=True)


def pipeline_cmd(text: str, name: str) -> str:
    text = text.replace('pushd "%~dp0.."', 'pushd "%~dp0..\\.."')
    text = text.replace("-s tests", "-s code\\tests")
    text = text.replace('if exist "..\\.venv\\Scripts\\python.exe" goto :havevenv\r\n', "")
    text = text.replace('if exist "..\\.venv\\Scripts\\python.exe" goto :havevenv\n', "")
    if name in ("02_download_pilot.cmd", "02b_download_cohort.cmd"):
        lines = text.replace("\r\n", "\n").split("\n")
        i = next(k for k, ln in enumerate(lines) if ln.startswith("pushd "))
        lines.insert(i + 1, CONFIRM.replace("\r\n", "\n").rstrip("\n"))
        text = "\n".join(lines)
    return text


def _norm(name: str) -> str:
    return re.sub(r"[-_.]+", "-", name).lower()


def project_packages(freeze_text: str):
    """Restrict a pip freeze to requirements.txt and the installed dependencies of those packages."""
    import importlib.metadata as md
    try:
        from packaging.requirements import Requirement
    except ImportError:
        try:
            from pip._vendor.packaging.requirements import Requirement
        except ImportError:
            Requirement = None
    roots = []
    for line in (PROJECT_ROOT / "release_template" / "requirements.txt").read_text(encoding="utf-8").splitlines():
        line = line.split("#")[0].strip()
        m = re.match(r"[A-Za-z0-9._-]+", line)
        if m:
            roots.append(_norm(m.group(0)))
    seen, stack = set(), list(roots)
    while stack:
        n = stack.pop()
        if n in seen:
            continue
        try:
            dist = md.distribution(n)
        except md.PackageNotFoundError:
            continue
        seen.add(n)
        for r in dist.requires or []:
            if Requirement is not None:
                try:
                    req = Requirement(r)
                    if req.marker is not None and not req.marker.evaluate({"extra": ""}):
                        continue
                    stack.append(_norm(req.name))
                    continue
                except Exception:
                    pass
            if "extra ==" in r:
                continue
            m = re.match(r"[A-Za-z0-9._-]+", r)
            if m:
                stack.append(_norm(m.group(0)))
    keep, dropped = [], 0
    for line in freeze_text.splitlines():
        m = re.match(r"([A-Za-z0-9._-]+)==", line.strip())
        if not m:
            continue
        if _norm(m.group(1)) in seen:
            keep.append(line.strip())
        else:
            dropped += 1
    return keep, dropped, sorted(set(roots) - seen)


DOCS_NOTE = ("> Public repository note: references below to analysis plans (`*_PLAN*.md`, `ANALYSIS_FREEZE*.md`), "
             "stage reports (`*_RESULT*.md`), decision records (`DECISOES*.md`), the execution log and run logs (`logs/`) "
             "point to the author's project records, which are not part of this repository. They are available from "
             "the author on request.\n\n")


def lock_file(freeze: bool) -> str:
    head = f"# Exact package versions. Written by hb.build_release on {dt.datetime.now():%Y-%m-%d %H:%M}.\n"
    if freeze:
        r = subprocess.run([sys.executable, "-m", "pip", "freeze", "--exclude-editable"], capture_output=True, text=True)
        if r.returncode == 0 and r.stdout.strip():
            keep, dropped, missing = project_packages(r.stdout)
            note = (f"# Missing from the build environment: {', '.join(missing)}\n" if missing else "")
            return (head + f"# Interpreter: Python {platform.python_version()} on {platform.platform()}\n"
                    "# pip freeze of the environment that ran the build, restricted to the packages in\n"
                    "# requirements.txt and their dependencies "
                    f"({len(keep)} kept, {dropped} unrelated packages of that environment left out).\n"
                    "# Install with the same Python minor version; scripts\\setup_windows.cmd checks this.\n"
                    + note + "".join(k + "\n" for k in keep))
    body = "".join(f"{n}=={v}\n" for n, v in LOGGED_VERSIONS if n != "python")
    body += "matplotlib>=3.8\ntqdm>=4.60\nh5py>=3.10\nscikit-learn>=1.3\n"
    return (head + "# INCOMPLETE: built outside the analysis environment (or pip freeze failed).\n"
            "# Versions below are those recorded in the analysis logs (Windows 11, Python "
            + dict(LOGGED_VERSIONS)["python"] + ").\n"
            "# Matplotlib, tqdm, h5py and scikit-learn versions were not logged. Rebuild on the analysis computer\n"
            "# with scripts\\19_build_release.cmd to replace this file by a full pip freeze.\n"
            "# Lines with == are pinned; the four lines with >= are lower bounds only.\n" + body)


def build(project: Path, freeze: bool, log) -> Path:
    dest = project / "submission_ready" / "latest"
    clean_dest(dest, project)
    b = Builder(project, dest, log)
    tmpl = project / "release_template"

    # templates (README, LICENSE, CITATION, .gitignore, requirements, config example, launchers, data README)
    for p in sorted(tmpl.rglob("*")):
        if p.is_file():
            rel = p.relative_to(tmpl).as_posix()
            b.put(p, rel, crlf=True if p.suffix == ".cmd" else None)

    # code and tests
    b.tree(project / "src" / "hb", "code/src/hb", keep=lambda p: p.suffix in {".py", ".md"})
    b.tree(project / "tests", "code/tests", keep=lambda p: p.suffix == ".py")

    # pipeline launchers
    for p in sorted((project / "scripts").glob("*.cmd")):
        if p.name in {"_env.cmd", "19_build_release.cmd", "20_publish_github.cmd"}:  # project tools, not analysis steps
            continue
        b.put(p, f"scripts/pipeline/{p.name}", text_transform=lambda t, n=p.name: pipeline_cmd(t, n), crlf=True)
    b.write("scripts/pipeline/_env.cmd", PIPELINE_ENV, crlf=True)

    # configuration
    cfg = (project / "config" / "download_config.json").read_text(encoding="utf-8")
    cfg = cfg.replace("data_manifest/ds008302_v1.0.0_manifest.tsv", "data_manifest/dataset_manifest.tsv")
    b.write("config/download_config.json", redact(cfg))

    # data manifest
    dm = project / "data_manifest"
    b.put(dm / "ds008302_v1.0.0_manifest.tsv", "data_manifest/dataset_manifest.tsv")
    hits = sorted(dm.glob("download_report_cohort_*.tsv"))  # the files of the cohort used, with their checksums
    if hits:
        b.put(hits[-1], "data_manifest/download_report_cohort.tsv")

    # result tables (every small text result of every stage, in the project folder structure)
    res = project / "results"
    b.tree(res, "results/tables",
           keep=lambda p: p.suffix in {".tsv", ".json", ".txt", ".md"}
           and not p.relative_to(res).as_posix().startswith(("figures/", "figures_three_verdicts/")))
    fig = res / "figures_three_verdicts"
    for p in sorted(fig.glob("figure_*_source_data.tsv")) + [fig / "figure_sources.tsv"]:
        b.put(p, f"results/figure_source_data/{p.name}")
    for ext in ("png", "pdf"):
        for p in sorted(fig.glob(f"figure_*.{ext}")):
            b.put(p, f"figures/{ext}/{p.name}")
    tv = res / "tables_three_verdicts"
    b.put(tv / "claims_check.tsv", "results/validation/claims_check.tsv")
    # the one log read by hb.make_tables_three_verdicts (share of partner blocks of the same length)
    for p in sorted((project / "logs").glob("symmetric_real_*.log"))[-1:]:
        b.put(p, f"results/validation/logs/{p.name}")

    # manuscript (current version only; earlier drafts stay in the project)
    md = project / "docs" / "manuscript"
    for name in ("manuscript_three_verdicts.md", "manuscript_three_verdicts.src.md", "references.tsv"):
        b.put(md / name, f"manuscript/{name}")

    # documents needed to understand and reproduce the analysis; decision records, execution logs,
    # analysis plans and stage reports stay in the project and are available from the author on request
    docs = project / "docs"
    for name in ("METHODS.md", "REPRODUCIBILITY.md", "RESULTS_INDEX.md"):
        b.put(docs / name, f"docs/{name}", text_transform=lambda t: DOCS_NOTE + t)

    b.write("requirements-lock.txt", lock_file(freeze))
    log(f"copied or written: {b.copied} files")
    return dest


def write_sums(dest: Path):
    lines = []
    for p in sorted(dest.rglob("*")):
        if p.is_dir() or ".git" in p.relative_to(dest).parts:
            continue
        rel = p.relative_to(dest).as_posix()
        if rel in VOLATILE:
            continue
        lines.append(f"{sha256(p)}  {rel}")
    (dest / "SHA256SUMS.txt").write_bytes(("\n".join(lines) + "\n").encode("utf-8"))
    return len(lines)


def status_doc(dest: Path, report: dict, freeze_used: bool) -> str:
    ok = report["passed"]
    rows = "\n".join(f"| {c['name']} | {'pass' if c['ok'] else 'FAIL'} | {c['detail']} |" for c in report["checks"])
    return f"""# RELEASE STATUS

Built {dt.datetime.now():%Y-%m-%d %H:%M} by `hb.build_release` with Python {platform.python_version()} on {platform.system()} {platform.release()}.

**Automated verification: {'PASSED' if ok else 'FAILED'}.** {'All checks below passed in the environment that built the release.' if ok else 'At least one check failed; see below and results/validation/verify_report_latest.md.'}

| Check | Result | Detail |
| --- | --- | --- |
{rows}

## What this status does and does not mean

- Verified here: integrity of the files, absence of local paths and private addresses, portable launchers, and exact regeneration of Table 3, the claim checks, the figure source data and the manuscript from the stored result tables.
- Not verified by the build: the unit tests (run `scripts\\run_tests.cmd`), and the reproduction of the analysis from raw data, which has not been attempted on any other computer.
- `requirements-lock.txt`: {'full pip freeze of the build environment' if freeze_used else 'INCOMPLETE, versions from the logs only; rebuild on the analysis computer'}.
- License: MIT for the code (`LICENSE`), CC BY 4.0 for the manuscript, documents, tables and figures (`LICENSE-CC-BY-4.0.md`); decision D-022.
- The release is therefore not yet labeled fully reproducible.
"""


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--no-freeze", action="store_true", help="do not run pip freeze (lock file from logged versions)")
    ap.add_argument("--skip-verify", action="store_true")
    ap.add_argument("--snapshot", help="copy the verified release to submission_ready/<NAME> (never overwrites)")
    a = ap.parse_args(argv)
    t0 = dt.datetime.now()
    log = lambda m: print(m, flush=True)
    dest = build(PROJECT_ROOT, not a.no_freeze, log)
    n = write_sums(dest)
    log(f"SHA256SUMS.txt: {n} files")
    report = {"passed": False, "checks": [{"name": "verification", "ok": False, "detail": "skipped (--skip-verify)"}]}
    if not a.skip_verify:
        from hb import verify_release
        report = verify_release.verify(dest, log=log)
    (dest / "docs" / "RELEASE_STATUS.md").write_bytes(status_doc(dest, report, not a.no_freeze).encode("utf-8"))
    log(f"verification: {'PASSED' if report['passed'] else 'FAILED'}")
    if a.snapshot:
        if not re.fullmatch(r"[A-Za-z0-9._-]+", a.snapshot) or a.snapshot == "latest":
            log("invalid snapshot name")
            return 2
        snap = PROJECT_ROOT / "submission_ready" / a.snapshot
        if snap.exists():
            log(f"snapshot {snap.name} exists; not overwritten")
            return 2
        if not report["passed"]:
            log("verification failed; snapshot not created")
            return 2
        shutil.copytree(dest, snap, ignore=shutil.ignore_patterns(".git"))
        log(f"snapshot written: submission_ready/{snap.name}")
    append_execution_log([
        f"## {dt.datetime.now():%Y-%m-%d %H:%M:%S} release built (hb.build_release)",
        f"- Command: python -m hb.build_release {' '.join(sys.argv[1:])}; Python {platform.python_version()} on {platform.system()}",
        f"- Output: submission_ready/latest ({n} files in SHA256SUMS.txt); verification {'PASSED' if report['passed'] else 'FAILED'}; "
        f"{(dt.datetime.now() - t0).total_seconds():.0f} s" + (f"; snapshot submission_ready/{a.snapshot}" if a.snapshot else ""),
    ])
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    sys.exit(main())
