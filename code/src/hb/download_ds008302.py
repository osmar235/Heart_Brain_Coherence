"""Resumable, checksum-verified downloader for NEMAR on008302 / OpenNeuro ds008302 v1.0.0.

Standard library only. Run through the .cmd wrappers in scripts/.

Sources
  - Small text files (json, tsv, README): GitHub mirror nemarDatasets/on008302 at a pinned commit.
  - Large files (.set, .tsv.gz): OpenNeuro public S3 bucket, pinned S3 version id.
Every file is checked against the size and SHA256 recorded in
data_manifest/ds008302_v1.0.0_manifest.tsv (built from the dataset's own git-annex keys).
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import io
import json
import os
import platform
import shutil
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
import os
import zipfile
from pathlib import Path

# HB_PROJECT_ROOT lets a relocated copy (for example the release package, where the code sits in code/src)
# name its own root; without it the root is two levels above this file, as before.
PROJECT_ROOT = Path(os.environ["HB_PROJECT_ROOT"]).resolve() if os.environ.get("HB_PROJECT_ROOT") else Path(__file__).resolve().parents[2]
CHUNK = 1024 * 1024
UA = "hb-ds008302-downloader/0.1 (python-urllib)"


class Log:
    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.path = path
        self.fh = open(path, "a", encoding="utf-8", newline="\n")

    def __call__(self, msg: str) -> None:
        line = f"{dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} {msg}"
        print(line, flush=True)
        self.fh.write(line + "\n")
        self.fh.flush()


def load_json(path: Path) -> dict:
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def resolve_data_root(cfg: dict) -> Path:
    lp = PROJECT_ROOT / "config" / "local_paths.json"
    root = cfg["default_raw_data_root"]
    if lp.exists():
        root = load_json(lp).get("raw_data_root", root)
    root = os.environ.get("HB_DATA_ROOT") or root  # the environment variable wins over both files
    p = Path(root)
    return p if p.is_absolute() else (PROJECT_ROOT / p)


def read_manifest(path: Path) -> list[dict]:
    with open(path, encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for r in rows:
        r["size_bytes"] = int(r["size_bytes"])
    return rows


def select(rows: list[dict], mode: str, cfg: dict, subjects=None, tasks=None, include_beh=None) -> list[dict]:
    if mode == "metadata":
        return [r for r in rows if r["source"] == "github_raw"]
    if mode == "cohort":
        c = cfg["cohort"]
        tasks = tasks or c["tasks"]
        include_beh = False if include_beh is None else include_beh
        want = {f"{r['subject']}" for r in rows if r["path"].endswith("_task-memory_eeg.set")}
        subjects = subjects or sorted(want)
    if mode == "pilot":
        p = cfg["pilot"]
        subjects = subjects or p["subjects"]
        tasks = tasks or p["tasks"]
        include_beh = p["include_beh_eyetracking"] if include_beh is None else include_beh
    if not subjects or not tasks:
        raise SystemExit("ERROR: subjects and tasks are required for this mode")
    known = {r["subject"] for r in rows}
    unknown = [s for s in subjects if s not in known]
    if unknown:
        raise SystemExit(f"ERROR: subjects not in manifest: {unknown}")
    out = []
    for r in rows:
        if r["subject"] == "root":
            out.append(r)
            continue
        if r["subject"] not in subjects:
            continue
        in_beh = "/beh/" in r["path"]
        if in_beh and not include_beh:
            continue
        if r["task"] != "na" and r["task"] not in tasks:
            continue
        out.append(r)
    return out


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(CHUNK)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def s3_urls(cfg: dict, row: dict) -> list[str]:
    base = cfg["openneuro_s3_base"].rstrip("/") + "/" + urllib.parse.quote(row["path"])
    urls = []
    if row.get("s3_version_id"):
        urls.append(base + "?versionId=" + urllib.parse.quote(row["s3_version_id"], safe=""))
    urls.append(base)
    return urls


def raw_url(cfg: dict, row: dict) -> str:
    return cfg["github_raw_base"].rstrip("/") + "/" + cfg["github_commit"] + "/" + urllib.parse.quote(row["path"])


def http_download(url: str, part: Path, expected: int, log: Log, timeout: int, label: str) -> None:
    """Download url into part, resuming from the current size of part."""
    part.parent.mkdir(parents=True, exist_ok=True)
    have = part.stat().st_size if part.exists() else 0
    if have > expected:
        part.unlink()
        have = 0
    if have == expected:
        return
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    if have:
        req.add_header("Range", f"bytes={have}-")
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        if have and resp.status != 206:
            log(f"    server ignored resume request, restarting {label}")
            have = 0
        mode = "ab" if have else "wb"
        t0 = last = time.time()
        got = have
        with open(part, mode) as fh:
            while True:
                b = resp.read(CHUNK)
                if not b:
                    break
                fh.write(b)
                got += len(b)
                now = time.time()
                if now - last >= 15:
                    rate = (got - have) / max(now - t0, 1e-6) / 1e6
                    log(f"    {label}: {got / 1e6:,.0f} of {expected / 1e6:,.0f} MB ({100 * got / max(expected, 1):.1f}%), {rate:.1f} MB/s")
                    last = now


def fetch_one(row: dict, urls: list[str], dest: Path, log: Log, cfg: dict, state: dict) -> str:
    """Return a status string. Final file exists only after size and SHA256 pass."""
    expected, sha = row["size_bytes"], row["sha256"]
    key = row["path"]
    if dest.exists():
        st = dest.stat()
        cached = state.get(key)
        if st.st_size == expected and cached and cached.get("sha256") == sha and cached.get("mtime_ns") == st.st_mtime_ns:
            return "ok_already_verified"
        if st.st_size == expected and sha256_file(dest) == sha:
            state[key] = {"sha256": sha, "mtime_ns": st.st_mtime_ns}
            return "ok_existing_verified"
        bad = dest.with_name(dest.name + ".mismatch")
        if bad.exists():
            bad.unlink()
        dest.rename(bad)
        log(f"    existing file failed verification, kept as {bad.name}")
    part = dest.with_name(dest.name + ".part")
    last_err = "no attempt"
    for url in urls:
        for attempt in range(1, cfg["retries"] + 1):
            try:
                http_download(url, part, expected, log, cfg["timeout_seconds"], Path(key).name)
                size = part.stat().st_size
                if size < expected:
                    raise IOError(f"incomplete: {size} of {expected} bytes")
                if size > expected:
                    part.unlink()
                    raise IOError(f"too large: {size} of {expected} bytes")
                got = sha256_file(part)
                if got != sha:
                    part.unlink()
                    raise IOError(f"SHA256 mismatch (got {got[:16]}...)")
                os.replace(part, dest)
                state[key] = {"sha256": sha, "mtime_ns": dest.stat().st_mtime_ns}
                return "ok_downloaded"
            except urllib.error.HTTPError as e:
                last_err = f"HTTP {e.code} {e.reason}"
                if e.code == 416 and part.exists():
                    part.unlink()
                if e.code in (400, 401, 403, 404):
                    log(f"    {last_err} from {url.split('?')[0]}; trying next source")
                    break
            except Exception as e:  # network drop, timeout, verification failure
                last_err = f"{type(e).__name__}: {e}"
            log(f"    attempt {attempt}/{cfg['retries']} failed: {last_err}")
            time.sleep(min(60, 3 * attempt))
    return "FAILED: " + last_err


def fetch_metadata_zip(sel: list[dict], root: Path, log: Log, cfg: dict, state: dict) -> dict:
    """Try one zip of the pinned commit; return {path: status} for files it satisfied."""
    done = {}
    todo = [r for r in sel if r["source"] == "github_raw" and not (root / r["path"]).exists()]
    if len(todo) < cfg["zip_threshold_files"]:
        return done
    url = cfg["github_zip_base"].rstrip("/") + "/" + cfg["github_commit"]
    log(f"  fetching pinned metadata archive ({len(todo)} files wanted): {url}")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=cfg["timeout_seconds"]) as resp:
            blob = resp.read()
        zf = zipfile.ZipFile(io.BytesIO(blob))
        prefix = zf.namelist()[0].split("/")[0] + "/"
        for r in todo:
            try:
                b = zf.read(prefix + r["path"])
            except KeyError:
                continue
            if len(b) == r["size_bytes"] and hashlib.sha256(b).hexdigest() == r["sha256"]:
                dest = root / r["path"]
                dest.parent.mkdir(parents=True, exist_ok=True)
                tmp = dest.with_name(dest.name + ".part")
                tmp.write_bytes(b)
                os.replace(tmp, dest)
                state[r["path"]] = {"sha256": r["sha256"], "mtime_ns": dest.stat().st_mtime_ns}
                done[r["path"]] = "ok_downloaded_zip"
        log(f"  archive supplied {len(done)} verified files; the rest go file by file")
    except Exception as e:
        log(f"  archive route unavailable ({type(e).__name__}: {e}); falling back to file by file")
    return done


def append_execution_log(lines: list[str]) -> None:
    if os.environ.get("HB_NO_EXECUTION_LOG"):  # verification reruns must not add entries
        return
    p = PROJECT_ROOT / "docs" / "EXECUTION_LOG.md"
    p.parent.mkdir(parents=True, exist_ok=True)
    new = not p.exists()
    with open(p, "a", encoding="utf-8", newline="\n") as fh:
        if new:
            fh.write("# EXECUTION LOG\n\nMechanical record of runs. Scientific decisions live in DECISOES.md.\n")
        fh.write("\n" + "\n".join(lines) + "\n")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mode", required=True, choices=["metadata", "pilot", "cohort", "custom"])
    ap.add_argument("--subjects", nargs="*", help="e.g. sub-124 sub-125 (custom mode, or override pilot list)")
    ap.add_argument("--tasks", nargs="*", choices=["memory", "rest", "video"])
    ap.add_argument("--include-beh", action="store_true", default=None, help="also fetch eye tracking files in beh/")
    ap.add_argument("--dry-run", action="store_true", help="list what would be fetched and stop")
    ap.add_argument("--verify-only", action="store_true", help="check files already on disk, download nothing")
    ap.add_argument("--allow-large", action="store_true", help="permit selections above the configured size cap")
    ap.add_argument("--config", default=str(PROJECT_ROOT / "config" / "download_config.json"))
    a = ap.parse_args(argv)

    cfg = load_json(Path(a.config))
    stamp = dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    log = Log(PROJECT_ROOT / "logs" / f"download_{a.mode}_{stamp}.log")
    root = resolve_data_root(cfg)
    manifest_path = PROJECT_ROOT / cfg["manifest"]
    rows = read_manifest(manifest_path)
    sel = select(rows, a.mode, cfg, a.subjects, a.tasks, a.include_beh)
    total = sum(r["size_bytes"] for r in sel)

    log(f"mode={a.mode} dry_run={a.dry_run} verify_only={a.verify_only}")
    log(f"python {platform.python_version()} on {platform.platform()}")
    log(f"project root: {PROJECT_ROOT}")
    log(f"raw data root: {root}")
    log(f"manifest: {manifest_path.name} sha256={sha256_file(manifest_path)}")
    log(f"dataset version {cfg['dataset_version']}, GitHub commit {cfg['github_commit']}")
    log(f"selected {len(sel)} files, {total / 1e9:.3f} GB")
    by_sub = {}
    for r in sel:
        by_sub[r["subject"]] = by_sub.get(r["subject"], 0) + r["size_bytes"]
    if len(by_sub) <= 12:
        for s, b in sorted(by_sub.items()):
            log(f"  {s}: {b / 1e9:.3f} GB")
    if "onedrive" in str(root).lower():
        log("WARNING: raw data root is inside OneDrive. Multi GB files will sync to the cloud. "
            "Consider pointing config/local_paths.json at a folder outside OneDrive.")

    if a.dry_run:
        for r in sel:
            if r["size_bytes"] > 50e6:
                log(f"  would fetch {r['path']} ({r['size_bytes'] / 1e9:.2f} GB)")
        log("dry run only, nothing downloaded")
        return 0

    cap = cfg["max_selection_gb"]
    if total / 1e9 > cap and not a.allow_large and not a.verify_only:
        log(f"ERROR: selection is {total / 1e9:.1f} GB, above the {cap} GB cap. Re-run with --allow-large if intended.")
        return 3

    root.mkdir(parents=True, exist_ok=True)
    state_path = root / ".download_state.json"
    state = load_json(state_path) if state_path.exists() else {}

    if not a.verify_only:
        missing = sum(r["size_bytes"] for r in sel if not (root / r["path"]).exists())
        free = shutil.disk_usage(root).free
        need = missing + int(cfg["free_space_margin_gb"] * 1e9)
        log(f"to download: {missing / 1e9:.3f} GB; free on target drive: {free / 1e9:.1f} GB")
        if free < need:
            log(f"ERROR: not enough free space (need {need / 1e9:.1f} GB including margin). Nothing downloaded.")
            return 3

    results = {}
    if not a.verify_only:
        results.update(fetch_metadata_zip(sel, root, log, cfg, state))
    n = len(sel)
    for i, r in enumerate(sel, 1):
        if r["path"] in results:
            continue
        dest = root / r["path"]
        big = r["size_bytes"] > 50e6
        if big:
            log(f"[{i}/{n}] {r['path']} ({r['size_bytes'] / 1e9:.2f} GB)")
        if a.verify_only:
            if not dest.exists():
                status = "MISSING"
            elif dest.stat().st_size == r["size_bytes"] and sha256_file(dest) == r["sha256"]:
                status = "ok_existing_verified"
            else:
                status = "FAILED: verification mismatch"
        else:
            urls = [raw_url(cfg, r)] if r["source"] == "github_raw" else s3_urls(cfg, r)
            status = fetch_one(r, urls, dest, log, cfg, state)
        results[r["path"]] = status
        if big or not status.startswith("ok"):
            log(f"    -> {status}")
        if big or i % 200 == 0:
            with open(state_path, "w", encoding="utf-8") as fh:
                json.dump(state, fh)
    with open(state_path, "w", encoding="utf-8") as fh:
        json.dump(state, fh)

    report = PROJECT_ROOT / "data_manifest" / f"download_report_{a.mode}_{stamp}.tsv"
    with open(report, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("path\tsize_bytes\tsha256\tstatus\n")
        for r in sel:
            fh.write(f"{r['path']}\t{r['size_bytes']}\t{r['sha256']}\t{results[r['path']]}\n")
    counts = {}
    for s in results.values():
        k = s.split(":")[0]
        counts[k] = counts.get(k, 0) + 1
    bad = [p for p, s in results.items() if not s.startswith("ok")]
    log(f"summary: {counts}")
    log(f"report: {report}")
    append_execution_log([
        f"## {dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')} download mode={a.mode}" + (" (verify only)" if a.verify_only else ""),
        f"- Command: python -m hb.download_ds008302 {' '.join(sys.argv[1:])}",
        f"- Python {platform.python_version()}, {platform.platform()}",
        f"- Dataset ds008302 v{cfg['dataset_version']}, GitHub commit {cfg['github_commit']}, manifest sha256 {sha256_file(manifest_path)}",
        f"- Raw data root: {root}",
        f"- Selected {len(sel)} files ({total / 1e9:.3f} GB); result counts: {counts}",
        f"- Report: data_manifest/{report.name}; log: logs/{log.path.name}",
        f"- Outcome: {'ALL VERIFIED (size and SHA256)' if not bad else str(len(bad)) + ' FILES NOT VERIFIED, see report'}",
    ])
    if bad:
        log(f"NOT COMPLETE: {len(bad)} files missing or failed. Re-run the same command to resume.")
        return 2
    log("ALL SELECTED FILES PRESENT AND VERIFIED (size and SHA256).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
