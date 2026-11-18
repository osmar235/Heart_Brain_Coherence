"""Environment and connectivity check. Standard library only. Downloads a few kilobytes at most."""
from __future__ import annotations

import datetime as dt
import hashlib
import platform
import shutil
import sys
import urllib.request
from pathlib import Path

from hb.download_ds008302 import PROJECT_ROOT, UA, load_json, raw_url, read_manifest, resolve_data_root, s3_urls, sha256_file


def main() -> int:
    out = []
    fails = 0

    def say(tag: str, msg: str) -> None:
        nonlocal fails
        if tag == "FAIL":
            fails += 1
        line = f"[{tag}] {msg}"
        print(line, flush=True)
        out.append(line)

    say("INFO", f"time {dt.datetime.now().isoformat(timespec='seconds')}")
    say("INFO", f"python {platform.python_version()} at {sys.executable}")
    say("INFO", f"platform {platform.platform()}")
    say("OK" if sys.version_info >= (3, 10) else "FAIL", "python version is 3.10 or newer")
    in_venv = sys.prefix != getattr(sys, "base_prefix", sys.prefix)
    say("OK" if in_venv else "WARN", f"running inside project .venv: {in_venv}")
    say("INFO", f"project root {PROJECT_ROOT}")

    cfg = load_json(PROJECT_ROOT / "config" / "download_config.json")
    mpath = PROJECT_ROOT / cfg["manifest"]
    rows = read_manifest(mpath)
    say("OK", f"manifest readable: {len(rows)} files, {sum(r['size_bytes'] for r in rows) / 1e9:.1f} GB in full dataset, sha256 {sha256_file(mpath)[:16]}...")

    root = resolve_data_root(cfg)
    probe = root
    while not probe.exists():
        probe = probe.parent
    free = shutil.disk_usage(probe).free / 1e9
    pilot = cfg["pilot"]
    need = sum(r["size_bytes"] for r in rows if r["subject"] in pilot["subjects"] and "/eeg/" in r["path"]) / 1e9
    say("INFO", f"raw data root {root}")
    say("OK" if free > need + cfg["free_space_margin_gb"] else "FAIL", f"free space {free:.1f} GB; pilot needs about {need:.1f} GB plus {cfg['free_space_margin_gb']} GB margin")
    if "onedrive" in str(root).lower():
        say("WARN", "raw data root is inside OneDrive; large files will sync. See config/local_paths.example.json")

    by = {r["path"]: r for r in rows}
    # 1. GitHub pinned text file
    r = by["participants.tsv"]
    try:
        b = urllib.request.urlopen(urllib.request.Request(raw_url(cfg, r), headers={"User-Agent": UA}), timeout=30).read()
        say("OK" if hashlib.sha256(b).hexdigest() == r["sha256"] else "FAIL", "GitHub pinned commit reachable and participants.tsv matches manifest SHA256")
    except Exception as e:
        say("FAIL", f"GitHub raw not reachable: {type(e).__name__}: {e}")
    # 2. smallest S3 object, full verification
    small = min((x for x in rows if x["source"] == "openneuro_s3"), key=lambda x: x["size_bytes"])
    ok = False
    for u in s3_urls(cfg, small):
        try:
            b = urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent": UA}), timeout=30).read()
            ok = hashlib.sha256(b).hexdigest() == small["sha256"]
            say("OK" if ok else "FAIL", f"OpenNeuro S3 object fetched ({'with' if 'versionId' in u else 'without'} version id), SHA256 {'matches' if ok else 'DIFFERS'}: {small['path']}")
            if ok:
                break
        except Exception as e:
            say("WARN", f"S3 fetch {'with' if 'versionId' in u else 'without'} version id failed: {type(e).__name__}: {e}")
    if not ok:
        say("FAIL", "could not fetch and verify a small file from OpenNeuro S3")
    # 3. range request on a pilot recording (first 1024 bytes only)
    big = by[f"{pilot['subjects'][0]}/eeg/{pilot['subjects'][0]}_task-rest_eeg.set"]
    try:
        req = urllib.request.Request(s3_urls(cfg, big)[0], headers={"User-Agent": UA, "Range": "bytes=0-1023"})
        with urllib.request.urlopen(req, timeout=30) as resp:
            cr = resp.headers.get("Content-Range", "")
            total = int(cr.split("/")[-1]) if "/" in cr else -1
            good = resp.status == 206 and total == big["size_bytes"]
            say("OK" if good else "FAIL", f"resume (Range) supported and remote size equals manifest for {big['path']}: status {resp.status}, remote {total}, manifest {big['size_bytes']}")
    except Exception as e:
        say("FAIL", f"range request on pilot file failed: {type(e).__name__}: {e}")

    logdir = PROJECT_ROOT / "logs"
    logdir.mkdir(exist_ok=True)
    lp = logdir / f"environment_{dt.datetime.now().strftime('%Y%m%d_%H%M%S')}.txt"
    lp.write_text("\n".join(out) + "\n", encoding="utf-8")
    print(f"\nSaved to {lp}")
    print("RESULT: " + ("READY" if fails == 0 else f"{fails} CHECK(S) FAILED, send the full output back before downloading"))
    return 0 if fails == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
