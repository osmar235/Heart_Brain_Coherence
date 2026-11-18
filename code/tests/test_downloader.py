"""Offline tests: a local HTTP server with Range support stands in for GitHub and S3."""
import hashlib
import http.server
import io
import json
import os
import sys
import tempfile
import threading
import unittest
import zipfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from hb import download_ds008302 as d  # noqa: E402


class RangeHandler(http.server.BaseHTTPRequestHandler):
    files = {}
    truncate = {}  # path -> bytes to send before dropping the connection (once)

    def log_message(self, *a):
        pass

    def do_GET(self):
        path = self.path.split("?")[0]
        if "versionId=bad" in self.path or path not in self.files:
            self.send_error(404)
            return
        data = self.files[path]
        start = 0
        rng = self.headers.get("Range")
        if rng:
            start = int(rng.split("=")[1].split("-")[0])
            if start >= len(data):
                self.send_error(416)
                return
            self.send_response(206)
            self.send_header("Content-Range", f"bytes {start}-{len(data) - 1}/{len(data)}")
        else:
            self.send_response(200)
        body = data[start:]
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        cut = self.truncate.pop(path, None)
        if cut is not None:
            self.wfile.write(body[:cut])
            self.wfile.flush()
            self.connection.close()
            return
        self.wfile.write(body)


class Quiet:
    def __call__(self, msg):
        pass


class DownloaderTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.srv = http.server.ThreadingHTTPServer(("127.0.0.1", 0), RangeHandler)
        threading.Thread(target=cls.srv.serve_forever, daemon=True).start()
        cls.base = f"http://127.0.0.1:{cls.srv.server_address[1]}"

    @classmethod
    def tearDownClass(cls):
        cls.srv.shutdown()

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        RangeHandler.files.clear()
        RangeHandler.truncate.clear()
        self.cfg = {"retries": 3, "timeout_seconds": 10, "openneuro_s3_base": self.base + "/s3",
                    "github_raw_base": self.base + "/raw", "github_commit": "abc", "github_zip_base": self.base + "/zip",
                    "zip_threshold_files": 2}
        d.time.sleep = lambda s: None

    def tearDown(self):
        self.tmp.cleanup()

    def row(self, path, data, source="openneuro_s3", vid="v1"):
        return {"path": path, "subject": path.split("/")[0], "task": "rest", "source": source,
                "size_bytes": len(data), "sha256": hashlib.sha256(data).hexdigest(), "s3_version_id": vid}

    def test_fresh_download_is_verified(self):
        data = os.urandom(3_000_000)
        r = self.row("sub-001/eeg/a.set", data)
        RangeHandler.files["/s3/" + r["path"]] = data
        st = d.fetch_one(r, d.s3_urls(self.cfg, r), self.root / r["path"], Quiet(), self.cfg, {})
        self.assertEqual(st, "ok_downloaded")
        self.assertEqual((self.root / r["path"]).read_bytes(), data)

    def test_resume_from_partial_file(self):
        data = os.urandom(2_500_000)
        r = self.row("sub-001/eeg/b.set", data)
        RangeHandler.files["/s3/" + r["path"]] = data
        dest = self.root / r["path"]
        dest.parent.mkdir(parents=True)
        dest.with_name(dest.name + ".part").write_bytes(data[:1_000_000])
        st = d.fetch_one(r, d.s3_urls(self.cfg, r), dest, Quiet(), self.cfg, {})
        self.assertEqual(st, "ok_downloaded")
        self.assertEqual(d.sha256_file(dest), r["sha256"])

    def test_dropped_connection_resumes(self):
        data = os.urandom(2_000_000)
        r = self.row("sub-001/eeg/c.set", data)
        RangeHandler.files["/s3/" + r["path"]] = data
        RangeHandler.truncate["/s3/" + r["path"]] = 700_000
        st = d.fetch_one(r, d.s3_urls(self.cfg, r), self.root / r["path"], Quiet(), self.cfg, {})
        self.assertEqual(st, "ok_downloaded")
        self.assertEqual(d.sha256_file(self.root / r["path"]), r["sha256"])

    def test_corrupt_remote_never_becomes_final_file(self):
        data = os.urandom(100_000)
        r = self.row("sub-001/eeg/d.set", data)
        RangeHandler.files["/s3/" + r["path"]] = bytes(len(data))
        dest = self.root / r["path"]
        st = d.fetch_one(r, d.s3_urls(self.cfg, r), dest, Quiet(), self.cfg, {})
        self.assertTrue(st.startswith("FAILED"))
        self.assertFalse(dest.exists())

    def test_bad_version_id_falls_back_to_plain_url(self):
        data = os.urandom(50_000)
        r = self.row("sub-001/eeg/e.set", data, vid="bad")
        RangeHandler.files["/s3/" + r["path"]] = data
        st = d.fetch_one(r, d.s3_urls(self.cfg, r), self.root / r["path"], Quiet(), self.cfg, {})
        self.assertEqual(st, "ok_downloaded")

    def test_existing_good_file_kept_and_bad_file_replaced(self):
        data = os.urandom(80_000)
        r = self.row("sub-001/eeg/f.set", data)
        RangeHandler.files["/s3/" + r["path"]] = data
        dest = self.root / r["path"]
        dest.parent.mkdir(parents=True)
        dest.write_bytes(data)
        self.assertEqual(d.fetch_one(r, d.s3_urls(self.cfg, r), dest, Quiet(), self.cfg, {}), "ok_existing_verified")
        dest.write_bytes(b"x" * len(data))
        self.assertEqual(d.fetch_one(r, d.s3_urls(self.cfg, r), dest, Quiet(), self.cfg, {}), "ok_downloaded")
        self.assertTrue(dest.with_name(dest.name + ".mismatch").exists())

    def test_metadata_zip_route_verifies_each_file(self):
        good, bad = b"participant_id\n", b"tampered\n"
        r1 = self.row("participants.tsv", good, source="github_raw")
        r2 = self.row("README.md", b"expected readme\n", source="github_raw")
        buf = io.BytesIO()
        with zipfile.ZipFile(buf, "w") as z:
            z.writestr("on008302-abc/participants.tsv", good)
            z.writestr("on008302-abc/README.md", bad)
        RangeHandler.files["/zip/abc"] = buf.getvalue()
        done = d.fetch_metadata_zip([r1, r2], self.root, Quiet(), self.cfg, {})
        self.assertEqual(list(done), ["participants.tsv"])
        self.assertFalse((self.root / "README.md").exists())


class ManifestTests(unittest.TestCase):
    def setUp(self):
        # The project root: HB_PROJECT_ROOT when set (release layout), else the folder above tests/.
        self.project = d.PROJECT_ROOT if os.environ.get("HB_PROJECT_ROOT") else Path(__file__).resolve().parents[1]
        self.cfg = json.loads((self.project / "config" / "download_config.json").read_text(encoding="utf-8"))
        self.rows = d.read_manifest(self.project / self.cfg["manifest"])

    def test_manifest_integrity(self):
        self.assertEqual(len(self.rows), 4506)
        self.assertEqual(sum(r["source"] == "openneuro_s3" for r in self.rows), 1754)
        for r in self.rows:
            self.assertEqual(len(r["sha256"]), 64)
            if r["source"] == "openneuro_s3":
                self.assertTrue(r["s3_version_id"])

    def test_pilot_selection(self):
        sel = d.select(self.rows, "pilot", self.cfg)
        subs = {r["subject"] for r in sel} - {"root"}
        self.assertEqual(subs, set(self.cfg["pilot"]["subjects"]))
        sets = [r for r in sel if r["path"].endswith("_eeg.set")]
        self.assertEqual(len(sets), 9)
        self.assertFalse(any("/beh/" in r["path"] for r in sel))
        self.assertLess(sum(r["size_bytes"] for r in sel) / 1e9, self.cfg["max_selection_gb"])

    def test_metadata_selection_has_no_recordings(self):
        sel = d.select(self.rows, "metadata", self.cfg)
        self.assertTrue(all(r["source"] == "github_raw" for r in sel))
        self.assertLess(sum(r["size_bytes"] for r in sel), 20e6)

    def test_unknown_subject_rejected(self):
        with self.assertRaises(SystemExit):
            d.select(self.rows, "custom", self.cfg, ["sub-999"], ["rest"])


if __name__ == "__main__":
    unittest.main()
