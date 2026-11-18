"""Tests of the cohort estimands (two variants) and of the frozen group statistics."""
import sys
import unittest
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
warnings.simplefilter("ignore")
from hb import cohort as C  # noqa: E402
from hb import cohort_stats as CS  # noqa: E402
import test_pilot_coherence as tp  # noqa: E402


class Variants(unittest.TestCase):
    def test_mask_covers_trial_start_to_prompt_only(self):
        _, _, blocks, _ = tp.session(0, False)
        b = blocks[0]
        n = int(round((b["t1"] - b["t0"]) * C.FS))
        m = C.masks(b, n)
        self.assertTrue(m["full"].all())
        expected = sum(p - s for s, p in zip(b["start"], b["probe"])) * C.FS
        self.assertAlmostEqual(m["encret"].sum(), expected, delta=12)
        self.assertLess(m["encret"].mean(), 0.9)

    def test_both_variants_detect_coupling_and_stay_near_zero_without(self):
        res = {}
        for coupled in (False, True):
            heart, brain, blocks, n = tp.session(1, coupled)
            rows = C.coherence_tables("sub-999", heart, {"alpha": brain}, {"alpha": np.zeros(n, bool)}, blocks, [0.0] * len(blocks))
            for v in C.VARIANTS:
                res[(coupled, v)] = np.mean([r[f"excess_{v}_LF"] for r in rows if not r["excluded"]])
        for v in C.VARIANTS:
            self.assertLess(abs(res[(False, v)]), 0.05)
            self.assertGreater(res[(True, v)], 0.07)


class EndToEnd(unittest.TestCase):
    def test_run_subject_from_a_raw_recording_into_empty_folders(self):
        """Regression test for the 2026-11-10 failure: output folders that do not exist yet."""
        import tempfile
        import test_pilot_qc as tq
        import test_slow_signals as ts
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            tq.make_recording(tmp / "raw", eeg_labels=ts.LABELS)
            old = C.PROJECT_ROOT
            C.PROJECT_ROOT = tmp / "project"
            try:
                status = C.run_subject(tmp / "raw", "sub-999", lambda m: None)
                again = C.run_subject(tmp / "raw", "sub-999", lambda m: None)
            finally:
                C.PROJECT_ROOT = old
            self.assertTrue(status.startswith("ok"), status)
            self.assertTrue(again.startswith("skipped"))
            self.assertTrue((tmp / "project" / "derivatives" / "cohort_coherence" / "sub-999_blocks.tsv").exists())
            self.assertTrue((tmp / "project" / "results" / "cohort_qc" / "slow" / "sub-999_slow_summary.tsv").exists())


class MatlabV73(unittest.TestCase):
    def test_hdf5_recording_loads_like_the_v5_one(self):
        """Twelve cohort recordings above 2 GB are MATLAB v7.3 (HDF5). Same content must come out."""
        import tempfile
        import h5py
        import test_pilot_qc as tq
        from hb import bidsio
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            tq.make_recording(tmp)
            p = bidsio.rec_paths(tmp, "sub-999", "memory")["set"]
            v5 = bidsio.load_set(p)
            h = tmp / "v73.set"
            with h5py.File(h, "w", userblock_size=512) as f:
                f["data"] = np.ascontiguousarray(v5["data"].T)
                f["srate"] = np.array([[v5["fs"]]])
                f["pnts"] = np.array([[float(v5["n_samples"])]])
                f["ref"] = np.array([[ord(c)] for c in "FCz"], dtype=np.uint16)
                refs = f.create_group("#refs#")

                def put(name, value):
                    if isinstance(value, str):
                        refs[name] = np.array([[ord(c)] for c in value], dtype=np.uint16)
                    else:
                        refs[name] = np.array([[value]], dtype=float)
                    return refs[name].ref
                rt = h5py.special_dtype(ref=h5py.Reference)
                g = f.create_group("chanlocs")
                g.create_dataset("labels", data=np.array([[put(f"c{i}", lab)] for i, lab in enumerate(v5["all_labels"])], dtype=rt))
                e = f.create_group("event")
                n = len(v5["events"])
                e.create_dataset("type", data=np.array([[put(f"t{i}", ev[0])] for i, ev in enumerate(v5["events"])], dtype=rt))
                e.create_dataset("latency", data=np.array([[put(f"l{i}", ev[1] * v5["fs"] + 1.0)] for i, ev in enumerate(v5["events"])], dtype=rt))
                e.create_dataset("duration", data=np.array([[put(f"d{i}", ev[2])] for i, ev in enumerate(v5["events"])], dtype=rt))
            with open(h, "r+b") as fh:                       # the header MATLAB writes in front of the HDF5 content
                head = b"MATLAB 7.3 MAT-file, Platform: PCWIN64, Created on: test HDF5 schema 1.00 ."
                fh.write(head.ljust(116, b" ") + b"\x00" * 8 + b"\x00\x02" + b"IM")
            with self.assertRaises(NotImplementedError):
                bidsio._load_v5(h)
            v73 = bidsio.load_set(h)
            self.assertEqual(v73["all_labels"], v5["all_labels"])
            self.assertEqual(v73["data"].shape, v5["data"].shape)
            np.testing.assert_array_equal(v73["data"], v5["data"])
            self.assertEqual([x[0] for x in v73["events"]], [x[0] for x in v5["events"]])
            np.testing.assert_allclose([x[1] for x in v73["events"]], [x[1] for x in v5["events"]], atol=1e-9)
            self.assertEqual(v73["ref"], "FCz")
            sub = bidsio.load_set(h, channels=["ECG", "PPG"])
            np.testing.assert_array_equal(sub["data"][0], v5["data"][-1])


def fake_tables(rng, n_sub=30, effect=0.0, cond=0.0):
    tables = {}
    for s in range(n_sub):
        rows, k = [], 0
        for L, nb in ((5, 2), (6, 7), (7, 6), (8, 1)):
            for _ in range(nb):
                k += 1
                for e in C.EEG:
                    r = {"sub": f"sub-{s:03d}", "eeg": e, "block": str(k), "length": str(L), "t0_s": str(300.0 * k), "excluded": ""}
                    for v in C.VARIANTS:
                        for band in C.BANDS:
                            mu = (effect + (cond if L == 7 else 0.0)) if (e == "alpha" and band == "LF") else 0.0
                            r[f"excess_{v}_{band}"] = str(mu + 0.04 * rng.standard_normal())
                    rows.append(r)
        tables[f"sub-{s:03d}"] = rows
    return tables


class GroupStats(unittest.TestCase):
    def test_span_rule(self):
        t = fake_tables(np.random.default_rng(0), n_sub=1)
        self.assertEqual(CS.span_of(next(iter(t.values()))), 6)
        rows = [r for r in next(iter(t.values())) if r["length"] != "7"]
        self.assertIsNone(CS.span_of(rows))

    def test_holm(self):
        np.testing.assert_allclose(CS.holm([0.01, 0.04, 0.03]), [0.03, 0.06, 0.06])

    def test_null_and_effect(self):
        rng = np.random.default_rng(1)
        null = CS.group_tests(CS.participant_table(fake_tables(rng)))
        self.assertGreater(min(r["p_holm"] for r in null), 0.05)
        eff = CS.group_tests(CS.participant_table(fake_tables(rng, effect=0.03, cond=0.04)))
        hit = [r for r in eff if r["eeg"] == "alpha" and r["band"] == "LF"]
        self.assertTrue(all(r["p_holm"] < 0.01 and r["mean"] > 0 for r in hit))
        others = [r for r in eff if not (r["eeg"] == "alpha" and r["band"] == "LF")]
        self.assertGreater(np.mean([r["p_holm"] > 0.05 for r in others]), 0.9)
        self.assertEqual({r["n"] for r in eff}, {30})


if __name__ == "__main__":
    unittest.main()
