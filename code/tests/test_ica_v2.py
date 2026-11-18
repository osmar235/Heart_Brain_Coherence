"""Tests of step 16a (ICA fit stage) on synthetic recordings."""
import json
import sys
import tempfile
import unittest
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
warnings.simplefilter("ignore")
from hb import ica_v2 as IV  # noqa: E402
import test_pilot_qc as tq  # noqa: E402
import test_slow_signals as ts  # noqa: E402


class BadChannels(unittest.TestCase):
    def test_criteria(self):
        rng = np.random.default_rng(0)
        fs, n = 250.0, 250 * 120
        common = 15 * rng.standard_normal(n)
        from scipy.signal import butter, sosfiltfilt
        common = sosfiltfilt(butter(4, 30.0, fs=fs, output="sos"), common)
        X = np.array([common + 3 * rng.standard_normal(n) for _ in range(20)])
        X = sosfiltfilt(butter(4, [1.0, 45.0], btype="band", fs=fs, output="sos"), X, axis=1)
        X[3] = 0.01 * rng.standard_normal(n)                      # flat
        X[7] *= 20                                                 # deviation
        hf = sosfiltfilt(butter(4, [60.0, 100.0], btype="band", fs=fs, output="sos"), 200 * rng.standard_normal(n))
        X[11] = X[11] + hf                                         # high frequency noise
        X[15] = sosfiltfilt(butter(4, [1.0, 45.0], btype="band", fs=fs, output="sos"), 15 * rng.standard_normal(n))  # uncorrelated
        names = [f"E{i}" for i in range(20)]
        bad = IV.find_bad_channels(X, fs, names, np.ones(n, bool))
        for k in ("E3", "E7", "E11", "E15"):
            self.assertIn(k, bad)
        self.assertLessEqual(len(bad), 5)
        IV.CRIT["max_bad_channels"], old = 2, IV.CRIT["max_bad_channels"]
        try:
            capped = IV.find_bad_channels(X, fs, names, np.ones(n, bool))
        finally:
            IV.CRIT["max_bad_channels"] = old
        self.assertEqual(set(capped), {"E3", "E7"})


class FitStep(unittest.TestCase):
    def test_fit_subject_runs_and_writes_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            raw = tmp / "raw"
            tq.make_recording(raw, eeg_labels=ts.LABELS, seconds=400, n_blocks=3, trials=6, common=60.0)
            old_out, old_root = IV.OUT, IV.bidsio.raw_root
            IV.OUT = lambda *p: tmp.joinpath("project", *p)
            try:
                sub, status = IV.fit_subject((str(raw), "sub-999", "picard"))
                self.assertTrue(status.startswith("ok"), status)
                d = tmp / "project" / "derivatives" / "ica_v2" / "picard"
                self.assertTrue((d / "sub-999-ica.fif").exists())
                qc = json.loads((d / "sub-999_qc.json").read_text())
                self.assertEqual(qc["n_components"], len(ts.LABELS) - qc["n_bad"] - 1)
                lines = (d / "sub-999_iclabel.tsv").read_text().strip().split("\n")
                self.assertEqual(len(lines) - 1, qc["n_components"])
                p = np.array([[float(v) for v in ln.split("\t")[1:8]] for ln in lines[1:]])
                np.testing.assert_allclose(p.sum(axis=1), 1.0, atol=1e-3)
                self.assertEqual(IV.fit_subject((str(raw), "sub-999", "picard"))[1], "skipped (already done)")
                msgs = []
                IV.summarize("picard", msgs.append)
                self.assertTrue(msgs and msgs[0].startswith("summary picard: 1 participants"))
            finally:
                IV.OUT = old_out


if __name__ == "__main__":
    unittest.main()
