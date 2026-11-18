"""Tests of the pilot coherence estimand on synthetic sessions with known truth."""
import sys
import unittest
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
warnings.simplefilter("ignore")
from hb import pilot_coherence as PC  # noqa: E402
from hb import synth as S  # noqa: E402
import test_null_benchmark_real as tn  # noqa: E402

FS = PC.FS


def session(seed, coupled):
    rng = np.random.default_rng(seed)
    blocks, dur = tn.fake_timing(rng, n_blocks=8)
    for i, b in enumerate(blocks):
        b.update(index=i + 1, n_trials=10, boundary=False)
    n = int(dur * FS)
    shared = 0.8 * S.band_noise(n, FS, (0.07, 0.13), rng) if coupled else None
    heart = tn.fake_signal(blocks, dur, rng, shared)
    brain = tn.fake_signal(blocks, dur, rng, shared)
    return heart, brain, blocks, n


class Estimand(unittest.TestCase):
    def run_one(self, seed, coupled, marked=None, interp_block=None):
        heart, brain, blocks, n = session(seed, coupled)
        interp = np.zeros(n, bool)
        if interp_block is not None:
            b = blocks[interp_block]
            interp[int(b["t0"] * FS):int(b["t1"] * FS)] = True
        marked = marked or [0.0] * len(blocks)
        return PC.run_participant("sub-999", heart, {"alpha": brain}, {"alpha": interp}, blocks, marked, lambda m: None)

    def test_no_coupling_gives_excess_near_zero(self):
        vals = []
        for seed in range(4):
            _, s = self.run_one(seed, False)
            vals.append(next(x for x in s if x["band"] == "LF")["mean_excess"])
        self.assertLess(abs(np.mean(vals)), 0.03)

    def test_coupling_gives_positive_excess_in_the_right_band(self):
        _, s = self.run_one(1, True)
        lf = next(x for x in s if x["band"] == "LF")
        hf = next(x for x in s if x["band"] == "HF")
        self.assertGreater(lf["mean_excess"], 0.08)
        self.assertGreater(lf["boot95_low"], 0.0)
        self.assertLess(abs(hf["mean_excess"]), 0.04)
        self.assertEqual(lf["blocks_retained"], 8)

    def test_exclusion_rules(self):
        marked = [0.0] * 8
        marked[2] = 0.2
        rows, s = self.run_one(1, True, marked=marked, interp_block=5)
        reasons = {r["block"]: r["excluded"] for r in rows}
        self.assertIn("NN intervals marked", reasons[3])
        self.assertIn("envelope samples interpolated", reasons[6])
        self.assertEqual(next(x for x in s if x["band"] == "LF")["blocks_retained"], 6)
        self.assertEqual(next(x for x in s if x["band"] == "LF")["pct_retained"], 75.0)

    def test_go_check_counts_participants(self):
        mk = lambda sub, e, band, pct: {"sub": sub, "eeg": e, "band": band, "pct_retained": pct, "mean_coi_share": 0.7}
        rows = [mk(sub, e, band, 100.0 if not (sub == "b" and e != "beta") else 50.0)
                for sub in "ab" for e in PC.EEG for band in PC.BANDS]
        lines, verdict = PC.go_check(rows)
        self.assertIn("passes: True", lines[0])
        self.assertIn("passes: False", lines[1])
        self.assertTrue(verdict.startswith("1 of 2"))


if __name__ == "__main__":
    unittest.main()
