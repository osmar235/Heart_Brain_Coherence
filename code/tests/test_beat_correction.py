"""Tests of the beat correction policy (proposal D-007) on synthetic beats with known damage."""
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
from hb import beat_correction, bidsio, cardiac  # noqa: E402
from hb import synth as S  # noqa: E402
import test_pilot_qc as tq  # noqa: E402

FS = 1000.0
N = int(900 * FS)


def beats_with_hrv(depth, seed=3):
    rng = np.random.default_rng(seed)
    mod = 0.7 * S.band_noise(N, FS, (0.04, 0.15), rng) + 0.7 * np.sin(2 * np.pi * 0.25 * np.arange(N) / FS)
    return S.ipfm_beats(mod, FS, hr_mean=1.2, depth=depth)


def damage(b, seed=5, k=12):
    rng = np.random.default_rng(seed)
    rr = np.diff(b)
    miss = rng.choice(np.arange(50, b.size - 50), k, replace=False)
    extra = rng.choice(np.arange(50, b.size - 50), k, replace=False)
    return np.sort(np.r_[np.delete(b, miss), b[extra] + rng.uniform(0.25, 0.45, k) * rr[extra]])


def tachogram_error_ms(b, times, nn, good):
    grid = np.arange(5, 890, 0.25)
    true = np.interp(grid, b[1:], np.diff(b))
    est = np.interp(grid, times[1:][good], nn[good])
    return 1000 * np.sqrt(np.mean((est - true) ** 2))


class Policy(unittest.TestCase):
    def test_clean_normal_variability_is_untouched(self):
        b = beats_with_hrv(0.04)
        res = cardiac.correct_beats(b)
        self.assertEqual(res["times"].size, b.size)
        self.assertEqual(int(res["artifact"].sum()), 0)

    def test_very_high_real_variability_is_mostly_kept(self):
        b = beats_with_hrv(0.16)                                  # RMSSD near 100 ms, like pilot sub-124
        rr = np.diff(b)
        self.assertGreater(1000 * np.sqrt(np.mean(np.diff(rr) ** 2)), 80)
        old_rule = (cardiac.rr_flags((b * FS).astype(int), FS)[1] > 0).mean()
        new_rule = (cardiac.classify_intervals(rr) > 0).mean()
        self.assertGreater(old_rule, 0.05)
        self.assertLess(new_rule, 0.03)

    def test_missed_and_extra_beats_are_repaired(self):
        b = beats_with_hrv(0.04)
        dam = damage(b)
        res = cardiac.correct_beats(dam)
        self.assertEqual(res["counts"]["extra_removed"], 12)
        self.assertEqual(res["counts"]["missed_inserted"], 12)
        raw_err = tachogram_error_ms(b, dam, np.diff(dam), np.ones(dam.size - 1, bool))
        fix_err = tachogram_error_ms(b, res["times"], res["nn"], ~res["artifact"])
        self.assertGreater(raw_err, 50)
        self.assertLess(fix_err, 6)

    def test_repair_under_high_variability_is_imperfect_but_helps(self):
        b = beats_with_hrv(0.16)
        dam = damage(b)
        res = cardiac.correct_beats(dam)
        raw_err = tachogram_error_ms(b, dam, np.diff(dam), np.ones(dam.size - 1, bool))
        fix_err = tachogram_error_ms(b, res["times"], res["nn"], ~res["artifact"])
        self.assertLess(fix_err, 0.5 * raw_err)
        self.assertLess(res["artifact"].mean(), 0.08)

    def test_template_screen_rejects_noise_detections(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            beats = tq.make_recording(root)
            rec = bidsio.load_set(bidsio.rec_paths(root, "sub-999", "memory")["set"], channels=["ECG"])
            det = cardiac.detect_rpeaks(rec["data"][0], FS)
            fake = np.array([int((beats[20] + 0.35) * FS), int((beats[60] + 0.40) * FS), int((beats[90] + 0.30) * FS)])
            pk = np.sort(np.r_[det["peaks"], fake])
            keep, corr = cardiac.template_screen(det["filtered"], pk, FS)
            self.assertTrue((~keep[np.isin(pk, fake)]).all())
            self.assertGreater(keep[~np.isin(pk, fake)].mean(), 0.995)

    def test_full_step_runs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tq.make_recording(root)
            row, blocks = beat_correction.process(root, "sub-999", "memory", root / "o", root / "b")
            self.assertLess(row["pct_beats_changed"], 1.0)
            self.assertGreater(row["ppg_corr_corrected"], 0.9)
            self.assertEqual(len(blocks), 2)
            self.assertTrue((root / "b" / "sub-999_task-memory_beats.tsv").exists())


if __name__ == "__main__":
    unittest.main()
