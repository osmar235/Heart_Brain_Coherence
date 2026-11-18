"""Synthetic benchmarks for the EEG envelope / cardiac wavelet coherence layer.

Signals with known coupling are generated with fixed seeds. No real data is used.
Expected values were established on 2026-11-09 (numpy 2.5, scipy 1.18, pycwt 0.4);
thresholds leave a margin around them.
"""
import hashlib
import sys
import unittest
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
warnings.simplefilter("ignore")
from hb import hbcoh as H  # noqa: E402
from hb import synth as S  # noqa: E402
from hb.vendor.emgsim_v3_3_0 import functions as emg  # noqa: E402

FS = 200.0
T = 600
N = int(T * FS)
VENDORED_SHA256 = "0aaf110c0ce9fc36929686ac0adb6c4fd3ea855162bf98a950b834da63c55f19"


def synthetic_pair(seed, coupled, lag_s=0.0, cardiac="tachogram"):
    """Heart beats and EEG whose slow modulations share (or do not share) a 0.07 to 0.13 Hz process."""
    rng = np.random.default_rng(seed)
    m_heart = S.band_noise(N, FS, (0.07, 0.13), rng)
    m_brain = m_heart if coupled else S.band_noise(N, FS, (0.07, 0.13), rng)
    bg_heart = S.band_noise(N, FS, (0.03, 0.45), rng)
    beats = S.ipfm_beats(0.7 * m_heart + 0.7 * bg_heart, FS)
    bg_brain = S.band_noise(N, FS, (0.03, 0.45), rng)
    eeg = S.am_eeg(0.7 * np.roll(m_brain, int(lag_s * FS)) + 0.7 * bg_brain, FS, rng)
    _, env = H.band_envelope(eeg, FS, (8, 12))
    _, card = (H.tachogram if cardiac == "tachogram" else H.beat_density)(beats, t1=T)
    ok = np.where(np.isfinite(card))[0]
    i0, i1 = ok[0], min(ok[-1], env.size - 1)
    return H.coherence(H.cwt_slow(card[i0:i1 + 1]), H.cwt_slow(env[i0:i1 + 1]))


def task_block(seed, intrinsic):
    """One 300 s block: heart and EEG both respond to trials every 25 s; optional shared intrinsic process."""
    rng = np.random.default_rng(seed)
    nb = int(300 * FS)
    ev_heart = S.trial_locked(nb, FS, 25.0, rng)
    ev_brain = S.trial_locked(nb, FS, 25.0, rng)
    shared = S.band_noise(nb, FS, (0.07, 0.13), rng) if intrinsic else 0.0
    beats = S.ipfm_beats(0.8 * ev_heart + 0.6 * S.band_noise(nb, FS, (0.03, 0.45), rng) + 0.7 * shared, FS)
    eeg = S.am_eeg(-0.8 * ev_brain + 0.6 * S.band_noise(nb, FS, (0.03, 0.45), rng) + 0.7 * shared, FS, rng)
    _, env = H.band_envelope(eeg, FS, (8, 12))
    _, card = H.tachogram(beats, t1=300)
    sl = slice(12, 1180)
    return H.cwt_slow(card[sl]), H.cwt_slow(env[sl])


def same_and_cross(blocks, band, n_shuffles=5):
    same = np.mean([H.band_summary(H.coherence(a, b), band)[0] for a, b in blocks])
    cross = []
    for rep in range(n_shuffles):
        p = emg.random_derangement(len(blocks), np.random.default_rng(rep))
        cross.append(np.mean([H.band_summary(H.coherence(blocks[i][0], blocks[p[i]][1]), band)[0] for i in range(len(blocks))]))
    return float(same), float(np.mean(cross))


class VendoredCode(unittest.TestCase):
    def test_vendored_file_is_unmodified(self):
        p = ROOT / "src" / "hb" / "vendor" / "emgsim_v3_3_0" / "functions.py"
        self.assertEqual(hashlib.sha256(p.read_bytes()).hexdigest(), VENDORED_SHA256)


class Calibration(unittest.TestCase):
    def test_frequency_axis(self):
        x = np.sin(2 * np.pi * 0.25 * np.arange(0, 600, 1 / H.FS_SLOW))
        c = H.cwt_slow(x)
        power = (np.abs(c["W"]) ** 2 / c["scales"][:, None]).mean(axis=1)
        peak = c["freqs_full"][power.argmax()]
        self.assertLess(abs(np.log2(peak / 0.25)), 1.5 * H.DJ)   # within 1.5 voices
        self.assertLess(c["freqs"].min(), 0.02)
        self.assertGreater(c["freqs"].max(), 0.40)

    def test_cone_of_influence_matches_analytic_rule(self):
        for seconds in (60, 120, 232):
            c = H.cwt_slow(np.random.default_rng(0).standard_normal(int(seconds * H.FS_SLOW)))
            coh = H.coherence(c, c)
            for f in (0.04, 0.10, 0.15, 0.40):
                r = int(np.argmin(np.abs(coh["freqs"] - f)))
                numeric = coh["inside"][r].sum() / H.FS_SLOW
                self.assertAlmostEqual(numeric, H.edge_safe_seconds(seconds * H.FS_SLOW, coh["freqs"][r]), delta=0.75)

    def test_sixty_second_record_cannot_support_low_lf(self):
        c = H.cwt_slow(np.random.default_rng(1).standard_normal(240))
        value, frac = H.band_summary(H.coherence(c, c), (0.04, 0.046))
        self.assertTrue(np.isnan(value))
        self.assertEqual(frac, 0.0)

    def test_nan_input_is_refused(self):
        x = np.random.default_rng(2).standard_normal(400)
        x[100] = np.nan
        with self.assertRaises(ValueError):
            H.cwt_slow(x)


class KnownCoupling(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.band = (0.07, 0.13)
        cls.null = np.array([H.band_summary(synthetic_pair(100 + s, False), cls.band)[0] for s in range(20)])

    def test_uncoupled_coherence_is_far_from_zero(self):
        # Smoothed wavelet coherence of independent signals is biased upward (about 0.38 here),
        # so a raw value means nothing without a null distribution.
        self.assertGreater(self.null.mean(), 0.25)
        self.assertLess(self.null.mean(), 0.50)

    def test_coupled_exceeds_null(self):
        for seed in range(3):
            v = H.band_summary(synthetic_pair(seed, True), self.band)[0]
            self.assertGreater(v, self.null.max())
            self.assertGreater(v, 0.55)

    def test_coupling_survives_a_three_second_lag(self):
        v = H.band_summary(synthetic_pair(0, True, lag_s=3.0), self.band)[0]
        self.assertGreater(v, self.null.max())

    def test_beat_density_representation_agrees_with_tachogram(self):
        a = H.band_summary(synthetic_pair(0, True, cardiac="tachogram"), self.band)[0]
        b = H.band_summary(synthetic_pair(0, True, cardiac="density"), self.band)[0]
        self.assertLess(abs(a - b), 0.05)

    def test_coupling_is_band_specific(self):
        coh = synthetic_pair(0, True)
        self.assertGreater(H.band_summary(coh, self.band)[0], H.band_summary(coh, H.BANDS["HF"])[0] + 0.2)


class TaskRhythmConfound(unittest.TestCase):
    def test_trial_locked_responses_create_lf_coherence_without_coupling(self):
        blocks = [task_block(500 + k, intrinsic=False) for k in range(8)]
        same, cross = same_and_cross(blocks, H.BANDS["LF"])
        self.assertGreater(same, 0.50)                # far above the uncoupled level of about 0.35
        self.assertLess(abs(same - cross), 0.06)      # and fully reproduced by pairing different blocks

    def test_cross_block_null_leaves_room_to_detect_intrinsic_coupling(self):
        blocks = [task_block(500 + k, intrinsic=True) for k in range(8)]
        same, cross = same_and_cross(blocks, (0.07, 0.13))
        self.assertGreater(same - cross, 0.04)


if __name__ == "__main__":
    unittest.main()
