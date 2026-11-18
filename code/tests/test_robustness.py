"""Tests of the four robustness checks (docs/ROBUSTNESS_PLAN_v1.md) on synthetic data with known truth."""
import csv
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
from hb import cohort as C  # noqa: E402
from hb import hbcoh as H  # noqa: E402
from hb import null_calibration as NC  # noqa: E402
from hb import robust_coherence as RC  # noqa: E402
from hb import robust_reprocess as RR  # noqa: E402
from hb import synth as S  # noqa: E402
import test_pilot_coherence as tp  # noqa: E402
import test_pilot_qc as tq  # noqa: E402
import test_slow_signals as ts  # noqa: E402


def artifact_pulse(fs):
    t = np.arange(int(0.6 * fs)) / fs - 0.1
    return 12 * np.exp(-0.5 * (t / 0.012) ** 2) - 4 * np.exp(-0.5 * ((t - 0.03) / 0.012) ** 2) + 3 * np.exp(-0.5 * ((t - 0.25) / 0.05) ** 2)


class CardiacFieldArtifact(unittest.TestCase):
    @staticmethod
    def one(seed, amp=2.0):
        fs, T = 200.0, 600
        n = int(T * fs)
        rng = np.random.default_rng(seed)
        mod = S.band_noise(n, fs, (0.07, 0.13), rng)
        beats = S.ipfm_beats(mod, fs, hr_mean=1.2, depth=0.12)
        beats = beats[(beats > 1) & (beats < T - 1)]
        t = np.arange(n) / fs
        eeg = 1.5 * np.sin(2 * np.pi * 10 * t + rng.uniform(0, 6)) + 2.0 * rng.standard_normal(n)
        pulse = amp * artifact_pulse(fs)
        dirty = eeg.copy()
        for bt in beats:
            i = int(round((bt - 0.1) * fs))
            dirty[i:i + pulse.size] += pulse
        _, nn = H.tachogram(beats, t1=T)
        clean, p0, p1 = RR.subtract_cfa(dirty, np.round(beats * fs).astype(int), fs)

        def lf(x):
            _, env = H.band_envelope(x, fs, (8, 12))
            ok = np.where(np.isfinite(nn))[0]
            sl = slice(ok[0] + 8, min(ok[-1], env.size) - 8)
            return H.band_summary(H.coherence(H.cwt_slow(nn[sl]), H.cwt_slow(env[sl])), (0.07, 0.13))[0]
        return lf(eeg), lf(dirty), lf(clean), p0, p1

    def test_artifact_creates_envelope_heart_coherence_and_subtraction_removes_it(self):
        """The mechanism behind check 1: artifact pulses whose density follows heart rate make the alpha
        envelope coherent with the heart series; template subtraction restores the artifact free level."""
        r = np.array([self.one(s) for s in range(3)])
        base, dirty, fixed = r[:, 0].mean(), r[:, 1].mean(), r[:, 2].mean()
        self.assertGreater(dirty, base + 0.1)
        self.assertLess(abs(fixed - base), 0.03)
        self.assertTrue((r[:, 4] < 0.15 * r[:, 3]).all())

    def test_clean_signal_is_left_almost_unchanged(self):
        fs = 200.0
        rng = np.random.default_rng(1)
        x = rng.standard_normal(int(300 * fs))
        beats = np.arange(1.0, 299.0, 0.8)
        y, _, _ = RR.subtract_cfa(x, np.round(beats * fs).astype(int), fs)
        self.assertLess(np.std(y - x) / np.std(x), 0.15)


class Respiration(unittest.TestCase):
    def test_edr_follows_r_amplitude_modulation(self):
        fs, T = 1000.0, 300
        rng = np.random.default_rng(2)
        beats = np.cumsum(np.full(int(T / 0.85) - 2, 0.85)) + 0.5
        resp = np.sin(2 * np.pi * 0.25 * beats)
        n = int(T * fs)
        ecg = 200 * rng.standard_normal(n)
        tt = np.arange(-40, 41) / fs
        for bt, r in zip(beats, resp):
            i = int(round(bt * fs))
            ecg[i - 40:i + 41] += (10000 * (1 + 0.15 * r)) * np.exp(-0.5 * (tt / 0.01) ** 2)
        edr, amp = RR.edr_series(ecg, fs, beats, int(T * 4))
        grid = np.arange(edr.size) / 4.0
        ref = np.sin(2 * np.pi * 0.25 * grid)
        sl = slice(40, -40)
        self.assertGreater(abs(np.corrcoef(edr[sl], ref[sl])[0, 1]), 0.9)

    def test_covariate_removal_cancels_shared_respiratory_drive(self):
        res = {}
        for use in (False, True):
            heart, brain, blocks, n = tp.session(4, False)
            rng = np.random.default_rng(9)
            resp = S.band_noise(n, 4.0, (0.2, 0.3), rng)
            heart = heart + 0.8 * resp
            brain = brain + 0.8 * np.roll(resp, 2)
            rows = C.coherence_tables("s", heart, {"alpha": brain}, {"alpha": np.zeros(n, bool)}, blocks, [0.0] * len(blocks),
                                      covariate=resp if use else None)
            res[use] = np.mean([r["excess_full_HF"] for r in rows if not r["excluded"]])
        self.assertGreater(res[False], 0.08)
        self.assertLess(abs(res[True]), 0.03)


class EndToEnd(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        tmp = Path(cls.tmp.name)
        cls.raw = tmp / "raw"
        tq.make_recording(cls.raw, eeg_labels=ts.LABELS, seconds=960, n_blocks=4, trials=10)
        cls.proj = tmp / "project"
        cls.mods = (C, RR, RC, NC)
        cls.old = [m.PROJECT_ROOT for m in cls.mods]
        for m in cls.mods:
            m.PROJECT_ROOT = cls.proj
        cls.status = C.run_subject(cls.raw, "sub-999", lambda m: None)

    @classmethod
    def tearDownClass(cls):
        for m, o in zip(cls.mods, cls.old):
            m.PROJECT_ROOT = o
        cls.tmp.cleanup()

    def test_1_frozen_run_ok(self):
        self.assertTrue(self.status.startswith("ok"), self.status)

    def test_2_reprocess_and_robust_coherence(self):
        r = RR.process(self.raw, "sub-999", lambda m: None)
        self.assertEqual(r["status"], "ok")
        z = np.load(self.proj / "derivatives" / "slow_cfa" / "sub-999_task-memory_slow.npz")
        f = np.load(self.proj / "derivatives" / "slow" / "sub-999_task-memory_slow.npz")
        np.testing.assert_array_equal(z["nn"], f["nn"])
        self.assertEqual(z["alpha"].shape, f["alpha"].shape)
        self.assertEqual(z["edr"].shape, f["nn"].shape)
        status, rsa = RC.run_subject(self.raw, "sub-999")
        self.assertEqual(status, "ok")
        for v in RC.ROBUST:
            with open(self.proj / "derivatives" / "robust" / f"sub-999_{v}_blocks.tsv") as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            self.assertEqual(len(rows), 12)
            self.assertTrue(any(not r["excluded"] for r in rows))
        self.assertEqual(RC.run_subject(self.raw, "sub-999")[0], "skipped (already done)")

    def test_3_null_calibration_subject(self):
        sub, status = NC.run_subject((str(self.raw), "sub-999", 2, 1))
        self.assertTrue(status.startswith("ok"), status)
        rows = NC.load_rows(self.proj / "derivatives" / "null_calibration" / "sub-999.tsv")
        self.assertEqual({r["scenario"] for r in rows}, {"none", "equal"})
        self.assertEqual(len([r for r in rows if r["scenario"] == "none"]), 2 * 3 * 2 * 2)
        none = [r["mean_excess"] for r in rows if r["scenario"] == "none" and r["variant"] == "full"]
        self.assertLess(abs(np.nanmean(none)), 0.05)


class SummaryInputs(unittest.TestCase):
    def test_observed_table_without_k_column_loads(self):
        """Regression test for the 2026-11-11 crash: the frozen participant table has no 'k' column."""
        with tempfile.TemporaryDirectory() as tmp:
            f = Path(tmp) / "participant_level.tsv"
            f.write_text("sub\teeg\tvariant\tband\tmean_excess\tin_H1\tin_H2\tn_blocks\tdiff_span1_minus_span\texcess_span\texcess_span1\n"
                         "sub-001\talpha\tfull\tLF\t0.01\t1\t0\t12\t\t\t\n")
            r = NC.load_rows(f)[0]
            self.assertEqual(r["k"], 0)
            self.assertAlmostEqual(r["mean_excess"], 0.01)


class ReviewSteps(unittest.TestCase):
    """Steps 13 and 14 (review response) on the synthetic participant of EndToEnd."""
    @classmethod
    def setUpClass(cls):
        EndToEnd.setUpClass()
        from hb import joint_calibration as JC, null_validation as NV
        cls.JC, cls.NV = JC, NV
        cls.old2 = (JC.PROJECT_ROOT, NV.PROJECT_ROOT)
        JC.PROJECT_ROOT = NV.PROJECT_ROOT = EndToEnd.proj
        RR.process(EndToEnd.raw, "sub-999", lambda m: None)

    @classmethod
    def tearDownClass(cls):
        cls.JC.PROJECT_ROOT, cls.NV.PROJECT_ROOT = cls.old2
        EndToEnd.tearDownClass()

    def test_shift_surrogate_keeps_heart_and_breaks_alignment(self):
        rng = np.random.default_rng(0)
        heart, brain, blocks, n = tp.session(1, True)
        b10 = [b for b in blocks if b["n_trials"] == 10]
        parts = self.JC.prepare_parts(heart, brain, b10)
        h, b = self.JC.shift_session(parts, b10, rng)
        z = self.JC.synthesize(n, b10, parts["heart_model"]) + parts["heart_resid"]
        np.testing.assert_allclose(h, z)
        real = C.coherence_tables("s", heart, {"alpha": brain}, {"alpha": np.zeros(n, bool)}, blocks, [0.0] * len(blocks))
        sh = C.coherence_tables("s", h, {"alpha": b}, {"alpha": np.zeros(n, bool)}, blocks, [0.0] * len(blocks))
        m = lambda rows: np.mean([r["same_full_LF"] for r in rows if not r["excluded"]])
        self.assertLess(m(sh), m(real) - 0.1)

    def test_joint_step_runs(self):
        sub, status = self.JC.run_subject((str(EndToEnd.raw), "sub-999", 2))
        self.assertTrue(status.startswith("ok"), status)
        rows = self.JC.load_rows(EndToEnd.proj / "derivatives" / "joint_calibration" / "sub-999.tsv")
        self.assertEqual({r["scenario"] for r in rows}, {"observed", "phase", "shift"})

    def test_null_validation_step_runs_and_is_unbiased_enough(self):
        sub, status = self.NV.run_subject((str(EndToEnd.raw), "sub-999", 2))
        self.assertTrue(status.startswith("ok"), status)
        with open(EndToEnd.proj / "derivatives" / "null_validation" / "sub-999.tsv") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        obs = [float(r["mean_excess"]) for r in rows if r["kind"] == "observed" and r["scenario"] == "independent" and r["mean_excess"] not in ("", "nan")]
        self.assertTrue(obs)
        self.assertLess(abs(np.mean(obs)), 0.1)

    def test_error_rate_monte_carlo(self):
        rng = np.random.default_rng(3)
        D = rng.standard_normal((400, 6)) * 0.01
        mu, half, res = self.NV.error_rates(D - D.mean(axis=0), 122, rng, reps=4000)
        self.assertAlmostEqual(float(res["at_zero_bias"][0].mean()), 0.05, delta=0.015)
        self.assertLessEqual(res["at_zero_bias"][1], 0.07)
        mu2, _, res2 = self.NV.error_rates(D - D.mean(axis=0) + 0.003, 122, rng, reps=4000)
        self.assertGreater(res2["at_estimated_bias"][1], 0.5)


if __name__ == "__main__":
    unittest.main()
