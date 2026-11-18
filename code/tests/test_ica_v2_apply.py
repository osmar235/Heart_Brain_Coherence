"""Tests of steps 16b and 16c on synthetic recordings."""
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
from hb import cohort as C  # noqa: E402
from hb import ica_v2 as IV  # noqa: E402
from hb import ica_v2_apply as IA  # noqa: E402
import test_pilot_qc as tq  # noqa: E402
import test_slow_signals as ts  # noqa: E402


class Units(unittest.TestCase):
    def test_components_to_remove(self):
        with tempfile.TemporaryDirectory() as d:
            p = Path(d) / "x.tsv"
            p.write_text("component\tbrain\tmuscle\teye\theart\tline_noise\tchannel_noise\tother\tlabel\tprobability\n"
                         "0\t0.05\t0.0\t0.95\t0\t0\t0\t0\teye\t0.95\n"
                         "1\t0.9\t0.1\t0\t0\t0\t0\t0\tbrain\t0.9\n"
                         "2\t0.25\t0.75\t0\t0\t0\t0\t0\tmuscle\t0.75\n"
                         "3\t0.1\t0.85\t0\t0\t0\t0\t0.05\tmuscle\t0.85\n"
                         "4\t0.1\t0\t0\t0\t0\t0\t0.9\tother\t0.9\n")
            self.assertEqual(IA.components_to_remove(p), [(0, "eye"), (3, "muscle")])

    def test_emg_adjust_removes_lagged_muscle_part(self):
        rng = np.random.default_rng(0)
        n = 4 * 600
        emg = rng.standard_normal(n)
        base = rng.standard_normal(n)
        env = 5 + base + 2.0 * np.r_[np.zeros(4), emg[:-4]]
        blocks = [{"t0": 0.0, "t1": 300.0}, {"t0": 300.0, "t1": 600.0}]
        y = IA.emg_adjust(env, emg, blocks)
        lagged = np.r_[np.zeros(4), emg[:-4]]
        self.assertGreater(abs(np.corrcoef(env, lagged)[0, 1]), 0.8)
        self.assertLess(abs(np.corrcoef(y, lagged)[0, 1]), 0.05)
        self.assertAlmostEqual(y[:1200].mean(), env[:1200].mean(), places=6)

    def test_clean_envelope_matches_v1_rule(self):
        x = np.ones(400)
        x[200] = 100.0
        y, bad = IA.clean_envelope(x)
        self.assertTrue(bad[196:205].all() and not bad[:190].any())
        self.assertAlmostEqual(y[200], 1.0)


class Steps(unittest.TestCase):
    def test_fit_apply_analyse(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            raw = tmp / "raw"
            proj = tmp / "project"
            tq.make_recording(raw, eeg_labels=ts.LABELS, seconds=960, n_blocks=4, trials=10, common=60.0)
            olds = (C.PROJECT_ROOT, IV.OUT, IA.OUT)
            C.PROJECT_ROOT = proj
            IV.OUT = IA.OUT = lambda *p: proj.joinpath(*p)
            try:
                C.run_subject(raw, "sub-999", lambda m: None)
                self.assertTrue(IV.fit_subject((str(raw), "sub-999", "picard"))[1].startswith("ok"))
                sub, status = IA.apply_subject((str(raw), "sub-999", "picard"))
                self.assertTrue(status.startswith("ok"), status)
                z = np.load(proj / "derivatives" / "slow_v2" / "picard" / "sub-999_task-memory_slow.npz")
                frozen = np.load(proj / "derivatives" / "slow" / "sub-999_task-memory_slow.npz")
                for e in ("theta", "alpha", "beta"):
                    for key in (f"{e}_ica", f"{e}_noica", f"emg_{e}"):
                        self.assertEqual(z[key].size, frozen["nn"].size)
                        self.assertTrue(np.all(np.isfinite(z[key])))
                    # without removal the v2 envelope tracks the v1 envelope closely (same bands, different reference)
                    r = np.corrcoef(z[f"{e}_noica"], frozen[e])[0, 1]
                    self.assertGreater(r, 0.5)
                z.close()  # Windows cannot delete the temporary folder while an .npz is open
                frozen.close()
                info = json.loads((proj / "derivatives" / "slow_v2" / "picard" / "sub-999_apply.json").read_text())
                self.assertIn("removed_by_class", info)
                sub, status = IA.analyse_subject((str(raw), "sub-999", 1, "picard"))
                self.assertTrue(status.startswith("ok"), status)
                (proj / "results" / "trial_shuffle").mkdir(parents=True)
                (proj / "results" / "trial_shuffle" / "gate.json").write_text(json.dumps(
                    {"passed": True, "bias_upper_edges": {f"{e} {b}": 0.0 for e in ("theta", "alpha", "beta") for b in ("LF", "HF")}}))
                msgs = []
                IA.summarize_analysis("picard", msgs.append)
                self.assertTrue(any(m.startswith("realism v2") for m in msgs))
                self.assertTrue((proj / "results" / "ica_v2" / "shuffle_calibrated_H1_picard.tsv").exists())
                self.assertEqual(IA.apply_subject((str(raw), "sub-999", "picard"))[1], "skipped (already done)")
            finally:
                C.PROJECT_ROOT, IV.OUT, IA.OUT = olds


if __name__ == "__main__":
    unittest.main()
