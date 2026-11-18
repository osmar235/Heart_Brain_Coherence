"""Test of step 06 on a synthetic recording with a known trial locked alpha response."""
import sys
import tempfile
import unittest
import warnings
from pathlib import Path

import numpy as np
import scipy.io as sio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
warnings.simplefilter("ignore")
from hb import beat_correction, bidsio, slow_signals  # noqa: E402
import test_pilot_qc as tq  # noqa: E402

LABELS = ("Fz", "AFz", "F1", "F2", "P1", "P2", "P3", "P4", "C3", "C1", "CP3", "Oz", "T7", "T8")


class SlowSignals(unittest.TestCase):
    def test_step_runs_and_recovers_trial_locked_alpha(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            tq.make_recording(root, eeg_labels=LABELS)
            p = bidsio.rec_paths(root, "sub-999", "memory")
            m = sio.loadmat(str(p["set"]))
            parsed = bidsio.parse_memory_events(bidsio.read_tsv(p["events"]), 200.0)
            t = np.arange(m["data"].shape[1]) / 1000.0
            gain = np.ones(t.size)
            for tr in parsed["trials"]:                       # alpha is suppressed from trial start to recall prompt
                gain[(t >= tr["t_start"]) & (t < tr["t_probe"])] = 0.2
            rng = np.random.default_rng(1)
            for i in (4, 5, 6, 7):                             # P1..P4
                m["data"][i] += (60 * gain * np.sin(2 * np.pi * 11 * t + rng.uniform(0, 6))).astype(np.float32)
            sio.savemat(str(p["set"]), {k: v for k, v in m.items() if not k.startswith("__")})
            beat_correction.process(root, "sub-999", "memory", root / "bc", root / "beats")
            rows = slow_signals.process(root, "sub-999", "memory", root, root / "slow", root / "beats", lambda s: None, min_trials=3)
            self.assertEqual(len(rows), 1)
            r = rows[0]
            self.assertEqual(r["length"], 4)
            self.assertEqual(r["alpha_channels"], "P1,P2,P3,P4")
            self.assertGreater(r["alpha_evoked_var_fraction"], 0.5)
            self.assertLess(r["theta_evoked_var_fraction"], r["alpha_evoked_var_fraction"])
            z = np.load(root / "slow" / "sub-999_task-memory_slow.npz")
            self.assertEqual(z["nn"].shape, z["alpha"].shape)
            z.close()  # Windows cannot delete the temporary folder while the .npz is open
            prof = np.load(root / "sub-999_trial_profiles.npz")
            tt, a = prof["L4_t"], prof["L4_alpha_mean"]
            self.assertLess(a[(tt > 3) & (tt < 10)].mean(), a[(tt > 15) & (tt < 22)].mean() - 1.0)
            prof.close()
            self.assertTrue((root / "sub-999_trial_profiles.png").exists())


if __name__ == "__main__":
    unittest.main()
