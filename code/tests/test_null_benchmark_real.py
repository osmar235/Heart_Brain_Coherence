"""Tests of the surrogate construction used by the rebuilt null benchmark."""
import sys
import unittest
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
warnings.simplefilter("ignore")
from hb import null_benchmark_real as NB  # noqa: E402
from hb import synth as S  # noqa: E402

FS = NB.FS


def fake_timing(rng, n_blocks=6, length=6):
    blocks, t = [], 10.0
    for _ in range(n_blocks):
        t0 = t
        start, probe = [], []
        for _ in range(10):
            start.append(t)
            probe.append(t + 2.5 + 2 * length + 2.0)
            t = probe[-1] + rng.uniform(4, 12)
        blocks.append(dict(t0=t0, t1=t, length=length, start=np.array(start), probe=np.array(probe),
                           end=np.array(start[1:] + [t])))
        t += 25.0
    return blocks, t + 10.0


def fake_signal(blocks, dur, rng, shared=None, amp=1.0):
    n = int(dur * FS)
    t = np.arange(n) / FS
    x = S.band_noise(n, FS, (0.01, 1.0), rng)
    for b in blocks:
        for k in range(10):
            x += amp * (1 + 0.3 * rng.standard_normal()) * ((t >= b["start"][k] + 2.5) & (t < b["probe"][k]))
    return x if shared is None else x + shared


class Surrogates(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        rng = np.random.default_rng(0)
        cls.blocks, cls.dur = fake_timing(rng)
        n = int(cls.dur * FS)
        cls.shared = 1.2 * S.band_noise(n, FS, (0.07, 0.13), rng)
        cls.heart = fake_signal(cls.blocks, cls.dur, rng, cls.shared)
        cls.brain = fake_signal(cls.blocks, cls.dur, rng, cls.shared)
        cls.parts = NB.prepare_parts(cls.heart, cls.brain, cls.blocks)

    def test_decomposition_adds_back_to_the_signal(self):
        z = NB.zscore_inblock(self.heart[:self.parts["n"]], self.blocks)
        rebuilt = NB.synthesize(self.parts["n"], self.blocks, self.parts["heart_model"]) + self.parts["heart_resid"]
        self.assertLess(np.abs(rebuilt - z).max(), 1e-9)
        self.assertGreater(self.parts["heart_evoked_var_fraction"], 0.1)

    def test_surrogate_destroys_real_coupling_but_keeps_spectrum(self):
        band = (0.07, 0.13)
        b0 = self.blocks[0]

        def coh(h, b):
            c = NB.H.coherence(NB.H.cwt_slow(NB.seg(h, b0)), NB.H.cwt_slow(NB.seg(b, b0)))
            return NB.H.band_summary(c, band)[0]
        zh = NB.zscore_inblock(self.heart[:self.parts["n"]], self.blocks)
        zb = NB.zscore_inblock(self.brain[:self.parts["n"]], self.blocks)
        real = coh(zh, zb)
        sur = [coh(*NB.surrogate_session(self.parts, self.blocks, "task_only", np.random.default_rng(s))) for s in range(6)]
        self.assertGreater(real, 0.75)                       # the fake "real" data are strongly coupled
        self.assertLess(np.mean(sur), real - 0.2)            # the surrogates are not
        h, _ = NB.surrogate_session(self.parts, self.blocks, "task_only", np.random.default_rng(1))
        s = NB.block_slices(self.blocks, self.parts["n"])[0]
        self.assertAlmostEqual(float(h[s].var()), float(zh[s].var()), delta=0.25 * float(zh[s].var()))

    def test_added_coupling_raises_excess_on_average(self):
        ex = {"task_only": [], "ongoing_HF": []}
        for seed in range(8):
            for scenario in ex:
                rng = np.random.default_rng([5, seed])
                h, b = NB.surrogate_session(self.parts, self.blocks, scenario, rng)
                ex[scenario].append(NB.evaluate(h, b, self.blocks, rng)["excess_resid_clock_HF"])
        self.assertGreater(np.mean(ex["ongoing_HF"]), np.mean(ex["task_only"]) + 0.01)
        self.assertLess(abs(np.mean(ex["task_only"])), 0.03)

if __name__ == "__main__":
    unittest.main()
