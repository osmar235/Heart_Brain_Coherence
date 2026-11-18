"""Tests of step 15 (docs/STEP15_PLAN_v1.md) on synthetic data with known truth."""
import csv
import json
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
from hb import null_validation as NV  # noqa: E402
from hb import trial_shuffle as TS  # noqa: E402
import test_pilot_coherence as tp  # noqa: E402
from test_robustness import EndToEnd  # noqa: E402

FS = 4.0


def blocks10(seed=0):
    _, _, blocks, n = tp.session(seed, False)
    return [b for b in blocks if b["n_trials"] == 10], n


class Model(unittest.TestCase):
    def test_default_reproduces_step14_independent(self):
        b10, n = blocks10()
        a = NV.simulate(b10, n, "independent", np.random.default_rng(5))
        b = TS.simulate(b10, n, np.random.default_rng(5))
        for k in a:
            np.testing.assert_allclose(a[k], b[k])

    def test_variants_change_only_their_part(self):
        b10, n = blocks10()
        base = TS.simulate(b10, n, np.random.default_rng(1))
        stat = TS.simulate(b10, n, np.random.default_rng(1), modulate=False)
        self.assertFalse(np.allclose(base["alpha"], stat["alpha"]))
        same = TS.simulate(b10, n, np.random.default_rng(1), recall_scale=True)
        np.testing.assert_allclose(base["beta"], same["beta"])

    def test_realism_metrics(self):
        b10, n = blocks10()
        mo = TS.simulate(b10, n, np.random.default_rng(2), share={k: 0.6 for k in TS.SIGNALS}, **TS.VARIANTS["mean_only"])
        lo = TS.simulate(b10, n, np.random.default_rng(2), share={k: 0.05 for k in TS.SIGNALS})
        r2_hi, _ = TS.realism_metrics(mo["alpha"], b10)
        r2_lo, vmi = TS.realism_metrics(lo["alpha"], b10)
        self.assertGreater(r2_hi, r2_lo)
        self.assertTrue(0 <= r2_lo < 0.5 and vmi > 0)


class Shuffle(unittest.TestCase):
    def test_derangement(self):
        rng = np.random.default_rng(0)
        for _ in range(50):
            p = TS.derangement(10, rng)
            self.assertEqual(sorted(p), list(range(10)))
            self.assertFalse(np.any(p == np.arange(10)))

    def test_trial_gains_are_permuted_and_outside_kept(self):
        b10, n = blocks10()
        t = np.arange(n) / FS
        x = np.zeros(n)
        gains = {}
        for bi, b in enumerate(b10):
            for k in range(10):
                g = 1.0 + bi * 10 + k
                gains[(bi, k)] = g
                x[(t >= b["start"][k]) & (t < b["probe"][k])] = g
        y = TS.trial_shuffle(x, b10, np.random.default_rng(3))
        for bi, b in enumerate(b10):
            got = []
            for k in range(10):
                m = (t >= b["start"][k] + 1) & (t < b["probe"][k] - 1)
                got.append(np.median(y[m]))
            self.assertEqual(sorted(got), sorted(gains[(bi, k)] for k in range(10)))
            self.assertTrue(all(g != gains[(bi, k)] for k, g in enumerate(got)))
            pre = (t >= b["t0"]) & (t < b["start"][0])
            np.testing.assert_allclose(y[pre], x[pre])

    def test_shuffle_removes_ongoing_coupling(self):
        heart, brain, blocks, n = tp.session(1, True)
        b10 = [b for b in blocks if b["n_trials"] == 10]
        z = np.zeros(n, bool)
        m = lambda h, b: np.mean([r["same_full_LF"] for r in C.coherence_tables("s", h, {"alpha": b}, {"alpha": z}, blocks, [0.0] * len(blocks)) if not r["excluded"]])
        real = m(heart, brain)
        sh = np.mean([m(heart, TS.trial_shuffle(brain, b10, np.random.default_rng(k))) for k in range(2)])
        self.assertLess(sh, real - 0.1)


def fake_rows(cohort, null, nsub, bias, rng, k=4):
    rows = []
    for s in range(nsub):
        sub = f"sub-{s:03d}"
        for e in TS.EEG:
            for v in ("full", "encret"):
                for band in ("LF", "HF"):
                    o = bias + 0.01 * rng.standard_normal()
                    rows.append({"cohort": cohort, "kind": "observed", "null": "", "k": "0", "sub": sub, "eeg": e, "variant": v,
                                 "band": band, "mean_excess": str(o), "in_H1": "1"})
                    for kk in range(k):
                        rows.append({"cohort": cohort, "kind": "surrogate", "null": null, "k": str(kk), "sub": sub, "eeg": e, "variant": v,
                                     "band": band, "mean_excess": str(0.01 * rng.standard_normal()), "in_H1": "1"})
    return rows


class Gate(unittest.TestCase):
    def test_summaries_and_gate_on_known_bias(self):
        import tempfile
        rng = np.random.default_rng(0)
        with tempfile.TemporaryDirectory() as d:
            old = TS.OUT
            TS.OUT = lambda *p: Path(d).joinpath(*p)
            try:
                rows = fake_rows("independent", "shuffle", 120, 0.0, rng) + fake_rows("matched", "shuffle", 120, 0.0, rng)
                self.assertTrue(TS.gate(TS.synthetic_summary(rows, "validate", [("independent", "shuffle", None), ("matched", "shuffle", None)], lambda m: None)))
                rows = fake_rows("independent", "shuffle", 120, 0.006, rng) + fake_rows("matched", "shuffle", 120, 0.0, rng)
                self.assertFalse(TS.gate(TS.synthetic_summary(rows, "validate", [("independent", "shuffle", None), ("matched", "shuffle", None)], lambda m: None)))
                real = fake_rows("real", "shuffle", 60, 0.012, rng, k=16)
                (Path(d) / "results" / "cohort_stats").mkdir(parents=True)
                with open(Path(d) / "results" / "cohort_stats" / "participant_level.tsv", "w", newline="") as fh:
                    w = csv.DictWriter(fh, ["sub", "eeg", "variant", "band", "mean_excess"], delimiter="\t")
                    w.writeheader()
                    for r in real:
                        if r["kind"] == "observed":
                            w.writerow({c: r[c] for c in ["sub", "eeg", "variant", "band", "mean_excess"]})
                msgs = []
                TS.real_summary(real, msgs.append, {f"{e} {b}": 0.0 for e, b in TS.CELLS})
                self.assertIn("max abs difference 0 ", msgs[-1])
                with open(Path(d) / "results" / "trial_shuffle" / "real_shuffle_calibrated_H1.tsv") as fh:
                    res = list(csv.DictReader(fh, delimiter="\t"))
                self.assertEqual(len(res), 24)
                self.assertTrue(all(float(r["p_holm"]) < 0.05 for r in res if r["k_used"] == "16" and r["variant"] == "full"))
                self.assertTrue(all(r["lower_bound_above_edge"] == "1" for r in res if r["k_used"] == "16" and r["variant"] == "full"))
            finally:
                TS.OUT = old

    def test_gate(self):
        row = lambda c, p, j: {"cohort": c, "null": "shuffle", "variant": "full", "eeg": TS.CELLS[j][0], "band": TS.CELLS[j][1], "p_bias_holm": p, "ci95_high": 0.001 * j}
        ok = [row(c, 0.5, j) for c in ("independent", "matched") for j in range(6)]
        self.assertTrue(TS.gate(ok))
        self.assertFalse(TS.gate(ok[:-1] + [row("matched", 0.01, 5)]))
        self.assertFalse(TS.gate(ok[:6]))
        self.assertAlmostEqual(TS.bias_edges(ok)[TS.CELLS[5]], 0.005)


class Steps(unittest.TestCase):
    """15a, 15b and 15c run end to end on the synthetic participant of test_robustness.EndToEnd."""
    @classmethod
    def setUpClass(cls):
        EndToEnd.setUpClass()
        cls.old = TS.PROJECT_ROOT
        TS.PROJECT_ROOT = EndToEnd.proj
        TS.OUT = lambda *p: EndToEnd.proj.joinpath(*p)

    @classmethod
    def tearDownClass(cls):
        TS.PROJECT_ROOT = cls.old
        TS.OUT = lambda *p: TS.PROJECT_ROOT.joinpath(*p)
        EndToEnd.tearDownClass()

    def test_all_steps(self):
        raw = str(EndToEnd.raw)
        sub, rows = TS.realism_subject((raw, "sub-999"))
        self.assertEqual(len(rows), 8)
        self.assertTrue(all(0 <= r["R2"] <= 1 for r in rows))
        TS.summarize_realism(rows, lambda m: None)
        matched = json.loads(TS.OUT("results", "trial_shuffle", "matched_shares.json").read_text())
        self.assertEqual(set(matched), set(TS.SIGNALS))
        self.assertTrue(TS.decompose_subject((raw, "sub-999", 1))[1].startswith("ok"))
        self.assertTrue(TS.validate_subject((raw, "sub-999", 1, matched))[1].startswith("ok"))
        self.assertTrue(TS.real_subject((raw, "sub-999", 2))[1].startswith("ok"))
        rows = TS.load("validate")
        self.assertEqual({(r["cohort"], r["kind"], r["null"]) for r in rows},
                         {("independent", "observed", ""), ("independent", "surrogate", "shuffle"),
                          ("matched", "observed", ""), ("matched", "surrogate", "shuffle"), ("matched", "surrogate", "phase")})
        self.assertEqual({(r["kind"], r["null"]) for r in TS.load("real")}, {("observed", ""), ("surrogate", "shuffle")})
        with open(TS.OUT("derivatives", "trial_shuffle", "real", "sub-999.tsv")) as fh:
            real = list(csv.DictReader(fh, delimiter="\t"))
        from hb.cohort_stats import participant_table
        with open(EndToEnd.proj / "derivatives" / "cohort_coherence" / "sub-999_blocks.tsv") as fh:
            frozen_rows = list(csv.DictReader(fh, delimiter="\t"))
        frozen = {(r["eeg"], r["variant"], r["band"]): r["mean_excess"] for r in participant_table({"sub-999": frozen_rows})}
        obs = [r for r in real if r["kind"] == "observed" and r["mean_excess"] not in ("", "nan")]
        self.assertTrue(obs)
        for r in obs:
            self.assertAlmostEqual(float(r["mean_excess"]), frozen[(r["eeg"], r["variant"], r["band"])], places=5)


if __name__ == "__main__":
    unittest.main()
