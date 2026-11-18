"""Tests of step 18 (docs/STEP18_PLAN_v1.md) on synthetic data with known truth."""
import sys
import unittest
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "tests"))
warnings.simplefilter("ignore")
from hb import symmetric_tests as ST  # noqa: E402
from hb import trial_shuffle as TS  # noqa: E402
import test_null_benchmark_real as tn  # noqa: E402

FS = 4.0


def session(seed, n_blocks=6, length=6):
    rng = np.random.default_rng(seed)
    blocks, dur = tn.fake_timing(rng, n_blocks=n_blocks, length=length)
    for i, b in enumerate(blocks):
        b.update(index=i + 1, n_trials=10, boundary=False)
    return blocks, int(dur * FS)


def trial_code(blocks, n):
    """Signal whose value inside trial k of any block is k (encoding) or k + 0.5 (recall), -1 outside trials."""
    t = np.arange(n) / FS
    x = np.full(n, -1.0)
    for b in blocks:
        for k in range(10):
            x[(t >= b["start"][k]) & (t < b["probe"][k])] = k
            x[(t >= b["probe"][k]) & (t < b["end"][k])] = k + 0.5
    return x


class Warp(unittest.TestCase):
    def test_trial_k_receives_partner_trial_k(self):
        bs, ns = session(1)
        bp, npart = session(2)
        y, exact = ST.cross_surrogate(np.zeros(ns), bs, trial_code(bp, npart), bp, np.random.default_rng(0))
        self.assertEqual(exact, len(bs))
        code_self = trial_code(bs, ns)
        t = np.arange(ns) / FS
        for b in bs:
            for k in range(10):
                mid_enc = (t > b["start"][k] + 1.0) & (t < b["probe"][k] - 1.0)
                mid_rec = (t > b["probe"][k] + 1.0) & (t < b["end"][k] - 1.0)
                np.testing.assert_allclose(y[mid_enc], code_self[mid_enc])
                if mid_rec.any():
                    np.testing.assert_allclose(y[mid_rec], code_self[mid_rec])

    def test_outside_blocks_unchanged(self):
        bs, ns = session(3)
        bp, npart = session(4)
        x = np.random.default_rng(5).standard_normal(ns)
        y, _ = ST.cross_surrogate(x, bs, np.random.default_rng(6).standard_normal(npart), bp, np.random.default_rng(0))
        inside = np.zeros(ns, bool)
        for b in bs:
            inside[int(round(b["t0"] * FS)):int(round(b["t1"] * FS))] = True
        np.testing.assert_array_equal(y[~inside], x[~inside])
        self.assertFalse(np.allclose(y[inside], x[inside]))

    def test_same_length_preferred(self):
        bs, _ = session(1, length=6)
        b_other, _ = session(2, length=8)
        b_same, _ = session(3, length=6)
        pool = b_other[:3] + b_same[:2]
        pick, exact = ST.choose_partner_blocks(bs, pool, np.random.default_rng(0))
        self.assertEqual(exact, len(bs))
        self.assertTrue(all(pool[j]["length"] == 6 for j in pick))
        pick, exact = ST.choose_partner_blocks(bs, b_other, np.random.default_rng(0))
        self.assertEqual(exact, 0)


class Plant(unittest.TestCase):
    def test_shared_share(self):
        bs, n = session(7, n_blocks=8)
        rng = np.random.default_rng(1)
        sig = {"heart": rng.standard_normal(n), "beta": rng.standard_normal(n)}
        for c in (0.0, 0.1, 0.4):
            p = ST.plant_all(sig, c, bs, np.random.default_rng(2))
            r = np.corrcoef(p["heart"], p["beta"])[0, 1]
            self.assertAlmostEqual(r, c, delta=0.05)


class EndToEnd(unittest.TestCase):
    def test_planted_coupling_exceeds_cross_surrogates_and_none_does_not(self):
        res = {}
        for c in (0.0, 0.15):
            vals = []
            for seed in range(3):
                bs, n = session(10 + seed, n_blocks=8)
                rng = np.random.default_rng(seed)
                base = {"heart": tn.fake_signal(bs, n / FS, rng), "beta": tn.fake_signal(bs, n / FS, rng)}
                sig = ST.plant_all(base, c, bs, np.random.default_rng(100 + seed))
                bp, npart = session(50 + seed, n_blocks=8)
                partner = ST.plant_all({"heart": tn.fake_signal(bp, npart / FS, rng), "beta": tn.fake_signal(bp, npart / FS, rng)},
                                       c, bp, np.random.default_rng(200 + seed))
                zeros = np.zeros(n, bool)
                obs = [r for r in ST.table_rows("sub-999", sig["heart"], "beta", sig["beta"], zeros, bs, [0.0] * len(bs))
                       if r["band"] == "LF" and r["variant"] == "full"][0]["mean_excess"]
                sur = []
                for kk in range(3):
                    b, _ = ST.cross_surrogate(sig["beta"], bs, partner["beta"], bp, np.random.default_rng(kk))
                    sur.append([r for r in ST.table_rows("sub-999", sig["heart"], "beta", b, zeros, bs, [0.0] * len(bs))
                                if r["band"] == "LF" and r["variant"] == "full"][0]["mean_excess"])
                vals.append(obs - np.mean(sur))
            res[c] = np.mean(vals)
        self.assertGreater(res[0.15], 0.03)
        self.assertLess(abs(res[0.0]), 0.03)
        self.assertGreater(res[0.15], res[0.0] + 0.03)


if __name__ == "__main__":
    unittest.main()


class Summaries(unittest.TestCase):
    def fake_rows(self, cohort, nulls, effect, k, subs, rng, bias=0.0):
        rows = []
        for s in subs:
            for e, b in ST.CELLS:
                true = effect if b == "LF" else 0.0
                rows.append({"cohort": cohort, "kind": "observed", "null": "", "k": "0", "sub": s, "eeg": e, "variant": "full",
                             "band": b, "mean_excess": str(true + bias + 0.01 * rng.standard_normal()), "in_H1": "1"})
                for null in nulls:
                    for kk in range(k):
                        rows.append({"cohort": cohort, "kind": "surrogate", "null": null, "k": str(kk), "sub": s, "eeg": e, "variant": "full",
                                     "band": b, "mean_excess": str(bias + 0.01 * rng.standard_normal()), "in_H1": ""})
        return rows

    def test_gate_detection_and_claim(self):
        import json
        import tempfile
        tmp = Path(tempfile.mkdtemp())
        old = ST.OUT
        ST.OUT = lambda *p: tmp.joinpath(*p)
        try:
            (tmp / "results" / "trial_shuffle").mkdir(parents=True)
            (tmp / "results" / "trial_shuffle" / "gate.json").write_text(json.dumps(
                {"passed": True, "bias_upper_edges": {f"{e} {b}": 0.003 for e, b in ST.CELLS}}))
            rng = np.random.default_rng(0)
            subs = [f"sub-{i:03d}" for i in range(40)]
            rows = self.fake_rows("independent", ["cross"], 0.0, 8, subs, rng) + self.fake_rows("matched", ["cross"], 0.0, 8, subs, rng)
            lv = {"low": 0.01, "mid": 0.02, "high": 0.04}
            for name, eff in (("low", 0.0), ("mid", 0.02), ("high", 0.05)):
                rows += self.fake_rows(f"plant_{name}", ["shuffle", "cross"], eff, 8, subs, rng)
            self.assertTrue(ST.summarize_synthetic(rows, lv, lambda m: None))
            import csv
            tab = list(csv.DictReader(open(tmp / "results" / "symmetric" / "synthetic_cells.tsv"), delimiter="\t"))
            det = {(r["cohort"], r["null"], r["band"]): r["detected_by_rule"] for r in tab if r["cohort"].startswith("plant") and r["eeg"] == "beta"}
            self.assertEqual(det[("plant_high", "shuffle", "LF")], "1")
            self.assertEqual(det[("plant_low", "shuffle", "LF")], "0")
            self.assertEqual(det[("plant_high", "cross", "HF")], "0")
            real = []
            for p in ST.PIPELINES:
                real += self.fake_rows(p, ["cross"], 0.03 if p == "v2" else 0.0, 16, subs, rng)
            lines = []
            ST.summarize_real(real, lines.append)
            claims = [m for m in lines if "CLAIM" in m]
            self.assertTrue(claims and all(" v2 " in m and " LF" in m for m in claims))
        finally:
            ST.OUT = old
