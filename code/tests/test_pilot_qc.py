"""Tests of the loader, event parser and R peak detector on a synthetic recording that mimics
the ds008302 file layout (top level EEGLAB variables, ECG last, large DC offset, inverted polarity)."""
import sys
import tempfile
import unittest
import warnings
from pathlib import Path

import numpy as np
import scipy.io as sio

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
warnings.simplefilter("ignore")
from hb import bidsio, cardiac, pilot_qc  # noqa: E402
from hb import synth as S  # noqa: E402

FS = 1000.0


def synthetic_ecg(beats, n, rng, invert=True, offset=405000.0):
    t = np.arange(n) / FS
    x = np.zeros(n)
    for b in beats:
        i0, i1 = int((b - 0.4) * FS), int((b + 0.5) * FS)
        if i0 < 0 or i1 >= n:
            continue
        tt = t[i0:i1] - b
        x[i0:i1] += (15000 * np.exp(-0.5 * (tt / 0.012) ** 2) - 3000 * np.exp(-0.5 * ((tt + 0.03) / 0.012) ** 2)
                     + 4000 * np.exp(-0.5 * ((tt - 0.25) / 0.05) ** 2) + 1500 * np.exp(-0.5 * ((tt + 0.18) / 0.03) ** 2))
    wander = 8000 * np.sin(2 * np.pi * 0.2 * t) + 3000 * S.band_noise(n, FS, (0.05, 0.5), rng)
    x = x + wander + 400 * rng.standard_normal(n) + 300 * np.sin(2 * np.pi * 50 * t)
    return (-x if invert else x) + offset


def make_recording(folder, sub="sub-999", task="memory", seconds=200, seed=0, eeg_labels=("Fz", "Cz", "Pz", "Oz", "T7", "T8"),
                   n_blocks=2, trials=3, common=0.0):
    rng = np.random.default_rng(seed)
    n = int(seconds * FS)
    mod = S.band_noise(n, FS, (0.04, 0.4), rng)
    beats = S.ipfm_beats(mod, FS, hr_mean=1.2, depth=0.04)
    beats = beats[(beats > 1.0) & (beats < seconds - 1.0)]
    labels = list(eeg_labels) + ["PPG", "ECG"]
    data = np.zeros((len(labels), n), np.float32)
    shared = common * rng.standard_normal(n) if common else 0.0
    for i in range(len(eeg_labels)):
        data[i] = 20 * rng.standard_normal(n) + shared
    ppg = np.zeros(n)
    for b in beats:
        i0 = int((b + 0.2) * FS)
        seg = np.exp(-0.5 * ((np.arange(-300, 500) / FS) / 0.09) ** 2)
        if i0 - 300 >= 0 and i0 + 500 < n:
            ppg[i0 - 300:i0 + 500] += 2000 * seg
    data[-2] = ppg + 131000 + 50 * rng.standard_normal(n)
    data[-1] = synthetic_ecg(beats, n, rng)
    # two blocks of three 4 digit trials
    ev, t = [], 5.0
    for blk in range(n_blocks):
        ev.append((t, "S 90")); t += 1.5
        for k in range(trials):
            if k:
                ev.append((t, "S170")); t += 0.5
            ev.append((t, "S104")); t += 2.0
            ok = "1" if (blk + k) % 2 == 0 else "0"
            for d in range(4):
                ev.append((t, f"{d + 1:02d}04{ok}")); t += 2.0
            t += 1.5
            ev.append((t, "S160")); t += 8.0
        ev.append((t, "S180")); t += 3.0
        ev.append((t, "S190")); t += 20.0
    eeg = folder / sub / "eeg"
    eeg.mkdir(parents=True)
    base = eeg / f"{sub}_task-{task}"
    event = np.zeros(len(ev) + 1, dtype=[("latency", "O"), ("duration", "O"), ("type", "O")])
    event[0] = (0.5, float("nan"), "boundary")
    for i, (o, v) in enumerate(ev, 1):
        event[i] = (o * FS + 1, 1.0, v)
    chanlocs = np.zeros(len(labels), dtype=[("labels", "O")])
    for i, lab in enumerate(labels):
        chanlocs[i] = (lab,)
    sio.savemat(str(base) + "_eeg.set", dict(data=data, srate=FS, pnts=float(n), nbchan=float(len(labels)),
                                                chanlocs=chanlocs, event=event, ref="FCz"))
    Path(str(base) + "_eeg.json").write_text('{"SamplingFrequency": 1000, "RecordingDuration": %s}' % seconds)
    Path(str(base) + "_channels.tsv").write_text("name\ttype\tunits\n" + "".join(f"{lab}\tEEG\tmicroV\n" for lab in labels))
    Path(str(base) + "_events.tsv").write_text("onset\tduration\tsample\tvalue\n-0.0005\t1\t-0.5\tboundary\n"
                                              + "".join(f"{o:.3f}\t0.001\t{int(o * FS)}\t{v}\n" for o, v in ev))
    return beats


class PilotQC(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        cls.root = Path(cls.tmp.name)
        cls.beats = make_recording(cls.root)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_loader_reads_layout_and_events(self):
        rec = bidsio.load_set(bidsio.rec_paths(self.root, "sub-999", "memory")["set"])
        self.assertEqual(rec["all_labels"][-2:], ["PPG", "ECG"])
        self.assertEqual(rec["data"].shape, (8, 200000))
        self.assertEqual(rec["events"][1][0], "S 90")
        self.assertAlmostEqual(rec["events"][1][1], 5.0, places=6)

    def test_rpeaks_found_despite_offset_inversion_and_wander(self):
        rec = bidsio.load_set(bidsio.rec_paths(self.root, "sub-999", "memory")["set"], channels=["ECG"])
        det = cardiac.detect_rpeaks(rec["data"][0], FS)
        self.assertEqual(det["polarity"], -1)
        found = det["peaks"] / FS
        err = np.array([np.min(np.abs(found - b)) for b in self.beats])
        self.assertGreater((err < 0.005).mean(), 0.995)                         # sensitivity
        extra = np.array([np.min(np.abs(self.beats - f)) for f in found])
        self.assertLess((extra > 0.05).mean(), 0.005)                           # false detections

    def test_missed_and_extra_beats_are_flagged(self):
        pk = (self.beats * FS).astype(int)
        damaged = np.sort(np.r_[np.delete(pk, 50), pk[100] + 250])
        rr, flag = cardiac.rr_flags(damaged, FS)
        self.assertGreater(flag[49], 0)
        self.assertTrue((flag[99:102] > 0).any())
        self.assertLess((flag > 0).sum(), 8)

    def test_memory_event_parser(self):
        rows = bidsio.read_tsv(bidsio.rec_paths(self.root, "sub-999", "memory")["events"])
        p = bidsio.parse_memory_events(rows, 200.0)
        self.assertEqual(len(p["blocks"]), 2)
        self.assertEqual([len(b["trials"]) for b in p["blocks"]], [3, 3])
        self.assertEqual([t["correct"] for t in p["trials"]], [True, False, True, False, True, False])
        self.assertTrue(all(t["length"] == 4 and len(t["digit_times"]) == 4 for t in p["trials"]))
        self.assertEqual(p["boundaries"], [])

    def test_full_qc_runs_and_reports(self):
        out = self.root / "out"
        row, blocks = pilot_qc.qc_recording(self.root, "sub-999", "memory", out, self.root / "rp", lambda m: None)
        self.assertEqual(row["status"], "ok")
        self.assertTrue(row["labels_match_channels_tsv"] and row["fs_matches_sidecar"] and row["event_types_match"])
        self.assertLess(row["max_event_time_diff_ms"], 1.0)
        self.assertGreater(row["pct_rr_ok"], 99.0)
        self.assertGreater(row["rr_ppg_interval_correlation"], 0.9)
        self.assertEqual(len(blocks), 2)
        self.assertTrue((out / "sub-999_task-memory_qc.png").exists())
        self.assertTrue((self.root / "rp" / "sub-999_task-memory_rpeaks.tsv").exists())


if __name__ == "__main__":
    unittest.main()
