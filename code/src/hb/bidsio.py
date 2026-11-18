"""Reading ds008302 recordings (EEGLAB .set, MATLAB v5, data embedded) and BIDS event tables."""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path

import numpy as np
import scipy.io as sio

from hb.download_ds008302 import PROJECT_ROOT, load_json, resolve_data_root


def raw_root() -> Path:
    return resolve_data_root(load_json(PROJECT_ROOT / "config" / "download_config.json"))


def rec_paths(root: Path, sub: str, task: str) -> dict:
    base = root / sub / "eeg" / f"{sub}_task-{task}"
    return {k: Path(str(base) + s) for k, s in
            (("set", "_eeg.set"), ("json", "_eeg.json"), ("channels", "_channels.tsv"), ("events", "_events.tsv"))}


def read_tsv(path: Path) -> list[dict]:
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _load_v5(path: Path) -> dict:
    m = sio.loadmat(str(path), squeeze_me=True, struct_as_record=False,
                    variable_names=["data", "srate", "pnts", "nbchan", "chanlocs", "event", "ref"])
    if isinstance(m["data"], str):
        raise ValueError(f"{path.name}: data stored in a separate file ({m['data']}), not supported")
    labels = [str(c.labels) for c in np.atleast_1d(m["chanlocs"])]
    fs = float(m["srate"])
    events = []
    for e in np.atleast_1d(m["event"]):
        dur = getattr(e, "duration", np.nan)
        dur = float(dur) if np.size(dur) == 1 and dur is not None else float("nan")
        events.append((str(e.type), (float(e.latency) - 1.0) / fs, dur))
    return dict(fs=fs, n_samples=int(m["pnts"]), labels=labels, data=m["data"], events=events, ref=str(m.get("ref", "")))


def _load_v73(path: Path) -> dict:
    """EEGLAB .set saved as MATLAB v7.3 (HDF5). Used by the recordings larger than 2 GB.
    The matrix is stored as samples x channels; a transposed view is returned, without copying."""
    try:
        import h5py
    except ImportError as e:
        raise ImportError("This recording is a MATLAB v7.3 file and needs the h5py package: "
                          "python -m pip install h5py") from e

    def text(f, ref):
        a = f[ref][()]
        return "".join(chr(int(c)) for c in np.asarray(a).ravel()) if a.dtype.kind in "ui" else ""

    def number(f, ref):
        a = np.asarray(f[ref][()], float).ravel()
        return float(a[0]) if a.size else float("nan")

    with h5py.File(str(path), "r") as f:
        fs = float(np.asarray(f["srate"][()]).ravel()[0])
        pnts = int(np.asarray(f["pnts"][()]).ravel()[0])
        labels = [text(f, r) for r in f["chanlocs"]["labels"][()].ravel()]
        ev = f["event"]
        types = ev["type"][()].ravel()
        lats = ev["latency"][()].ravel()
        durs = ev["duration"][()].ravel() if "duration" in ev else [None] * len(types)
        events = [(text(f, t), (number(f, la) - 1.0) / fs, number(f, d) if d is not None else float("nan"))
                  for t, la, d in zip(types, lats, durs)]
        ref = "".join(chr(int(c)) for c in np.asarray(f["ref"][()]).ravel()) if "ref" in f else ""
        d = f["data"]
        arr = np.empty(d.shape, d.dtype)
        d.read_direct(arr)
    data = arr.T if arr.shape[1] == len(labels) else arr
    return dict(fs=fs, n_samples=pnts, labels=labels, data=data, events=events, ref=ref)


def load_set(path: Path, channels=None) -> dict:
    """Load an EEGLAB .set file (MATLAB v5, or v7.3 for recordings above 2 GB).
    channels=None keeps all rows; otherwise a list of labels to keep.

    Returns dict(fs, n_samples, labels, data (float32, channels x samples), events [(type, latency_s, duration_samples)]).
    The whole data matrix is read into memory once (1.5 GB for a typical task recording, up to 4.5 GB).
    """
    try:
        r = _load_v5(path)
    except NotImplementedError:
        r = _load_v73(path)
    labels, data = r["labels"], r["data"]
    if data.shape[0] != len(labels) or data.shape[1] != r["n_samples"]:
        raise ValueError(f"{path.name}: data shape {data.shape} disagrees with nbchan/pnts")
    out_labels = labels
    if channels is not None:
        idx = [labels.index(c) for c in channels]
        data = np.ascontiguousarray(data[idx])
        out_labels = list(channels)
    return dict(fs=r["fs"], n_samples=r["n_samples"], labels=out_labels, all_labels=labels, data=data,
                events=r["events"], ref=r["ref"])


def parse_memory_events(rows: list[dict], duration_s: float) -> dict:
    """Trials, blocks and discontinuities from a task-memory events.tsv.

    Digit codes are [position][length][flag], flag 1 correct, 0 incorrect, x unknown.
    A trial is correct when every digit flag is 1.
    """
    ev = [(float(r["onset"]), r["value"]) for r in rows]
    trials, blocks, cur, start = [], [], None, None
    for o, v in ev:
        if v == "S 90":
            if not blocks or o - blocks[-1]["t0"] > 1.0:
                blocks.append({"t0": o, "t1": None, "trials": []})
            start = o
        elif v == "S170":
            start = o
        m = re.fullmatch(r"S1(0\d|1[0-3])", v)
        if m:
            cur = {"length": int(m.group(1)), "t_start": start if start is not None else o, "t_baseline": o,
                   "digit_times": [], "flags": [], "t_probe": None, "block": len(blocks) - 1}
        d = re.fullmatch(r"(\d\d)(\d\d)([01x])", v)
        if d and cur is not None:
            cur["digit_times"].append(o)
            cur["flags"].append(d.group(3))
        if v == "S160" and cur is not None:
            cur["t_probe"] = o
            cur["correct"] = None if "x" in cur["flags"] else all(f == "1" for f in cur["flags"])
            trials.append(cur)
            if blocks:
                blocks[-1]["trials"].append(cur)
            cur = None
        if v in ("S180", "S190", "S200") and blocks and blocks[-1]["t1"] is None:
            blocks[-1]["t1"] = o
    for b in blocks:
        b["length"] = b["trials"][0]["length"] if b["trials"] else None
        if b["t1"] is None and b["trials"]:
            b["t1"] = b["trials"][-1]["t_probe"]
    boundaries = sorted({round(o, 3) for o, v in ev if v == "boundary" and 1.0 < o < duration_s - 1.0})
    return dict(trials=trials, blocks=blocks, boundaries=boundaries,
                overflow=[o for o, v in ev if v == "Buffer Overflow"])


def sidecar(path: Path) -> dict:
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)
