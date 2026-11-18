"""R peak detection and beat quality flags. No external ECG package; method is fully visible here."""
from __future__ import annotations

import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.signal import butter, find_peaks, sosfiltfilt


def detect_rpeaks(ecg, fs):
    """Detect R peaks. Returns dict(peaks (sample index), polarity (+1 or -1), filtered (polarity corrected)).

    Steps: 5 to 25 Hz zero phase band pass; polarity from the larger of the two extreme tails;
    slope energy (squared derivative, 120 ms moving mean); candidates at least 300 ms apart;
    keep candidates above 30% of the local 75th percentile of candidate energies (20 s window);
    place each peak at the maximum of the filtered signal within 60 ms.
    """
    ecg = np.asarray(ecg, float)
    x = sosfiltfilt(butter(3, [5, 25], btype="band", fs=fs, output="sos"), ecg - np.median(ecg))
    hi, lo = np.percentile(x, 99.9), -np.percentile(x, 0.1)
    polarity = 1 if hi >= lo else -1
    x = x * polarity
    energy = uniform_filter1d(np.gradient(x) ** 2, int(0.12 * fs))
    cand, _ = find_peaks(energy, distance=int(0.3 * fs))
    h = energy[cand]
    half = int(10 * fs)
    left = np.searchsorted(cand, cand - half)
    right = np.searchsorted(cand, cand + half)
    ref = np.array([np.percentile(h[a:b], 75) for a, b in zip(left, right)])
    cand = cand[h > 0.3 * ref]
    w = int(0.06 * fs)
    peaks = np.array([max(c - w, 0) + int(np.argmax(x[max(c - w, 0):c + w + 1])) for c in cand], dtype=int)
    peaks = np.unique(peaks)
    keep = np.r_[True, np.diff(peaks) > int(0.25 * fs)]
    return dict(peaks=peaks[keep], polarity=polarity, filtered=x)


def rr_flags(peaks, fs, lo=0.3, hi=2.0, rel=0.20, k=11):
    """RR intervals (s) ending at each beat from the second on, with flags.

    flag 0 ok; 1 outside [lo, hi] s; 2 deviates more than rel from the local median of k intervals.
    """
    rr = np.diff(peaks) / fs
    flag = np.zeros(rr.size, int)
    pad = k // 2
    padded = np.pad(rr, pad, mode="edge")
    med = np.array([np.median(padded[i:i + k]) for i in range(rr.size)])
    flag[np.abs(rr - med) > rel * med] = 2
    flag[(rr < lo) | (rr > hi)] = 1
    return rr, flag


def ppg_peaks(ppg, fs):
    """Pulse peaks from PPG: 0.5 to 8 Hz band pass, peaks at least 300 ms apart with prominence above 0.5 SD."""
    x = sosfiltfilt(butter(2, [0.5, 8], btype="band", fs=fs, output="sos"), np.asarray(ppg, float) - np.median(ppg))
    p, _ = find_peaks(x, distance=int(0.3 * fs), prominence=0.5 * x.std())
    return p


# ----------------------------------------------------------------------------- correction policy (proposal D-007)
def template_screen(x, peaks, fs, half_ms=50, min_corr=0.5):
    """Correlation of each detected beat with the recording's median QRS shape.

    x is the filtered, polarity corrected ECG from detect_rpeaks. Returns (keep mask, correlations).
    Detections that do not look like this participant's QRS (noise bursts) get keep=False.
    """
    w = int(half_ms * fs / 1000)
    ok = (peaks > w) & (peaks < x.size - w)
    seg = np.stack([x[p - w:p + w + 1] for p in peaks[ok]])
    tpl = np.median(seg, axis=0)
    tz = (tpl - tpl.mean()) / (tpl.std() + 1e-12)
    sz = seg - seg.mean(axis=1, keepdims=True)
    corr = (sz @ tz) / (np.linalg.norm(sz, axis=1) * np.sqrt(tz.size) + 1e-12)
    out = np.zeros(peaks.size)
    out[ok] = corr
    keep = out >= min_corr
    keep[~ok] = True          # beats cut by the record edge are not judged here
    return keep, out


def _rolling_qd(v, k=91):
    """Quartile deviation of |v| in a centred window of k beats."""
    a = np.abs(v)
    pad = k // 2
    p = np.pad(a, pad, mode="reflect") if a.size > pad else np.pad(a, pad, mode="edge")
    win = np.lib.stride_tricks.sliding_window_view(p, k)
    q1, q3 = np.percentile(win, [25, 75], axis=1)
    return (q3 - q1) / 2.0


def _rolling_median(v, k=11):
    pad = k // 2
    win = np.lib.stride_tricks.sliding_window_view(np.pad(v, pad, mode="edge"), k)
    return np.median(win, axis=1)


def classify_intervals(rr, alpha=5.2, c1=0.13, c2=0.17):
    """Interval artifact classes with time varying thresholds, after Lipponen and Tarvainen (2019,
    J Med Eng Technol 43:173). Thresholds adapt to each person's own variability, so large but
    regular variability is not flagged.

    Returns an int array per interval: 0 normal, 1 ectopic pattern, 2 missed beat (interval too long
    by about one beat), 3 extra beat (this interval plus the next make one normal interval),
    4 other long or short interval.
    """
    rr = np.asarray(rr, float)
    n = rr.size
    d = np.r_[0.0, np.diff(rr)]
    th1 = alpha * _rolling_qd(d) + 1e-9
    dn = d / th1
    med = _rolling_median(rr)
    m = rr - med
    m[m < 0] *= 2.0
    th2 = alpha * _rolling_qd(m) + 1e-9
    mn = m / th2
    cls = np.zeros(n, int)
    dpad = np.r_[0.0, dn, 0.0, 0.0]            # dn[i] is dpad[i + 1]
    for i in range(n):
        s11 = dn[i]
        nb = (dpad[i], dpad[i + 2])
        s12 = max(nb) if s11 > 0 else min(nb)
        if (s11 > 1 and s12 < -c1 * s11 - c2) or (s11 < -1 and s12 > -c1 * s11 + c2):
            cls[i] = 1
            continue
        nx = (dpad[i + 2], dpad[i + 3])
        s22 = min(nx) if s11 >= 0 else max(nx)
        if (s11 > 1 and s22 < -1) or (s11 < -1 and s22 > 1) or abs(mn[i]) > 3:
            if abs(rr[i] / 2.0 - med[i]) < th2[i]:
                cls[i] = 2
            elif i + 1 < n and abs(rr[i] + rr[i + 1] - med[i]) < th2[i]:
                cls[i] = 3
            else:
                cls[i] = 4
    return cls


def correct_beats(r_times, max_passes=3):
    """Apply the correction policy to beat times (s).

    Extra beats are removed, single missed beats are inserted at the midpoint, and intervals with an
    ectopic or other abnormal pattern are kept in time but marked so that they are interpolated, never
    used as data. Returns dict(times, nn, artifact (bool per interval), counts).
    """
    t = np.asarray(r_times, float).copy()
    inserted = np.zeros(t.size, bool)
    counts = {"extra_removed": 0, "missed_inserted": 0, "ectopic_marked": 0, "other_marked": 0}
    for _ in range(max_passes):
        cls = classify_intervals(np.diff(t))
        extra = np.where(cls == 3)[0]
        missed = np.where(cls == 2)[0]
        if extra.size == 0 and missed.size == 0:
            break
        drop = extra + 1                                     # the beat that ends an "extra" interval
        mid = (t[missed] + t[missed + 1]) / 2.0
        keep = np.ones(t.size, bool)
        keep[drop] = False
        t2 = np.r_[t[keep], mid]
        ins2 = np.r_[inserted[keep], np.ones(mid.size, bool)]
        order = np.argsort(t2)
        t, inserted = t2[order], ins2[order]
        counts["extra_removed"] += int(extra.size)
        counts["missed_inserted"] += int(missed.size)
    nn = np.diff(t)
    cls = classify_intervals(nn)
    artifact = cls > 0
    artifact |= inserted[1:] | inserted[:-1]                 # intervals touching an inserted beat are not data
    counts["ectopic_marked"] = int((cls == 1).sum())
    counts["other_marked"] = int(((cls == 4) | (cls == 2) | (cls == 3)).sum())
    return dict(times=t, nn=nn, artifact=artifact, inserted=inserted, counts=counts)
