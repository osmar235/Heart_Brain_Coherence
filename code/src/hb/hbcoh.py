"""EEG envelope and cardiac series wavelet coherence on a common slow timebase.

Thin layer over the vendored EMG simulator wavelet functions
(src/hb/vendor/emgsim_v3_3_0/functions.py, unmodified). Nothing here changes the
transform; this module only prepares slow signals, fixes the scale grid for the
0.04 to 0.40 Hz range, and summarises coherence inside bands, time windows and
the cone of influence.
"""
from __future__ import annotations

import numpy as np
from scipy.interpolate import interp1d
from scipy.signal import butter, hilbert, resample_poly, sosfiltfilt

from hb.vendor.emgsim_v3_3_0 import functions as emg

FS_SLOW = 4.0                      # Hz, common timebase for envelopes and cardiac series
BANDS = {"LF": (0.04, 0.15), "HF": (0.15, 0.40)}
DJ = 1.0 / 12.0                    # 12 voices per octave
W0 = 6


# ----------------------------------------------------------------- cardiac series
def tachogram(r_times, fs_out=FS_SLOW, t0=0.0, t1=None, kind="cubic"):
    """Evenly sampled NN interval series (seconds) from R peak times (seconds).

    Each interval is placed at the time of the beat that ends it. Samples before the
    second beat or after the last beat are NaN (never extrapolated).
    """
    r = np.asarray(r_times, float)
    if r.size < 4 or np.any(np.diff(r) <= 0):
        raise ValueError("r_times must be increasing with at least 4 beats")
    t1 = r[-1] if t1 is None else t1
    t = np.arange(t0, t1, 1.0 / fs_out)
    f = interp1d(r[1:], np.diff(r), kind=kind, bounds_error=False, fill_value=np.nan)
    return t, f(t)


def beat_density(r_times, fs_out=FS_SLOW, t0=0.0, t1=None, fs_fine=20.0):
    """Low passed heartbeat point process (beats per second), sampled at fs_out.

    Same idea as the binary beat series of von Jakitsch et al. 2024, brought to the
    slow timebase with an anti aliasing filter instead of being kept at 1000 Hz.
    """
    r = np.asarray(r_times, float)
    t1 = r[-1] if t1 is None else t1
    n_fine = int(np.floor((t1 - t0) * fs_fine))
    counts = np.bincount(np.floor((r[(r >= t0) & (r < t0 + n_fine / fs_fine)] - t0) * fs_fine).astype(int), minlength=n_fine)
    q = int(round(fs_fine / fs_out))
    y = resample_poly(counts.astype(float) * fs_fine, 1, q)
    return t0 + np.arange(y.size) / fs_out, y


# ----------------------------------------------------------------- EEG envelope
def band_envelope(x, fs, band, fs_out=FS_SLOW, order=4):
    """Band limited Hilbert amplitude envelope, anti alias filtered and resampled to fs_out."""
    x = np.asarray(x, float)
    sos = butter(order, [band[0], band[1]], btype="band", fs=fs, output="sos")
    env = np.abs(hilbert(sosfiltfilt(sos, x)))
    q = fs / fs_out
    if abs(q - round(q)) > 1e-9:
        raise ValueError("fs must be an integer multiple of fs_out")
    y = env
    q = int(round(q))
    while q > 1:                               # staged decimation keeps the FIR short
        step = next(k for k in (10, 8, 5, 4, 2, q) if q % k == 0)
        y = resample_poly(y, 1, step)
        q //= step
    return np.arange(y.size) / fs_out, y


# ----------------------------------------------------------------- coherence
def prepare(x):
    """Linear detrend and z score. NaN is not allowed: cut the segment first."""
    x = np.asarray(x, float)
    if not np.all(np.isfinite(x)):
        raise ValueError("signal contains NaN or inf; select a clean contiguous segment")
    t = np.arange(x.size)
    x = x - np.polyval(np.polyfit(t, x, 1), t)
    return emg.zscore(x)


def cwt_slow(x, fs=FS_SLOW, dj=DJ, fmax=0.5):
    """Cached Morlet CWT on the fixed slow grid: periods from about 2/fs s up to about 132 s."""
    s0 = 2.0 / fs
    J = int(round(8.0 / dj))
    return emg.compute_cwt_for_coherence(prepare(x), fs, fmax=fmax, dj=dj, s0=s0, J=J, w0=W0)


def coherence(cx, cy):
    """Squared wavelet coherence from two cached CWTs. Returns dict with freqs, period, coi, Rsq, inside."""
    freqs, period, coi, Rsq = emg.coherence_from_cached_cwts(cx, cy)
    inside = period[:, None] <= coi[None, :]
    return dict(freqs=freqs, period=period, coi=coi, Rsq=Rsq, inside=inside, fs=cx["fs"])


def band_summary(coh, band, t_start=None, t_stop=None, min_valid_fraction=0.5):
    """Mean squared coherence in a frequency band, inside the COI, inside [t_start, t_stop) seconds.

    Returns (value, valid_fraction). Value is NaN when fewer than min_valid_fraction of
    the time frequency samples of the window lie inside the cone of influence.
    """
    f = coh["freqs"]
    rows = (f >= band[0]) & (f < band[1])
    n = coh["Rsq"].shape[1]
    t = np.arange(n) / coh["fs"]
    cols = np.ones(n, bool)
    if t_start is not None:
        cols &= t >= t_start
    if t_stop is not None:
        cols &= t < t_stop
    R = coh["Rsq"][np.ix_(rows, cols)]
    ok = coh["inside"][np.ix_(rows, cols)] & np.isfinite(R)
    frac = float(ok.mean()) if ok.size else 0.0
    if frac < min_valid_fraction:
        return float("nan"), frac
    return float(R[ok].mean()), frac


def edge_safe_seconds(n_samples, freq, fs=FS_SLOW):
    """Seconds of a record of n_samples that lie inside the COI at a given frequency (analytic, Morlet w0=6)."""
    fourier_factor = 4 * np.pi / (W0 + np.sqrt(2 + W0 ** 2))
    e_fold = np.sqrt(2) * (1.0 / freq) / fourier_factor
    return max(0.0, n_samples / fs - 2 * e_fold)
