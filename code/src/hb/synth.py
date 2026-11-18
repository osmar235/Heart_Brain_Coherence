"""Synthetic signals with known coupling, for tests and benchmarks only."""
from __future__ import annotations

import numpy as np
from scipy.signal import butter, sosfiltfilt


def band_noise(n, fs, band, rng):
    """Unit variance Gaussian noise limited to band (Hz)."""
    sos = butter(3, band, btype="band", fs=fs, output="sos")
    x = sosfiltfilt(sos, rng.standard_normal(n + int(40 * fs)))[int(20 * fs):int(20 * fs) + n]
    return x / x.std()


def ipfm_beats(mod, fs, hr_mean=1.1, depth=0.06):
    """Integral pulse frequency modulation: beat times for heart rate hr_mean*(1+depth*mod) in Hz."""
    rate = hr_mean * (1.0 + depth * np.asarray(mod, float))
    phase = np.cumsum(rate) / fs
    k = np.arange(1, int(np.floor(phase[-1])) + 1)
    return np.interp(k, phase, np.arange(phase.size) / fs)


def am_eeg(mod, fs, rng, carrier_hz=10.0, depth=0.5, noise_sd=0.5):
    """Oscillation at carrier_hz whose amplitude follows 1+depth*mod, plus white noise."""
    t = np.arange(len(mod)) / fs
    amp = np.clip(1.0 + depth * np.asarray(mod, float), 0.05, None)
    return amp * np.sin(2 * np.pi * carrier_hz * t + rng.uniform(0, 2 * np.pi)) + noise_sd * rng.standard_normal(len(mod))


def trial_locked(n, fs, period_s, rng, width_s=6.0, jitter=0.3):
    """Sum of smooth bumps, one per trial, amplitudes jittered. Models a task evoked slow response."""
    t = np.arange(n) / fs
    y = np.zeros(n)
    for onset in np.arange(2.0, n / fs, period_s):
        a = 1.0 + jitter * rng.standard_normal()
        y += a * np.exp(-0.5 * ((t - onset - width_s) / (width_s / 2.5)) ** 2)
    return (y - y.mean()) / y.std()
