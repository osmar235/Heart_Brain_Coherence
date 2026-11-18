"""Wavelet coherence between heart period and an EEG band envelope, on synthetic signals.

Run from the repository folder after installing the package:  python examples/coherence_example.py
"""
import numpy as np
from hb import hbcoh

# Synthetic example: R peak times (s) and one EEG channel sampled at 500 Hz,
# both modulated at 0.1 Hz so that their slow fluctuations are coherent in LF.
rng = np.random.default_rng(1)
rr = 0.85 + 0.05 * np.sin(2 * np.pi * 0.1 * np.arange(350) * 0.85) + 0.01 * rng.standard_normal(350)
r_times = np.cumsum(rr)
fs_eeg = 500.0
t = np.arange(int(300 * fs_eeg)) / fs_eeg
eeg = (1 + 0.5 * np.sin(2 * np.pi * 0.1 * t)) * np.sin(2 * np.pi * 10 * t) + rng.standard_normal(t.size)

t0, t1 = np.ceil(r_times[1]), np.floor(min(r_times[-1], t[-1]))            # common window, seconds
_, nn = hbcoh.tachogram(r_times, t0=t0, t1=t1)                                # heart period at 4 Hz
_, env = hbcoh.band_envelope(eeg[int(t0 * fs_eeg):int(t1 * fs_eeg)], fs_eeg, (8.0, 13.0))  # alpha envelope at 4 Hz
n = min(nn.size, env.size)
coh = hbcoh.coherence(hbcoh.cwt_slow(nn[:n]), hbcoh.cwt_slow(env[:n]))
lf, valid = hbcoh.band_summary(coh, hbcoh.BANDS["LF"])
print(f"LF coherence {lf:.3f} (share of samples inside the cone of influence {valid:.2f})")
