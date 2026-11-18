# Provenance of the vendored wavelet code

- File: `functions.py`, byte identical copy. Do not edit it here.
- Source on the authors' PC: `<USER_HOME>\projects\Version3.3.1_EMGSim\functions.py`
- Internal version string: 3.3.0 (folder name says 3.3.1)
- SHA256: 0aaf110c0ce9fc36929686ac0adb6c4fd3ea855162bf98a950b834da63c55f19 (matches that folder's SHA256SUMS.txt)
- Copied: 2026-11-09
- Authors per CITATION.cff: Pinto Neto O, Pinho TOR, Wang Y, Balbinot G, Kennedy DM
- Software archive DOI listed in CITATION.cff: 10.5281/zenodo.19703513 (the file itself says to update it if V3 gets a new DOI)
- Associated article as given by O. Pinto Neto: "A computational framework for EMG simulation and coherence-based biomarker analysis on neural crosstalk in bimanual coordination", Computer Methods and Programs in Biomedicine 2027;288:109657. https://doi.org/10.1016/j.cmpb.2026.109657 (volume 288, January 2027; details supplied by the author).
- License: none stated in the source folder copied on 2026-11-09. On 2026-11-17 O. Pinto Neto stated that the EMG simulation framework is public on his GitHub under the MIT license; this copy is redistributed under those terms (decision D-022). Checked on 2026-11-17: https://github.com/osmar235/EMG_Simulator_WaveletCoherence has an MIT license (Copyright (c) 2025 osmar235); on that day its main branch held functions.py version 2.0.0. Later on 2026-11-17 Osmar published the Version3.3.1_EMGSim folder there: tags v3.3.0 and v3.3.1 (commit 71b6719, GitHub release v3.3.1) hold a functions.py with the same SHA256 as this copy (0aaf110c...). Zenodo version DOI of that release: to be added.

Functions used by this project: `compute_cwt_for_coherence`, `coherence_from_cached_cwts`, `compute_wavelet_coherence`, `smoothwavelet` (indirectly), `random_derangement`, `phase_randomize`.

Relation to the 2024 Biomedical Engineering Letters paper: that analysis ran in MATLAB. This Python implementation (pycwt Morlet transform plus Torrence and Webster / Grinsted smoothing) is a documented SUBSTITUTION, not the same code. See DECISOES proposal D-005.
