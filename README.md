# Heart and brain coherence: code and results of "One Dataset, Three Verdicts"

Code, result tables and figures of the study "One Dataset, Three Verdicts: Why Positive and Negative Controls Matter for Brain and Heart Coherence" (Osmar Pinto Neto; manuscript in preparation, not peer reviewed).

The study reanalyzes electrocardiography and EEG from 122 adults performing an adaptive digit span task (OpenNeuro ds008302). It asks whether slow wavelet coherence between heart period and EEG amplitude envelopes survives negative controls (synthetic participants without coupling, surrogate nulls, artifact controls), and whether those controls could have detected a coupling of the observed size (positive controls).

## Use the coherence functions

Install the package (Python 3.10 or newer):

```
pip install git+https://github.com/osmar235/Heart_Brain_Coherence
```

Then, from a folder with this repository, `python examples/coherence_example.py` runs the whole chain on synthetic signals:

```python
from hb import hbcoh

_, nn = hbcoh.tachogram(r_times, t0=t0, t1=t1)            # heart period from R peak times, 4 Hz
_, env = hbcoh.band_envelope(eeg_segment, fs_eeg, (8, 13))  # alpha amplitude envelope, 4 Hz
coh = hbcoh.coherence(hbcoh.cwt_slow(nn), hbcoh.cwt_slow(env))
lf, valid = hbcoh.band_summary(coh, hbcoh.BANDS["LF"])     # mean coherence in 0.04 to 0.15 Hz
```

Main modules in `code/src/hb/`:

| Module | What it does |
| --- | --- |
| `hbcoh` | Heart period and beat density series, EEG band envelopes, Morlet wavelet coherence on a common 4 Hz timebase, band summaries inside the cone of influence |
| `trial_shuffle`, `symmetric_tests` | Surrogate nulls (trial shuffle, cross participant) and planted coupling (positive controls) |
| `null_benchmark`, `null_calibration` | Synthetic participants without coupling, on the real trial timing, and bias limits |
| `vendor/emgsim_v3_3_0` | The wavelet functions (version 3.3.0) of our EMG simulation framework, https://github.com/osmar235/EMG_Simulator_WaveletCoherence, copied unmodified |

## Reproduce the results (Windows)

From a CMD window in this folder:

1. `scripts\setup_windows.cmd` creates `.venv` and installs the packages (exact versions from `requirements-lock.txt` when your Python has the same minor version as the build; otherwise the lower bounds of `requirements.txt`). To use an environment you already have, run `set HB_PYTHON=<path to its python.exe>` and skip this step.
2. `scripts\run_tests.cmd` runs the automated tests on synthetic signals (about 5 minutes).
3. `scripts\verify_release.cmd` checks the files against `SHA256SUMS.txt` and regenerates Table 3, the figure source data and the manuscript from the stored result tables, byte for byte.
4. To rerun the analysis from the raw recordings (about 186 GB): set the data location with `set HB_DATA_ROOT=D:\data\ds008302` and run the steps listed by `scripts\run_analysis.cmd`, in order. Downloads start only when you run the download steps, and every file is checked against its SHA256 in `data_manifest/dataset_manifest.tsv`.

`docs/METHODS.md` describes every step, `docs/REPRODUCIBILITY.md` lists versions, seeds and tolerances, and `docs/RESULTS_INDEX.md` connects each number of the manuscript to its result file.

## Contents

| Folder | Contents |
| --- | --- |
| `code/` | The analysis package `hb` and its tests |
| `examples/` | A short example of the coherence functions |
| `scripts/` | Windows launchers for setup, tests, verification and each analysis step |
| `config/`, `data_manifest/` | Dataset configuration and the manifest of the recordings, with SHA256 of every file |
| `results/tables/` | The numerical results of every analysis stage |
| `results/figure_source_data/` | The numbers plotted in each figure |
| `figures/` | Figures 1 to 3 (PNG and PDF) |
| `manuscript/` | The manuscript, its source and the reference list |
| `docs/` | Methods, reproducibility notes, results index, verification status |

## Data

OpenNeuro ds008302 version 1.0.0, CC0: doi:10.18112/openneuro.ds008302.v1.0.0 (NEMAR mirror on008302: doi:10.82901/nemar.on008302). The recordings are not included here.

## How to cite

If you use this code or these results, please cite the article (reference to be added on publication) and this repository (see `CITATION.cff`; Zenodo DOI to be added). If you use the wavelet functions, please also cite the EMG simulation framework (https://github.com/osmar235/EMG_Simulator_WaveletCoherence):

Pinto Neto O, Pinho TOR, Wang Y, Balbinot G, Kennedy DM. A computational framework for EMG simulation and coherence-based biomarker analysis on neural crosstalk in bimanual coordination. Computer Methods and Programs in Biomedicine. 2027;288:109657. doi:10.1016/j.cmpb.2026.109657

## License

Code (`code/`, `scripts/`, `examples/` and the configuration files): MIT license, see `LICENSE`. The vendored wavelet functions in `code/src/hb/vendor/emgsim_v3_3_0/` are distributed under the MIT license of https://github.com/osmar235/EMG_Simulator_WaveletCoherence. Manuscript, documents, result tables and figures: CC BY 4.0, see `LICENSE-CC-BY-4.0.md`.

Contact: Osmar Pinto Neto, Department of Kinesiology, California State University San Marcos (osmar@csusm.edu).
