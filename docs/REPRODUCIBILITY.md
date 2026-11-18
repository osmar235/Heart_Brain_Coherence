> Public repository note: references below to analysis plans (`*_PLAN*.md`, `ANALYSIS_FREEZE*.md`), stage reports (`*_RESULT*.md`), decision records (`DECISOES*.md`), the execution log and run logs (`logs/`) point to the author's project records, which are not part of this repository. They are available from the author on request.

# REPRODUCIBILITY

## What can be reproduced, and how

| Level | What is reproduced | Needs raw data | Time | Command (release folder) |
| --- | --- | --- | --- | --- |
| A. Release integrity | File checksums, required files, privacy scan, portable paths | No | seconds | `scripts\verify_release.cmd` |
| B. Tables, figures and manuscript from the result tables | Table 3, claim checks, figure source data, Figures 1 to 3, the assembled manuscript | No | about a minute | `scripts\run_figures.cmd`, then `scripts\verify_release.cmd` |
| C. Unit tests | 85 tests of the pipeline on synthetic signals | No | about 5 to 10 minutes | `scripts\run_tests.cmd` |
| D. Synthetic analyses | Steps 14, 15b (synthetic part), 18a: use only the public event files | Event files only (small) | hours | `scripts\run_analysis.cmd <step>` |
| E. Full cohort | Every step from the raw recordings | Yes, 185.9 GB | about four days of computing in total | `scripts\run_analysis.cmd <step>` in order |

Level B compares regenerated files byte for byte with the stored ones. Figure images (PNG, PDF) are not compared byte for byte because they depend on the Matplotlib version and fonts; their numerical content is the source data, which is compared.

## Software

Cohort analyses (steps 1 to 18) ran on Windows 11 with Python 3.14.3 and, as recorded in the logs: NumPy 2.4.3, SciPy 1.17.1, MNE 1.13.2, mne-icalabel 0.10.0, python-picard 0.8.2. onnxruntime 1.31.1 was requested by `16a_ica_fit.cmd`. The exact list of installed packages is `requirements-lock.txt` in the release; when the release was built outside the analysis environment that file says so and lists only the versions recorded in the logs.

Some steps were run by Claude in a Linux workspace with Python 3.10 (NumPy 2.2.6, SciPy 1.15.3, Matplotlib 3.10): the synthetic pilot of step 18 and the figure, table and manuscript scripts. The execution log names the environment of each run.

## Random seeds

Every random step uses a seed built from a fixed constant, the participant number and the surrogate index (`numpy.random.default_rng([constant, participant, ...])`).

| Step | Constant | Module |
| --- | --- | --- |
| Pilot coherence | 20261111 | `pilot_coherence.SEED` |
| Phase calibration (step 12) | 20261112 | `null_calibration.SEED` |
| Joint calibration (step 13) | 20261113 | `joint_calibration.SEED` |
| Synthetic validation (step 14) | 20261114 | `null_validation.SEED` |
| Trial shuffle validation, real, decomposition (step 15) | 20261115, 20261116, 20261117 | `trial_shuffle.SEED_VALID`, `SEED_REAL`, `SEED_DECOMP` |
| ICA fit (step 16a) | 20261118 | `ica_v2.SEED` |
| ICA analysis (step 16c) | 20261119 | `ica_v2_apply.SEED` |
| Symmetric tests (step 18) | 20261120 | `symmetric_tests.SEED18` |

## Determinism and tolerances

Within one environment the pipeline is deterministic: step 15c recomputed the observed values of step 1 and found a maximum absolute difference of 0 over 1,464 values (`logs/trial_shuffle_real_20261112_210753.log`). Across environments, floating point results can differ in the last digits (different BLAS libraries, NumPy versions, and the ICA solver's iterations). No cross platform rerun of the cohort has been made, so no cross platform tolerance has been measured. Result tables are stored rounded to five decimals; a difference below 1e-5 in a regenerated table is within rounding. The Table 3 and figure scripts reproduced the stored files exactly in the Linux workspace (Python 3.10) from tables made on Windows (Python 3.14), see `EXECUTION_LOG.md`, 2026-11-17.

## Data

OpenNeuro ds008302, version 1.0.0, CC0 (doi:10.18112/openneuro.ds008302.v1.0.0; NEMAR mirror on008302, doi:10.82901/nemar.on008302). Metadata pinned to GitHub commit `af3998a993adb108dc29a58b70080bcc72a42098` of nemarDatasets/on008302. The manifest `data_manifest/dataset_manifest.tsv` (release) lists every file of the dataset with its size and SHA256 (4,506 files, 216.9 GB); `data_manifest/download_report_cohort.tsv` lists the 742 files of the cohort selection used here (185.9 GB) and their verification status. The raw data are not part of the release.

## Known limitations of reproducibility

- Steps 1 to 18 have been run once, on one computer. Their reproduction elsewhere has not been tested.
- The release was assembled after the analyses; the code in it is the code at release time. Modules changed after their run were changed only to add the portable path options (environment variables with the old behavior as default); this is checked by the unit tests and by the identical regeneration of figure source data.
- `requirements.txt` gives lower bounds only. Exact pins require `requirements-lock.txt` generated on the analysis computer.
- The vendored wavelet code has no license statement; see `LICENSE`.
