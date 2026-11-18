> Public repository note: references below to analysis plans (`*_PLAN*.md`, `ANALYSIS_FREEZE*.md`), stage reports (`*_RESULT*.md`), decision records (`DECISOES*.md`), the execution log and run logs (`logs/`) point to the author's project records, which are not part of this repository. They are available from the author on request.

# METHODS: the analysis pipeline, step by step

Technical companion to the Methods of the manuscript. Each step lists its module in `src/hb/` (release: `code/src/hb/`), its Windows launcher in `scripts/` (release: `scripts/pipeline/`), its inputs and outputs, and the plan that fixed its rules. Times are those of the authors' runs (Windows 11, 6 worker processes); they scale roughly with the number of workers.

## 0. Environment and data

| Step | Module | Launcher | What it does |
| --- | --- | --- | --- |
| 00 | `check_environment` | `00_check_environment.cmd` | Python version, project root, manifest, free space, reachability of the pinned GitHub commit and of OpenNeuro S3 |
| 00b | unittest | `00b_run_tests.cmd` | Automated tests on synthetic signals (85 tests) |
| 01, 02, 02b | `download_ds008302` | `01_download_metadata.cmd`, `02_download_pilot.cmd`, `02b_download_cohort.cmd` | Downloads metadata, the three pilot participants, then the cohort selection (memory task, 742 files, 185.9 GB), verifying every file against the SHA256 of `data_manifest/ds008302_v1.0.0_manifest.tsv`. Resumable. Downloads only when the user runs these launchers. |

Raw data location: `HB_DATA_ROOT` environment variable, else `config/local_paths.json` (`raw_data_root`), else `data/raw/ds008302` under the project root.

## 1. Pilot and benchmarks (2026-11-09 and 10)

| Step | Module | Plan or document |
| --- | --- | --- |
| 03 pilot inventory | `pilot_qc` | `docs/PILOT_INVENTORY.md` |
| 03b beat correction (pilot) | `beat_correction` | D-007 |
| 05 wavelet core test | unittest `test_wavelet_core` | D-005 |
| 05b, 05c null benchmarks on real timing | `null_benchmark`, `null_benchmark_real` | `docs/NULL_BENCHMARK.md` |
| 06 slow signals (pilot) | `slow_signals` | `docs/PILOT_COHERENCE_PLAN.md` |
| 07 pilot coherence | `pilot_coherence` | `docs/PILOT_COHERENCE_PLAN.md`, `docs/PILOT_COHERENCE_RESULT.md` |

## 2. Primary analysis (stage 1 of Table 1)

`08_run_cohort.cmd` (`hb.cohort`) per participant: ECG R peaks (5 to 25 Hz band pass, slope energy, 300 ms refractory period, adaptive threshold; QRS template screening), interval classification after Lipponen and Tarvainen, marked intervals bridged in a 4 Hz NN series (`cardiac`, `beat_correction`); EEG decimated to 200 Hz, median reference, region envelopes (fourth order Butterworth band pass, Hilbert amplitude, region mean, 4 Hz, robust outlier interpolation) (`slow_signals`); trial locked averages removed; within block coherence minus between block reference (`pilot_coherence`, `hbcoh`, vendored wavelet code). `09_cohort_stats.cmd` (`hb.cohort_stats`): participant means, one sample t tests, Holm over six cells, Wilcoxon, span contrast. Plan: `docs/ANALYSIS_FREEZE_v1.md`. Outputs: `derivatives/` (per participant), `results/cohort_stats/`.

## 3. Sensitivity checks (stage 2)

`10_robust_reprocess.cmd` and `11_robust_coherence.cmd`: cardiac field template subtraction and ECG derived respiration regression. `12_null_calibration.cmd`: phase randomized surrogates (8 per participant). Plan: `docs/ROBUSTNESS_PLAN_v1.md`. Outputs: `results/robustness/`, `results/null_calibration/`.

## 4. Joint pipeline and synthetic validation (stage 3)

`13_joint_calibration.cmd` (all adjustments together, circular shift null) and `14_null_validation.cmd` (synthetic participants without coupling on each participant's real timing; independent and shared trial scenarios). Plan: `docs/REVIEW_RESPONSE_PLAN_v1.md`. Outputs: `results/joint_calibration/`, `results/null_validation/`.

## 5. Trial shuffle null (stage 4)

`15a_realism_decomposition.cmd`, `15b_shuffle_validation.cmd`, `15c_shuffle_real.cmd` (`hb.trial_shuffle`): realism metrics and matched shares, decomposition of the synthetic bias, validation of the trial shuffle null on two synthetic cohorts (gate written to `results/trial_shuffle/gate.json` with the upper bias limits), then the real data (16 surrogates). Plan: `docs/STEP15_PLAN_v1.md`.

## 6. ICA cleaning (stage 5)

`16a_ica_fit.cmd` (`hb.ica_v2 fit`, about 5 h): 250 Hz, notch 50 and 100 Hz, 1 to 100 Hz, bad channels (at most 12), average reference, extended infomax ICA with the Picard solver, ICLabel. `16b_ica_apply.cmd` and `16c_ica_analyse.cmd` (`hb.ica_v2_apply`): removal of components most probably eye, muscle, heart, line noise or channel noise (probability at least 0.8), envelopes as in step 2, 55 to 95 Hz muscle proxy, three variants (ICA, no removal, muscle proxy regressed), trial shuffle null (16 surrogates). Requires MNE, mne-icalabel, python-picard and onnxruntime. Plan: `docs/ANALYSIS_FREEZE_v2_EEG.md` (amendment 1).

## 7. Symmetric tests (stage 6)

`18a_symmetric_synthetic.cmd` (`hb.symmetric_tests synthetic`, about 9 h): cross participant null validated on the independent and matched synthetic cohorts; planted coupling at three levels with the trial shuffle and cross participant nulls. `18b_symmetric_real.cmd` (`real`, then `inject`, about 8 h): cross participant null on the real data in three pipelines (16 surrogates); the middle planted level injected into the real ICA data (8 surrogates, paired seeds). Plan: `docs/STEP18_PLAN_v1.md`; levels from a synthetic pilot (`results/symmetric/levels.json`).

## 8. Tables, figures and manuscript

`hb.make_tables_three_verdicts` (Table 3, claim checks), `hb.make_figures_three_verdicts` (Figures 1 to 3 with source data and source checksums), `hb.build_manuscript` (citations numbered in order of first use, Table 3 inserted). These read result tables only. `hb.build_release` builds the portable folder `submission_ready/latest/`; `hb.verify_release` checks it.

## Vendored code

`src/hb/vendor/emgsim_v3_3_0/functions.py` is a byte identical copy of the wavelet and surrogate functions of the authors' EMG simulation framework (version string 3.3.0; SHA256 and provenance in `src/hb/vendor/emgsim_v3_3_0/PROVENANCE.md`). It replaced the MATLAB implementation of the 2024 study as a documented substitution (D-005).
