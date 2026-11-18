> Public repository note: references below to analysis plans (`*_PLAN*.md`, `ANALYSIS_FREEZE*.md`), stage reports (`*_RESULT*.md`), decision records (`DECISOES*.md`), the execution log and run logs (`logs/`) point to the author's project records, which are not part of this repository. They are available from the author on request.

# RESULTS INDEX

Connects each manuscript table, figure and important numerical claim to the result file that holds it and to the code that wrote that file. Paths are relative to the project root; in the release folder the same files are under `results/tables/` (for example `results/cohort_stats/group_tests.tsv` becomes `results/tables/cohort_stats/group_tests.tsv`). Cells are identified by `eeg` (theta, alpha, beta) and `band` (LF, HF). "Primary variant" means `variant = full`.

## Code that produced each result folder

| Result folder | Module | Windows launcher | Plan written before the run | Execution log |
| --- | --- | --- | --- | --- |
| `results/cohort_stats/` | `hb.cohort`, `hb.cohort_stats` | `08_run_cohort.cmd`, `09_cohort_stats.cmd` | `docs/ANALYSIS_FREEZE_v1.md` | `logs/cohort_*.log`, `logs/cohort_stats_*.log` |
| `results/robustness/` | `hb.robust_reprocess`, `hb.robust_coherence` | `10_robust_reprocess.cmd`, `11_robust_coherence.cmd` | `docs/ROBUSTNESS_PLAN_v1.md` | `logs/robust_*.log` |
| `results/null_calibration/` | `hb.null_calibration` | `12_null_calibration.cmd` | `docs/ROBUSTNESS_PLAN_v1.md` | `logs/null_calibration_*.log` |
| `results/joint_calibration/` | `hb.joint_calibration` | `13_joint_calibration.cmd` | `docs/REVIEW_RESPONSE_PLAN_v1.md` | `logs/joint_calibration_*.log` |
| `results/null_validation/` | `hb.null_validation` | `14_null_validation.cmd` | `docs/REVIEW_RESPONSE_PLAN_v1.md` | `logs/null_validation_*.log` |
| `results/trial_shuffle/` | `hb.trial_shuffle` | `15a`, `15b`, `15c` | `docs/STEP15_PLAN_v1.md` | `logs/trial_shuffle_*.log` |
| `results/ica_v2/` | `hb.ica_v2`, `hb.ica_v2_apply` | `16a`, `16b`, `16c` | `docs/ANALYSIS_FREEZE_v2_EEG.md` | `logs/ica_v2_*.log` |
| `results/symmetric/` | `hb.symmetric_tests` | `18a`, `18b` | `docs/STEP18_PLAN_v1.md` | `logs/symmetric_*.log` |
| `results/tables_three_verdicts/` | `hb.make_tables_three_verdicts` | `19_build_release.cmd` (also run directly) | not an analysis | `docs/EXECUTION_LOG.md` |
| `results/figures_three_verdicts/` | `hb.make_figures_three_verdicts` | `19_build_release.cmd` (also run directly) | not an analysis | `docs/EXECUTION_LOG.md` |
| `results/figures/` (historical, version 2) | `hb.make_figures` | `17_figures.cmd` | not an analysis | `docs/EXECUTION_LOG.md` |

## Tables

| Item | Source | How it was made |
| --- | --- | --- |
| Table 1 (chronology) | Plan documents listed above (dates in their headers) and `docs/DECISOES_PROPOSED_2026-11-09.md` | Written by hand from the documents |
| Table 2 (what each test can exclude) | Methods of each plan; outcomes from the files listed under Verdicts 1 to 3 below | Written by hand |
| Table 3 (final status of the six cells) | `results/tables_three_verdicts/table_3.tsv` and `.md` | `hb.make_tables_three_verdicts` from `cohort_stats/group_tests.tsv`, `ica_v2/shuffle_calibrated_H1_picard.tsv`, `symmetric/real_cross_calibrated_H1.tsv`, `symmetric/synthetic_cells.tsv`, `trial_shuffle/gate.json`, `symmetric/gate_cross.json` |
| Numerical claim checks | `results/tables_three_verdicts/claims_check.tsv` | Same script: margins of lower limits against limits, encoding and retention values, partner block share |

## Figures

| Figure | Source data stored with the figure | Source tables (with SHA256 prefixes) | Script |
| --- | --- | --- | --- |
| 1 | `results/figures_three_verdicts/figure_1_source_data.tsv` | `results/figures_three_verdicts/figure_sources.tsv` | `hb.make_figures_three_verdicts` |
| 2 | `figure_2_source_data.tsv` | same | same |
| 3 | `figure_3_source_data.tsv` | same | same |

## Numerical claims by section

### Methods and sample

| Claim | Source |
| --- | --- |
| 124 participants, 122 with a working memory recording; 99 women, 23 men; age 20.5 ± 2.3 (18 to 32) | `docs/DATA_AUDIT.md` (from the dataset's participants.tsv) |
| Block median 258 s; median trial start interval 25 s | `docs/DATA_AUDIT.md`, `docs/COHORT_RESULT_v1.md` |
| Median heart rate 76.2 bpm (IQR 71.1 to 84.6); 1.6% of intervals marked (maximum 9.9%) | `docs/COHORT_RESULT_v1.md` (from `results/cohort_qc/beats/*_beat_summary.tsv`) |
| Median 18 blocks (11 to 30); 89% to 91% of block and signal pairs retained | `results/cohort_stats/block_exclusions.tsv`, `docs/COHORT_RESULT_v1.md` |
| 70% of cross participant partner blocks with the same sequence length | `logs/symmetric_real_20261116_084503.log` (218,926 of 313,056); `claims_check.tsv` |
| Planted levels 0.032, 0.045, 0.064 | `results/symmetric/levels.json`; pilot `results/symmetric/pilot_levels.tsv` |
| Seeds | `docs/REPRODUCIBILITY.md` |

### Verdict 1

| Claim | Source (row keys) |
| --- | --- |
| Six of six cells; Holm p 3 × 10⁻⁹ to 0.044 | `cohort_stats/group_tests.tsv`: `hypothesis = H1_existence`, `variant = full`, `family = primary`, column `p_holm` |
| Beta LF +0.0107 (+0.0076 to +0.0139), d 0.61, 68% above zero; theta HF +0.0018 | same rows, columns `mean`, `ci95_low`, `ci95_high`, `cohen_d`, `share_positive` |
| Wilcoxon agreement | same rows, `p_wilcoxon` |
| Unchanged without pilot participants | `cohort_stats/group_tests_without_pilot_participants.tsv` |
| Span contrast, all Holm p ≥ 0.58 | `cohort_stats/group_tests.tsv`, H2 rows; simulation in `docs/ROBUSTNESS_RESULT_v1.md` |

### Verdict 2

| Claim | Source (row keys) |
| --- | --- |
| Synthetic phase calibrated LF bias theta +0.0045 (+0.0016 to +0.0074), alpha +0.0027, beta +0.0059; HF +0.0001 to +0.0009; p values | `null_validation/null_validation_cells.tsv`, `scenario = independent`, `variant = full`, calibrated statistic |
| Shared trial scenario +0.0136 to +0.0159 | same file, `scenario = shared_trial` |
| Observed phase calibrated theta LF +0.0043, alpha LF +0.0059, beta LF +0.0099 | `null_calibration/check2_H1_calibrated.tsv`, `calibrated_excess` |
| Realism: R² real 0.07 to 0.14 against 0.14 to 0.17; variance modulation about 1.7 times | `trial_shuffle/realism_summary.tsv` |
| Matched model phase bias theta +0.0026, alpha +0.0043, beta +0.0036 | `trial_shuffle/validate_cells.tsv`, `cohort = matched`, `null = phase` |
| Decomposition: no gain variability beta LF +0.0095; identical responses and stationary noise +0.0061 (+0.0027 to +0.0096) | `trial_shuffle/decompose_cells.tsv`, `cohort = no_gain`, `mean_only` |
| Estimator bias +0.0007 to +0.0011; 4 of 16 surrogate cohorts | `null_calibration/check4_surrogate_cohort_tests.tsv` |
| Circular shift null means −0.0007 to +0.0003; calibrated +0.0062 to +0.0111 | `joint_calibration/joint_calibrated_H1.tsv`, `null = shift` |
| Cardiac field 97% to 109% remaining | `robustness/cfa_vs_frozen_H1.tsv`, `share_of_frozen_effect_remaining` |
| ECG derived respiration coherence median 0.46 | `robustness/rsa_check.tsv` |
| Trial shuffle synthetic bias, smallest Holm p 0.57 and 0.43; upper bias limits | `trial_shuffle/validate_cells.tsv` (`null = shuffle`); limits in `trial_shuffle/gate.json` |
| First pipeline under trial shuffle; beta LF +0.0099 (+0.0066 to +0.0131) | `trial_shuffle/real_shuffle_calibrated_H1.tsv`, `k_used = 16`, `variant = full` |
| ICA removed median 10.5 of 61 components; 7 muscle (0 to 23); 2 eye (2 to 6) | `logs/ica_v2_apply_picard_20261113_100203.log`; `docs/ICA_V2_RESULT_v1.md` |
| After ICA and after muscle proxy: all six cells; beta LF +0.0108 (+0.0076 to +0.0140) d 0.60; +0.0059 (+0.0031 to +0.0088) Holm p 0.0003 | `ica_v2/shuffle_calibrated_H1_picard.tsv`, `variant_eeg = v2` or `v2_emgadj`, `variant = full`, `k_used = 16` |
| Encoding and retention only (beta LF unchanged) | same file, `variant = encret`; `claims_check.tsv` |

### Verdict 3

| Claim | Source (row keys) |
| --- | --- |
| Planted effects about 0.0066, 0.013, 0.027; detection per cell | `symmetric/synthetic_cells.tsv`, `cohort = plant_low`, `plant_mid`, `plant_high`, columns `mean`, `ci95_low`, `detected_by_rule` |
| Low level missed in alpha LF (0.0041 against 0.0047) and beta LF (0.0041 against 0.0046) under trial shuffle | same file, `null = shuffle` |
| Low level detected in beta LF under cross participant by less than 0.0001 | same file, `null = cross`; `claims_check.tsv` (0.000014) |
| Cross participant synthetic bias, smallest Holm p 0.76 and 0.61; limits | `symmetric/synthetic_cells.tsv`, `cohort = independent` and `matched`; `symmetric/gate_cross.json` |
| Real data cross participant: beta LF +0.0102, +0.0115, +0.0072 (+0.0045 to +0.0100) | `symmetric/real_cross_calibrated_H1.tsv`, `cohort = v1`, `v2`, `v2_emgadj`, `k_used = all` |
| Margins 0.0003 (K = 16) and 0.0002 (K = 8); trial shuffle shortfall 0.0015 | `claims_check.tsv` (0.000334, 0.000244, −0.001492) |
| Injection: added over real, trial shuffle, theta LF +0.0138, alpha LF +0.0075, beta LF +0.0040 (+0.0013 to +0.0067); HF +0.0137 to +0.0172 | `symmetric/inject_cells.tsv`, columns `delta_vs_real`, `delta_ci95_low`, `delta_ci95_high` |
| Detection after injection in all cells | same file, `detected_by_rule` |
