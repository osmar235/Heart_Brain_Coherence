# One Dataset, Three Verdicts: Why Positive and Negative Controls Matter for Brain–Heart Coherence

**Author:** Osmar Pinto Neto¹,*

¹ Department of Kinesiology, California State University San Marcos, San Marcos, CA 92096, USA.

*Corresponding author:* Osmar Pinto Neto (osmar@csusm.edu).

**Draft status:** version 3.2, 2026-11-17, prepared for PNAS Nexus (decision D-024). Framework approved by O.P.N. (decision D-019); authorship decision D-021. All numbers come from the result tables listed in the results index of the code release.

## Abstract

**Background.** Coherence between heart period and electroencephalographic (EEG) amplitude envelopes is used as evidence of brain and heart coupling. A statistically significant coherence is not sufficient evidence of physiological coupling: negative controls must challenge plausible alternative explanations, and positive controls must establish what the procedure can detect. The purpose of this study was to determine whether an apparently robust slow coherence between heart period and EEG envelopes during working memory survives negative controls designed to explain it without coupling, and whether those controls could have detected a coupling of the observed size.

**Methods.** We reanalyzed electrocardiography and 63 channel EEG from 122 adults performing an adaptive digit span task (OpenNeuro ds008302). Squared wavelet coherence between heart period and frontal theta, posterior alpha and left centroparietal beta envelopes was estimated in the low frequency (LF, 0.04 to 0.15 Hz) and high frequency (HF, 0.15 to 0.40 Hz) ranges. Negative controls were synthetic participants without coupling built on each participant's real trial timing, three surrogate nulls, a cross participant null, and artifact controls (cardiac field, respiration, independent component analysis, ICA, and a muscle proxy). Positive controls were planted couplings of three sizes in synthetic and in real data. All stages after the first were post hoc.

**Results.** The initial analysis found excess coherence in all six cells (Holm adjusted p from 3 × 10⁻⁹ to 0.044). Synthetic participants without coupling reproduced LF effects of similar size, and usual controls did not remove them. Under a trial shuffle null and after ICA and muscle proxy adjustment, no cell met all criteria. The positive controls then showed that a planted effect of about 0.013 was detected in every cell, whereas one of about 0.007 was missed in two LF cells. Frontal theta and HF effects were absent where such effects were detectable. Left centroparietal beta in LF exceeded the cross participant null in every pipeline (+0.0072, 95% CI +0.0045 to +0.0100 after muscle proxy adjustment) but not the trial shuffle criterion after that adjustment.

**Conclusions.** The same data supported three successive readings: coupling in every cell, no supported coupling, and one unresolved candidate. Negative controls removed most of the initial evidence, and positive controls showed where the absence of evidence was informative and where it was not.

## Introduction

Cardiac autonomic regulation and cortical activity are functionally related. The neurovisceral integration framework proposes that prefrontal networks supporting executive control contribute to the regulation of cardiac vagal activity \[1, 2\]. Heart period and cortical rhythms also share slow regulators, including arterial pressure rhythms near 0.1 Hz and fluctuations of arousal, so their covariation over tens of seconds could index how effort is regulated during cognitive work. Coherence between heart period and band limited EEG amplitude envelopes is one approach to measuring such covariation in the low frequency (LF, 0.04 to 0.15 Hz) and high frequency (HF, 0.15 to 0.40 Hz) ranges of heart rate variability (HRV) \[3, 4, 5\]. Relations between neural and cardiac rhythms have been reported at rest \[6\], and associations between HRV and EEG activity have been reported during working memory \[7\]. Coherence does not by itself identify a direct interaction between the two systems: it can also arise from a shared driver, from the organization of the task, or from artifacts.

Structured cognitive tasks make a positive result easy to obtain. Heart period and EEG envelopes both respond to events that repeat every 16 to 34 s (about 0.029 to 0.063 Hz), which overlaps the lower LF range. Subtracting the average response does not remove what varies from trial to trial \[8\], and task activations inflate connectivity estimates in functional MRI by a related mechanism \[9\]. Phase randomized surrogates, the usual null for coherence, may not reproduce the nonstationarity of task data \[10, 11\]. Blinks, eye movements and muscle activity are locked to task events and reach frontal, central and temporal electrodes. Earlier work from our group showed that a signal transformation can change coherence based inferences for the electromyogram, a different setting that illustrates why analysis choices need validation \[12, 13\].

Two kinds of control address this problem, and they answer different questions. A negative control asks whether the result can be produced without the effect of interest \[14\]. A positive control asks whether the procedure would have detected the effect had it been present. A claim has been tested only to the extent that the test could have contradicted it \[15, 16, 17\], and the same holds for a refutation: a nonsignificant result after a control is evidence of absence only for effects the procedure can detect \[18\]. Analysis decisions taken after seeing results add a further risk in either direction \[19\].

We examined these issues in a public dataset of verbal working memory with simultaneous electrocardiography (ECG) and EEG \[20, 21\]. The accompanying preprint recorded the ECG but did not analyze it \[20\]. Working memory load modulates frontal theta, posterior alpha and centroparietal beta activity in this paradigm \[20, 22, 23\]. The study began as a search for brain and heart coupling during the task. The initial analysis found effects in every cell tested, and the attempt to challenge that result, and then to challenge the challenge, became the subject of this paper. Each stage followed a plan written before that stage was run, but every stage after the first was designed after earlier results had been seen (Table 1); none was preregistered. The purpose of this study was to determine whether an apparently robust slow coherence between heart period and EEG envelopes during verbal working memory survives negative controls designed to explain it without coupling, and whether those controls could have detected a coupling of the observed size.

## Methods

### Dataset and task

We analyzed version 1.0.0 of the Cognitive overload and rest dataset (OpenNeuro ds008302; NEMAR on008302), released under CC0 \[21\]. The original study was approved by the ethics committee of Ural Federal University and all participants gave written informed consent \[20\]. Of 124 participants, 122 had a working memory recording (99 women, 23 men; age 20.5 ± 2.3 years, range 18 to 32). In each trial a cue (0.5 s) and a fixation period (2 s) preceded digits shown for 0.5 s every 2 s; the recall prompt appeared 4 s after the last digit and recall was typed without time limit \[20\]. Trials came in blocks of ten of the same length, adapted to performance until at least six blocks had been completed at each of two adjacent lengths \[20\]. Blocks lasted a median of 258 s, and the median time between trial starts was 25 s. EEG came from 63 scalp channels (actiCHamp, Brain Products) at 1000 Hz, with ECG on an auxiliary input of the same amplifier \[20\]. An earlier public dataset from the same group \[24\] is a different recording set (ds003838) and was not used.

### Heart period and EEG envelopes

R peaks were detected in ECG filtered from 5 to 25 Hz (slope energy, 300 ms refractory period, adaptive threshold), screened against each participant's median QRS, and intervals classified after Lipponen and Tarvainen \[25\]; abnormal intervals were marked and bridged in a 4 Hz NN series. Median heart rate was 76.2 bpm (interquartile range 71.1 to 84.6) and a median of 1.6% of intervals were marked (maximum 9.9%).

Regions and bands followed the original study \[20\] and define six cells (three regions by two heart period bands): frontal midline theta (4 to 8 Hz; AFz, Fz, F1, F2), posterior alpha (9 to 14 Hz; P1 to P6, PO3, PO4) and left centroparietal beta (16 to 22 Hz; C3, C1, CP3, CP5, C5, CP1). In the first pipeline, EEG was decimated to 200 Hz and referenced to the median of all channels, without independent component analysis. Each channel was band pass filtered (fourth order Butterworth, zero phase), its Hilbert amplitude taken, averaged over the region and decimated to 4 Hz; samples above the median plus 8 robust standard deviations were interpolated. The second pipeline is described under Artifact controls.

### Wavelet coherence and the estimand

We used the wavelet coherence functions of our EMG simulation framework \[13\]: Morlet transform (ω₀ = 6) and squared coherence after Grinsted et al. \[26\], smoothed after Torrence and Webster \[27, 28\], with 12 voices per octave from 0.5 s. Within each block both series were detrended and z scored; values outside the cone of influence were excluded, and a band value was the unweighted mean over the time and scale samples of the band inside the cone, required to cover at least half of them.

Trial locked responses were removed first. For each sequence length, an encoding average aligned on trial start and a recall average aligned on the recall prompt (up to 12 s, ending at the next trial) were formed across all trials of that length, including the block analyzed, and subtracted from every trial. For each block we then computed the coherence between its own heart and EEG residuals and subtracted the mean coherence between its heart residual and the EEG residuals of every other block of the same length, each pair truncated to the shorter block. This excess was averaged over the whole block (primary) or only from trial start to recall prompt (secondary). Blocks were excluded for fewer than ten trials, a recording discontinuity, more than 10% of marked heart intervals or interpolated envelope samples, insufficient coverage inside the cone, or no partner block; no exclusion depended on a coherence value.

### Negative controls

#### Synthetic participants without coupling

For each real participant we simulated one synthetic participant on that participant's real block and trial timing, taken from the public events files. Heart period and three envelope like signals had no ongoing coupling. Each responded to every trial with an encoding response and a recall response lasting the real self paced recall period and scaled with its duration, shaped by a gamma kernel (time constants 1 to 3 s), with trial to trial gain variability (SD 0.4) and latency jitter (SD 1.0 s) \[8\], over background noise whose variance changed over minutes. The trial locked share of variance was 0.30 for heart and 0.20 for each envelope. In the independent scenario heart and EEG trial responses varied independently; in a shared trial scenario their gains were correlated (0.7) and latencies shared. A matched model multiplied each signal's evoked share by the ratio of the real to the synthetic median share of variance removed by the trial averages; no coherence estimate entered this step. Here, “without coupling” means no direct ongoing interaction between the simulated heart period and EEG envelope processes; it does not exclude common task driven responses, particularly in the shared trial scenario. To locate the source of the bias, model components were switched off one at a time in 61 participants: gain variability, latency jitter, noise nonstationarity, recall duration scaling, and all four. Each synthetic participant was analyzed exactly as the real data.

#### Null models

All nulls were applied per participant and calibrated the same way: the participant's observed excess minus the mean excess of its surrogates, tested across participants.

1. **Phase randomized.** In each block the trial locked model was kept (EEG trial gains shuffled) and the residuals of heart period and EEG were phase randomized independently by the discrete Fourier transform (magnitudes kept, no amplitude adjustment) \[10\]; stationary within a block.
2. **Circular shift.** Heart period unchanged; the EEG residual rotated within each block by at least 20 s from either end, which keeps each signal's own nonstationarity \[11\].
3. **Trial shuffle.** Heart period unchanged; within each block the EEG envelope of every trial (start to next start) was replaced by that of another trial of the same block, chosen by a random derangement and warped piecewise linearly so that trial start, recall prompt and end matched. This keeps each signal's responses to every trial and their timing, and breaks the pairing of trial specific fluctuations. It counts covariation of trial responses as signal.
4. **Cross participant.** Heart period unchanged; for each surrogate a partner participant was drawn at random, and every complete block received the EEG envelope of a partner block with the same sequence length (the nearest length when the partner had none; 70% of blocks had an exact match in the real data), warped piecewise linearly so that block start, every trial start, every recall prompt and block end matched. This keeps the participant's task timing and replaces the EEG by that of another person doing the same task. It removes any relation specific to the person, including covariation of a person's own trial responses, which the trial shuffle null removes only in part.

The phase randomized null was evaluated in the independent synthetic cohort; the trial shuffle and cross participant nulls were evaluated in the independent and matched cohorts (8 surrogates per participant). A null was retained if no cell of its calibrated synthetic bias differed from zero after Holm correction over six cells. This criterion tests bias within the simulated scenarios; failure to reject zero bias is not proof of calibration under other data generating mechanisms. For each retained null, the upper bias limit of a cell was the larger upper 95% limit of the synthetic bias in the two cohorts. A Monte Carlo error rate estimate from one synthetic cohort of the size of the real one was abandoned as a criterion, because the uncertainty of the estimated bias was comparable with the standard error of the observed mean. The circular shift null was applied to the real data only and was not evaluated on synthetic data.

#### Artifact controls

- **Cardiac field artifact** \[29\]: an R locked moving average template (80 beats, −0.10 to +0.50 s) subtracted from each channel.
- **ECG derived respiration** \[30\]: R wave amplitude, 0.08 to 0.6 Hz, regressed out at lags of −2 to +2 s within each block.
- **Independent component analysis** (second pipeline, following the original study \[20\] with automatic criteria): 250 Hz, notch at 50 and 100 Hz, 1 to 100 Hz, bad channels marked by fixed criteria and interpolated, common average reference, extended infomax ICA (Picard solver) fitted on the task blocks; components classified with ICLabel \[31\], and every component most probably eye, muscle, heart, line noise or channel noise with probability of at least 0.8 removed. The same pipeline without removal was analyzed as a control.
- **Muscle proxy**: the 55 to 95 Hz amplitude envelope of the same channels after ICA, regressed out of each region's envelope within each block at lags of −2 to +2 s.

### Positive controls

The positive controls measured the detection sensitivity of the final procedure for one specified coupling pattern. They do not establish sensitivity to every possible pattern of physiological coupling.

**Planted coupling in synthetic data.** A shared ongoing process with equal power in the LF and HF ranges was mixed in phase into the heart period and the three envelopes of the matched synthetic model, carrying a share c of each signal's variance. Each participant had an independent process. A pilot in 13 synthetic participants (every tenth) related c to the added excess coherence; three levels were then fixed (c = 0.032, 0.045 and 0.064) to target added LF effects of 0.005, 0.010 and 0.020, the range of the observed effects. The rule for choosing the levels was amended once, after the pilot and before any other part was run, because the pilot grid was too coarse for the original rule. The three planted cohorts (122 participants each) were calibrated with the trial shuffle and the cross participant nulls (8 surrogates).

**Planted coupling in real data.** The middle level was mixed into the real heart period and envelopes after ICA, with an independent process per participant. Trial shuffle and cross participant surrogates used the seeds, and for the cross participant null the partners, of the corresponding analyses of the real data, so the comparison with the values before injection is paired (8 surrogates).

Whether the muscle proxy regression would also remove a true neural coupling cannot be tested by injection: a process uncorrelated with the proxy survives the regression by construction, and one added to the proxy is removed by construction.

### Statistics and decision rules

The participant was the unit of inference. For each of six cells (three regions by two bands), the mean excess across participants was evaluated by a two sided one sample t test with Holm correction across the six cells \[32\]. We report means, 95% confidence intervals (CIs), and Cohen's d. The following rules were written before the corresponding stage was run.

- **Separation rule.** Under the trial shuffle or the cross participant null, a cell exceeded the null only if it was Holm significant with 16 and with 8 surrogates and its lower 95% confidence limit exceeded the upper bias limit of that cell and null at both surrogate counts.
- **Muscle rule.** After ICA, centroparietal beta in LF was to be called possibly muscular if it met the separation rule under the trial shuffle null but not after the muscle proxy was regressed out.
- **Detection rule.** A planted cell counted as detected if it was Holm significant and its lower 95% confidence limit exceeded the upper bias limit of the null used.
- **Reading of the positive controls.** If a planted level was detected in a cell, the absence of an observed effect in that cell counts as evidence against a coupling of that size and pattern; if it was not detected, the absence is uninformative at that size.
- **Reading of the cross participant null.** Exceeding it shows covariation beyond what the task timing produces in any participant. It does not by itself show ongoing coupling, because covariation of a person's own trial responses also exceeds it.

A contrast between individual span and one digit beyond was specified in the first plan and is reported descriptively.

### Chronology of the analysis

Table 1 lists the analysis stages. Each plan was written before its stage was run. The first was an internal analysis plan, not a public preregistration. Every later stage was prompted by the results of earlier stages and is post hoc with respect to the original question.

| Stage | Plan written | What it tested | Main finding | Rules fixed before its results | Prompted by earlier results |
| --- | --- | --- | --- | --- | --- |
| 1. Primary analysis | 10 Nov 2026 | Excess coherence against a between block reference; span contrast | Excess in six of six cells | Yes | No |
| 2. Sensitivity checks | 11 Nov 2026 | Cardiac field, phase randomized calibration, respiration | Effects persisted | Yes | Yes |
| 3. Joint pipeline and synthetic validation | 11 Nov 2026 | All adjustments together, circular shift null, synthetic participants without coupling | Synthetic data without coupling reproduced the LF effects | Yes | Yes |
| 4. Trial shuffle null | 12 Nov 2026 | Realism and decomposition of the synthetic model, trial shuffle null | No detectable synthetic bias; only beta LF exceeded its limit | Yes | Yes |
| 5. ICA cleaning | 12 Nov 2026 (amended 13 Nov) | ICA, control without removal, muscle proxy | Theta effects absent; beta LF missed the rule after the muscle proxy | Yes | Yes |
| 6. Symmetric tests | 14 Nov 2026 | Cross participant null, planted coupling in synthetic and real data | Sensitivity limited near 0.007 in two LF cells; beta LF exceeded the cross participant null | Yes, with one amendment of the planted levels after a synthetic pilot | Yes |

*Table 1. Chronology of the analysis stages. The dated analysis plans and decision records are available from the author on request.*

### Reproducibility

The cohort analyses ran in Python 3.14.3 (NumPy 2.4.3, SciPy 1.17.1, MNE 1.13.2, mne-icalabel 0.10.0, python-picard 0.8.2) on Windows 11. The wavelet functions were vendored unmodified from the EMG framework with their checksum. Every random step used a fixed seed per participant and surrogate. The pipeline has 85 automated tests on synthetic signals. Each step is a numbered script that writes its tables and a dated log, and every run is recorded in an execution log.

### Use of generative AI

Claude (Anthropic; Claude Opus 5.5 **\[and earlier Claude Opus versions used in the project, to confirm]**) was used to write and test the analysis code, design the simulations, write the figure scripts, and draft and edit the manuscript text. ChatGPT (OpenAI; **\[model and version to confirm]**) was used to review the manuscript. All uses were directed by the author, who checked the code outputs, analyses, references and claims and takes responsibility for the content. The figures were drawn by the analysis scripts from the result tables; no AI image generator was used.

## Results

### Sample

All 122 recordings were processed. Participants completed a median of 18 blocks (range 11 to 30), and 89% (theta) to 91% (alpha, beta) of block and signal pairs were retained. Every participant contributed at least six blocks for alpha and beta, and 121 for theta.

### Verdict 1: the initial analysis found excess coherence in every cell

Within block coherence exceeded the between block reference in all six cells, each significant after Holm correction (adjusted p from 3 × 10⁻⁹ to 0.044; Figure 1A). The effect was largest for centroparietal beta in LF (+0.0107 of squared coherence, 95% CI +0.0076 to +0.0139; d = 0.61; 68% of participants above zero) and smallest for frontal theta in HF (+0.0018). Wilcoxon tests agreed, and excluding the three participants used to develop the pipeline did not change the result. No cell differed between span and one digit beyond (all Holm adjusted p ≥ 0.58); simulations showed that this contrast can produce spurious differences, so it is not evidence of equivalence.

> **Box 1. The paper we could have written.** *This box is a counterfactual. It shows the interpretation that the initial evidence could have supported had the analysis stopped at this point. It is not a finding of this study.* "Slow fluctuations of heart period were coherent with frontal theta, posterior alpha and centroparietal beta envelopes throughout a demanding working memory task, in a large sample and after removal of trial locked responses. The result agrees with the neurovisceral integration model, in which prefrontal networks that support executive control also regulate cardiac vagal output \[1, 2\]. Frontal midline theta, a marker of working memory demand \[23\], covaried with heart period in both frequency ranges, as expected if control networks and autonomic outflow are engaged together. The findings extend reports of coupling between neural and cardiac rhythms at rest \[6\] and of associations between HRV and EEG during working memory \[7\], and suggest that slow brain and heart coherence could serve as a marker of cognitive effort." Each sentence of this paragraph is compatible with the initial result and cites real literature. The remaining sections show why the paragraph is not supported.

### Verdict 2: negative controls explained most of the evidence without coupling

**Task timing reproduced the effects.** In synthetic participants without direct coupling, built on the real trial timing, the same analysis followed by phase randomized calibration left a positive bias in every LF cell: theta +0.0045 (95% CI +0.0016 to +0.0074), alpha +0.0027 and beta +0.0059, against +0.0001 to +0.0009 in HF (Figure 2A). Each synthetic LF bias differed from zero (theta p = 0.0025, alpha p = 0.029, beta p = 1 × 10⁻⁵). The observed phase calibrated values were of the same size (theta LF +0.0043, alpha LF +0.0059, beta LF +0.0099). With heart period and EEG trial gains correlated through a common task response, still without ongoing coupling, the calibrated LF excess reached +0.0136 to +0.0159. The real trial locked share of variance was smaller than in the model (median R² 0.07 to 0.14 against 0.14 to 0.17), and the real residual variance was more unevenly spread within trials (about 1.7 times the model's modulation); in a model matched to the real shares, phase randomized calibration remained biased in all three LF cells (theta +0.0026, alpha +0.0043, beta +0.0036, all intervals above zero).

**The bias came from self paced trial timing.** Switching model components off (Figure 2B) showed that the bias did not require trial to trial variability: without gain variability it grew (beta LF +0.0095), and with identical responses on every trial and stationary noise it remained (beta LF +0.0061, 95% CI +0.0027 to +0.0096). Self paced recall makes trials differ in length while trial averages have a fixed length, so even identical responses leave task locked residuals in both signals.

**Usual controls did not remove it.** The between block reference is part of the estimand and did not remove the bias. Phase randomized calibration removed a small estimator bias (+0.0007 to +0.0011) but not the timing bias. The circular shift null gave null means near zero (−0.0007 to +0.0003) and left all three LF cells significant (calibrated +0.0062 to +0.0111). Cardiac field subtraction left 97% to 109% of each effect. ECG derived respiration was coherent with heart period in HF (median 0.46), but regressing it out did not remove the HF effects.

**The trial shuffle null showed no detectable bias in the simulated cohorts.** In both synthetic cohorts the trial shuffle calibrated bias did not differ from zero in any cell (smallest Holm adjusted p = 0.57 in the independent model and 0.43 in the matched model; Figure 2A). The upper bias limits were +0.0028 (theta LF), +0.0047 (alpha LF), +0.0046 (beta LF) and +0.0020 to +0.0021 in HF. Under this null, all six cells of the first pipeline remained significant, but only beta LF exceeded its limit (+0.0099, 95% CI +0.0066 to +0.0131).

**Artifact removal eliminated the frontal theta effects and reduced beta LF.** ICA removed a median of 10.5 of 61 components per participant: 7 muscle (range 0 to 23) and 2 eye (2 to 6). After ICA (Figure 1B, Table 3), the frontal theta effects under the trial shuffle null were close to zero (LF from +0.0050 in the first pipeline to +0.0006; HF from +0.0017 to +0.0003). ICA removed several classes of components, so this reduction cannot be attributed to eye components alone. Alpha LF stayed near +0.005 and did not exceed its limit, and alpha HF and beta HF were no longer significant. Beta LF was unchanged after ICA (+0.0108, 95% CI +0.0076 to +0.0140; d = 0.60) and met the separation rule. After muscle proxy regression it declined to +0.0059 (95% CI +0.0031 to +0.0088; Holm adjusted p = 0.0003): the estimate remained above the upper bias limit (+0.0046), but its lower confidence limit did not, so it no longer met the rule and was classified as possibly muscular under the rule fixed in advance.

At this point no cell met every criterion. This result does not show that coupling is absent; whether it is evidence of absence depends on what the procedure could have detected.

### Verdict 3: positive controls showed where the negative result was informative

**A planted effect of about 0.013 was detected in every cell; one of about 0.007 was missed in two LF cells.** In synthetic data the three planted levels produced calibrated LF effects of about 0.0066, 0.013 and 0.027, somewhat above the targets (Figure 3). The middle and high levels were detected in all six cells under both nulls. The low level was detected under the trial shuffle null in theta LF and in the three HF cells (planted HF effects of 0.0042 to 0.0049), but not in alpha LF (lower confidence limit 0.0041 against a limit of 0.0047) or beta LF (0.0041 against 0.0046). Under the cross participant null the low level was detected in all six cells, in beta LF by a margin smaller than 0.0001.

**The cross participant null showed no detectable bias in the simulated cohorts, and in the real data only beta LF exceeded it.** No cell of the cross participant synthetic bias differed from zero (smallest Holm adjusted p = 0.76 and 0.61). Its upper bias limits were +0.0034 (theta LF), +0.0041 (alpha LF), +0.0042 (beta LF) and +0.0017 to +0.0019 in HF. In the real data (Figure 1C, Table 3), beta LF met the separation rule in the first pipeline (+0.0102, 95% CI +0.0071 to +0.0133), after ICA (+0.0115, +0.0084 to +0.0146) and after muscle proxy regression (+0.0072, +0.0045 to +0.0100). In the last case the lower confidence limit exceeded the limit by 0.0003 with 16 surrogates and by 0.0002 with 8, a margin within the uncertainty of the limit itself. No other cell met the rule in any pipeline. Values under the cross participant null were close to those under the trial shuffle null (beta LF after ICA 0.0115 against 0.0108; after the muscle proxy 0.0072 against 0.0059).

**Injection into the real data was detected in every cell, with an unexpected dependence on the existing effect.** After injection of the middle level into the real data, all six cells were detected under both nulls. In HF the added amount over the value before injection was similar across regions (+0.0137 to +0.0172). In LF it was smaller where the observed effect was larger: +0.0138 in theta (95% CI +0.0109 to +0.0166), +0.0075 in alpha (+0.0050 to +0.0101) and +0.0040 in beta (+0.0013 to +0.0067) under the trial shuffle null, with the same pattern under the cross participant null. Because beta LF met the rule before injection, its detection after injection is not informative by itself.

**Final status of the six cells.** Table 3 combines the negative and positive controls. Frontal theta LF and the three HF cells showed no supported effect in cells where the low planted level was detected. Alpha LF showed an effect near 0.005 that differed from zero and never exceeded a limit, at a size the trial shuffle rule did not detect in that cell. Beta LF remained an unresolved candidate.

| Cell | Initial analysis | Trial shuffle, ICA | Trial shuffle, ICA and muscle proxy | Limit, trial shuffle | Cross participant, ICA | Cross participant, ICA and muscle proxy | Limit, cross participant | Smallest planted effect detected | Status |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Frontal theta LF | +0.0054 (+0.0028 to +0.0079) | +0.0006 (−0.0020 to +0.0032) | +0.0007 (−0.0019 to +0.0033) | +0.0028 | +0.0012 (−0.0014 to +0.0037) | +0.0004 (−0.0021 to +0.0029) | +0.0034 | 0.0066 (low) | No supported effect; absence informative at the low planted size |
| Frontal theta HF | +0.0018 (+0.0001 to +0.0036) | +0.0003 (−0.0013 to +0.0020) | +0.0004 (−0.0014 to +0.0021) | +0.0020 | +0.0003 (−0.0012 to +0.0017) | +0.0001 (−0.0014 to +0.0016) | +0.0017 | 0.0042 (low) | No supported effect; absence informative at the low planted size |
| Posterior alpha LF | +0.0071 (+0.0044 to +0.0097) | +0.0046 (+0.0016 to +0.0075) | +0.0038 (+0.0008 to +0.0068) | +0.0047 | +0.0053 (+0.0025 to +0.0082) | +0.0045 (+0.0017 to +0.0074) | +0.0041 | 0.0145 (mid) | Neither supported nor excluded (below the detected size) |
| Posterior alpha HF | +0.0028 (+0.0014 to +0.0043) | +0.0011 (−0.0004 to +0.0027) | +0.0007 (−0.0009 to +0.0023) | +0.0020 | +0.0014 (−0.0001 to +0.0030) | +0.0013 (−0.0003 to +0.0029) | +0.0017 | 0.0042 (low) | No supported effect; absence informative at the low planted size |
| Centroparietal beta LF | +0.0107 (+0.0076 to +0.0139) | +0.0108 (+0.0076 to +0.0140) | +0.0059 (+0.0031 to +0.0088) | +0.0046 | +0.0115 (+0.0084 to +0.0146) | +0.0072 (+0.0045 to +0.0100) | +0.0042 | 0.0131 (mid) | Unresolved candidate: shuffle ICA, cross ICA, cross ICA+proxy met the rule |
| Centroparietal beta HF | +0.0027 (+0.0012 to +0.0042) | +0.0017 (+0.0003 to +0.0031) | +0.0009 (−0.0006 to +0.0023) | +0.0021 | +0.0021 (+0.0006 to +0.0036) | +0.0011 (−0.0004 to +0.0027) | +0.0019 | 0.0049 (low) | No supported effect; absence informative at the low planted size |

*Table 3. Evidence status of the six cells. Values are mean calibrated excess squared coherence with 95% CI (16 surrogates per participant; n = 121 for theta, 122 for alpha and beta). Limit: upper bias limit of the null in synthetic participants without coupling. "Smallest planted effect detected" refers to the trial shuffle detection rule in synthetic data, for the planted pattern tested. Generated from the result tables by the table script of the code release.*

| Test | Question | What it can exclude | What it cannot exclude | Outcome |
| --- | --- | --- | --- | --- |
| Between block reference | Does coherence exceed that of blocks paired at random? | Estimator bias common to all blocks | Structure locked to each block's own trial timing | Six of six cells positive |
| Phase randomized surrogates | Does coherence exceed that of signals with the same spectra and no phase relation? | Bias of smoothed coherence for stationary signals | Shared nonstationarity locked to the task | LF bias of +0.003 to +0.006 in synthetic data without coupling |
| Circular shift | Does coherence exceed that of the same signals misaligned in time? | Relations that need exact alignment | Not evaluated on synthetic data; moves trial locked structure away from the trials | Three LF cells positive |
| Synthetic participants without coupling | Can task timing alone produce the effect? | That a positive result implies coupling | That real signals behave as the model | Reproduced the LF effects |
| Trial shuffle null | Does coherence exceed that of the participant's own trials paired at random? | Timing bias (no detectable bias in two simulated scenarios) | Covariation of trial responses; bias in untested scenarios; effects near 0.007 in alpha LF and beta LF | Only beta LF exceeded its limit; not after the muscle proxy |
| Cross participant null | Does coherence exceed that with another person's EEG on the same task timing? | Timing bias and anything common to participants | Covariation of a person's own trial responses; bias in untested scenarios | Only beta LF exceeded its limit, in all three pipelines |
| Cardiac field subtraction | Is the effect the cardiac electric field in the EEG? | A cardiac field of the modeled shape | Residual field of other shape | 97% to 109% of each effect remained |
| ECG derived respiration | Is the HF effect respiratory? | Linear respiratory contribution visible in R wave amplitude | Respiration not captured by the ECG | HF effects remained; pattern unstable across analyses |
| ICA | Do the effects depend on artifact components? | Contributions of components classified as artifact | Artifact not isolated in a component; which class mattered | Theta effects absent; beta LF unchanged |
| Muscle proxy regression | Does beta LF depend on 55 to 95 Hz activity? | Contributions linearly related to the proxy | Whether the removed part is muscular or broadband neural | Beta LF fell from +0.0108 to +0.0059 |
| Planted coupling, synthetic | Would the procedure detect a coupling of known size? | Insensitivity at about 0.013 for the planted pattern | Sensitivity to other phase relations, bands or intermittent coupling | 0.013 detected in all cells; 0.007 missed in alpha LF and beta LF (trial shuffle) |
| Planted coupling, real data | Is a planted coupling detected with real noise and trial structure? | Loss of sensitivity due to real signal properties, in cells without an existing effect | Sensitivity where an effect already exists | Detected in all cells; added LF amount smaller where the observed effect was larger |

*Table 2. The tests, the question each one asks, what it can and cannot exclude, and its outcome.*

## Discussion

The purpose of this study was to determine whether an apparently robust slow coherence between heart period and EEG envelopes during verbal working memory survives negative controls designed to explain it without coupling, and whether those controls could have detected a coupling of the observed size. We analyzed 122 adults, challenged the initial result with synthetic participants without coupling, surrogate nulls and artifact controls, and then challenged those controls with planted couplings of known size. There were three main findings: (1) the initial analysis found significant excess coherence in all six cells; (2) negative controls showed that task timing and EEG artifacts could produce most of this evidence without coupling, and under a trial shuffle null with artifact controls no cell met all criteria; and (3) positive controls showed that this negative result was informative for frontal theta and for the HF cells, uninformative for posterior alpha in LF, and left one unresolved candidate, left centroparietal beta in LF. These three verdicts are successive levels of evaluation of the same evidence, not three equally valid conclusions.

### Negative controls: what explained the initial evidence

Removing the average trial response only removes what is identical across trials. Trial to trial variability leaves task locked structure in the residuals and inflates connectivity estimates \[8, 9\]. Our decomposition points to a simpler source in self paced tasks: even identical responses leave residuals when trials differ in length, because the averages have a fixed length. Both signals then carry residual variance at the same moments of each trial. A between block reference pairs blocks with different timing, and phase randomized surrogates spread residual variance evenly in time, so neither reproduces this structure. Nulls that keep each signal's trial structure did so in the scenarios we simulated. Smoothed wavelet coherence is also never zero for unrelated signals \[26\], which adds a small estimator bias that calibration removes.

The frontal theta effects disappeared after ICA. Blinks and eye movements around recall are plausible contributors, because their timing coincides with heart period responses, but several classes of components were removed together and the design does not isolate the eye components. HF results changed between analyses and are not interpreted; HF heart period is strongly influenced by respiratory sinus arrhythmia \[33\], and respiration was only inferred from the ECG.

### Positive controls: what the negative result means

A nonsignificant result after a control is evidence against an effect only at sizes the procedure can detect \[18\]. The planted couplings gave this information for one pattern, an in phase broadband process. For that pattern, the procedure detected effects of about 0.013 everywhere and effects of about 0.004 to 0.007 in theta LF and in HF. The absence of frontal theta and HF effects after ICA is therefore evidence against couplings of those sizes and that pattern in these cells. In alpha LF and beta LF, a planted effect of about 0.007 was missed under the trial shuffle rule. The alpha LF result (about 0.005 in every analysis) is thus neither supported nor excluded.

Stopping after the second verdict would have replaced one unsupported conclusion with another. The negative controls were added after a positive result and none would have been added after a negative one; the positive controls corrected this asymmetry only at the end. The injection result adds a caution of its own. In LF the injected process added less where the observed effect was larger. One possible explanation is partial cancellation in the cross spectrum between the in phase injected process and an existing relation with a different, consistent phase. We did not test this explanation, and the positive controls say nothing about couplings with other phase relations.

### The beta LF candidate

Left centroparietal beta in LF is an unresolved candidate, not a confirmed coupling. In its favor, it exceeded two nulls that showed no detectable bias in simulation, was unchanged by ICA, and had similar values under the trial shuffle and cross participant nulls, which suggests that covariation of trial responses contributes little to it. Against it, the muscle proxy halved it, and after that adjustment the two nulls disagreed by narrow margins: it exceeded the cross participant limit by 0.0002 to 0.0003 and fell short of the trial shuffle limit by 0.0015. The residual of 0.006 to 0.007 is at the size the trial shuffle rule did not detect in this cell, so that failure is uninformative. The muscle proxy cannot distinguish residual muscle activity from broadband neural activity shared by the beta and 55 to 95 Hz envelopes. The cell was one of six and was not predicted in advance. Finally, the region is contralateral to the right hand, and recall was typed; whether shared motor preparation of sensorimotor beta and heart period could produce the effect is a hypothesis we did not test. Restricting the analysis to encoding and retention left the first pipeline value unchanged, which argues against typing during recall as the whole explanation, although LF coherence is smoothed over tens of seconds. Our earlier study of auditory stimulation found lower coherence between LF heart rate variability and frontal and temporal beta envelopes during high pitch tones than during silence or low pitch tones \[4\]; the present reanalysis does not determine whether those contrasts were affected by the problems described here, and that study used continuous conditions without the self paced trial structure examined here.

### Methodological implications

Studies of brain and heart coherence during structured tasks should:

1. simulate signals without coupling on the real task timing and pass them through the complete analysis;
2. use nulls that keep each signal's trial structure, such as trial shuffling or pairing across participants on aligned timing, and report their bias in simulation with its uncertainty;
3. report effects against the bias that remains, not only against zero;
4. plant couplings of known size in synthetic and real data, report the smallest size detected in each cell, and interpret nonsignificant cells only against that size;
5. remove eye and muscle components, and measure muscle activity directly for regions near the scalp edge;
6. state which decisions were taken after seeing results.

### Limitations

- This is one dataset and one task. The study illustrates the use of controls; it does not estimate how often published coherence results would change under them.
- All stages after the first were post hoc. Plans written before each stage limit, but do not remove, the influence of earlier results on later choices, and they are not a preregistration.
- The nulls showed no detectable bias in two simulated scenarios with parameters largely fixed in advance. This does not establish their validity for other task structures or forms of dependence, and the upper bias limits depend on the model.
- The positive controls used one coupling pattern at three sizes. Sensitivity to lagged, band specific or intermittent coupling was not measured.
- The trial shuffle null counts covariation of trial responses as signal, and the cross participant null counts all person specific covariation as signal. An effect exceeding either could reflect a shared driver or task related responses, not a direct interaction between heart and brain.
- The first pipeline had no ICA; ICA used automatic classification without visual inspection; respiration and muscle activity were inferred, not measured; rest recordings in the dataset were not analyzed; sequence lengths visited once were dropped; the sample was young and mostly female.

### A proposed experiment to resolve the open questions

The remaining questions can be answered by a study designed for them:

1. **One cell specified in advance.** Left centroparietal beta in LF only, with the trial shuffle and cross participant nulls and the separation rule registered before access to the data.
2. **Recall without right hand typing.** Spoken or delayed recall, or typing with the left hand in half of the blocks. An effect that moves with the responding hand reflects motor preparation; one that stays does not.
3. **Measured muscle and respiration.** EMG over temporal and neck muscles and a respiration belt, so that both contributions are measured.
4. **Fixed trial durations in part of the session,** which removes the main source of timing bias and measures its contribution.
5. **Power for the remaining size.** Sample size and number of surrogates chosen so that an effect of 0.006 is detectable under the trial shuffle rule, using the planted couplings of this study for the power analysis.
6. **Phase as well as magnitude,** to test the cancellation suggested by the injection result and to plant couplings with the observed phase.

Two resources exist without new recordings: the rest recordings of this dataset, which have no trial timing, and an earlier dataset from the same group with a different design \[24\]. We reserve both for analyses with plans written before their coupling outcomes are examined.

### Conclusion

One dataset supported three successive readings of slow coherence between heart period and EEG envelopes during working memory. The initial analysis indicated coupling in every cell. Negative controls showed that task timing and artifacts could account for most of that evidence. Positive controls showed that the resulting absence of evidence was informative for frontal theta and for the HF range, uninformative for posterior alpha in LF, and that left centroparietal beta in LF remains a candidate that these data neither confirm nor exclude. A significant coherence was not sufficient evidence of physiological coupling, and a nonsignificant one after controls was not sufficient evidence of its absence.

## Data and code availability

The recordings are public under CC0 as OpenNeuro ds008302, version 1.0.0 (doi:10.18112/openneuro.ds008302.v1.0.0), with a NEMAR mirror (on008302; doi:10.82901/nemar.on008302). Analysis code, synthetic data generators, result tables, figure source data and figure scripts are available at https://github.com/osmar235/Heart_Brain_Coherence (code under the MIT license; documents, tables and figures under CC BY 4.0). The dated analysis plans, decision records and execution logs are available from the author on request. The release is archived on Zenodo (version 3.1.1, doi:10.5281/zenodo.24328792). The wavelet functions were copied without modification from the published EMG framework \[13\] (release v3.3.1 of https://github.com/osmar235/EMG_Simulator_WaveletCoherence; doi:10.5281/zenodo.24328375), and the release records their checksum.

## Ethics

This study reanalyzed public deidentified data. The original researchers reported approval by the ethics committee of Ural Federal University and written informed consent \[20\]. Institutional review board approval was not required for this secondary analysis of public deidentified data.

## Funding

This research received no external funding.

## Competing interests

Osmar Pinto Neto is a cofounder of Neurometra, a company that works on physiological measurement. Neurometra provided no funding or other support for this study.

## Author contributions

O.P.N. is the sole author and was responsible for conceptualization, methodology, software, validation, formal analysis, investigation, data curation, visualization, writing (original draft), writing (review and editing) and project administration.

## References

1. Thayer JF, Lane RD. A model of neurovisceral integration in emotion regulation and dysregulation. J Affect Disord. 2000;61(3):201-216. doi:10.1016/S0165-0327(00)00338-4.
2. Thayer JF, Hansen AL, Saus-Rose E, Johnsen BH. Heart rate variability, prefrontal neural function, and cognitive performance: the neurovisceral integration perspective on self-regulation, adaptation, and health. Ann Behav Med. 2009;37(2):141-153. doi:10.1007/s12160-009-9101-z.
3. Task Force of the European Society of Cardiology and the North American Society of Pacing and Electrophysiology. Heart rate variability: standards of measurement, physiological interpretation and clinical use. Circulation. 1996;93(5):1043-1065. doi:10.1161/01.CIR.93.5.1043.
4. von Jakitsch CB, Pinto Neto O, Pinho TOR, Ribeiro W, Pereira R, Baltatu OC, Osório RAL. High and low pitch sound stimuli effects on heart-brain coupling. Biomed Eng Lett. 2024;14(2):331-339. doi:10.1007/s13534-023-00340-5.
5. Pinto Neto O, Pinheiro AO, Pereira VL Jr, Pereira R, Baltatu OC, Campos LA. Morlet wavelet transforms of heart rate variability for autonomic nervous system activity. Appl Comput Harmon Anal. 2016;40(1):200-206. doi:10.1016/j.acha.2015.07.002.
6. Sargent KS, Martinez EL, Reed AC, Guha A, Bartholomew ME, Diehl CK, Chang CS, Salama S, Popov T, Thayer JF, Miller GA, Yee CM. Oscillatory coupling between neural and cardiac rhythms. Psychol Sci. 2024;35(5):517-528. doi:10.1177/09567976241235932.
7. Lee D, Kwon W, Heo J, Park JY. Associations between heart rate variability and brain activity during a working memory task: a preliminary electroencephalogram study on depression and anxiety disorder. Brain Sci. 2022;12(2):172. doi:10.3390/brainsci12020172.
8. Truccolo WA, Ding M, Knuth KH, Nakamura R, Bressler SL. Trial-to-trial variability of cortical evoked responses: implications for the analysis of functional connectivity. Clin Neurophysiol. 2002;113:206-226. doi:10.1016/S1388-2457(01)00739-8.
9. Cole MW, Ito T, Schultz D, Mill R, Chen R, Cocuzza C. Task activations produce spurious but systematic inflation of task functional connectivity estimates. NeuroImage. 2019;189:1-18. doi:10.1016/j.neuroimage.2018.12.054.
10. Theiler J, Eubank S, Longtin A, Galdrikian B, Farmer JD. Testing for nonlinearity in time series: the method of surrogate data. Physica D. 1992;58(1-4):77-94. doi:10.1016/0167-2789(92)90102-S.
11. Chavez M, Cazelles B. Detecting dynamic spatial correlation patterns with generalized wavelet coherence and non-stationary surrogate data. Sci Rep. 2019;9:7389. doi:10.1038/s41598-019-43571-2.
12. Pinto Neto O, Christou EA. Rectification of the EMG signal impairs the identification of oscillatory input to the muscle. J Neurophysiol. 2010;103(2):1093-1103. doi:10.1152/jn.00792.2009.
13. Pinto Neto O, Pinho TOR, Wang Y, Balbinot G, Kennedy DM. A computational framework for EMG simulation and coherence-based biomarker analysis on neural crosstalk in bimanual coordination. Comput Methods Programs Biomed. 2027;288:109657. doi:10.1016/j.cmpb.2026.109657.
14. Lipsitch M, Tchetgen Tchetgen E, Cohen T. Negative controls: a tool for detecting confounding and bias in observational studies. Epidemiology. 2010;21(3):383-388. doi:10.1097/EDE.0b013e3181d61eeb.
15. Popper KR. The Logic of Scientific Discovery. London: Hutchinson; 1959.
16. Platt JR. Strong inference. Science. 1964;146(3642):347-353. doi:10.1126/science.146.3642.347.
17. Mayo DG. Statistical Inference as Severe Testing: How to Get Beyond the Statistics Wars. Cambridge: Cambridge University Press; 2018. doi:10.1017/9781107286184.
18. Lakens D. Equivalence tests: a practical primer for t tests, correlations, and meta-analyses. Soc Psychol Personal Sci. 2017;8(4):355-362. doi:10.1177/1948550617697177.
19. Gelman A, Loken E. The statistical crisis in science. Am Sci. 2014;102(6):460-465. doi:10.1511/2014.111.460.
20. Kotyusov AI, Kasanov DA, Zaitseva AI, Kotchoubey B, Pavlov YG. One digit too many: brain rhythms and autonomic responses at the limits of working memory capacity. bioRxiv. 2026. doi:10.64898/2026.11.01.763481.
21. Pavlov YG. Cognitive overload and rest dataset. Simultaneous EEG + PPG + ECG + eye-tracking (version 1.0.0). OpenNeuro ds008302. 2026. doi:10.18112/openneuro.ds008302.v1.0.0. NEMAR on008302, doi:10.82901/nemar.on008302.
22. Kosachenko AI, Kasanov D, Kotyusov AI, Pavlov YG. EEG and pupillometric signatures of working memory overload. Psychophysiology. 2023;60(6):e14275. doi:10.1111/psyp.14275.
23. Pavlov YG, Kotchoubey B. Oscillatory brain activity and maintenance of verbal and visual working memory: a systematic review. Psychophysiology. 2022;59(5):e13735. doi:10.1111/psyp.13735.
24. Pavlov YG, Kasanov D, Kosachenko AI, Kotyusov AI, Busch NA. Pupillometry and electroencephalography in the digit span task. Sci Data. 2022;9:325. doi:10.1038/s41597-022-01414-2.
25. Lipponen JA, Tarvainen MP. A robust algorithm for heart rate variability time series artefact correction using novel beat classification. J Med Eng Technol. 2019;43(3):173-181. doi:10.1080/03091902.2019.1640306.
26. Grinsted A, Moore JC, Jevrejeva S. Application of the cross wavelet transform and wavelet coherence to geophysical time series. Nonlinear Process Geophys. 2004;11(5/6):561-566. doi:10.5194/npg-11-561-2004.
27. Torrence C, Compo GP. A practical guide to wavelet analysis. Bull Am Meteorol Soc. 1998;79(1):61-78. doi:10.1175/1520-0477(1998)079\<0061:APGTWA\>2.0.CO;2.
28. Torrence C, Webster PJ. Interdecadal changes in the ENSO-monsoon system. J Clim. 1999;12(8):2679-2690. doi:10.1175/1520-0442(1999)012\<2679:ICITEM\>2.0.CO;2.
29. Dirlich G, Vogl L, Plaschke M, Strian F. Cardiac field effects on the EEG. Electroencephalogr Clin Neurophysiol. 1997;102(4):307-315. doi:10.1016/S0013-4694(96)96506-2.
30. Moody GB, Mark RG, Zoccola A, Mantero S. Derivation of respiratory signals from multi-lead ECGs. In: Computers in Cardiology 1985. Vol 12. Washington, DC: IEEE Computer Society Press; 1985:113-116.
31. Pion-Tonachini L, Kreutz-Delgado K, Makeig S. ICLabel: an automated electroencephalographic independent component classifier, dataset, and website. NeuroImage. 2019;198:181-197. doi:10.1016/j.neuroimage.2019.05.026.
32. Holm S. A simple sequentially rejective multiple test procedure. Scand J Stat. 1979;6(2):65-70.
33. Grossman P, Taylor EW. Toward understanding respiratory sinus arrhythmia: relations to cardiac vagal tone, evolution and biobehavioral functions. Biol Psychol. 2007;74(2):263-285. doi:10.1016/j.biopsycho.2005.11.014.

## Figures

**Figure 1.** The three verdicts on one axis. A: initial analysis, first pipeline, excess coherence against the between block reference. B: trial shuffle calibrated excess after ICA and after ICA with the muscle proxy regressed out. C: cross participant calibrated excess for the same two pipelines. Points are means across participants with 95% CIs; dashed segments are the upper bias limits of each null in synthetic participants without coupling. File: figures/png/figure_1.png.

**Figure 2.** Negative control with synthetic participants. A: calibrated excess coherence in synthetic participants without coupling (122 per cohort, real trial timing) under phase randomized and trial shuffle calibration, with the observed phase calibrated values (open diamonds). B: phase calibrated synthetic bias in LF with model components switched off (61 participants, 4 surrogates each, same random draws in every version). File: figures/png/figure_2.png.

**Figure 3.** Positive controls. Calibrated excess coherence of couplings planted at three sizes in synthetic participants (122 per level, 8 surrogates) under the trial shuffle null (A) and the cross participant null (B), with the upper bias limits (dashed) and the observed values after ICA and after the muscle proxy. File: figures/png/figure_3.png.
