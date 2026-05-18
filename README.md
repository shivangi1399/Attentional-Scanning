# Attentional Scanning

Analysis code for studying **attentional sampling** using electrophysiological recordings (LFP and MUA) from macaque visual cortex area V4. The project investigates whether pre-stimulus oscillatory phase at critical time predicts behavioral outcomes (hit vs. miss) and neural response amplitudes, and how phase, amplitude, coherence, and reaction time relate to each other.

## Datasets

Two macaque monkeys performing an attentional-sampling task, recorded with 64-channel electrode arrays in area V4. Data types: LFP (Local Field Potential), MUA (Multi-Unit Activity), eye position, and pupil size.

## Project Structure

### `preprocessing/`
Semi-automated artifact rejection for LFP, MUA, and eye data (blink detection), followed by z-scoring across sessions for cross-session averaging.

### `LFP_MUA_data/`
Average ERP analysis (hit vs. miss condition differences) with permutation-based statistics and pixel-based multiple comparison correction. Time-frequency analysis (TFR) using multi-taper convolution and power spectrum analysis with FOOOF to isolate oscillatory peaks from the 1/f aperiodic component.

### `Phase_analysis/`
Pre-stimulus phase estimation using autoregressive (AR) models across 40 log-spaced frequencies (2-80 Hz). Includes inter-trial coherence (ITC) computation and critical time estimation (when the stimulus signal reaches the brain).

### `Correlation_analysis/`
Circular-linear correlation between pre-stimulus oscillatory phase and LFP/MUA ERP amplitude, computed per channel with permutation testing (1000 permutations). Analyses run separately for hits and misses across difficulty levels.

### `Coherence_analysis/`
LFP-MUA coherence computed per session using spectral analysis framework with permutation-based statistics.

### `Phase_coherence/`
Phase coherence across channels and travelling wave analysis. Generates visualizations of per-channel phase patterns across stimulus locations and frequency bands.

### `multiple_linear_reg/`
Multiple linear regression relating predictors (phase at critical time, MUA baseline, amplitude at critical time) to dependent variables (reaction time, MUA ERP, LFP ERP, hit/miss). Uses logistic regression for hit/miss classification. Permutation testing (1000 permutations) with max-stat thresholds for family-wise error rate control. Also includes reaction time computation via Engbert2003 saccade detection and pupil data processing.

### `sampling_compare/`
Cross-pipeline / cross-hypothesis framework that asks **at which levels (positions, difficulty bins, channels) the preferred phase is shared vs. varies**. Defines four hypotheses (H1, H2, H3, H1+H4) as "where in the pooling hierarchy is `abs()` taken" and tests them across the three independent metrics from `Phase_coherence/`, `Correlation_analysis/`, and `multiple_linear_reg/` (R² and R_phase).
See `sampling_compare/README.md` for the full framework, the Way-1/Way-2 pooling rules, and the matched-permutation logic.

### `scanning/`
Phase-progression analysis: per-channel test of whether preferred phase changes **systematically** with stimulus position (peripheral → foveal). Computes the Mardia/Jupp circular-linear correlation `r_cl` between position index and the per-position preferred phase. `r_cl` is a narrow test — it only catches approximately-monotonic progression. For the broader "do positions disagree on preferred phase in any way?" question, see `sampling_compare/compare_hypotheses_per_chan.m` (paired H2 − H1).

### `RF_Mapping/`
Receptive field mapping for V4 electrodes using high-gamma LFP power back-projected onto bar stimulus geometry.

## Dependencies

- [FieldTrip](https://www.fieldtriptoolbox.org/) - Electrophysiology analysis toolbox
- [CircStat2012a](https://github.com/circstat/circstat-matlab) - Circular statistics
- [FOOOF](https://fooof-tools.github.io/fooof/) - Parameterization of neural power spectra
- [slurmfun](https://github.com/esi-neuroscience/slurmfun) - SLURM-based parallel computation
