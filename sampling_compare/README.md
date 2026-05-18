# Sampling Comparison: H1, H2, H3, H4

## The core idea

We measure the relationship between pre-stimulus oscillatory **phase** and a
behavioural / neural outcome (DV: MUA amplitude, LFP amplitude, RT, hit/miss).
The quantity of interest is the **preferred phase** — the LFP phase angle at
which the DV is systematically largest.

A trial is described by several **levels** that could each be a source of
preferred-phase variation:

- **stimulus position** (which target was attended)
- **difficulty level** (the dE00 staircase value, binned into quantiles)
- **channel** (cortical recording site within the array)
- **animal**

Each hypothesis is a claim about **which levels share the same optimal phase**:

| Hypothesis | "Optimal phase is shared across…" | Folder              |
|------------|-----------------------------------|---------------------|
| H1         | …everything (positions, difficulties, channels, animals) | `complex/` |
| H2         | …everything **except positions**  | `abs_per_pos/`      |
| H3         | …everything **except positions AND difficulty levels** | `abs_per_pos_diff/` |
| H4         | …everything **except channels**   | `abs_per_chan/` (combined with H1) |

H1 ⊂ H2 ⊂ H3 form a nested chain along the **trial** axis (progressively
allowing more sources of variation across trial conditions).

**Why H1+H4, not just H4?**
H4 is not a standalone hypothesis — it only describes what happens at the
*channel* level ("channels don't need to agree on a preferred phase"). It says
nothing about how trials are pooled *within* a channel. That within-channel
recipe is still governed by H1, H2, or H3. So the full label is always
a combination:

| Label   | Within-channel trial recipe | Across-channel recipe |
|---------|-----------------------------|-----------------------|
| H1      | pool all trials (Way 1)     | pool all channels (Way 1) |
| H1+H4   | pool all trials (Way 1)     | abs per channel (Way 2)   |
| H2+H4   | abs per position (Way 2)    | abs per channel (Way 2)   |
| H3+H4   | abs per cell (Way 2)        | abs per channel (Way 2)   |

H4 sits on an **orthogonal** axis to H1/H2/H3 and can in principle be
combined with any of them. Only H1+H4 is currently implemented (`abs_per_chan/`)
because it is the most direct channel-level counterpart to H1 — same trial
pooling, only the channel aggregation relaxed.

H2+H4 and H3+H4 are not implemented because they are not needed for the core
diagnostic. The trial-axis comparisons (H1 vs H2, H2 vs H3) and the
channel-axis comparison (H1 vs H1+H4) are already orthogonal — each isolates
one source of variation independently. H2+H4 would only answer the
second-order question "given that positions already differ, do channels also
differ within each position?", which is not a primary hypothesis here.
Relaxing more levels also makes significance easier to reach, so H2+H4 /
H3+H4 would be weaker, less informative tests.

### How "relaxing a level" is done in code

Within each channel, we compute a complex resultant `mean(DV · exp(i·phase))`
and decide **at which level of the pooling hierarchy to take `abs()`**. Once
`abs()` is taken, direction is discarded. Anything aggregated *above* that
point is asking only "is there *some* phase-DV relationship here?", not "do
these things agree on a preferred phase?".

- **H1 (`complex/`)**: take `abs()` only at the very end, after pooling all
  trials, all channels, all animals. The pool is large only if all phasors
  point in the same direction → strong claim that everything shares an
  optimal phase.
- **H2 (`abs_per_pos/`)**: take `abs()` per stimulus position before
  averaging. Positions don't have to agree on a preferred phase; we only
  require each position to individually show a phase-DV relationship.
- **H3 (`abs_per_pos_diff/`)**: take `abs()` per (position × difficulty bin)
  cell. Now neither positions nor difficulty levels need to agree.
- **H1+H4 (`abs_per_chan/`)**: pool all trials in complex space within a
  channel (the H1 trial-level recipe), then take `abs()` per channel before
  averaging across channels. Channels don't have to agree on a preferred
  phase.

### Reading a result

A single hypothesis tells you whether a phase-DV relationship exists *under
its assumptions*. The diagnostic power comes from comparing pairs:

| Pattern                   | Interpretation                                  |
|---------------------------|-------------------------------------------------|
| H2 sig, H1 not            | Preferred phase varies across **positions**     |
| H3 sig, H2 not            | Preferred phase varies across **difficulty levels** |
| H1+H4 sig, H1 not         | Preferred phase varies across **channels** (no single cortical rhythm) |

> **Visual pattern ≠ formal test.** The table above describes what you
> might *see* in the H1-H4 overlay. It is **not** a valid significance
> test on its own: relaxing a level is structurally biased upward (Jensen
> for magnitudes, fitting flexibility for R²), so a per-hypothesis
> threshold being crossed at H_n but not at H_{n-1} can happen under the
> null. The proper test is the **paired permutation test** on the
> difference H_n − H_{n-1}, described in *Two different statistical
> questions* below.

> **Correlation caveat:** the H1 vs H1+H4 channel comparison is only
> meaningful for **coherence and regression**. `circ_corrcl` returns a
> non-negative scalar magnitude by construction — it does not produce a
> complex phasor that can be vector-summed across channels. This means the
> correlation pipeline is always doing Way 2 (magnitude per channel, then
> mean) regardless of which folder the script lives in, and
> `Correlation_analysis/complex/` and `Correlation_analysis/abs_per_chan/`
> produce numerically identical results. The `abs_per_chan` correlation
> scripts exist for structural symmetry but add no new information.

---

## What each script is doing, conceptually

There are three independent ways of measuring "phase predicts DV", run in
parallel as three pipelines:

| Pipeline                  | Measure                                             | What it captures |
|---------------------------|-----------------------------------------------------|------------------|
| `Phase_coherence/`        | `abs(mean(DV · exp(i·phase)))`                      | DV-amplitude-weighted resultant of unit phasors. Sensitive to monotonic phase-tuning. |
| `Correlation_analysis/`   | `circ_corrcl(phase, DV)` (or POS / ITC for hit-miss)| Circular-linear correlation. Catches non-sinusoidal tuning that coherence may miss. |
| `multiple_linear_reg/`    | DV ~ pup + MUA_baseline + amp + sin(φ) + cos(φ)     | Tests phase contribution **after partialling out** confounds (pupil, MUA baseline, amplitude). |

For each pipeline there are four folders — one per hypothesis — containing
analysis scripts that implement the recipe described above:

```
Phase_coherence/{complex, abs_per_pos, abs_per_pos_diff, abs_per_chan}/Coh_*.m
Correlation_analysis/{complex, abs_per_pos, abs_per_pos_diff, abs_per_chan}/Corr_*.m
multiple_linear_reg/{complex, abs_per_pos, abs_per_pos_diff, abs_per_chan}/regress_stats_R2_*.m
```

Each analysis script does the same five things:

1. **Loads pre-computed phase data** (`ph_all_sess.mat`) from
   `multi_lin_reg/cp10_till_100/` — this is the shared input across all
   hypotheses and pipelines.
2. **Computes the observed metric per channel** following the hypothesis's
   recipe (where to take `abs()`, what to stratify on).
3. **Runs a permutation null** by shuffling the DV / labels and recomputing
   the metric. Same shuffle is shared across channels so the cross-channel
   mean of perm distributions is a valid null. Permutation work is
   parallelised over channels via `slurmfun`.
4. **Aggregates** to channel-average (per animal) and monkey-average
   (across animals) using the same arithmetic mean of magnitudes that the
   recipe prescribes.
5. **Saves** per-channel and aggregate results under
   `results_{animal}/<pipeline>/<hypothesis>/cp10_till_100/...` and
   `results_combined/<pipeline>/<hypothesis>/cp10_till_100/...`.

### Comparison scripts (`sampling_compare/`)

These read the saved results and assemble cross-pipeline / cross-hypothesis
figures. None of them do new statistics — they only consume what the analysis
scripts wrote.

| Script                                               | What it produces |
|------------------------------------------------------|------------------|
| `<hypothesis>/compare_all_measures*.m`               | One figure per hypothesis: coherence + correlation + regression panels for all 4 DVs (per-animal and monkey-average), with permutation thresholds shaded. |
| `<hypothesis>/compare_preferred_phase*.m`            | Heatmaps of preferred phase (channels × frequency), polar histograms at top frequencies, pairwise circular correlation across DVs (per-animal and monkey-average). For H2 specifically, recomputes the preferred phase from `ph_all_sess.mat` because the saved H2 magnitudes don't carry direction. |
| `compare_hypotheses.m` (top level)                   | Overlays H1, H2, H3 magnitude curves on the same axes per pipeline / DV / animal, plus the monkey-average. This is the figure that directly answers "where does the preferred phase vary?" by comparing significance across the nested hypotheses. **Monkey-average scope only.** |
| `compare_hypotheses_per_chan.m` (top level)          | Same paired Jensen-corrected H2−H1 test as `compare_hypotheses.m`, but applied at **per-animal channel-average** and **per-channel within each animal** levels. Produces a 4 DV × 4 metric channel-avg grid, a 4×4 grid of "# sig channels vs freq", and an 8×8 per-channel grid PDF per (pipeline × DV). |
| `aggregate_regression_per_channel_nulls.m` (top level) | One-time aggregator. Consolidates the 1000 perm shards under `multi_lin_reg/<hyp>/cp10_till_100/perm_R[_pos]/<dv>/<ch>/perm_NNNN.mat` into a single `<ch>/per_channel_null.mat` per channel with fields `null_R2_phase` and `null_R_phase_mag`, so the per-channel regression paired test can do one load per channel instead of 1000. Idempotent (skips channels whose aggregate already exists). |

---

## How the three pipelines aggregate β-style information

All three pipelines (coherence, correlation, regression) end up with a per-
unit "arrow" — a complex resultant that has both a length (effect strength)
and a direction (preferred phase). Coherence's arrow is
`mean(DV·exp(i·phase))`; regression's is `β_complex = β_cos + i·β_sin`
recovered from the model fit. Combining arrows across a level (positions,
cells, channels, animals) can be done two ways:

- **Way 1 — vector-sum, then measure.** Add the arrows tip-to-tail; measure
  the final arrow's length / angle. Arrows pointing in **different**
  directions cancel → smaller resultant. Sensitive to whether the level
  agrees on direction.
- **Way 2 — measure each, then mean.** Measure each arrow's length on its
  own; arithmetic-mean the lengths (and circular-mean the angles). Direction
  doesn't matter once measured. Arrows can disagree freely.

The mapping onto hypotheses is exact: **at every level the hypothesis says
is allowed to disagree on optimal phase, the recipe uses Way 2.** Wherever
the hypothesis says "this level should agree", it uses Way 1.

| Level                          | H1     | H2     | H3     | H1+H4  |
|--------------------------------|--------|--------|--------|--------|
| Across trials (within stratum) | Way 1  | Way 1  | Way 1  | Way 1  |
| Across positions (per channel) | —      | **Way 2** | (subsumed by cells) | — |
| Across cells (per channel)     | —      | —      | **Way 2** | —    |
| Across channels                | **Way 1** | Way 2  | Way 2  | Way 2  |
| Across animals                 | **Way 1** | Way 2  | Way 2  | Way 2  |

All three pipelines now follow this table. For coherence and correlation
this is automatic — the metric *is* the magnitude of a complex resultant, so
the "Way 1 vs Way 2" choice is the "complex pool vs magnitude pool" choice
the README has been describing. For regression, the magnitude metric is R²
(which has no direction, so Way 2 is the only option) and the directional
metric is `R_phase = sqrt(β_sin² + β_cos²)` together with
`phi_pref = atan2(β_sin, β_cos)`. The four regression scripts now apply Way 1
or Way 2 at each level to match the table:

- `multi_lin_reg/abs_per_pos/regress_stats_R2_abs_per_pos.m` — Way 2 across
  positions for R_phase / phi_pref (magnitude per position, then mean of
  magnitudes; angle per position, then circular mean of angles). R² is
  also Way 2 across positions, which it already was.
- `multi_lin_reg/abs_per_pos_diff/regress_stats_R2_abs_per_pos_diff.m` —
  same but across (position × difficulty) cells.
- `multi_lin_reg/complex/regress_stats_R2.m` — Way 1 across channels for
  R_phase / phi_pref: vector-sum the per-channel β's, then take magnitude
  / angle. Way 1 across animals for the monkey-average.
- `multi_lin_reg/abs_per_chan/regress_stats_R2_abs_per_chan.m` — Way 2
  across channels: arithmetic mean of per-channel R_phase magnitudes;
  circular mean of per-channel phi_pref angles. Way 2 across animals.

R² recipes are unchanged at every level. Only `R_phase` and `phi_pref` are
affected.

---

## Two different statistical questions, two different tests

Two questions get asked of these analyses, and they need different machinery.

### Question A — "Is phase predictive of DV at all, under this hypothesis?"

This is the **per-hypothesis** test. Each pipeline already does it with the
standard recipe:

1. Compute the metric on real data (a single curve over frequencies).
2. Shuffle the DV; recompute; repeat ~1000 times under the same hypothesis.
3. Take the max-stat across frequencies of each shuffle to control for
   multiple comparisons; the 95th percentile is the threshold.
4. The observed curve is significant at any frequency where it exceeds the
   threshold.

This works well because the **null and the observed are on the same scale**
(same recipe, just shuffled labels). A high value is evidence that the data
diverges from chance under that recipe.

For regression specifically, the partial R² is the natural metric here:
it directly answers "after controlling for MUA and Amp, do sin/cos
predictors explain non-trivial variance?" The threshold is built from the
same partial R² statistic under shuffles → directly comparable.

### Question B — "Does relaxing the level reveal more structure than chance?"

This is the H1 → H2 → H3 → H1+H4 comparison. It's not the same as Question
A and **cannot be answered with per-hypothesis thresholds alone**.

The reason: relaxing a level is structurally biased upward, even under the
null.

- For magnitudes (coherence, correlation, regression `R_phase`): by
  Jensen's inequality, `mean(|x_i|) ≥ |mean(x_i)|`. Taking `abs()` earlier
  (Way 2) cannot return less than taking it later (Way 1), regardless of
  signal — a pure recipe artefact.
- For partial R²: each extra stratum is its own regression with its own
  betas. Smaller groups and more parameters per data point → more in-sample
  noise gets fit → mean R² across strata grows even when there is no
  signal.

So observed `H_n − H_{n-1} > 0` on its own means nothing. Some of it is
always there — the question is whether the gap exceeds what flexibility
alone gives you.

### The paired permutation test (`compare_hypotheses.m`)

Fix: build the null for the **difference**, not for each curve in
isolation. For each permutation index `i`, with the same shuffled DV:

```
H1_perm(i, f), H2_perm(i, f), H3_perm(i, f), H4_perm(i, f)
```

are all computed from the **same shuffle** (every analysis script seeds
with `rng(2025)` before generating perm indices, so the shuffles are
matched across hypotheses). The null distribution of

```
null_diff(perm, f) = H_n_perm(perm, f) − H_{n-1}_perm(perm, f)
```

captures **exactly** the Jensen / flexibility advantage that relaxing the
level provides under chance. The threshold is the 95th percentile of the
max-stat across frequencies; the observed `H_n − H_{n-1}` is significant
only where it exceeds that threshold.

A significant difference now has a clean interpretation:

| Significant pair | Interpretation                                        |
|------------------|-------------------------------------------------------|
| H2 − H1          | Positions disagree on preferred phase                 |
| H3 − H2          | Difficulty levels disagree (within position)          |
| H1+H4 − H1       | Channels disagree on preferred phase                  |

### What each pipeline persists for the paired test

The paired test needs `[nPerm × nFreq]` per-perm null curves at the
monkey-average level for each hypothesis. They are saved into
`monkey_avg_results.mat` under `perm_monkey_avg`:

| Pipeline   | Field used                       | What it is |
|------------|----------------------------------|-----------|
| Coherence  | `perm_monkey_avg`                | per-perm magnitude curve |
| Correlation| `perm_monkey_avg`                | per-perm magnitude curve |
| Regression | `perm_monkey_avg.phase`          | per-perm partial R² for sin+cos |
| Regression | `perm_monkey_avg.R_phase`        | per-perm `|complex β|` (built from full-model sin/cos betas saved by `regress_perm_R*.m`) |

For regression, the R_phase null curve is constructed using the same
Way-1 / Way-2 recipe as the observed curve at each hypothesis level — so
H1 uses complex mean across channels and animals, while H2/H3/H1+H4 use
mean of magnitudes.

`compare_hypotheses.m` loads these, plugs them into `perm{p,h,d}`, and
fires the paired test automatically. Old result files lacking the field
fall back to "no paired test" gracefully.

### Why R² alone wasn't enough

When the only question was Question A — "is phase predictive at all here?"
— per-hypothesis R² thresholds did the job because the observed and null
are on the same scale and built from the same recipe. The H1-H4 chart
adds Question B, which compares **across recipes**, and the Jensen /
flexibility bias means same-scale per-hypothesis thresholds aren't a fair
yardstick: a "more significant" R² at H2 vs H1 doesn't tell you positions
disagree, only that splitting was allowed. Paired differences are the
only way to ask "did relaxing this level reveal signal beyond what
relaxing alone buys you?".

### Per-animal and per-channel variants of the paired test

`compare_hypotheses.m` runs the paired Jensen-corrected H2 − H1 test
**only at the monkey-average level** — it averages the per-perm null
across channels and animals before testing. That's the right scope for
"across the cortical population, do positions disagree on preferred
phase?", but it can hide effects that exist only at a finer scope:

- **Dilution.** If position-dependence is real in, say, 8 of 64 channels
  and absent in the rest, the channel-mean H2 − H1 shrinks by ~8/64.
  The Jensen-corrected null shrinks too, but a sparse signal can still
  be buried.
- **Channel cancellation.** H1 vector-sums per-channel resultants
  (Way 1 across channels), so when channels prefer different phases the
  pooled H1 collapses. H2 averages per-position magnitudes (Way 2
  across channels) and doesn't cancel. Heterogeneity across channels
  inflates the H2 − H1 gap at the population level even when no
  individual channel has a clean position-disagreement.

`compare_hypotheses_per_chan.m` applies the same paired test at two
finer scopes, using only existing on-disk data — no new permutations:

1. **Per-animal channel-average.** Reads each animal's existing
   `channel_avg_results.mat` for H1 and H2 and runs the paired test on
   the channel-mean nulls. Cheaper than re-aggregating across animals,
   and lets you check whether one animal carries the population-level
   signal disproportionately.
2. **Per-channel within each animal.** Loads each channel's own H1 and
   H2 per-perm null files (matched by the shared `rng(2025)` seed) and
   runs the paired test channel by channel. The 8 × 8 grid output shows
   which channels (if any) reject "shared phase across positions" at
   the FWER-corrected level.

The matched-perm requirement still holds at both scopes — the upstream
pipelines all seed with `rng(2025)` before generating their perm
indices, so perm `i` at H1 is the same trial shuffle as perm `i` at H2
at every channel of every animal. The paired difference at any scope
isolates real position-disagreement from the Jensen / flexibility
advantage.

---

## File organisation

    sampling_compare/
      complex/                                      H1 comparison scripts
      abs_per_pos/                                  H2 (H2+H4) comparison scripts
      abs_per_pos_diff/                             H3 comparison scripts
      abs_per_chan/                                 H4 (H1+H4) comparison scripts
      compare_hypotheses.m                          Monkey-average overlay + paired test
      compare_hypotheses_per_chan.m                 Per-animal + per-channel paired test
      aggregate_regression_per_channel_nulls.m      One-time aggregator for regression
                                                    per-channel paired test
      README.md                                     this file

Per-analysis results:

    results_{animal}/
      phase_coherence/{complex, abs_per_pos, abs_per_pos_diff, abs_per_chan}/cp10_till_100/...
      phase_correlation/{complex, abs_per_pos, abs_per_pos_diff, abs_per_chan}/cp10_till_100/...
      multi_lin_reg/{complex, abs_per_pos, abs_per_pos_diff, abs_per_chan}/cp10_till_100/...
        + perm_R[_pos]/<dv>/<ch>/per_channel_null.mat   produced by
                                                        aggregate_regression_per_channel_nulls.m
      multi_lin_reg/cp10_till_100/                  shared input data
                                                    (ph_all_sess.mat,
                                                     frequency.mat) —
                                                    all compare scripts
                                                    load these from here

    results_combined/
      phase_coherence/{complex, abs_per_pos, abs_per_pos_diff, abs_per_chan}/cp10_till_100/...   monkey averages
      phase_correlation/{complex, abs_per_pos, abs_per_pos_diff, abs_per_chan}/cp10_till_100/... monkey averages
      multi_lin_reg/{complex, abs_per_pos, abs_per_pos_diff, abs_per_chan}/cp10_till_100/...     monkey averages

    Plots/
      phase_coherence/{...}/{hermes,klecks,monkey_avg}/cp10_till_100/...
      phase_correlation/{...}/{hermes,klecks,monkey_avg}/cp10_till_100/...
      multi_lin_reg/{...}/{hermes,klecks,monkey_avg}/cp10_till_100/...
      sampling_compare/                             per-hypothesis preferred-phase figures
      sampling_compare/hypotheses/monkey_avg/       monkey-avg H1/H2/H3 overlays
                                                    from compare_hypotheses.m
      sampling_compare/hypotheses/<animal>/         per-animal channel-avg + per-channel
                                                    paired H2−H1 outputs from
                                                    compare_hypotheses_per_chan.m
