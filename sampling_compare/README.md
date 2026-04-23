# Sampling Comparison: H1, H2, H3, and H4

## What these hypotheses are really about

This project measures the relationship between pre-stimulus oscillatory phase and
behavioral or neural outcomes (MUA amplitude, LFP amplitude, reaction time, hit/miss).
The core quantity is the **preferred phase** — the LFP phase angle at which a given
dependent variable (DV) is systematically largest.

H1, H2 and H3 are not just different averaging recipes. They correspond to
**progressively weaker claims about where a shared optimal phase lives**:

| Hypothesis | Scientific claim |
| H1 | There is **one** optimal preferred phase, shared across all trials (i.e. across all stimulus positions and difficulty levels, and, after averaging, across channels and animals). |
| H2 | There is a phase-DV relationship **within each stimulus position**, but the preferred phase is **allowed to differ across positions**. |
| H3 | There is a phase-DV relationship **within each (position × difficulty level) cell**, with preferred phase allowed to differ across positions AND difficulty levels. |
| H4 | There is a phase-DV relationship **within each channel**, but the preferred phase is **allowed to differ across channels**. This is orthogonal to H1/H2/H3 (it stratifies on channel rather than on trial condition) and can be combined with any of them. |

Each relaxation of the hypothesis is implemented by deciding **at which level of the
pooling hierarchy to collapse complex numbers into magnitudes** (`abs()`). Once `abs()`
is taken, direction is discarded — everything averaged above that level is treated as
"is there SOME phase-DV relationship here?" without requiring directions to align.

---

## H1 — One shared optimal phase (`complex/`)

**What it does**

All trials across positions and difficulty levels are pooled together. For each channel
and frequency, the coherence is computed as the DV-amplitude-weighted mean phasor across
all trials:

    coh_complex = mean( DV .* exp(i * phase) )   [over all trials]

`abs()` is only taken at the very end, after all pooling. Channel and monkey averages
are also done in the complex domain (`mean(coh_complex)` across channels / animals)
before `abs()` or `angle()` is applied.

**The underlying hypothesis**

"**There is a single preferred phase that is the same across all stimulus positions,
difficulty levels, channels, and animals.**"

The pooled complex mean is large only if trials from different conditions all have
phasors pointing in roughly the same direction. If preferred phases differ across
conditions, the phasors cancel and the pooled `|coh_complex|` shrinks.

**What a result means**

- Significant H1 → there is a preferred phase that is consistent across everything
  pooled. Strong evidence for a shared optimal phase.
- Non-significant H1 → **ambiguous**: either (a) no phase-DV relationship exists, OR
  (b) relationships exist but the preferred phase is inconsistent across conditions
  so phasors cancel when pooled. H1 cannot distinguish these two cases.

---

## H2 — Per-position coherence, position-agnostic average (`abs_per_pos/`)

**What it does**

Trials are grouped by stimulus position (trialinfo column 16). Within each position, the
complex coherence is computed and then the **absolute value is taken immediately**:

    coh_pos(p) = abs( mean( DV_p .* exp(i * phase_p) ) )   [per position p]
    coh        = mean( coh_pos )                             [arithmetic mean of magnitudes]

Taking `abs()` per position collapses the direction information at that level — so when
we then average across positions, we're asking about the **strength** of each position's
coherence, not whether positions agree on a preferred phase. Channel and monkey averages
are then plain arithmetic means of these magnitudes.

**The underlying hypothesis**

"**Each stimulus position has its own phase-DV relationship. The preferred phase does
not need to be the same across positions; we only require that, on average, positions
individually show coherence.**"

We have relaxed the H1 claim: we no longer require positions to share a preferred phase,
only that each one produces a phase-DV relationship.

**What a result means**

- Significant H2 → each position on average shows a phase-DV relationship, but we have
  made **no claim** about whether their preferred phases align.
- Significant H2 with non-significant H1 → strong hint that preferred phases differ
  across positions: the per-position relationship is real, but pooling destroys it.
- Significant H2 and significant H1 → preferred phase is both present and consistent
  across conditions.

**Limitation: phase is lost**

Because `abs()` is applied before averaging across positions, the saved `coh` is a
real-valued magnitude. The channel and monkey averages carry no directional information.
Recovering a preferred phase requires going back to the raw data and using a different
strategy (see below).

---

## H2 preferred-phase variant (in `compare_preferred_phase_abs_per_pos.m`)

Because H2's saved output contains no phase, the compare script recomputes a preferred
phase directly from `ph_all_sess.mat` by **taking `angle()` per position instead of
`abs()`**, then circular-averaging directions across positions:

    coh_complex_p = mean( DV_p .* exp(i * phase_p) )   [complex per position]
    phi_p         = angle( coh_complex_p )               [preferred phase per position]
    phi_ch        = angle( mean( exp(i * phi_p) ) )      [circular mean across positions]

Conceptually this is still "H2-flavoured" — each position contributes independently,
with equal weight, regardless of that position's coherence magnitude. It answers a
different question than the H1 preferred phase: it asks **"if each position has its own
preferred phase, what's the average preferred direction across positions?"**, rather
than "what is the single preferred direction of the pooled data?"

---

## H3 — Per-(position × difficulty) coherence (proposed)

**Motivation**

H2 corrects for unequal sampling across stimulus positions but ignores the fact that
difficulty level (dE00, trialinfo column 18) varies within positions and is not
uniformly sampled. Difficulty level directly modulates hit rate and reaction time, so
if it is unevenly sampled within positions, H2's per-position estimates are still
contaminated by a difficulty bias — and importantly, preferred phase might also vary
with difficulty level.

**What it would do**

Trials are grouped by the joint condition (position × difficulty level). Within each
cell, complex coherence is computed, then `abs()` is taken (for magnitude, H3 analogue
of H2) or `angle()` is taken (for preferred phase):

    coh_cell(p,d) = abs( mean( DV_{p,d} .* exp(i * phase_{p,d}) ) )   [H3 magnitude]
    phi_cell(p,d) = angle( mean( DV_{p,d} .* exp(i * phase_{p,d}) ) ) [H3 preferred phase]

    coh    = mean over all (p,d) cells of coh_cell
    phi_ch = angle( mean( exp(i * phi_cell) ) )  over all (p,d) cells

Channel and monkey averages follow at successive levels.

**The underlying hypothesis**

"**Each (position × difficulty level) combination has its own phase-DV relationship. The
preferred phase is allowed to differ across positions AND difficulty levels; we only
require that each cell individually shows coherence.**"

This is the weakest of the three claims about shared optimal phase.

**What a result means**

- Significant H3 → every condition cell, on average, has a phase-DV relationship.
- The gap between H3 and H2 significance tests whether difficulty level is itself a
  source of preferred-phase variation.
- The gap between H2 and H1 significance tests whether position is a source.

**Trade-off**

H3 is the most stringent bias correction, but cells with very few trials (< 2) must be
excluded. If the difficulty staircase produces highly unequal cell counts, sparse cells
reduce the effective number of independent observations entering the average.

---

## H4 — Per-channel coherence (proposed)

**Motivation**

H1/H2/H3 all stratify along the **trial** dimension (stimulus conditions). They are
silent on a separate implicit assumption baked into H1: **that channels within an
animal share a single preferred phase**. In H1, the channel average is a complex
mean — `angle(mean(coh_complex across channels))` — which only produces a large
resultant if channels agree on their preferred direction. This amounts to assuming a
"single cortical sampling rhythm" across the electrode array.

This assumption is often implausible. Different channels may record from different
cortical layers, columns, or slightly different functional regions, each of which
could plausibly sample at its own preferred phase. H4 tests whether the data support
the per-channel relationship without requiring cross-channel alignment.

H4 is orthogonal to H1/H2/H3: it can be combined with any of them (e.g. "H2 at the
trial level, H4 at the channel level").

**What it would do**

Compute complex coherence per channel using any of the trial-level recipes (H1, H2, or
H3), then take `abs()` per channel before averaging across channels:

    coh_complex_ch = mean( DV .* exp(i * phase) )     [per channel, via H1/H2/H3 recipe]
    coh_ch         = abs( coh_complex_ch )              [magnitude per channel]
    coh_animal     = mean( coh_ch )                     [arithmetic mean over channels]

(For the per-channel preferred phase, the complex value `coh_complex_ch` is kept and
its angle extracted: `phi_ch = angle(coh_complex_ch)`. The animal-level preferred
phase is then a circular mean: `angle(mean(exp(i * phi_ch)))`.)

Note: the current H2 implementation already effectively does this at the channel level,
because `abs()` is taken so early that channels are never averaged in complex space.
So "H2 + H4" is what the existing `abs_per_pos/` scripts compute. What is new is
"H1 + H4" — pool all trials in complex space within a channel (the H1 recipe), then
take `abs()` per channel before averaging across channels. This isolates the
cross-channel claim from the cross-condition claim.

**The underlying hypothesis**

"**Each channel has its own phase-DV relationship. The preferred phase is allowed to
differ across channels; we only require that each channel individually shows
coherence.**"

This relaxes the implicit H1 assumption that the same oscillatory rhythm gates neural
sampling identically across the electrode array.

**What a result means**

- Significant H1+H4 with non-significant H1 → channels each have a phase-DV
  relationship but their preferred phases do not align. The "single cortical
  sampling rhythm" interpretation of H1 is not supported — the signal is there, but
  it's channel-specific.
- Significant H1 and significant H1+H4 → channels share a preferred phase (H1 is
  well-supported) AND each channel individually shows it. The strongest possible
  result.
- Non-significant H1+H4 → channels individually don't show phase-DV coherence; the
  apparent effect at the pooled-channel level is either absent or only emerges from
  cross-channel pooling (rare, usually indicates a statistical artifact).

**Trade-off**

H4 has lower sensitivity than H1 when the true preferred phase really is shared across
channels, because each channel's null distribution is wider (fewer trials per unit).
It also gives up the ability to talk about a "cortical" preferred phase — only channel-
specific ones remain.

---

## Summary

H1/H2/H3 form a nested hierarchy along the **trial** axis, progressively allowing the
preferred phase to vary across stimulus conditions:

    H1 ⊂ H2 ⊂ H3
    (strongest claim)      (weakest claim)
    one optimal phase      one per position      one per (position × difficulty)

H4 is on an **orthogonal axis** — it allows preferred phase to vary across channels —
and can be combined with any of H1/H2/H3.

| Hypothesis | Stratifies on | Where `abs()` is taken | Allows preferred phase to vary across... |
|---|---|---|---|
| H1 | Nothing (all pooled) | After pooling all trials & channels | Nothing — preferred phase assumed shared everywhere |
| H2 | Stimulus position | Within each position | Positions |
| H3 | Position × difficulty | Within each (position × difficulty) cell | Positions AND difficulty levels |
| H4 | Channel | Within each channel | Channels (orthogonal: combine with H1/H2/H3) |

A common workflow is to run multiple variants and use the pattern of significance to
diagnose *where* the preferred phase is stable vs. variable:

- H2 significant, H1 not → preferred phase varies by **position**
- H3 significant, H2 not → preferred phase varies by **difficulty level**
- H1+H4 significant, H1 not → preferred phase varies by **channel** (no single cortical rhythm)

---

## File organisation

    sampling_compare/
      complex/                  H1 comparison script
      abs_per_pos/              H2 comparison script
      README.md                 this file

Per-analysis results:

    results_{animal}/
      phase_coherence/complex/          H1 coherence
      phase_coherence/abs_per_pos/      H2 coherence
      phase_correlation/complex/        H1 correlation (ITC / circ-lin)
      phase_correlation/abs_per_pos/    H2 correlation
      multi_lin_reg/complex/            H1 regression
      multi_lin_reg/abs_per_pos/        H2 regression
      multi_lin_reg/cp10_till_100/      shared input data (ph_all_sess.mat)

    results_combined/
      phase_coherence/{complex,abs_per_pos}/         monkey-average coherence
      phase_correlation/{complex,abs_per_pos}/       monkey-average correlation
      multi_lin_reg/{complex,abs_per_pos}/           monkey-average regression

    Plots/
      phase_coherence/{complex,abs_per_pos}/{hermes,klecks,monkey_avg}/
      phase_correlation/{complex,abs_per_pos}/{hermes,klecks,monkey_avg}/
      multi_lin_reg/{complex,abs_per_pos}/{hermes,klecks,monkey_avg}/
      sampling_compare/                 preferred phase comparison figures
