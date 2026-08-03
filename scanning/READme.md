# `scanning/` — traveling-wave analyses of V4 LFP

**Hypothesis.** Attention scans the stimulus locations, so the optimal (preferred) phase should shift systematically across cortical space — a traveling wave whose origin moves with the stimulated retinotopic patch.

Two animals (`hermes`, `klecks`), 64-channel 8×8 V4 arrays, `cp10_till_100`, `dv = 'lfp'` unless stated. Animals are **never pooled at the channel level**; they are combined by replication (and, where implemented, by a pooled standardised test).

---

## Flow

```
  ph_all_sess.mat  (trial-level: phase_all, LFP/MUA amplitudes, trialinfo)
        │
        ├──────────────► phase_progression.m ──► phase_progression.mat
        │                 (SLURM: 1 job/channel)      pref_phase (ch × freq × pos)
        │                                             coh_mag    (ch × freq × pos)
        │                                                 │
        │                                                 ├──► cortical_planar_wave_PGD.m
        │                                                 ├──► cortical_wave_type_classification.m
        │                                                 ├──► cortical_planar_wave_derotation.m
        │                                                 └──► stimulus_loc_traveling_wave.m
        │
        └──────────────► functions/trial_position_sums_chan.m ──► trial_position_sums/
                          (SLURM: 1 job/channel)                  S(ch,pos,freq), W(ch)
                                                                  │
                                    used by the 'trial' ESTIMATOR in the two
                                    de-rotation scripts (cache shared by both)
```

`phase_progression.m` is the only script that touches raw trials in the main path. Everything else consumes `phase_progression.mat` — except the `trial` estimator, which goes back to the trials via the cached per-location sums.

Also read (not produced here): `phase_coherence/complex/.../coherence.mat` + `coh_perm_complex.mat`, used to build the per-channel coherence-significance mask that gates which channels enter each analysis.

---

## The scripts, and what each one answered

### 1. `phase_progression.m` — the producer, and test #1

**Question.** For each channel, does the preferred phase depend *systematically* on stimulus position?

**Method.** Per (channel × frequency): complex coherence per position, `c_p = mean over trials of y·exp(i·phase)`; `angle(c_p)` = preferred phase, `|c_p|` = coherence magnitude. Statistic = circular–linear correlation `r_cl` between position index (peripheral→foveal) and the preferred-phase vector. Null = shuffle position labels across trials.

**Critical design point.** The permutation null is **synchronised** — permutation *k* uses `rng(perm_seed_base + k)`, the same trial→position relabelling in every channel (`functions/phase_progression_chan.m`). Array LFP is spatially correlated; independent per-channel shuffles made the channel-average null far too tight and the channel-average significance anti-conservative. The old broadband "significance" was that artifact.

**Answer: WEAK.** After the fix, the channel-average is significant at only ~**2/35 frequencies (hermes)** and ~**3/35 (klecks)**, matching the flat per-channel picture. MUA / RT / hit_miss are weaker still than LFP.

Also produces `pref_phase` and `coh_mag`, which every other script below runs on.

### 2. `cortical_planar_wave_PGD.m` — test #2, existence of a planar wave

**Question.** Is there a planar traveling wave across the 8×8 array, and at which frequencies?

**Method.** Build the preferred-phase map on the array using only coherence-significant channels. `PGD = |mean(grad φ)| / mean(|grad φ|)` — 1 = gradient arrows all aligned (planar wave), 0 = random. Null = shuffle phases across electrode locations, then **cluster-corrected across frequency** so it reports a band, not isolated bins. Three views: `COLLAPSED` (positions merged), `PER-POSITION` (includes the wave-**origin** check against RF-driven channels, i.e. test #3), and `CONSENSUS` (collapsed-sig AND cross-position direction agreement AND sig in ≥ N positions).

**Answer: POSITIVE — and it replicates.** A planar wave exists in **both animals** at **theta (~4–6 Hz)** and **beta (~13–25 Hz)**. Beta speed ~28–51 cm/s (physically plausible); theta ~5 cm/s (borderline / near-synchronous). Phase tilt is gentle (~19–86° across the array). Propagation **direction differs between animals** (hermes ~10°, klecks ~166°).

**Answer to test #3 (origin), computed here: NULL.** The wave focus is not at the RF-driven patch (p ≈ 0.36–0.78) and does not track that patch across positions (p ≈ 0.42–0.90), in either animal.

> ⚠️ Look at the per-frequency speeds it reports: theta ≈ 5 cm/s, beta ≈ 28–51 cm/s. **Speed rising with frequency is the constant-wavenumber signature** — see §4. PGD is scale-free and evaluated one frequency at a time, so it cannot distinguish a propagating wave from a fixed spatial phase gradient.

### 3. `cortical_wave_type_classification.m` — what kind of pattern is it?

**Question.** Phase can be spatially organised in several ways. Which one is this?

**Method.** On the same collapsed map: **planar** (PGD), **radial** source/sink (net divergence = mean Laplacian), **rotating/spiral** (net curl + count of phase singularities: 2×2 loops whose wrapped phase sums to ±2π), **synchronous/none** (mean |grad| not above null). Works on **wrapped** phase with local neighbour differences and deliberately does *not* 2-D unwrap — unwrapping would destroy the very singularities a spiral is defined by.

**Answer: PLANAR.** PGD beats its null; radial divergence, rotational curl and singularity counts all sit at or below theirs. Not radial, not rotational, not spiral.

### 4. `cortical_planar_wave_derotation.m` — is the cortical wave dispersion-free?

**Question.** PGD says the gradients are aligned. Does the phase ramp also *scale with frequency*, as a real wave must?

**Method.** Rotate-to-overlap in cortical space. Each electrode sits at `r_c` mm on the array; a plane wave of speed `v`, direction `θ` predicts a lag `k·d_c(θ)`, `k = 2πf/v`, `d_c(θ) = x_c cosθ + y_c sinθ`. De-rotate and measure the resultant, swept over **frequency × speed × direction** (30 speeds 1–300 cm/s × 24 directions). Coherent across electrodes. Null = shuffle electrode ↔ array coordinate.

```
   a real wave:  one v explains all f  ->  k ∝ f  ->  HORIZONTAL band
   constant k:   fixed, frequency-independent phase offset  ->  DIAGONAL ridge
```

**Answer: NEGATIVE — diagonal ridge, no replication.** Absolute R never clears its threshold (R max ≈ thr, 0 significant cells in either animal). Gain is significant in each animal separately (23 cells hermes, 86 klecks) but **0 cells replicate**. The best-fit speed rises with frequency rather than staying flat.

> **PGD-positive and de-rotation-negative is not a contradiction.** PGD checks whether the gradient arrows point the same way, one frequency at a time. De-rotation additionally checks the *size* of the ramp and forces one speed to explain all frequencies. There is a spatial phase tilt; it is not dispersion-consistent.

Also note the baselines: `R0` ≈ **0.947 (hermes)** and **0.992 (klecks)**. Klecks' entire gain ceiling is `1 − 0.992 = 0.008`, so its 86 "significant" cells are statistically real but physiologically negligible — the electrodes are already ~99% phase-aligned.

### 5. `stimulus_loc_traveling_wave.m` — does the wave track the scanned stimulus?

**Question.** The actual scanning hypothesis: as the stimulus moves, does the preferred phase ramp with distance?

**Method.** Same rotate-to-overlap, but the distance axis is **stimulus eccentricity in degrees**, so speed is in deg/s. Three modes (full detail in [the appendix below](#appendix--the-three-modes-of-stimulus_loc_traveling_wavem)):

| mode | de-rotates by | says |
|---|---|---|
| `visual` | `k·d_p` | phase depends on where the stimulus is |
| `visual_coherent` | `k·(d_p + D_c)` | + a per-electrode delay fixed across positions (wave sweeping out from the fovea) |
| `visual_arrival` | `k·(d_p + a(c,p))` | the wave starts *where the stimulus lands* |

**Answer: NEGATIVE in all three modes.** Significance in hermes only, confined to **3.33–4.44 Hz** (the 3 lowest usable frequency bins), gains 1.0–2.7× threshold, with the significant region running to the top edge of the speed range. Klecks shows nothing anywhere. **0 cells replicate in any mode.**

The decisive number is the baseline: `R0` ≈ **0.80 (hermes)**, **0.93 (klecks)**. The per-location preferred phases are already 80–93 % aligned before any de-rotation, so there is almost no position-dependent phase for a wave model to organise. That is test #1's weak result, restated as a single number.

### 6. `traveling_wave_H2_H1.m` — deprecated

Descriptive 9-page PDF (phase heatmaps, phase-change maps, circular-variance maps, per-channel phase trajectories, Moran's I spatial autocorrelation, 2-D plane fit). Uses H2-H1 significance thresholds instead of coherence significance; its own header points you to the two `cortical_*` scripts. Untracked in git. Still useful for visualisation.

---

## Overall answer

> There **is** an intrinsic planar traveling wave across V4 cortex (theta + beta, replicated in both animals) — but it is **not locked to the scanned stimulus**, and it is not dispersion-consistent. Phase progression across positions is weak, the wave origin does not track the driven patch, and every de-rotation test produces a constant-wavenumber diagonal rather than a constant-speed band, with zero cross-animal replication.
>
> **Ongoing planar wave: yes. Scanning wave: unsupported.**

---

## The two estimators — two levels of "coherence"

Both de-rotation scripts take `ESTIMATOR = 'phase' | 'trial'`. The three modes say **what** is de-rotated; the estimator says **at which level the alignment is measured**. This distinction matters because the word "coherence" is used for two different things in this pipeline.

### Level 1 — coherence across TRIALS (the `phase_coherence/` sense)

Computed once, upstream, in `phase_progression.m`:

```
   c(channel, freq, position) = mean over TRIALS of  y_t · e^(i·phi_t)

       |c|      = coherence magnitude   -> saved as coh_mag
       angle(c) = preferred phase       -> saved as pref_phase
```

One complex number per (channel, frequency, position). This is the same construction as the `phase_coherence/` folder.

### Level 2 — alignment of those preferred phases ACROSS LOCATIONS

What the `phase` estimator does with them:

```
                | Σ_p  coh_mag · e^( i ( pref_phase − k·d_p ) ) |
   R_c(f,v)  =  ──────────────────────────────────────────────────
                            Σ_p  coh_mag
```

The level-1 coherence enters only as a **weight**; what is summed is the *preferred phases*. So R asks "do the 16 per-location preferred phases line up after de-rotation?"

```
   trials ──(level 1: coherence)──► one complex number per (ch, freq, position)
                                     │
                                     └──(level 2: resultant over positions)──► R
```

### The `trial` estimator collapses the two levels into one

Every trial is rotated by the model's prediction for the location *that trial* had, and one coherence is taken over the whole trial set:

```
   c(f,v) = ( 1/W ) · Σ over ALL TRIALS  y_t · e^( i ( phi_t − k·d_p(t) ) ) ,
   W      = Σ_t |y_t| ,        R = |c|
```

At `k = 0` this is **exactly the phase coherence of the `phase_coherence/` pipeline** with all locations pooled.

### What R, R0 and the gain mean under each

| | `phase` | `trial` |
|---|---|---|
| what R measures | how similar the per-location **preferred phases** are after de-rotation | the **phase coherence** over all trials after de-rotation |
| what R0 is (`k = 0`) | how similar they already were | the pooled phase coherence, no rotation |
| typical R0 | **0.80–0.93** (stimulus script), **0.947–0.992** (cortical) | an ordinary coherence, order 0.05–0.2 |
| gain ceiling `1 − R0` | tiny — 0.07 for klecks, 0.008 in the cortical script | large, but limited by how much position-dependence actually exists |
| units comparable across estimators? | **no** — never compare an R from one with an R from the other | |

`gain = R − R0` is computed identically for both. Nothing about the gain readout was removed or changed.

### Why `trial` is the better-behaved estimator

**1. No binning loss.** `phase` estimates 16 separate coherences from ~1/16 of the trials each, then fits to those noisy intermediates. `trial` fits the unbinned data directly.

**2. It fixes a weighting flaw.** `phase` weights location *p* by `|c_p|` — but under noise `|c_p| ~ 1/√n_p`, so a location with **fewer** trials gets a **larger** weight. Backwards. The two estimators differ by exactly this:

```
   phase :  uses  S(p,f) / n_p     (the per-location MEAN, then weighted by |c_p|)
   trial :  uses  S(p,f)           (the raw SUM — every trial counts once)
```

**3. An interpretable baseline.** `R0` becomes a coherence directly comparable with the `phase_coherence/` results, instead of a phase-similarity index pinned near 0.9 by construction.

**It will not manufacture an effect.** If preferred phase barely depends on location, both estimators say so. `trial` has more power to see a small one and a cleaner baseline — that is all.

### Implementation

The de-rotation depends on a trial only through its location, so the estimator collapses onto per-location complex **sums**:

```
   S(c,p,f) = Σ over trials at location p of  y_t · e^(i·phi_tf)
   c(f,v)   = ( 1/W ) · Σ_p  S(c,p,f) · e^(−i·k·d_p)
```

`functions/trial_position_sums_chan.m` computes `S` for the observed labels and for all `nPerm` shuffles — **one SLURM job per channel**, same `slurmfun` pattern and the same synchronised `rng(perm_seed_base + k)` seeding as `phase_progression_chan.m`, so a permutation is the same relabelling in every channel. All three modes then reduce to re-weightings of `S`, done in-script.

```
   cache: results_<animal>/scanning/trial_position_sums/cp10_till_100/<dv>/<ch>/

   stimulus_loc_traveling_wave.m   needs S_obs AND S_perm   (~286 MB/animal)
                                   — its null shuffles trial->location labels
   cortical_planar_wave_derotation.m  needs S_obs only      (~0.3 MB/animal)
                                   — its null shuffles electrode->array coordinate
```

The cache is shared: whichever script runs first pays for the jobs. `RECOMPUTE_TRIAL_SUMS = true` forces resubmission (do this after changing `dv`, `nPerm`, or the worker).

Combining animals (replication vs pooled) and how to read a frequency × speed grid are covered in the appendix below.

## Gotchas worth knowing

- **Coordinate frames.** `trialinfo` col 16/17 are **fixation-centred** pixels; `RF_Center_X/Y` in the RF summary tables are **screen** pixels. Differencing them without subtracting (840, 525) is off by the screen centre. `elec_rf_deg` handles this.
- **`ppd` is not in the repo.** `PIX_PER_DEG = struct('hermes',53.24,'klecks',50.56)`, from rig calibration. `sessInfo.ppd` does not exist in the saved `.RF` files despite being referenced elsewhere.
- **Extrapolated RF centres.** 18/64 (hermes) and 21/64 (klecks) rows are `Status = Extrapolated` — the Gaussian fit failed and the centre was filled from the median of 8-connected array neighbours, or a plane fit. `RF_VALID_ONLY` excludes them. For `visual_coherent` this is a *bias* concern (a smoothed `D_c` becomes a proxy for array position, which already carries the intrinsic wave); for `visual_arrival` it is mostly a *noise* concern.
- **Channel indexing is 1:1.** Both animals are 64-channel 8×8. Klecks' V4 is labelled V4-64…V4-127 in the raw RF data but is still 64 channels. No remap is needed between phase-progression index, RF-fit index, and array layout: `ch_col = ceil((1:64)'/8); ch_row = 8 - mod((1:64)'-1,8);`
- **NaN in permutation nulls.** `mean`/`std` do **not** omit NaN by default; one NaN permutation poisons that cell for all permutations, and NaN then spreads through the cross-animal sum. Both de-rotation scripts use `'omitnan'` and warn rather than silently producing `thr = NaN` or `-inf` (which would mark every cell significant).

---

## Where things land

```
  results_<animal>/scanning/phase_progression/cp10_till_100/<dv>/   phase_progression.mat
  results_<animal>/scanning/trial_position_sums/cp10_till_100/<dv>/ S(ch,pos,freq) cache
  results_combined/scanning/<analysis>/cp10_till_100/<dv>/          *.mat
  Plots/scanning/<analysis>/cp10_till_100/<dv>/                     *.pdf
```

Output names carry `[_trial]` for the estimator and `[_validRF]` for the RF filter, so all combinations coexist rather than overwriting.

---

# Appendix — the three modes of `stimulus_loc_traveling_wave.m`

All three ask the same question — *does the preferred phase carry a delay proportional to distance?* — and differ **only in which distance they de-rotate by**. Null, thresholds, figures and animal-combination logic are shared.

### The shared estimator

Each (channel, location) contributes a weighted phasor `w · e^(iφ)`. A wave of speed `v` predicts a lag `k · distance`, with `k = 2πf/v`. De-rotate by that predicted lag and measure how well the phasors line up:

```
   R  =  | SUM  w · e^( i ( phi  -  k * distance ) ) |  /  SUM w        0 <= R <= 1

   wrong v                        right (f, v)
   ↖  ↑  ↗  →                     ↗  ↗  ↗  ↗
   fanned out,  R small           aligned,  R -> 1
```

Sweep frequency × speed. **A real wave gives a horizontal band** (one speed explains all frequencies). A diagonal ridge means constant wavenumber = a fixed frequency-independent phase offset = not a wave.

---

### What the three modes de-rotate by

All three include `k · d_p` (stimulus eccentricity). They differ in what they add on top.

Illustration: 3 electrodes, 6 stimulus locations at eccentricities 0.5 … 4.2 deg.

```
MODE 1   visual              extra term:  none

               p=1    2      3      4      5      6
   ch 12       -      -      -      -      -      -
   ch 40       -      -      -      -      -      -        no per-electrode term
   ch 57       -      -      -      -      -      -


MODE 2   visual_coherent     extra term:  k * D_c      (D_c = RF eccentricity)

               p=1    2      3      4      5      6
   ch 12      2.1    2.1    2.1    2.1    2.1    2.1       ROWS ARE FLAT
   ch 40      0.9    0.9    0.9    0.9    0.9    0.9       no p in the formula
   ch 57      3.5    3.5    3.5    3.5    3.5    3.5


MODE 3   visual_arrival      extra term:  k * a(c,p)   (RF-to-stimulus distance)

               p=1    2      3      4      5      6
   ch 12      1.6    0.9    0.0    0.7    1.4    2.1       dips at p=3
   ch 40      0.4    0.3    1.2    1.9    2.6    3.3       dips at p=1
   ch 57      3.0    2.3    1.4    0.7    0.0    0.7       dips at p=5
                            ^                    ^
                     stimulus is ON that electrode's RF -> arrives first
```

(1-D numbers for clarity; the code uses the 2-D `hypot` of the RF-to-stimulus offset.)

---

### The physical picture behind Modes 2 and 3

```
MODE 2 — the wave always starts at the FOVEA and sweeps outward.
         Moving the stimulus does NOT change any electrode's delay.

      fovea
       (o) ~~~~~~~~~~> ~~~~~~~~~~> ~~~~~~~~~~>
                ch40         ch12         ch57
               D=0.9        D=2.1        D=3.5
        delay is a fixed property of the electrode


MODE 3 — the wave starts WHERE THE STIMULUS LANDS and spreads from there.
         Moving the stimulus changes EVERY electrode's delay.

    stimulus at p=3                    stimulus at p=5
          (*)                                    (*)
      <~~~ | ~~~>                            <~~~ | ~~~>
    ch40  ch12  ch57                      ch40  ch12  ch57
     far  NEAR   far                       far   far  NEAR
```

Mode 2 cannot express Mode 3's idea: a per-electrode constant structurally cannot depend on where the stimulus is.

---

### How channels are combined — and why

```
INCOHERENT  (visual, visual_arrival)     COHERENT  (visual_coherent)

  ch1 ↗↗↗↗  -> |.| = 0.8 ┐                ch1 ↗↗↗↗ ┐
  ch2 ↘↘↘↘  -> |.| = 0.8 ├ mean -> 0.8    ch2 ↘↘↘↘ ├ summed -> cancels -> 0.1
  ch3 →→→→  -> |.| = 0.8 ┘                ch3 →→→→ ┘

  magnitude FIRST, then average           sum FIRST, magnitude LAST
  each channel judged on its own          channels must AGREE with each other
  immune to per-channel phase offsets     needs a common phase reference
```

The rule that decides which one a mode needs:

```
   de-rotation(c,p)  =  [ row mean for channel c ]  +  [ variation across p ]
                          flat within the channel      dips / rises in the row
                          DELETED by per-channel |.|   SURVIVES per-channel |.|
```

- Mode 2's term is **all** row-mean → a per-channel `|.|` would erase it completely (you'd just recompute Mode 1) → **must** sum coherently.
- Mode 3's term **varies within the row** → survives the per-channel `|.|` → can stay incoherent, and keeps that robustness.

---

### Summary

| | `visual` | `visual_coherent` | `visual_arrival` |
|---|---|---|---|
| De-rotates by | `k·d_p` | `k·(d_p + D_c)` | `k·(d_p + a(c,p))` |
| Extra term shape | — | flat per channel | dips at that channel's RF |
| Channel combination | incoherent | **coherent** | incoherent |
| Uses RF centres | no | yes (`D_c`) | yes (`a(c,p)`) |
| Affected by `RF_VALID_ONLY` | no | yes | yes |
| Hypothesis | phase depends on stimulus position | one wave sweeping out from the fovea, stimulus-independent | the wave starts where the stimulus lands |
| Blind to | any per-electrode delay | — | the flat part of the retinotopic offset |
| Noise floor of R | ~0.22 (16 locations) | ~0.03 (all c×p terms) | ~0.22 |
| Interpreting a null | meaningful | **ambiguous** (offsets only lower R) | meaningful |
| Interpreting a hit | limited in meaning | trustworthy | trustworthy |

R values are **not comparable across modes** — the noise floors differ. Compare each mode only against its own permutation threshold.

---

### How this relates to `cortical_planar_wave_derotation.m`

Both scripts run the **same estimator** — de-rotate by `k · distance`, measure alignment, sweep frequency × speed — on **orthogonal axes of the same `φ(channel × position)` matrix**:

```
                    stimulus position p ──────►
                 ┌──────────────────────────┐
    channel c    │                          │
        │        │        φ(c, p)           │
        ▼        └──────────────────────────┘

   cortical script:  collapse →, fit the ramp ↓   (across electrodes)
   stimulus script:  fit the ramp →              (along the stimulus axis)
```

That much is a clean mirror. But "where the electrode is" can be measured in **two different spaces**, and that turns the pair into a 2×2:

- **cortical mm** — the electrode's physical position on the 8×8 array
- **visual degrees** — the eccentricity of its receptive field

```
                          electrode position measured in
                     VISUAL degrees        |   CORTICAL mm
                     (RF eccentricity)     |   (array coords)
   ──────────────────────────────────────────────────────────────────────
   delay is a        Mode B                |   cortical_planar_wave_
   per-electrode     visual_coherent       |   derotation.m
   CONSTANT                                |
   ──────────────────────────────────────────────────────────────────────
   delay depends     Mode C                |   ← MISSING
   on (electrode,    visual_arrival        |     (cortical arrival mode)
   stimulus)                               |
```

(Mode A sits outside the table — it has no electrode term at all.)

So **Mode B and the cortical de-rotation script are the same model with different rulers**: both de-rotate each electrode by a per-electrode constant and sum coherently. One uses `k·d_c(θ)` in mm with the propagation direction fitted; the other uses `k·D_c` in degrees of RF eccentricity.

#### Consequences

- **Don't port Modes A or B into the cortical script.** A would re-implement the stimulus script (the cortical script has no position axis); B already exists there, in its better cortical-ruler form with direction swept.
- **The empty cell is the only genuinely new test.** The cortical script *cannot* express the arrival model, because it discards the required axis on its first line: `z = sum(W .* exp(1i*PHI), 2)`. Once positions are collapsed, "the delay depends on which stimulus" is unsayable.
- **It would go in the stimulus script**, as a 4th mode with distances in cortical mm — that script already has the position axis, the RF centres in the right frame, and per-mode units (`Vsets` / `Vunit` are already per-mode cell arrays, so `V_CORTICAL` in `cm/s` drops straight in).
- **Dependency:** converting RF-to-stimulus separation into cortical mm needs a log-polar CMF calibration (`E2_DEG`, `MM_PER_LOGDEG`). These are **not in the repo** — they'd have to come from the V4 literature. A log-polar mapping is preferable to a "nearest driven patch" mapping because it doesn't depend on which channels are driven, so it sidesteps the saturation problem (hermes positions 10–16 drive zero channels).

#### One asymmetry that is easy to miss

The two scripts **collapse the other axis differently**, and it's a deliberate choice, not an oversight:

```
   cortical script:  collapses positions with a COMPLEX SUM      -> keeps phase
   stimulus Mode A:  collapses channels with MAGNITUDE-then-mean -> discards phase
```

Electrodes are assumed to share a clock, so their relative phase is meaningful and worth preserving. Mode A declines to assume that — which buys robustness to per-channel offsets at the cost of being blind to per-electrode delays.

#### Would the missing cell change anything?

Probably not decisive, for a reason that is ruler-independent:

```
    gain ceiling = 1 − R0 ,   and R0 is the k = 0 baseline,
                              which involves NO distance at all
```

`R0` is identical under any distance metric. A different ruler can only *redistribute* where a given phase pattern lands in the (f, v) grid — it can rescue an effect smeared across speeds by cortical magnification, but it cannot create phase variance the data doesn't contain. With `R0` ≈ 0.80 (hermes) and 0.93 (klecks), there is very little position-dependent phase to explain in either ruler.

---

### Combining animals: two criteria, side by side

Channels are never pooled across animals. Instead the last two columns of each figure report two different combination rules:

| | **replication** (col 3) | **pooled standardised z** (col 4) |
|---|---|---|
| What is shown | mean of the two animals' grids | mean of the two animals' *z-scored* grids |
| What the contour means | `sig(hermes) AND sig(klecks)` | pooled z above its own max-stat threshold |
| Is there a test on it? | no — the AND of two per-animal tests | yes — paired-permutation null on the pooled grid |
| Can one animal drive it? | **no** | **yes** |
| Power | low | higher (aggregates weak evidence) |
| Role | **primary** | secondary, always report which animal drives it |

Why standardise before pooling: the animals sit on very different scales (`R0` ≈ 0.80 vs 0.93, thresholds differ 4–6×), so a raw average is dominated by the animal with the bigger numbers. Each cell is converted to `z = (obs − null mean) / null sd` against that animal's own permutation distribution first.

The pooled null pairs permutation *b* of one animal with permutation *b* of the other — legitimate because the animals are independent, so the pairing is arbitrary and samples the product null.

```
   Zobs(f,v)    = mean over animals of z_a(f,v)
   Znull(f,v,b) = mean over animals of z_a_null(f,v,b)
   thr_pool     = 95th pctile of  max over grid of Znull
```

**With n = 2 neither rule supports a population-level inference** — between-animal variance isn't estimable. Replication is simply the more conservative descriptive rule. Note also that averaging roughly halves a one-animal effect, so a strong single-animal result can easily fail the pooled test too.

### Reading the output

Both PDFs are laid out **rows = modes, columns = hermes / klecks / mean+replication / pooled z**:

```
   Plots/scanning/phase_alignment_wave/cp10_till_100/lfp/
       phase_alignment_grids_validRF.pdf     absolute R
       phase_alignment_gain_validRF.pdf      gain dR = R - R0, with best-speed line

   results_combined/scanning/phase_alignment_wave/cp10_till_100/lfp/
       phase_alignment_validRF.mat
```

```
   speed                          speed
     ^                              ^
     |   ############               |              ###
     |                              |        ###
     |                              |   ###
     |                              |###
     +---------------> freq         +---------------> freq

     HORIZONTAL BAND                DIAGONAL RIDGE
     one speed fits all freqs       best speed rises with freq
     = REAL TRAVELLING WAVE         = constant wavenumber, NOT a wave
```

Support = horizontal band, plausible speed, **replicated in both animals**. The `_validRF` tag means RF-based modes used `Valid_Gaussian` channels only (~46/64 hermes, ~43/64 klecks); without the tag, interpolated `Extrapolated` centres were included too.

---

### Implementation pointers

| what | where |
|---|---|
| Mode list | `stimulus_loc_traveling_wave.m` → `metrics = {'visual','visual_coherent','visual_arrival'}` |
| Mode 1 estimator | `align_grid` |
| Mode 2 estimator | `align_grid_coherent` |
| Mode 3 estimator | `align_grid_arrival` |
| `a(c,p)` table | main loop: `a_arr = hypot(rf_deg(:,1) - stim_deg(:,1).', rf_deg(:,2) - stim_deg(:,2).')` |
| RF centres → degrees | `elec_rf_deg` (subtracts the screen centre 840/525 before `/ppd`) |
| Valid-only toggle | `RF_VALID_ONLY` |

Speed is shared between both terms in Mode 3 (`Vs = v`), so no extra free parameter is introduced — the grid stays frequency × speed.
