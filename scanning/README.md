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
                                    used by the 'coherence' ESTIMATOR in the two
                                    de-rotation scripts (cache shared by both)
```

`phase_progression.m` is the only script that touches raw trials in the main path. Everything else consumes `phase_progression.mat` — except the `coherence` estimator, which goes back to the trials via the cached per-location sums.

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

**Answer: NEGATIVE — diagonal ridge, no replication.** Run of 2026-08-04, both estimators, `nPerm = 1000`, max-stat over the 35 × 30 × 24 grid.

| estimator | animal | R0 (median) | max R | thr | sig R | thr gain | sig gain |
|---|---|---|---|---|---|---|---|
| `phase` | hermes | 0.947 | 0.9910 | 0.9915 | **0** | 0.0147 | 23 |
| `phase` | klecks | 0.9916 | 0.9971 | 0.9972 | **0** | 0.0029 | 86 |
| `coherence` | hermes | 0.953 | 0.9942 | 0.9942 | 3 | 0.0112 | 35 |
| `coherence` | klecks | 0.9916 | 0.9971 | 0.9973 | **0** | 0.0030 | 87 |

**Replicated cells: 0 for R under both estimators; 0 gain cells under `phase`, 1/1050 under `coherence`.** Pooled: 256 (`phase`) and 343 (`coherence`) cells — sitting on the diagonal.

Three things make this a clean negative:

**1. The two animals are significant in disjoint frequency bands.** That is *why* nothing replicates — it is not a threshold accident:

```
   hermes  sig gain at  4-6 / 35 frequencies,  25.5 - 80.0 Hz
   klecks  sig gain at   12 / 35 frequencies,   5.6 - 25.5 Hz
```

Neither band matches the other, and neither matches PGD's theta + beta (§2).

**2. The best-fit speed rises with frequency, in proportion — the constant-wavenumber signature, quantified:**

```
   klecks   f  5.6 -> 25.5 Hz   (x4.6)      v  34.5 -> 202.4 cm/s   (x5.9)
   hermes   f 25.5 -> 80.0 Hz   (x3.1)      v  62.2 -> 166.3 cm/s   (x2.7)
```

`v ∝ f` means `k = 2πf/v` is **constant**. One fixed phase offset fits every frequency. A real wave is the opposite: one fixed *speed* fits every frequency, and `k` grows with `f`.

**3. The ceiling is severe, and here it cannot be engineered away.** `R0` medians are **0.947 (hermes)** and **0.9916 (klecks)** — electrodes are already ~95–99 % phase-aligned with nothing de-rotated. Klecks' entire gain ceiling is `1 − 0.9916 = 0.008`, so its 86–87 "significant" cells are statistically real and physiologically negligible.

> Note the contrast with §5: there, switching to the `coherence` estimator dropped `R0` from ~0.9 to ~0.07 and removed the ceiling objection. **It does not do that here** — 0.947 → 0.953 and 0.9916 → 0.9916. The reason is structural: this script's statistic is always a resultant *across electrodes*, `R = |Σ_c w_c e^(i(φ_c − k·d_c))| / Σ_c w_c`, and that normaliser divides `|Zmap|` back out. The estimator changes the weights, not the scale of `R`. The ceiling is a property of the analysis, not of the estimator.

> **PGD-positive and de-rotation-negative is not a contradiction** — see [the reconciliation below](#why-pgd-is-positive-but-every-de-rotation-test-is-negative).

### 5. `stimulus_loc_traveling_wave.m` — does the wave track the scanned stimulus?

**Question.** The actual scanning hypothesis: as the stimulus moves, does the preferred phase ramp with distance?

**Method.** Same rotate-to-overlap, but the distance axis is **stimulus eccentricity in degrees**, so speed is in deg/s. Three modes (full detail in [the appendix below](#appendix--the-three-modes-of-stimulus_loc_traveling_wavem)):

| mode | de-rotates by | says |
|---|---|---|
| `visual` | `k·d_p` | phase depends on where the stimulus is |
| `visual_coherent` | `k·(d_p + D_c)` | + a per-electrode delay fixed across positions (wave sweeping out from the fovea) |
| `visual_arrival` | `k·(d_p + a(c,p))` | the wave starts *where the stimulus lands* |

**Answer: NEGATIVE in all three modes, under BOTH estimators.** Run of 2026-08-04, `nPerm = 1000`, `alpha = 0.05`, max-stat corrected over the whole 35 × 30 = 1050-cell grid. Numbers below are `all RF centres` (the `_validRF` variant is in the note at the end of this section).

**Replication — the primary criterion — is ZERO everywhere.**

| estimator | mode | replicated R | replicated gain | pooled sig |
|---|---|---|---|---|
| `phase` | `visual` | **0** | **0** | 2 (thr 2.99) |
| `phase` | `visual_coherent` | **0** | **0** | 1 (thr 3.27) |
| `phase` | `visual_arrival` | **0** | **0** | 87 (thr 2.69) |
| `coherence` | `visual` | **0** | **0** | 40 (thr 3.45) |
| `coherence` | `visual_coherent` | **0** | **0** | 1 (thr 3.27) |
| `coherence` | `visual_arrival` | **0** | **0** | 395 (thr 2.80) |

Per animal (hermes / klecks), the gain never clears threshold in klecks at all, and in hermes only in the `phase` estimator:

| estimator | mode | sig-gain hermes | sig-gain klecks |
|---|---|---|---|
| `phase` | `visual` | 10 / 1050 | 0 |
| `phase` | `visual_coherent` | 14 / 1050 | 0 |
| `phase` | `visual_arrival` | 21 / 1050 | 0 |
| `coherence` | all three modes | **0** | **0** |

#### Why this is a strong negative and not a power problem

**1. The `phase` baseline is nearly saturated.** `R0` medians are **0.80 (hermes)** and **0.93 (klecks)** — the per-location preferred phases are already 80–93 % aligned before any de-rotation. Gain ceilings of 0.20 and 0.07. That is test #1's weak result restated as one number.

**2. The `coherence` estimator removes that objection, and still finds nothing.** Its `R0` is an ordinary phase coherence — medians **0.063 (hermes)**, **0.078 (klecks)** — leaving ample headroom. It does not use it:

| mode | max R0 | max R (best cell in the whole grid) |
|---|---|---|
| `visual` | 0.081 / 0.096 | 0.080 / 0.095 |
| `visual_coherent` | 0.077 / 0.096 | 0.076 / 0.084 |
| `visual_arrival` | 0.081 / 0.096 | 0.079 / 0.086 |

(hermes / klecks). **`max R ≤ max R0` in every mode and both animals** — the best de-rotation anywhere in the grid recovers no more coherence than doing nothing at all. Gain thresholds are 0.0003–0.0007, i.e. the null itself barely moves, and the observed gain does not reach even that.

**3. The pooled hits are the diagonal, not a wave.** `visual_arrival` produces the largest pooled counts (87 and 395 cells) and the largest z values (up to 7.4). Look at where they sit on the grid: a **diagonal ridge**, speed rising with frequency = constant wavenumber = a fixed, frequency-independent phase offset. A wave requires a horizontal band. This is the same signature as PGD's theta-5 / beta-28–51 cm/s split (§2) and the cortical de-rotation result (§4). Pooled significance confirms the ridge is not noise; it does not make it a wave. With `n = 2` a pooled cell can also be ~100 % one animal — the script prints the per-animal z at the peak cell for exactly this reason.

> ⚠️ **`_validRF` caveat.** In the 2026-08-04 run the `phase` estimator's `visual_coherent` and `visual_arrival` modes were invalid under `RF_VALID_ONLY = true`: dropping the Extrapolated RF centres leaves NaN geometry, which an unrestricted `randperm` shuffled into the used channel set, making every permutation NaN and collapsing the threshold to `-inf` (hence a spurious 1050/1050 "significant"). Fixed — those two helpers now permute only among finite-geometry channels, and all six helpers call `null_guard`. **Re-run needed before quoting `_validRF` numbers for those two modes.** Unaffected: everything in the tables above (all-RF), and the `coherence` estimator in both variants — its `_validRF` results agree with the all-RF ones (0 replicated, 0 sig-gain, 40 / 0 / 477 pooled).

### 6. `traveling_wave_H2_H1.m` — deprecated

Descriptive 9-page PDF (phase heatmaps, phase-change maps, circular-variance maps, per-channel phase trajectories, Moran's I spatial autocorrelation, 2-D plane fit). Uses H2-H1 significance thresholds instead of coherence significance; its own header points you to the two `cortical_*` scripts. Untracked in git. Still useful for visualisation.

---

## Overall answer

> There **is** a robust, replicated, planar spatial phase gradient across V4 (theta + beta, both animals). But it is **stationary, not propagating**, and it is **not organised by the scanned stimulus**. Every de-rotation test produces a constant-wavenumber diagonal rather than a constant-speed band, with essentially zero cross-animal replication.
>
> **Planar phase structure: yes. Traveling wave: no. Scanning wave: unsupported.**

Five tests, each asking a stricter question than the last:

| test | question | result |
|---|---|---|
| #1 phase progression | does preferred phase depend on stimulus position? | weak — 2/35 and 3/35 frequencies |
| #2 PGD | are the phase-gradient arrows aligned across cortex? | **positive, replicated** — theta + beta |
| #3 origin | is the focus at the RF-driven patch, and does it follow it? | null (p ≈ 0.36–0.90) |
| #4 cortical de-rotation | does the ramp scale with frequency (one `v` for all `f`)? | negative — 0 replicated R, 1/1050 gain cell |
| #5 stimulus de-rotation | does the phase ramp with stimulus distance? | negative — **0 replicated, 3 modes × 2 estimators** |

**Why #2 passes and #4 fails is not a contradiction** — PGD is scale-free and single-frequency, so a fixed phase offset scores perfectly on it. [Full reconciliation below.](#why-pgd-is-positive-but-every-de-rotation-test-is-negative)

**The ceiling objection is closed for #5, and structurally unclosable for #4.** The obvious complaint was that `R0` sat at 0.80–0.99, leaving nothing for a wave to win:

- **#5 (stimulus):** the `coherence` estimator drops `R0` to an ordinary phase coherence (0.06–0.08), so the headroom is real. De-rotation still buys nothing — `max R ≤ max R0` in every mode and both animals. **The negative result is not a baseline artifact.**
- **#4 (cortical):** the same switch does *not* lower `R0` (0.947 → 0.953, 0.9916 → 0.9916), because that statistic is always a resultant across electrodes and the normaliser divides `|Zmap|` back out. Its gain ceiling of 0.008 (klecks) is a property of the analysis. So #4's "significant" gain cells are real but negligible in size, and #4 rests on **replication and the diagonal geometry**, not on the gain magnitude.

---

## Why PGD is positive but every de-rotation test is negative

This is the central puzzle of the whole folder, and it has a clean answer: **the three tests ask three different questions, and only PGD's is weak enough to pass.**

```
   #2 PGD        at ONE frequency, do the phase-gradient arrows point
                 the same way across the array?                          -> YES

   #4 cortical   does the phase ramp GROW with frequency, as k = 2πf/v
      de-rotation with ONE v for all f?                                  -> NO

   #5 stimulus   does the preferred phase ramp with STIMULUS distance?   -> NO
      de-rotation
```

### PGD vs cortical de-rotation: direction is not propagation

PGD is `|mean(grad φ)| / mean(|grad φ|)`. Two properties make it permissive:

- **It is scale-free.** The normaliser divides the gradient magnitude straight out, so PGD sees only the *direction* of the arrows, never how big the phase steps are.
- **It is evaluated one frequency at a time.** It never compares frequencies, so it cannot notice that the ramp fails to scale.

A **fixed, frequency-independent spatial phase offset** therefore scores a perfect PGD at every frequency — the arrows are identical and perfectly aligned — while being the opposite of a traveling wave.

The discriminator is what stays constant across frequency:

```
   constant TIME delay τ      ->  Δφ = 2πf·τ  ->  k ∝ f     ->  v fixed
                                                              HORIZONTAL band
                                                              = TRAVELING WAVE

   constant PHASE offset Δφ   ->  k independent of f  ->  v = 2πf/k ∝ f
                                                              DIAGONAL ridge
                                                              = STATIONARY pattern
```

Your data is unambiguously the second. Three independent readings agree:

| evidence | reading |
|---|---|
| PGD's own per-frequency speeds: theta ≈ 5, beta ≈ 28–51 cm/s | `v` rising with `f` |
| §4 klecks: `f` ×4.6 → `v` ×5.9; hermes: `f` ×3.1 → `v` ×2.7 | `v ∝ f`, so `k` const |
| §4/§5 grids: diagonal ridge, never a horizontal band | `k` const |

PGD was never wrong. It correctly detected an **aligned spatial phase gradient**. De-rotation adds the requirement that the gradient *scale with frequency*, and that is what fails.

**This also rules out a whole class of explanations.** Anything based on signals propagating — axonal conduction, attention sweeping across cortex, any finite-speed mechanism — produces a constant *time* delay and therefore `k ∝ f`. The observation of constant `k` is incompatible with all of them. What remains are stationary explanations: a fixed anatomical/laminar phase offset across the array, reference or volume-conduction structure, or a genuine standing (non-propagating) phase pattern.

### PGD vs the stimulus test: different axes entirely

These two are not even measured along the same dimension.

```
   PGD / cortical de-rotation :  distance = mm across the ARRAY
   stimulus de-rotation       :  distance = degrees in the VISUAL FIELD
```

A phase gradient can exist across cortex without being organised by where the stimulus is. The scanning hypothesis needs the *stimulus* axis, and three separate tests say it is not there:

1. **Test #1** — preferred phase barely depends on stimulus position: significant at 2/35 (hermes) and 3/35 (klecks) frequencies.
2. **Test #3** — the wave focus is not at the RF-driven patch (p ≈ 0.36–0.78) and does not follow it across positions (p ≈ 0.42–0.90).
3. **Test #5** — `R0` = 0.80/0.93 under `phase`: the per-location preferred phases are already 80–93 % aligned *before* any de-rotation, so there is almost no position-dependent phase for any wave model to organise. Under `coherence`, where that ceiling is gone, `max R ≤ max R0` in every mode and both animals.

Test #1 and test #5's `R0` are the same fact stated twice — one as a correlation, one as a baseline.

### Putting it together

```
   there IS a spatial phase structure across V4        (PGD, replicated)
        it is PLANAR, not radial or spiral             (test #3 classification)
        it is STATIONARY, not propagating              (constant k, both scripts)
        it is NOT organised by stimulus position       (tests #1, #3-origin, #5)
```

So the honest summary is not "we found nothing". It is: **there is a robust, replicated, planar spatial phase gradient across V4 — and it is a fixed pattern, not a wave, and not tied to the scanned stimulus.**

## The two estimators — `phase` (PHASE ALIGNMENT) vs `coherence` (PHASE COHERENCE)

Both de-rotation scripts set `ESTIMATORS = {'phase','coherence'}` and compute **both in a single run**. The three modes say **what** is de-rotated; the estimator says **what is being vector-summed**, and therefore what the number `R` actually means.

> Naming note: the `coherence` estimator was called `trial` in earlier versions. The SLURM worker and its on-disk cache keep the name `trial_position_sums` — that is literally what they hold (per-**trial** sums grouped by position). Only the estimator was renamed.

### In one line

**The estimator swaps what sits in each grid cell: an alignment index, or an actual phase coherence. Nothing else changes.**

```
   cell (f, v) of the grid contains:

   'phase'      vector-sums the nPos PER-LOCATION PREFERRED-PHASE vectors
                R = resultant length / sum of weights
                -> "how concentrated are these nPos angles"
                -> NOT a phase coherence

   'coherence'  vector-sums ALL TRIALS directly
                R = |Σ_t y_t·e^(iφ_t)| / Σ_t|y_t|
                -> a genuine phase coherence (the phase_coherence/ measure)
```

**Identical under both** — this is worth being blunt about, because it is easy to assume the gain is what distinguishes them:

| shared | |
|---|---|
| the de-rotation | `k = 2πf/v`, spin by `k·d` |
| the three modes | `visual` / `visual_coherent` / `visual_arrival` |
| the baseline | `R0` = the same statistic at `k = 0` |
| **the gain** | **`dR = R − R0` — computed identically. `R − R0` is NOT the difference between the estimators.** |
| the null | synchronised label shuffle |
| thresholds, replication, pooling | |

### Two things that are easy to get backwards

- **`R − R0` is a subtraction of two magnitudes** (resultant lengths / coherences), not of phases. The phase subtraction `φ_p − k·d_p` happens *inside*, before the vector sum, in both estimators.
- **`phase` is not a coherence.** Its `R` is a weighted resultant length over locations — a circular-concentration measure. The trial-level coherence enters only as the **weight** `coh_mag`; the trials were averaged away upstream. This is why its `R0` sits at 0.80–0.99: a handful of angles are easily concentrated.

### What each run produces

Three figures per RF setting, both estimators in one run:

```
   stimulus_loc_grids[_validRF].pdf
       THE DE-ROTATION R. 6 rows = 3 modes x 2 estimators, estimator named in
       every panel title. 4 cols = animals + mean/replication + pooled z.
       Colour limits are set PER ESTIMATOR BLOCK — a fixed [0,1] scale renders
       the phase-coherence block (R ~ 0.08) as flat blue.

   stimulus_loc_gain_phase[_validRF].pdf        gain dR, PHASE ALIGNMENT
   stimulus_loc_gain_coherence[_validRF].pdf    gain dR, PHASE COHERENCE
       3 rows x 3 cols: animals + mean/replication. No pooled column (below).

   stimulus_loc_wave[_validRF].mat
       BOTH estimators: results.G.(estimator).(mode)(animal), results.C.(estimator)
```

The cortical script mirrors this: `derotation_R_grids.pdf` (2 rows = the two estimators), plus `derotation_gain_grids_<est>.pdf` and `derotation_direction_speed_<est>.pdf` per estimator.

**Why the gain figures carry no pooled column.** Because `gain = R − R0(f)` and `R0` depends on neither speed nor the shuffle:

```
   muG = muR − R0 ,   sdG = sdR
   zG  = (R − R0 − muG)/sdG = (R − muR)/sdR = zR
```

The pooled gain panel would be a pixel-for-pixel copy of the pooled R panel (verified: `max|zR − zG| = 1.75e-13`, identical thresholds), so it is drawn once, on the grids figure. Both scripts now assert this at runtime and warn if it ever fails. The **per-animal** R and gain panels are genuinely different and both are kept.

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

### The `coherence` estimator collapses the two levels into one

Every trial is rotated by the model's prediction for the location *that trial* had, and one coherence is taken over the whole trial set:

```
   c(f,v) = ( 1/W ) · Σ over ALL TRIALS  y_t · e^( i ( phi_t − k·d_p(t) ) ) ,
   W      = Σ_t |y_t| ,        R = |c|
```

At `k = 0` this is **exactly the phase coherence of the `phase_coherence/` pipeline** with all locations pooled.

**Order of operations.** De-rotation happens on the *trials*, before any coherence is computed. There is exactly **one** coherence per `(f, v)` cell — not one per position.

```
   trial   location   d_p    y_t    phi_t      phi_t − k·d_p   <- de-rotated
   -----   --------   ---    ---    ------     -------------
     1        p1       0     2.1     10°            10°
     2        p3       4     1.4    170°            10°        (k·d = 160°)
     3        p2       2     0.8     90°            10°        (k·d =  80°)
     4        p4       6     2.2    250°            10°        (k·d = 240°)

   then ONE coherence over all of them:
       R(f,v)  = | Σ_t y_t · e^(i(phi_t − k·d_p(t))) | / W      phases shifted
       R0(f)   = | Σ_t y_t · e^(i·phi_t)             | / W      phases untouched
       gain    = R − R0                                          de-rotated MINUS actual
```

The raw phase column (10°, 170°, 90°, 250°) is scattered → low coherence. The de-rotated column is all 10° → high coherence. That difference *is* the gain.

**`S(p,f)` is a shortcut, not a comparison.** Every trial at the same location gets the same rotation `e^(−i·k·d_p)`, so the trials at a location can be pre-added:

```
   Σ_t  y_t·e^(i(phi_t − k·d_p(t)))  =  Σ_p  e^(−i·k·d_p) · [ Σ_{t at p} y_t·e^(i·phi_t) ]
                                                             \_______________________/
                                                                       S(p,f)
```

`S(p,f)` is cached once, then every `(f,v)` cell is a cheap re-weighting of ~16 numbers. It is never divided by anything or read on its own — per-location coherences are never computed, and never compared with each other.

**Three ways to misread this:**

- *"De-rotation adds coherence to a position."* No — `|S_p·e^(−i·k·d_p)| = |S_p|`. Rotating never changes a vector's length. Only the **sum across locations** changes, because rotation changes directions and therefore how much the terms cancel.
- *"The wave predicts position 3 is more coherent than position 1."* No — it predicts a **phase** relation, `phi_3 − phi_1 = k·(d_3 − d_1)`. Same coherence, later phase. (Only distance *differences* matter, which is why `d_vis = ecc − min(ecc)` is safe: a common offset rotates every term equally and leaves `|Σ|` unchanged.)
- *"gain = actual − de-rotated."* Backwards. `gain = R − R0`, de-rotated minus actual, so a real wave gives a **positive** gain.

### What R, R0 and the gain mean under each

| | `phase` (PHASE ALIGNMENT) | `coherence` (PHASE COHERENCE) |
|---|---|---|
| what R measures | how similar the per-location **preferred phases** are after de-rotation | the **phase coherence** over all trials after de-rotation |
| what R0 is (`k = 0`) | how similar they already were | the pooled phase coherence, no rotation |
| typical R0 (measured, `visual`, validRF) | median **0.80** hermes / **0.93** klecks (cortical script: 0.947 / 0.992) | median **0.063** hermes / **0.078** klecks |
| gain ceiling `1 − R0` | tiny — 0.07 for klecks, 0.008 in the cortical script | large in principle; bounded in practice by how much position-dependence exists |
| gain threshold actually obtained | 0.031 hermes / 0.008 klecks | 0.0006 hermes / 0.0007 klecks |
| units comparable across estimators? | **no** — never compare an R from one with an R from the other | |

`gain = R − R0` is computed identically for both. Nothing about the gain readout was removed or changed.

A **lower `R0` under `coherence` is good news**, not a weaker result — it means headroom for a gain to appear. Compare each estimator against its own threshold, never against the other's.

### Why `coherence` is the better-behaved estimator

**1. No binning loss.** `phase` estimates 16 separate coherences from ~1/16 of the trials each, then fits to those noisy intermediates. `coherence` fits the unbinned data directly.

**2. It fixes a weighting flaw.** `phase` weights location *p* by `|c_p|` — but under noise `|c_p| ~ 1/√n_p`, so a location with **fewer** trials gets a **larger** weight. Backwards. The two estimators differ by exactly this:

```
   phase     :  uses  S(p,f) / n_p   (the per-location MEAN, then weighted by |c_p|)
   coherence :  uses  S(p,f)         (the raw SUM — every trial counts once)
```

**3. An interpretable baseline.** `R0` becomes a coherence directly comparable with the `phase_coherence/` results, instead of a phase-similarity index pinned near 0.9 by construction.

**It will not manufacture an effect.** If preferred phase barely depends on location, both estimators say so. `coherence` has more power to see a small one and a cleaner baseline — that is all.

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

Both estimators run together. The R grids figure holds both; the gain (and, in the cortical script, direction/speed) figures carry `_phase` / `_coherence`, and `[_validRF]` marks the RF filter, so all combinations coexist rather than overwriting.

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

The R-grid PDF is laid out **rows = modes x estimators, columns = hermes / klecks / mean+replication / pooled z**; the gain PDFs are **rows = modes, columns = hermes / klecks / mean+replication** (no pooled column — see the estimator section):

```
   Plots/scanning/stimulus_loc_wave/cp10_till_100/lfp/
       stimulus_loc_grids_validRF.pdf            absolute R, BOTH estimators
       stimulus_loc_gain_phase_validRF.pdf       gain dR, PHASE ALIGNMENT
       stimulus_loc_gain_coherence_validRF.pdf   gain dR, PHASE COHERENCE

   results_combined/scanning/stimulus_loc_wave/cp10_till_100/lfp/
       stimulus_loc_wave_validRF.mat             both estimators in one file
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
