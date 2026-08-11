# `scanning/` — traveling-wave analyses of V4 LFP

**Hypothesis under test.** Attention scans the stimulus locations, so the optimal (preferred) phase should shift systematically across cortical space — a traveling wave whose origin moves with the stimulated retinotopic patch.

Two animals (`hermes`, `klecks`), 64-channel 8×8 V4 arrays, `cp10_till_100`, `dv = 'lfp'` unless stated. Animals are **never pooled at the channel level**; they are combined by replication and by a pooled standardised test.

> **The two animals do not have the same number of stimulus locations** — hermes has **16**, klecks **9** (`positions` in `phase_progression.mat`). Wherever this file says "16 locations" it is quoting hermes; substitute 9 for klecks. Nothing in the code hard-codes 16, but every per-location noise floor is `~1/√nPos`, so klecks' floors are correspondingly higher.

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
        ├──────────────► functions/trial_position_sums_chan.m ──► trial_position_sums/
        │                 (SLURM: 1 job/channel)                  S(ch,pos,freq), W(ch)
        │                                                         │
        │                           used by the 'coherence' ESTIMATOR in the two
        │                           de-rotation scripts (cache shared by both)
        │
        └──────────────► erp_latency_wave.m   (goes back to clean_lfp.mat —
                          the only script that uses TIME-domain LFP)
```

`phase_progression.m` is the only script that touches raw trials in the main path. Everything else consumes `phase_progression.mat` — except the `coherence` estimator, which goes back to the trials via the cached per-location sums, and `erp_latency_wave.m`, which needs the time-domain signal.

Also read (not produced here): `phase_coherence/complex/cp10_till_100/<dv>/all_loc_difflev/<ch>/coherence.mat` + `coh_perm_complex.mat`, used to build the per-channel coherence-significance mask (`load_coh_sig_mask`).

That mask gates **§2, §4 and §5 only** (`CH_FILTER = 'significant'` in each). §3 uses its own reliability quantile instead (`RELIABLE_Q = 0.5`, drop the bottom half of `coh_mag` per frequency), and §1 and §6 apply no such gate at all — so the channel sets are not identical across the six scripts.

---

## The quantities, in one place

Every script below is built from the same handful of numbers. This section defines them once.

### The phase map

```
   c(channel, freq, position) = mean over TRIALS of  y_t · e^(i·phi_t)

       |c|      = coherence magnitude   -> saved as coh_mag
       angle(c) = preferred phase       -> saved as pref_phase
```

One complex number per (channel, frequency, position): **how reliably** that channel locks to that frequency, and **at what phase**. Everything downstream is a statement about `pref_phase` weighted by `coh_mag`.

### Wave quantities

| symbol | meaning | units |
|---|---|---|
| `f` | frequency | Hz |
| `v` | propagation speed | cm/s (cortex) or deg/s (visual field) |
| `k` | wavenumber = phase change per unit distance, `k = 2πf/v` | rad/mm or rad/deg |
| `d` | distance along the propagation axis | mm or deg |
| `τ` | time delay between two sites | s |
| `Δφ` | phase difference between two sites, `Δφ = k·d = 2πf·τ` | rad |

**The one relation the whole folder turns on:**

```
   a constant TIME delay τ   ->  Δφ = 2πf·τ  grows with frequency  ->  k ∝ f
   a constant PHASE offset   ->  Δφ fixed                          ->  k constant
```

Propagation of any kind — axonal conduction, attention sweeping across cortex, any finite-speed mechanism — takes *time*, so it must give `k ∝ f`. A pattern frozen in space gives constant `k`. **This is the discriminator every test below is trying to reach.**

### PGD

```
   PGD = |mean(grad phi)| / mean(|grad phi|)          0 <= PGD <= 1
```

The phase gradient is computed at every electrode; PGD asks whether those arrows **point the same way**. 1 = perfectly aligned (a plane), 0 = random.

Two properties decide how it must be read:

- **Scale-free** — the denominator divides the gradient magnitude straight out, so PGD sees only arrow *direction*, never how big the phase steps are.
- **Single-frequency** — it never compares frequencies, so it cannot see whether the ramp scales with `f`.

**A significant PGD therefore says the phase map is a plane. It does not say the plane is moving.** A frozen phase offset scores a perfect PGD at every frequency.

The direction of the tilt comes from the same complex number: with `g = gx + i·gy` (real = column, imag = row), `PGD = |mean(g)|/mean(|g|)` and `dir = angle(mean(g))` are its **modulus and argument**.

### The bias floor — why PGD's speeds cannot carry a wave claim

`GMAG = |mean(grad φ)|` is a vector mean over only ~56 finite-difference sites, so random scatter does **not** cancel completely. `GMAG_null` — the same statistic on a *shuffled* map of the same phases — measures what survives from nothing, and it is not small:

```
   GMAG / GMAG_null   ~1.0  ->  the measured tilt IS the finite-array artifact
                      ~1.2  ->  five-sixths of it could be artifact
```

Over the significant bands the observed ratios run **0.89–1.67 (hermes)** and **1.53–2.91 (klecks)**, so `k = GMAG/spacing` is inflated, and since `v = 2πf/k`, every speed from this script is correspondingly deflated. The floor also varies with frequency, so the bias is not even a constant offset — it changes across the spectrum.

The band report therefore prints `GMAG/GMAG_null` next to every band, and `gradient_magnitude_per_animal.pdf` draws the floor as a dashed line under the measured tilt. **Where the curves meet, the tilt and the speed are not interpretable.**

Treat every speed from this script as descriptive. Whether the ramp *scales with frequency* — the question that separates a wave from a frozen offset — is answered by de-rotation (§4), which fits the wave model to all electrodes at once and so never forms this local derivative.

### `R`, `R0` and the gain (de-rotation)

```
   R  =  | SUM  w · e^( i ( phi  -  k * distance ) ) |  /  SUM w        0 <= R <= 1

   wrong v                        right (f, v)
   ↖  ↑  ↗  →                     ↗  ↗  ↗  ↗
   fanned out,  R small           aligned,  R -> 1
```

| quantity | definition | meaning |
|---|---|---|
| `R(f,v)` | resultant after de-rotating by `k·distance` | how well the phasors line up **under this wave model** |
| `R0(f)` | the same statistic at `k = 0` | how well they lined up **already**, before any model |
| gain `dR` | `R − R0` | what the wave model **bought**. De-rotated minus actual, so a real wave gives a **positive** gain |
| ceiling | `1 − R0` | the most any model could possibly gain |

`R0` is the baseline that makes a gain interpretable: `R = 0.99` is impressive only if `R0` was not already 0.98. **The ceiling is the first thing to check before reading any gain.**

### How to read a frequency × speed grid

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

This is **the** discriminator between propagation and a frozen offset, and this script is where it is decided. Support for a wave = horizontal band, at a plausible speed, **replicated in both animals**.

### Statistics vocabulary

- **Synchronised permutations.** Permutation *k* uses `rng(perm_seed_base + k)` — the *same* trial→label relabelling in every channel. Array LFP is spatially correlated; independent per-channel shuffles make the channel-average null far too tight and the channel-average significance anti-conservative.
- **Max-statistic correction.** The threshold is the 95th percentile of the per-permutation *maximum over the whole grid*, so it controls the family-wise error across all cells at once (e.g. 35 × 30 = 1050 cells, or 35 × 30 × 24 = 25 200 with direction swept). One threshold for the entire grid — which is strict, and deliberately so.
- **Cluster correction across frequency.** Used where the axis is frequency alone: cluster mass = `Σ(stat − thr)` over each supra-threshold run, compared with the null's largest cluster mass. Reports a *band*, not isolated bins.

---

## The scripts, and what each one asks

### 1. `phase_progression.m` — the producer, and test #1

**Question.** For each channel, does the preferred phase depend *systematically* on stimulus position?

**Method.** Per (channel × frequency): complex coherence per position, `c_p = mean over trials of y·exp(i·phase)`; `angle(c_p)` = preferred phase, `|c_p|` = coherence magnitude. Statistic = circular–linear correlation `r_cl` between position index (peripheral→foveal) and the preferred-phase vector. Null = shuffle position labels across trials, **synchronised** across channels (`functions/phase_progression_chan.m`).

Also produces `pref_phase` and `coh_mag`, which every other script runs on.

### 2. `cortical_planar_wave_PGD.m` — test #2, existence of a planar wave

**Question.** Is there a planar traveling wave across the 8×8 array, and at which frequencies?

**Method.** Build the preferred-phase map on the array using only coherence-significant channels. Compute PGD; null = shuffle phases across electrode locations; **cluster-corrected across frequency**. Three views:

| view | what it does |
|---|---|
| `COLLAPSED` | positions merged (reliability-weighted circular mean) — a wave survives only if consistent across positions |
| `PER-POSITION` | wave fitted separately per stimulus position, plus the wave-**origin** check against RF-driven channels (test #3) and cross-position direction agreement |
| `CONSENSUS` | collapsed-sig **AND** cross-position direction agreement (Rayleigh) **AND** sig in ≥ N positions |

Each significant band prints its `GMAG/GMAG_null` margin. **This script does not test propagation** — it says the phase map is a plane, not that the plane moves. For that, see §4.

> **The origin check has no null.** `origin_align(f,p)` is the angle between the fitted propagation axis and the driven-patch→array-centre axis, folded to 0–90°, and the script prints its mean — there is no permutation distribution and no p-value, so it cannot currently return a verdict either way. 45° is what a random axis gives; treat departures from it as descriptive until a null is added.

Outputs: `pgd_existence_per_animal.pdf`, `pgd_existence_combined.pdf`, `gradient_magnitude_per_animal.pdf`, `per_position_existence.pdf`, `wave_origin_<animal>.pdf`, `coverage_per_position.pdf`, `consensus_wave.pdf`.

### 3. `cortical_wave_type_classification.m` — what kind of pattern is it?

**Question.** Phase can be spatially organised in several ways. Which one is this?

**Method.** On the collapsed map, test four patterns against their own nulls:

| pattern | statistic |
|---|---|
| planar | PGD |
| radial (source/sink) | net divergence = mean Laplacian |
| rotating/spiral | net curl + count of phase singularities (2×2 loops whose wrapped phase sums to ±2π) |
| synchronous/none | mean \|grad\| not above null |

Works on **wrapped** phase with local neighbour differences and deliberately does *not* 2-D unwrap — unwrapping would destroy the very singularities a spiral is defined by.

> **Its channel filter is not §2's.** This script has no `CH_FILTER`; it drops the bottom `RELIABLE_Q = 0.5` quantile of `coh_mag` per frequency. So its map is built from a different channel set than the PGD map it is compared against — close, but not the same electrodes.

Outputs: `wave_type_vs_freq.pdf`, `wave_type_strip.pdf`, `phase_maps_focus.pdf`. **No `.mat` is saved** — this is the one script whose results exist only as figures.

### 4. `cortical_planar_wave_derotation.m` — is the cortical wave dispersion-free?

**Question.** PGD says the gradients are aligned. Does the phase ramp also *scale with frequency*, as a real wave must?

**Method.** Rotate-to-overlap in cortical space. Each electrode sits at `r_c` mm on the array; a plane wave of speed `v`, direction `θ` predicts a lag `k·d_c(θ)` with `d_c(θ) = x_c cosθ + y_c sinθ`. De-rotate, measure the resultant, sweep **frequency × speed × direction** (30 speeds 1–300 cm/s × 24 directions). Summed coherently across electrodes. Null = shuffle electrode ↔ array coordinate.

> **Structural note on the ceiling here.** This script's statistic is always a resultant *across electrodes*, `R = |Σ_c w_c e^(i(φ_c − k·d_c))| / Σ_c w_c`, and that normaliser divides `|Zmap|` back out. Switching estimator changes the *weights*, not the scale of `R` — so unlike §5, the `coherence` estimator does **not** lower `R0` here. The ceiling is a property of the analysis, not of the estimator.

Outputs: `derotation_R_grids.pdf` (2 rows = the two estimators), plus `derotation_gain_grids_<est>.pdf` and `derotation_direction_speed_<est>.pdf` per estimator.

### 5. `stimulus_loc_traveling_wave.m` — does the wave track the scanned stimulus?

**Question.** The actual scanning hypothesis: as the stimulus moves, does the preferred phase ramp with distance?

**Method.** Same rotate-to-overlap, but the distance axis is **stimulus eccentricity in degrees**, so speed is in deg/s (30 speeds, 1–200 deg/s; no direction sweep — the axis is already 1-D). Three modes (full detail in [the appendix](#appendix--the-three-modes-of-stimulus_loc_traveling_wavem)):

| mode | de-rotates by | says |
|---|---|---|
| `visual` | `k·d_p` | phase depends on where the stimulus is |
| `visual_coherent` | `k·(d_p + D_c)` | + a per-electrode delay fixed across positions (wave sweeping out from the fovea) |
| `visual_arrival` | `k·(d_p + a(c,p))` | the wave starts *where the stimulus lands* |

Both estimators run in a single execution. Outputs: `stimulus_loc_grids[_validRF].pdf`, `stimulus_loc_gain_phase[_validRF].pdf`, `stimulus_loc_gain_coherence[_validRF].pdf`, `stimulus_loc_wave[_validRF].mat`.

### 6. `erp_latency_wave.m` — the time-domain cross-check

**Question.** Everything above measures phase. This measures **time** directly: does the LFP response arrive later at some electrodes than others, and by how much per millimetre?

**Method.** Pool trials across sessions → stimulus-locked ERP per channel → zero-phase low-pass (`filtfilt`; a *causal* filter would inject its own group delay, i.e. fabricate the quantity being measured) → per-channel latency by cross-correlation against the array-mean ERP → regress latency on distance. Slope⁻¹ = speed. Null = channel shuffle.

Two axes, answering different questions:

```
   'cortical'      latency vs mm along the PGD-fitted propagation direction
                   -> tests cortical propagation (§2/§4)

   'retinotopic'   latency vs RF eccentricity, foveal -> peripheral
                   -> tests the scanning hypothesis (§5)
```

Then the comparison that decides it, against PGD's saved wavenumbers:

```
   k_pred(f) = 2*pi*f * (dtau/dd)      from the ERP latency slope
   k_obs(f)  = k_corr                  from planar_wave_existence.mat

   k_obs scales with f AND matches |k_pred|  ->  the tilt IS a time delay
   k_obs flat while |k_pred| rises           ->  fixed phase offset
```

> ⚠️ **Compare MAGNITUDES, not signed values.** `k_obs = |mean(grad φ)|/spacing` is a magnitude and can never be negative, whereas `k_pred` carries the sign of the latency slope. A negative `k_pred` therefore says nothing on its own — the meaningful quantities are `|k_pred|` versus `k_obs`, and separately whether latency rises or falls along the PGD direction.

**Two independent things the comparison tests, and only the first is decisive:**

```
   1. SCALING   does k_obs grow in proportion to f, as a fixed time delay must?
                a flat k_obs while |k_pred| rises = fixed offset, whatever the speeds
   2. MAGNITUDE do the two speeds agree numerically?
                they can scale alike and still disagree by a constant factor
```

Scaling is the strong test because it is calibration-free. Magnitude agreement is worth reporting but a mismatch has mundane explanations — the evoked broadband transient and the ongoing narrowband phase gradient need not be the same phenomenon.

**Channel mapping.** V4 labels are mapped onto canonical 1..64 slots **by label number** (`V4-n` → slot `n` for hermes, `n − 64` for klecks), never by sorted position: no session contains all 64 channels and *which* are missing varies per session, so positional selection would put a different electrode in row *i* in different sessions. Absent channels stay NaN and contribute no trials. Verified against `phase_progression.mat`'s zero-trial channels — see `v4_channel_slots`.

**Assumptions, status after the first run:**

| # | assumption | status |
|---|---|---|
| 1 | `clean_data.time` zeroed on the stimulus event | **discharged** — ERPs are flat pre-0, response onset ~60–70 ms, peak ~100–150 ms |
| 2 | channel slot matches `phase_progression.mat` | **discharged** — absent-slot sets reproduce its zero-trial channels exactly in both animals |
| 3 | pooling sessions across electrode drift | open — adds noise, no bias in the slope |

Outputs: `erp_latency.pdf` (per animal: all channel ERPs, then latency-vs-distance on each axis) and `erp_latency.mat` (`results.L(ia).cortical` / `.retinotopic` with `slope`, `speed`, `r`, `p`, `nullslope`; plus `k_obs`, `k_pred`, `fHz`, `pgd_sig`).

### 7. `traveling_wave_H2_H1.m` — deprecated

Descriptive 9-page PDF (phase heatmaps, phase-change maps, circular-variance maps, per-channel phase trajectories, Moran's I spatial autocorrelation, 2-D plane fit). Uses H2-H1 significance thresholds instead of coherence significance. Untracked in git. Still useful for visualisation.

---

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

**Identical under both** — worth being blunt about, because it is easy to assume the gain is what distinguishes them:

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
- **`phase` is not a coherence.** Its `R` is a weighted resultant length over locations — a circular-concentration measure. The trial-level coherence enters only as the **weight** `coh_mag`; the trials were averaged away upstream. This is why its `R0` sits high: a handful of angles are easily concentrated.

### Level 1 — coherence across TRIALS (the `phase_coherence/` sense)

Computed once, upstream, in `phase_progression.m` — the `c(channel, freq, position)` defined at the top of this file. One complex number per (channel, frequency, position), the same construction as the `phase_coherence/` folder.

### Level 2 — alignment of those preferred phases ACROSS LOCATIONS

What the `phase` estimator does with them:

```
                | Σ_p  coh_mag · e^( i ( pref_phase − k·d_p ) ) |
   R_c(f,v)  =  ──────────────────────────────────────────────────
                            Σ_p  coh_mag
```

The level-1 coherence enters only as a **weight**; what is summed is the *preferred phases*. So R asks "do the `nPos` per-location preferred phases (16 hermes, 9 klecks) line up after de-rotation?"

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

`S(p,f)` is cached once, then every `(f,v)` cell is a cheap re-weighting of `nPos` numbers (16 hermes, 9 klecks). It is never divided by anything or read on its own — per-location coherences are never computed, and never compared with each other.

**Three ways to misread this:**

- *"De-rotation adds coherence to a position."* No — `|S_p·e^(−i·k·d_p)| = |S_p|`. Rotating never changes a vector's length. Only the **sum across locations** changes, because rotation changes directions and therefore how much the terms cancel.
- *"The wave predicts position 3 is more coherent than position 1."* No — it predicts a **phase** relation, `phi_3 − phi_1 = k·(d_3 − d_1)`. Same coherence, later phase. (Only distance *differences* matter, which is why `d_vis = ecc − min(ecc)` is safe: a common offset rotates every term equally and leaves `|Σ|` unchanged.)
- *"gain = actual − de-rotated."* Backwards. `gain = R − R0`, de-rotated minus actual, so a real wave gives a **positive** gain.

### What R, R0 and the gain mean under each

| | `phase` (PHASE ALIGNMENT) | `coherence` (PHASE COHERENCE) |
|---|---|---|
| what R measures | how similar the per-location **preferred phases** are after de-rotation | the **phase coherence** over all trials after de-rotation |
| what R0 is (`k = 0`) | how similar they already were | the pooled phase coherence, no rotation |
| gain ceiling `1 − R0` | small — the baseline is high by construction | large in principle; bounded in practice by how much position-dependence exists |
| units comparable across estimators? | **no** — never compare an R from one with an R from the other | |

`gain = R − R0` is computed identically for both.

A **lower `R0` under `coherence` is good news**, not a weaker result — it means headroom for a gain to appear. Compare each estimator against its own threshold, never against the other's.

### Why `coherence` is the better-behaved estimator

**1. No binning loss.** `phase` estimates `nPos` separate coherences from ~1/`nPos` of the trials each (1/16 hermes, 1/9 klecks), then fits to those noisy intermediates. `coherence` fits the unbinned data directly.

**2. It fixes a weighting flaw.** `phase` weights location *p* by `|c_p|` — but under noise `|c_p| ~ 1/√n_p`, so a location with **fewer** trials gets a **larger** weight. Backwards. The two estimators differ by exactly this:

```
   phase     :  uses  S(p,f) / n_p   (the per-location MEAN, then weighted by |c_p|)
   coherence :  uses  S(p,f)         (the raw SUM — every trial counts once)
```

**3. An interpretable baseline.** `R0` becomes a coherence directly comparable with the `phase_coherence/` results, instead of a phase-similarity index pinned near 0.9 by construction.

**It will not manufacture an effect.** If preferred phase barely depends on location, both estimators say so. `coherence` has more power to see a small one and a cleaner baseline — that is all.

### Why the gain figures carry no pooled column

Because `gain = R − R0(f)` and `R0` depends on neither speed nor the shuffle:

```
   muG = muR − R0 ,   sdG = sdR
   zG  = (R − R0 − muG)/sdG = (R − muR)/sdR = zR
```

The pooled gain panel would be a pixel-for-pixel copy of the pooled R panel, so it is drawn once, on the grids figure. Both scripts **assert this at runtime** and warn if it ever fails. The **per-animal** R and gain panels are genuinely different and both are kept.

### Implementation

The de-rotation depends on a trial only through its location, so the estimator collapses onto per-location complex **sums**:

```
   S(c,p,f) = Σ over trials at location p of  y_t · e^(i·phi_tf)
   c(f,v)   = ( 1/W ) · Σ_p  S(c,p,f) · e^(−i·k·d_p)
```

`functions/trial_position_sums_chan.m` computes `S` for the observed labels and for all `nPerm` shuffles — **one SLURM job per channel**, same `slurmfun` pattern and the same synchronised `rng(perm_seed_base + k)` seeding as `phase_progression_chan.m`. All three modes then reduce to re-weightings of `S`, done in-script.

```
   cache: results_<animal>/scanning/trial_position_sums/cp10_till_100/<dv>/<ch>/

   stimulus_loc_traveling_wave.m   needs S_obs AND S_perm
                                   — its null shuffles trial->location labels
                                   on disk: 248 MB hermes, 142 MB klecks
                                   (the size scales with nPos: 16 vs 9)
   cortical_planar_wave_derotation.m  needs S_obs only      (~0.3 MB/animal)
                                   — its null shuffles electrode->array coordinate
```

The cache is shared: whichever script runs first pays for the jobs. `RECOMPUTE_TRIAL_SUMS = true` forces resubmission (do this after changing `dv`, `nPerm`, or the worker).

---

## Combining animals: two criteria, side by side

Channels are never pooled across animals. Instead the last two columns of each figure report two different combination rules:

| | **replication** (col 3) | **pooled standardised z** (col 4) |
|---|---|---|
| What is shown | mean of the two animals' grids | mean of the two animals' *z-scored* grids |
| What the contour means | `sig(hermes) AND sig(klecks)` | pooled z above its own max-stat threshold |
| Is there a test on it? | no — the AND of two per-animal tests | yes — paired-permutation null on the pooled grid |
| Can one animal drive it? | **no** | **yes** |
| Power | low | higher (aggregates weak evidence) |
| Role | **primary** | secondary, always report which animal drives it |

Why standardise before pooling: the animals sit on very different scales (their `R0` and thresholds differ several-fold), so a raw average is dominated by the animal with the bigger numbers. Each cell is converted to `z = (obs − null mean) / null sd` against that animal's own permutation distribution first.

The pooled null pairs permutation *b* of one animal with permutation *b* of the other — legitimate because the animals are independent, so the pairing is arbitrary and samples the product null.

```
   Zobs(f,v)    = mean over animals of z_a(f,v)
   Znull(f,v,b) = mean over animals of z_a_null(f,v,b)
   thr_pool     = 95th pctile of  max over grid of Znull
```

**With n = 2 neither rule supports a population-level inference** — between-animal variance isn't estimable. Replication is simply the more conservative descriptive rule. Note also that averaging roughly halves a one-animal effect, so a strong single-animal result can easily fail the pooled test too.

---

## Gotchas worth knowing

- **Coordinate frames.** `trialinfo` col 16/17 are **fixation-centred** pixels; `RF_Center_X/Y` in the RF summary tables are **screen** pixels. Differencing them without subtracting (840, 525) is off by the screen centre. `functions/elec_rf_deg.m` handles this, and exists in exactly one place for that reason.
- **`ppd` is not in the repo.** `PIX_PER_DEG = struct('hermes',53.24,'klecks',50.56)`, from rig calibration. `sessInfo.ppd` does not exist in the saved `.RF` files despite being referenced elsewhere.
- **Extrapolated RF centres.** 18/64 (hermes) and 21/64 (klecks) rows are `Status = Extrapolated` — the Gaussian fit failed and the centre was filled from the median of 8-connected array neighbours, or a plane fit. `RF_VALID_ONLY` excludes them. For `visual_coherent` this is a *bias* concern (a smoothed `D_c` becomes a proxy for array position, which already carries the intrinsic wave); for `visual_arrival` it is mostly a *noise* concern.
- **Channel indexing is 1:1.** Both animals are 64-channel 8×8. Klecks' V4 is labelled V4-64…V4-127 in the raw RF data but is still 64 channels. No remap is needed between phase-progression index, RF-fit index, and array layout: `ch_col = ceil((1:64)'/8); ch_row = 8 - mod((1:64)'-1,8);`
- **NaN in permutation nulls.** `mean`/`std` do **not** omit NaN by default; one NaN permutation poisons that cell for all permutations, and NaN then spreads through the cross-animal sum. Both de-rotation scripts use `'omitnan'` and call `null_guard`, warning rather than silently producing `thr = NaN` or `-inf` — which would mark **every** cell significant. Permutations must also shuffle only among channels with finite geometry; shuffling a NaN-geometry channel into a used slot is what caused exactly that failure once.
- **The CMF calibration is missing.** Converting visual degrees to cortical mm needs `E2_DEG` and `A_MM` (`M(E) = M0/(E+E2)`, `D(E) = A·ln(1+E/E2)`). These are **not in the repo** and would have to come from the V4 literature.

---

## Where things land

```
  results_<animal>/scanning/phase_progression/cp10_till_100/<dv>/   phase_progression.mat
  results_<animal>/scanning/trial_position_sums/cp10_till_100/<dv>/ S(ch,pos,freq) cache
  results_combined/scanning/<analysis>/cp10_till_100/<dv>/          *.mat
  Plots/scanning/<analysis>/cp10_till_100/<dv>/                     *.pdf
```

`<analysis>` is `phase_progression`, `planar_wave_existence` (§2), `wave_type` (§3), `planar_wave_derotation` (§4), `stimulus_loc_wave` (§5), `erp_latency` (§6). §3 writes PDFs only — there is no `wave_type` folder under `results_combined`.

Both estimators run together. The R grids figure holds both; the gain (and, in the cortical script, direction/speed) figures carry `_phase` / `_coherence`, and `[_validRF]` marks the RF filter, so all combinations coexist rather than overwriting.

The R-grid PDF is laid out **rows = modes × estimators, columns = hermes / klecks / mean+replication / pooled z**; the gain PDFs are **rows = modes, columns = hermes / klecks / mean+replication**:

```
   Plots/scanning/stimulus_loc_wave/cp10_till_100/lfp/
       stimulus_loc_grids_validRF.pdf            absolute R, BOTH estimators
       stimulus_loc_gain_phase_validRF.pdf       gain dR, PHASE ALIGNMENT
       stimulus_loc_gain_coherence_validRF.pdf   gain dR, PHASE COHERENCE

   results_combined/scanning/stimulus_loc_wave/cp10_till_100/lfp/
       stimulus_loc_wave_validRF.mat             both estimators in one file
```

Colour limits on the grids figure are set **per estimator block** — a fixed `[0,1]` scale renders the phase-coherence block (R ~ 0.08) as flat blue.

The `_validRF` tag means RF-based modes used `Valid_Gaussian` channels only (~46/64 hermes, ~43/64 klecks); without the tag, interpolated `Extrapolated` centres were included too.

---

# Appendix — the three modes of `stimulus_loc_traveling_wave.m`

All three ask the same question — *does the preferred phase carry a delay proportional to distance?* — and differ **only in which distance they de-rotate by**. Null, thresholds, figures and animal-combination logic are shared.

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
| Noise floor of R | ~1/√nPos: ~0.25 hermes (16 loc), ~0.33 klecks (9 loc) | ~0.03 (all c×p terms) | same as `visual` |
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
- **Dependency:** the log-polar CMF calibration, which is not in the repo (see Gotchas). A log-polar mapping is preferable to a "nearest driven patch" mapping because it doesn't depend on which channels are driven, so it sidesteps the saturation problem (hermes positions 10–16 drive zero channels).

#### One asymmetry that is easy to miss

The two scripts **collapse the other axis differently**, and it's a deliberate choice, not an oversight:

```
   cortical script:  collapses positions with a COMPLEX SUM      -> keeps phase
   stimulus Mode A:  collapses channels with MAGNITUDE-then-mean -> discards phase
```

Electrodes are assumed to share a clock, so their relative phase is meaningful and worth preserving. Mode A declines to assume that — which buys robustness to per-channel offsets at the cost of being blind to per-electrode delays.

#### Why a different ruler cannot rescue a null

```
    gain ceiling = 1 − R0 ,   and R0 is the k = 0 baseline,
                              which involves NO distance at all
```

`R0` is identical under any distance metric. A different ruler can only *redistribute* where a given phase pattern lands in the (f, v) grid — it can rescue an effect smeared across speeds by cortical magnification, but it **cannot create phase variance the data doesn't contain**.

---

### Implementation pointers

| what | where |
|---|---|
| Mode list | `stimulus_loc_traveling_wave.m` → `metrics = {'visual','visual_coherent','visual_arrival'}` |
| Mode 1 estimator | `align_grid` (phase) / `coh_grid` (coherence) |
| Mode 2 estimator | `align_grid_coherent` / `coh_grid_coherent` |
| Mode 3 estimator | `align_grid_arrival` / `coh_grid_arrival` |
| `a(c,p)` table | main loop: `a_arr = hypot(rf_deg(:,1) - stim_deg(:,1).', rf_deg(:,2) - stim_deg(:,2).')` |
| RF centres → degrees | `functions/elec_rf_deg.m` (subtracts the screen centre 840/525 before `/ppd`) |
| Valid-only toggle | `RF_VALID_ONLY` |
| Gradient bias floor | `cortical_planar_wave_PGD.m` → `GMAG_null` (shuffle null), printed as `GMAG/GMAG_null` per band |
| Frequency-scaling test | `cortical_planar_wave_derotation.m` → `vbest` ridge (§4); **not** in the PGD script |
| NaN-safe permutation guard | `null_guard`, called by all six `*_grid*` helpers |
