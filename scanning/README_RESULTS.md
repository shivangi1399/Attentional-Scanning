# `scanning/` — traveling-wave analyses of V4 LFP 

**This file is what the analyses found and what it means.**
For how the code works and what every quantity means, see **[`README.md`](README.md)** — it defines `PGD`, `R`, `R0`, the gain, the slope test, the two estimators, the three stimulus modes, and the combination rules used throughout this file. Nothing here re-explains method.

**Hypothesis.** Attention scans the stimulus locations, so the optimal (preferred) phase should shift systematically across cortical space — a traveling wave whose origin moves with the stimulated retinotopic patch.

Two animals (`hermes`, `klecks`), 64-channel 8×8 V4 arrays, `cp10_till_100`, `dv = 'lfp'`. Animals are never pooled at the channel level.

---

## Overall answer

> There **is** a robust, replicated, planar spatial phase gradient across V4 (theta + beta, both animals). In almost every band its wavenumber is **constant across frequency** — a fixed phase offset, not a propagating signal — and nothing about it is **organised by the scanned stimulus**. Cross-animal replication of any wave-consistent effect is essentially zero.
>
> **Planar phase structure: yes. Traveling wave: not supported, with one unresolved exception (hermes beta). Scanning wave: unsupported — now in the time domain as well as in phase.**

Six tests, each asking a stricter question than the last:

| test | question | result |
|---|---|---|
| #1 phase progression | does preferred phase depend on stimulus position? | weak — 2/35 and 3/35 frequencies |
| #2 PGD | are the phase-gradient arrows aligned across cortex? | **positive, replicated** — theta + beta |
| #3 origin | is the focus at the RF-driven patch, and does it follow it? | null (p ≈ 0.36–0.90) |
| #4 cortical de-rotation | does the ramp scale with frequency (one `v` for all `f`)? | negative — 0 replicated R, 1/1050 gain cell |
| #5 stimulus de-rotation | does the phase ramp with stimulus distance? | negative — **0 replicated, 3 modes × 2 estimators** |
| #6 ERP latency | is there a real TIME delay across cortex / with eccentricity? | delays exist and are significant, but `k_obs` does not scale with them (klecks, hermes theta); retinotopic slopes have **opposite signs** across animals |

**Test #6 is the independent one.** Tests #1–#5 all read the same `pref_phase` matrix; #6 goes back to the time-domain LFP and measures delay directly. It agrees: the phase gradients are not conduction delays, and nothing sweeps with eccentricity consistently.

**Why #2 passes and #4 fails is not a contradiction** — PGD is scale-free and single-frequency, so a fixed phase offset scores perfectly on it. [Full reconciliation below.](#why-pgd-is-positive-but-every-de-rotation-test-is-negative)

**The ceiling objection is closed for #5, and structurally unclosable for #4.** The obvious complaint was that `R0` sat at 0.80–0.99, leaving nothing for a wave to win:

- **#5 (stimulus):** the `coherence` estimator drops `R0` to an ordinary phase coherence (0.06–0.08), so the headroom is real. De-rotation still buys nothing — `max R ≤ max R0` in every mode and both animals. **The negative result is not a baseline artifact.**
- **#4 (cortical):** the same switch does *not* lower `R0` (0.947 → 0.953, 0.9916 → 0.9916), for the structural reason given in [`README.md` §4](README.md). Its gain ceiling of 0.008 (klecks) is a property of the analysis. So #4's "significant" gain cells are real but negligible in size, and #4 rests on **replication and the diagonal geometry**, not on the gain magnitude.

---

## Results per test

### 1. `phase_progression.m` — does preferred phase depend on stimulus position?

**WEAK.** The channel-average is significant at only ~**2/35 frequencies (hermes)** and ~**3/35 (klecks)**, matching the flat per-channel picture. MUA / RT / hit_miss are weaker still than LFP.

> These are the numbers *after* the synchronised-permutation fix. The earlier broadband "significance" was an artifact of independent per-channel shuffles, which made the channel-average null far too tight.

### 2. `cortical_planar_wave_PGD.m` — is there a planar wave?

**POSITIVE — and it replicates.** A planar wave exists in **both animals** at **theta (~4–6 Hz)** and **beta (~13–25 Hz)**.

| readout | hermes | klecks |
|---|---|---|
| beta speed | ~28–51 cm/s (physically plausible) | ~28–51 cm/s |
| theta speed | ~5 cm/s (borderline / near-synchronous) | ~5 cm/s |
| phase tilt across array | 19–86° | 19–86° |
| propagation direction | ~10° | ~166° |

**Test #3 (origin), computed here: NULL.** The wave focus is not at the RF-driven patch (p ≈ 0.36–0.78) and does not track that patch across positions (p ≈ 0.42–0.90), in either animal.

> ⚠️ **Speed rises with frequency** — theta ≈ 5 cm/s, beta ≈ 28–51 cm/s. That is the constant-wavenumber signature, and PGD cannot see it (scale-free, one frequency at a time). Quantified by the [slope test](#the-slope-test-results), which also finds the one band that behaves differently.

### 3. `cortical_wave_type_classification.m` — what kind of pattern?

**PLANAR.** PGD beats its null; radial divergence, rotational curl and singularity counts all sit at or below theirs. Not radial, not rotational, not spiral.

### 4. `cortical_planar_wave_derotation.m` — does the cortical ramp scale with frequency?

**NEGATIVE — diagonal ridge, no replication.** Run of 2026-08-04, both estimators, `nPerm = 1000`, max-stat over the 35 × 30 × 24 grid.

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

**2. The best-fit speed rises with frequency, in proportion:**

```
   klecks   f  5.6 -> 25.5 Hz   (x4.6)      v  34.5 -> 202.4 cm/s   (x5.9)
   hermes   f 25.5 -> 80.0 Hz   (x3.1)      v  62.2 -> 166.3 cm/s   (x2.7)
```

`v ∝ f` means `k = 2πf/v` is **constant**. One fixed phase offset fits every frequency. A real wave is the opposite: one fixed *speed* fits every frequency, and `k` grows with `f`.

(This holds for the bands where *de-rotation* is significant. Run over PGD's own significant bands instead, the same test isolates one exception — see the [slope test](#the-slope-test-results).)

**3. The ceiling is severe, and here it cannot be engineered away.** `R0` medians are **0.947 (hermes)** and **0.9916 (klecks)** — electrodes are already ~95–99 % phase-aligned with nothing de-rotated. Klecks' entire gain ceiling is `1 − 0.9916 = 0.008`, so its 86–87 "significant" cells are statistically real and physiologically negligible.

> Contrast with §5: there, switching to the `coherence` estimator dropped `R0` from ~0.9 to ~0.07 and removed the ceiling objection. Here it does not — 0.947 → 0.953 and 0.9916 → 0.9916 — for the structural reason in [`README.md` §4](README.md). The ceiling is a property of the analysis, not of the estimator.

### 5. `stimulus_loc_traveling_wave.m` — does the wave track the scanned stimulus?

**NEGATIVE in all three modes, under BOTH estimators.** Run of 2026-08-04, `nPerm = 1000`, `alpha = 0.05`, max-stat corrected over the whole 35 × 30 = 1050-cell grid. Numbers are `all RF centres`; the `_validRF` variant is in the caveat at the end.

**Replication — the primary criterion — is ZERO everywhere.**

| estimator | mode | replicated R | replicated gain | pooled sig |
|---|---|---|---|---|
| `phase` | `visual` | **0** | **0** | 2 (thr 2.99) |
| `phase` | `visual_coherent` | **0** | **0** | 1 (thr 3.27) |
| `phase` | `visual_arrival` | **0** | **0** | 87 (thr 2.69) |
| `coherence` | `visual` | **0** | **0** | 40 (thr 3.45) |
| `coherence` | `visual_coherent` | **0** | **0** | 1 (thr 3.27) |
| `coherence` | `visual_arrival` | **0** | **0** | 395 (thr 2.80) |

Per animal, the gain never clears threshold in klecks at all, and in hermes only under `phase`:

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

**3. The pooled hits are the diagonal, not a wave.** `visual_arrival` produces the largest pooled counts (87 and 395 cells) and the largest z values (up to 7.4). Look at where they sit on the grid: a **diagonal ridge**, speed rising with frequency = constant wavenumber. A wave requires a horizontal band. Same signature as PGD's theta/beta split (§2) and the cortical de-rotation result (§4). Pooled significance confirms the ridge is not noise; it does not make it a wave. With `n = 2` a pooled cell can also be ~100 % one animal — the script prints the per-animal z at the peak cell for exactly this reason.

**Measured gain thresholds, for reference:**

| | `phase` | `coherence` |
|---|---|---|
| typical `R0` (visual, validRF) | 0.80 hermes / 0.93 klecks | 0.063 / 0.078 |
| gain threshold obtained | 0.031 / 0.008 | 0.0006 / 0.0007 |

> ⚠️ **`_validRF` caveat — re-run needed.** In the 2026-08-04 run the `phase` estimator's `visual_coherent` and `visual_arrival` modes were invalid under `RF_VALID_ONLY = true`: dropping the Extrapolated RF centres leaves NaN geometry, which an unrestricted `randperm` shuffled into the used channel set, making every permutation NaN and collapsing the threshold to `-inf` (hence a spurious 1050/1050 "significant"). **Fixed in code** — those two helpers now permute only among finite-geometry channels, and all six helpers call `null_guard`. **Do not quote `_validRF` numbers for those two modes until re-run.** Unaffected: everything in the tables above (all-RF), and the `coherence` estimator in both variants — its `_validRF` results agree with the all-RF ones (0 replicated, 0 sig-gain, 40 / 0 / 477 pooled).

### 6. `erp_latency_wave.m` — the time-domain cross-check

**Run 2026-08-05.** Trials pooled across all sessions (hermes 31, klecks 25; ~23 000 and ~16 000 trials). ERPs are clean — flat baseline, response onset ~60–70 ms, peak ~100–150 ms — so `t = 0` is the stimulus event and the latencies are meaningful.

**There ARE significant latency gradients, on every axis in both animals:**

| animal | axis | slope | implied speed | r | p |
|---|---|---|---|---|---|
| hermes | cortical | −1.47 ms/mm | 68 cm/s | −0.26 | 0.029 |
| klecks | cortical | −4.16 ms/mm | 24 cm/s | −0.48 | <0.001 |
| hermes | retinotopic | **+3.15 ms/deg** | 318 deg/s | +0.31 | 0.036 |
| klecks | retinotopic | **−4.13 ms/deg** | 242 deg/s | −0.46 | 0.001 |

Both cortical speeds are physiologically plausible, and both cortical slopes are negative — latency *falls* along each animal's own PGD direction, consistently across animals.

> `k_obs` is a magnitude and can never be negative; `k_pred` carries the latency slope's sign. Only `|k_pred|` vs `k_obs` is meaningful — see [`README.md` §6](README.md).

#### What it says about the phase gradient

The decisive comparison is **scaling**: a fixed time delay forces `k ∝ f`.

```
   band                    f       k_obs              |k_pred|        k_obs/|k_pred|
   hermes  2.2- 6.1 Hz   x2.75   x0.82 (FALLS)        x2.75           29.4 -> 8.7
   hermes 13.3-25.5 Hz   x1.92   x2.05 (rises ~ f)    x1.92            2.4 -> 2.5
   klecks  4.4-23.3 Hz   x5.25   x1.09 (FLAT)         x5.25            0.8 -> 0.2
```

**1. Klecks is settled: the gradient is not a conduction delay.** Its `k_obs` is flat (0.094 → 0.102 rad/mm) across a 5.25× frequency range, while a real time delay of the measured size would have driven `|k_pred|` up by the same 5.25×. This is the textbook fixed-offset signature, now measured directly in the time domain rather than inferred from the de-rotation grid.

**2. Hermes theta is settled the same way** — `k_obs` actually *falls* while `|k_pred|` rises 2.75×.

**3. Hermes beta is not confirmed, and not cleanly refuted.** It is the one band where `k_obs` scales like a wave (×2.05 for a ×1.92 frequency change). But the magnitudes disagree by a consistent factor of ~2.4:

```
   from the phase gradient   28 cm/s   (3.6 ms/mm)   <- the prediction this test was built for
   from the ERP latency      68 cm/s   (1.5 ms/mm)
```

So the ERP does find a significant, plausibly-scaled timing gradient along that axis — just not the same speed. Two readings are open: the evoked broadband transient and the ongoing beta phase gradient may simply be different phenomena, or one of the two speed estimates is biased. The measurement is not sharp enough to choose: channel latency sd is 5.1 ms (hermes) / 7.9 ms (klecks) against a ~12 ms predicted spread across the 3.5 mm array, and `r` = −0.26 explains only 7 % of the variance.

**4. The scanning hypothesis fails again, now in the time domain.** The retinotopic slopes have **opposite signs** in the two animals: hermes' response arrives later at greater eccentricity (+3.15 ms/deg), klecks' arrives earlier (−4.13 ms/deg). Each is individually significant; together they describe no consistent foveal→peripheral sweep. This is an independent test of §5's conclusion, using time rather than phase, and it agrees.

> ⚠️ Caveats: the ±30 ms cross-correlation search bound was reached on 1 hermes channel (visible as the outlier in the figure). Correlations are modest throughout (|r| = 0.26–0.48). Sessions are pooled, which adds noise but no slope bias.

### 7. `traveling_wave_H2_H1.m` — deprecated

Descriptive only; superseded by the two `cortical_*` scripts. Still useful for visualisation.

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

### Direction is not propagation

PGD is scale-free and single-frequency ([why](README.md#pgd)), so a **fixed, frequency-independent spatial phase offset** scores a perfect PGD at every frequency while being the opposite of a traveling wave.

Your data is unambiguously that second case. Three independent readings agree:

| evidence | reading |
|---|---|
| PGD's own per-frequency speeds: theta ≈ 5, beta ≈ 28–51 cm/s | `v` rising with `f` |
| §4 klecks: `f` ×4.6 → `v` ×5.9; hermes: `f` ×3.1 → `v` ×2.7 | `v ∝ f`, so `k` const |
| §4/§5 grids: diagonal ridge, never a horizontal band | `k` const |

PGD was never wrong. It correctly detected an **aligned spatial phase gradient**. De-rotation adds the requirement that the gradient *scale with frequency*, and that is what fails.

**This also rules out a whole class of explanations.** Anything based on signals propagating — axonal conduction, attention sweeping across cortex, any finite-speed mechanism — produces a constant *time* delay and therefore `k ∝ f`. Constant `k` is incompatible with all of them. What remains are stationary explanations: a fixed anatomical/laminar phase offset across the array, reference or volume-conduction structure, or a genuine standing (non-propagating) phase pattern.

### The slope test: results

Run over PGD's own cluster-significant bands ([how it works](README.md#the-slope-test--what-turns-pgd-into-a-wave-claim); figure `slope_test_per_animal.pdf`):

| animal | band | v across the band | slope_v | slope_k | verdict |
|---|---|---|---|---|---|
| hermes | 2.2–6.1 Hz | 2.3 → 7.8 cm/s | **+1.25** | −0.25 | fixed phase offset |
| hermes | **13.3–25.5 Hz** | **28.8 → 26.9 cm/s** | **−0.10** | **+1.10** | **constant speed — wave-like** |
| klecks | 4.4–23.3 Hz | 29.8 → 143.8 cm/s | **+1.42** | −0.42 | fixed offset (script labels it *steeper than constant-k*) |

> ⚠️ **Hermes' beta band is the exception, and it must not be flattened into the null.** It holds ~28 cm/s across a *doubling* of frequency, with `k` rising in proportion (0.291 → 0.596 rad/mm), and 28 cm/s is a physiologically plausible cortical conduction speed. That is the horizontal-band signature, in one animal, in one band.
>
> What stops it being a positive result:
> - **it does not replicate** — klecks over the same range gives `slope_v = +1.42`, the opposite;
> - **de-rotation does not confirm it** — hermes' significant gain sits at 25.5–80 Hz, not 13–25 Hz. Most likely the ceiling (`R0` = 0.947 → gain ceiling 0.05, threshold 0.0147) plus max-stat correction over 25 200 cells;
> - **the ERP latency test does not confirm it either** (§6) — it finds a real, significant timing gradient along that axis, but at **68 cm/s against the 28 cm/s this band predicts**, a consistent factor-2.4 mismatch;
> - **n = 1 animal, 1 band, 8 frequency bins.**
>
> Treat it as the one lead worth chasing, not as evidence of a scanning wave — it is in *cortical* coordinates and says nothing about stimulus position. §6 was built to settle it and did not: the band survives the scaling test in both measurements, and fails on magnitude agreement.

**Direction is also frequency-dependent, which constrains the explanation:**

```
   hermes    2.2- 4.4 Hz   dir ~ 350°
             5.0- 6.1 Hz   dir ~ 183°      <- ~170° flip (sign inversion?)
            13.3-25.5 Hz   dir ~   6-19°

   klecks    4.4 Hz  216°  ->  10.0 Hz  152°  ->  23.3 Hz  160°   (smooth ~75° rotation)
```

A rigid geometric artifact — array tilt, systematic depth gradient, reference structure — would give **one direction at every frequency**, because geometry does not know about frequency. It does not. So "it is just electrode geometry" is weakened, and different bands are doing genuinely different things. Something like two spatially fixed generators whose relative contribution shifts with frequency would produce both the constant `k` and the rotating direction.

**Honest status: the theta gradients are not explained.** Constant `k` rules out propagation; frequency-dependent direction rules out a single fixed geometric offset.

### PGD vs the stimulus test: different axes entirely

These two are not even measured along the same dimension.

```
   PGD / cortical de-rotation :  distance = mm across the ARRAY
   stimulus de-rotation       :  distance = degrees in the VISUAL FIELD
```

A phase gradient can exist across cortex without being organised by where the stimulus is. The scanning hypothesis needs the *stimulus* axis, and three separate tests say it is not there:

1. **Test #1** — preferred phase barely depends on stimulus position: 2/35 (hermes) and 3/35 (klecks) frequencies.
2. **Test #3** — the wave focus is not at the RF-driven patch (p ≈ 0.36–0.78) and does not follow it across positions (p ≈ 0.42–0.90).
3. **Test #5** — `R0` = 0.80/0.93 under `phase`: the per-location preferred phases are already 80–93 % aligned *before* any de-rotation, so there is almost no position-dependent phase for any wave model to organise. Under `coherence`, where that ceiling is gone, `max R ≤ max R0` in every mode and both animals.
4. **Test #6** — the retinotopic ERP latency slopes have **opposite signs** in the two animals (+3.15 vs −4.13 ms/deg). Each is significant on its own; together they describe no consistent sweep. This one uses *time*, not phase, so it is not a restatement of the others.

Test #1 and test #5's `R0` are the same fact stated twice — one as a correlation, one as a baseline. Test #6 is genuinely independent of both.

A different distance ruler cannot rescue this: `R0` is the `k = 0` baseline and involves no distance at all, so [no re-mapping can create phase variance the data doesn't contain](README.md#why-a-different-ruler-cannot-rescue-a-null).

### Putting it together

```
   there IS a spatial phase structure across V4        (PGD, replicated)
        it is PLANAR, not radial or spiral             (test #3 classification)
        it is STATIONARY, not propagating              (constant k, both scripts,
                                                        and no matching TIME delay, #6)
        it is NOT organised by stimulus position       (tests #1, #3-origin, #5,
                                                        and #6's retinotopic axis)

   there ARE real conduction delays across the array   (#6: 24-68 cm/s, both animals)
        but they do NOT account for the phase gradient  (k_obs does not scale with them)
```

So the honest summary is not "we found nothing". It is: **there is a robust, replicated, planar spatial phase gradient across V4, and there are real evoked conduction delays across the same array — but they are not the same thing, the gradient is a fixed pattern rather than a wave, and neither is tied to the scanned stimulus.**

---

## Open items

| item | status |
|---|---|
| Re-run `stimulus_loc_traveling_wave.m` with `RF_VALID_ONLY = true` | needed — `phase`/`visual_coherent` and `phase`/`visual_arrival` `_validRF` numbers are invalid until then (bug fixed in code) |
| Hermes beta, 13.3–25.5 Hz | **the one open lead.** Wave-like scaling in both the phase gradient and the ERP latency, but the speeds disagree 2.4× (28 vs 68 cm/s) and it does not replicate in klecks |
| The theta gradients | unexplained — constant `k` rules out propagation, the ERP latency test (§6) independently rules out a conduction delay, and frequency-dependent direction rules out fixed geometry |
| Why the two cortical speeds differ 3× | hermes 68 vs klecks 24 cm/s from the same measurement. Could be genuine, could be the modest fit quality (\|r\| = 0.26 vs 0.48) |
| Cortical arrival mode (the missing 2×2 cell) | blocked on a log-polar CMF calibration (`E2_DEG`, `A_MM`) that is not in the repo |
| `RECOMPUTE_TRIAL_SUMS = false` in both de-rotation scripts | optional — the worker is unchanged and the cache is intact, so this skips 128 SLURM jobs |
