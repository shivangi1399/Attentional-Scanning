% =====================================================================
% Traveling-wave detection by PHASE-VECTOR ALIGNMENT across stimulus
% locations, swept over frequency × speed ("rotate-to-overlap"
% test).
%
% Idea : each stimulus location p has a preferred-phase vector e^{iφ(p)}. 
% If processing scans the locations as a TRAVELING WAVE of speed v, then 
% the phase at each location is delayed relative to a reference by 
% Δt(p) = d(p)/v, i.e. rotated by
%       θ(p) = 2π f d(p) / v .
% De-rotating each vector by θ(p) should make all locations OVERLAP
% (resultant R → 1) — but only at the wave's true (f, v). So we sweep a
% grid of (frequency × speed) and read off the parameter region where the
% de-rotated vectors align. A GENUINE wave shows a band of high R at a
% (near-)constant speed ACROSS frequencies (speed is frequency-independent),
% which is the key discriminator from a chance per-frequency alignment.
%
% Everything is in VISUAL degrees / deg-per-second. d(p) = stimulus
% eccentricity (deg). Three ways to combine channels are computed:
%
%   VISUAL (incoherent) — the wave is across STIMULUS LOCATIONS only:
%     per channel:  R_c(f,v) = | Σ_p w_p e^{i(φ_p - θ_p)} | / Σ_p w_p ,
%                   θ_p = 2πf d_p / v  (w_p = coherence magnitude). The |·|
%                   is taken PER CHANNEL, so each channel's own constant phase
%                   offset cancels and channels can be averaged — but that
%                   also makes any per-electrode delay INVISIBLE.
%     grid       :  R(f,v) = mean over coherence-significant channels.
%
%   VISUAL_COHERENT — ONE cortical wave across BOTH stimulus locations AND
%     electrodes. Each channel also has an RF-centre eccentricity Dc_c (deg,
%     from the Gaussian RF fit), and the wave delays the electrode by k·Dc_c:
%       R(f,v) = | Σ_{c,p} w_cp e^{i(φ_cp − k d_p − k Dc_c)} | / Σ w_cp ,
%     summed COHERENTLY across channels (per-channel phase kept), so the
%     per-electrode delay is testable. Assumes phases share a common clock;
%     unmodelled offsets can only LOWER R, so a band is trustworthy but a
%     null is ambiguous. Kept alongside the incoherent grid as a comparison.
%
%   VISUAL_ARRIVAL (incoherent) — the ARRIVAL-TIME model: the stimulus lands
%     on its retinotopic patch and the information spreads outward at the SAME
%     speed v, so electrode c hears about location p only after a(c,p)/v, where
%     a(c,p) = visual-field separation between c's RF centre and location p:
%       R_c(f,v) = | Σ_p w_cp e^{i(φ_cp − k d_p − k a(c,p))} | / Σ_p w_cp ,
%       R(f,v)   = mean over coherence-significant channels of R_c .
%     Unlike Dc_c, a(c,p) VARIES across locations within a channel (it dips
%     where the stimulus falls on that channel's RF), so it survives the
%     per-channel |·| and can be fitted INCOHERENTLY — the only mode that
%     expresses "the wave starts where the stimulus is" while keeping
%     robustness to per-channel offsets and an interpretable null. Because
%     |·|_c deletes any channel constant, this tests the position-VARYING part
%     of the arrival delay (the row mean of a is projected out).
%
% Also computed for each mode: the de-rotation GAIN dR = R(f,v) − R0(f),
% where R0 is the un-rotated (k=0) coherence — the "actual" coherence with no
% wave assumed; dR isolates the wave-specific phase ramp.
%
% =====================================================================
% TWO ESTIMATORS — 'phase' (PHASE ALIGNMENT) vs 'coherence' (PHASE COHERENCE)
% =====================================================================
% Setting: ESTIMATORS = {'phase','coherence'} — BOTH are computed in one run.
%
% The three modes above describe WHAT is de-rotated (which distance goes into
% k*d). The ESTIMATOR describes WHAT IS BEING VECTOR-SUMMED, and therefore what
% the number R actually means. This is the single most important distinction in
% this script, so it is spelled out here in full.
%
% -------------------------------------------------------------------------
%  BOTH estimators share, identically:
%       the de-rotation          k = 2*pi*f/v, spin by k*d
%       the three modes          visual / visual_coherent / visual_arrival
%       the baseline             R0 = the same statistic at k = 0
%       the gain                 dR = R - R0        <-- COMMON TO BOTH.
%                                                       dR is NOT what
%                                                       distinguishes them.
%       the null                 synchronised location-label shuffle
%       thresholds, replication, pooling, figures
%  They differ in ONE thing: what the vectors being summed are.
% -------------------------------------------------------------------------
%
%   ESTIMATOR 'phase'  —  PHASE ALIGNMENT  (two-stage; the original estimator)
%     Stage 1 (upstream, phase_progression.m): per-location coherence across
%       trials,  c(ch,f,p) = mean over trials at location p of y*exp(i*phase)
%       -> pref_phase = angle(c),  coh_mag = |c|.
%     Stage 2 (here): vector-sum those nPos PER-LOCATION PREFERRED-PHASE
%       vectors, one per location:
%         R_c(f,v) = | SUM_p coh_mag * exp(i(pref_phase - k*d_p)) | / SUM_p coh_mag
%
%     *** R IS NOT A PHASE COHERENCE. *** It is a weighted RESULTANT LENGTH
%     over locations — a circular-concentration measure answering "how
%     similar are these nPos angles after de-rotation?". The trial-level
%     coherence enters only as the WEIGHT coh_mag; the trials themselves have
%     already been averaged away in stage 1.
%     Consequence: R0 = "how similar were the preferred phases already", which
%     runs ~0.80-0.93 in this data. High by construction (a handful of angles
%     are easily concentrated), leaving a gain ceiling of only 0.07-0.20.
%
%   ESTIMATOR 'coherence'  —  PHASE COHERENCE  (one-stage, trial-level)
%     No per-location stage at all. Every TRIAL is rotated by the model's
%     prediction for the location THAT TRIAL had, and ONE coherence is taken
%     over the whole trial set:
%       c(f,v) = (1/W) * SUM over ALL TRIALS y_t*exp(i(phi_tf - k*d_p(t))) ,
%       W = SUM_t |y_t| ,   R = |c| .
%
%     *** R IS A GENUINE PHASE COHERENCE *** — at k = 0 it is EXACTLY the
%     measure computed by Phase_coherence/functions/phase_coherence.m, with all
%     stimulus locations pooled (only the normaliser differs: /W instead of
%     /nTrials, so R is bounded in [0,1]).
%     Consequence: R0 is a coherence in the familiar sense (~0.06-0.08 here),
%     so there is real headroom for a gain to appear.
%
%     NOTE: no per-location coherences are ever computed or compared. Locations
%     enter only through which rotation each trial receives.
%
%   Side by side, at the level of the sum:
%       'phase'      SUM over nPos LOCATION VECTORS  (trials pre-averaged)
%       'coherence'  SUM over ALL TRIALS             (no pre-averaging)
%
%   And the one-line algebraic difference:
%       'phase'      uses  S(p,f)/n_p   the per-location MEAN, weighted by
%                                       |c_p| ~ 1/sqrt(n_p)  -> a location with
%                                       FEWER trials gets a LARGER weight
%       'coherence'  uses  S(p,f)       the raw SUM -> every trial counts once
%
%   R VALUES ARE NOT COMPARABLE BETWEEN THE TWO. Judge each against its own
%   threshold. A lower R0 under 'coherence' is good news (headroom), not a
%   weaker result.
%
%   Implementation of 'coherence': the de-rotation factor depends on a trial
%   only through its location, so the estimator collapses onto per-location
%   complex SUMS
%       S(p,f) = SUM over trials at p of y_t*exp(i*phi_tf)
%       c(f,v) = (1/W) * SUM_p S(p,f)*exp(-i*k*d_p) ,
%   and all three modes reduce to re-weightings of S. One SLURM job per channel
%   (functions/trial_position_sums_chan.m) produces S for the observed labels
%   and for all nPerm shuffles; everything else is done here from S.
%   (The worker and its on-disk cache keep the name 'trial_position_sums'
%   because that is literally what they hold — per-TRIAL sums grouped by
%   position. Only the ESTIMATOR was renamed 'trial' -> 'coherence'.)
%
% Pipeline, PER ANIMAL (never pooled): grid -> animals -> mean grid +
%   "significant in both" replication.
% Significance: permute the location<->distance assignment (plus electrode<->Dc
%   for the coherent mode, or electrode<->a(c,:) rows for the arrival mode)
%   with a SINGLE synchronised shuffle shared
%   across the grid, recompute, take its max -> max-statistic threshold that
%   corrects across the entire (f,v) grid. R0 is shuffle-invariant, giving a
%   clean max-stat null on the gain dR too.
%
% COORDINATE FRAME (resolved from code/RF_Mapping/chan_loc_mua.m):
%   The stimulus positions and the electrode RF centres are stored with
%   DIFFERENT origins in the raw files, so both are converted to ONE common
%   frame here — FIXATION-CENTRED DEGREES, with (0,0) = fovea:
%     stimulus   trialinfo col-16/17 are ALREADY fixation-centred pixels
%                (chan_loc_mua.m:212 does x_screen = x_target_pix + 840),
%                so:   xy_deg = [col16, col17] / ppd.
%     RF centre  the *_channel_target_summary.txt stores SCREEN pixels
%                (the mapping code wrote xRF_screen = ... + 840), so the
%                screen centre must be SUBTRACTED first:
%                      xy_deg = ([RF_Center_X, RF_Center_Y] - [840 525]) / ppd.
%   Both then yield eccentricity = hypot(x,y) deg measured from the same
%   origin, so any electrode-minus-stimulus difference is frame-safe.
%   ppd = PIX_PER_DEG, the PER-ANIMAL rig calibration (the raw .RF sessInfo
%   did not store ppd). VISUAL distance d(p) = stimulus eccentricity (deg).
%
% RF_VALID_ONLY toggle (affects VISUAL_COHERENT and VISUAL_ARRIVAL — the two
%   modes that use RF centres; plain VISUAL never does and is unchanged by it):
%   keep only channels whose RF centre is a real
%   2-D Gaussian fit (Status == 'Valid_Gaussian'). Channels marked
%   'Extrapolated' had a FAILED fit and their centre was filled in from the
%   median of their 8-connected array neighbours, or from a plane fit across
%   the array (RF_Mapping/mapping_lfp.m:513) — i.e. inferred under a
%   retinotopic-smoothness assumption rather than measured. Output file names
%   and figure titles carry a '_validRF' tag when the toggle is on, so both
%   variants can coexist on disk.
%
% OUTPUT — THREE figures per run. Names carry the '_validRF' tag when
% RF_VALID_ONLY is on; the two gain figures also carry the estimator name.
%   Plots/scanning/stimulus_loc_wave/cp10_till_100/<dv>/
%     stimulus_loc_grids[_validRF].pdf
%         THE DE-ROTATION R. 6 rows = 3 modes x 2 estimators (estimator named
%         in every panel title, because R means different things in the two
%         blocks), 4 cols = animals + mean/replication + pooled z.
%         The pooled z column lives HERE ONLY — see note below.
%     stimulus_loc_gain_phase[_validRF].pdf
%         Gain dR = R - R0 for the PHASE-ALIGNMENT estimator. 3 rows x 3 cols.
%     stimulus_loc_gain_coherence[_validRF].pdf
%         Gain dR = R - R0 for the PHASE-COHERENCE estimator. 3 rows x 3 cols.
%   results_combined/scanning/stimulus_loc_wave/cp10_till_100/<dv>/
%       stimulus_loc_wave[_validRF].mat (BOTH estimators in one file:
%                                        results.G.(estimator).(mode)(animal))
%
% WHY THE GAIN FIGURES HAVE NO POOLED COLUMN. Because gain = R - R0(f) and R0
% depends on neither speed nor the shuffle, the pooled standardised statistics
% are algebraically IDENTICAL:
%       muG = muR - R0 ,  sdG = sdR
%       zG  = (R - R0 - muG)/sdG = (R - muR)/sdR = zR
% (verified on saved output: max|zR - zG| = 1.75e-13, identical thresholds).
% A pooled gain panel would therefore be a pixel-for-pixel copy of the pooled R
% panel, so it is drawn once, on the grids figure. The PER-ANIMAL R and gain
% panels are genuinely different and both are kept.
% =====================================================================

clearvars; close all; clc

%% ─── Dependencies (only the 'coherence' estimator needs slurmfun) ────
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/scanning/functions

%% ─── Settings ────────────────────────────────────────────────────────
animals    = {'hermes','klecks'};
dv         = 'lfp';
base       = '/mnt/hpc/projects/MWSampling/4Shivangi';
nPerm      = 1000;
alpha      = 0.05;
CH_FILTER  = 'significant';        % coherence-sig channels only (as in planar wave)
COH_SIG_ALPHA = 0.05;
MIN_LOC    = 4;                    % min valid stimulus locations per channel
FREQ_RANGE = [2 100];               % Hz swept (low-freq band the hypothesis lives in)
NV         = 30;                   % # speed steps
V_VISUAL   = logspace(log10(1),   log10(200), NV);   % deg/s
rng(2025);

% Experiment calibration
PIX_PER_DEG   = struct('hermes', 53.24, 'klecks', 50.56);   % pixels per degree (per-animal rig calibration)

% Electrode RF centres (for the COHERENT mode). Each channel's RF centre comes
% from the Gaussian-fit centre already written by RF_Mapping/mapping_lfp.m into
% the per-channel summary table (the SAME file cortical_planar_wave_PGD.m
% reads). Its RF_Center_X/Y are SCREEN pixels, so elec_rf_deg subtracts the
% screen centre (840, 525) before /ppd to reach the fixation-centred degrees
% the stimulus positions already use — see the COORDINATE FRAME note above.
% Both animals are 8x8 = 64-channel arrays and the summary has 64 rows 1:1
% with the phase-progression channel index, so no channel remap is needed.
RF_DATE   = '20170829';   % RF session both animals share
SCREEN_XY = [1680 1050];  % screen pixels; fixation/fovea at the centre

% Use ONLY channels whose RF centre is a real Gaussian fit (Status ==
% 'Valid_Gaussian')? 'Extrapolated' centres were interpolated from array
% neighbours after a failed fit, not measured (see header). Affects
% visual_coherent (via Dc) and visual_arrival (via a(c,p)); plain visual never
% uses RF centres and is unchanged either way.
RF_VALID_ONLY = true;

% WHICH ESTIMATORS TO RUN — see the TWO ESTIMATORS block in the header.
% Both are computed in a SINGLE run; the grids figure shows them side by side
% and each gets its own gain figure.
%   'phase'     PHASE ALIGNMENT.  R = resultant length of the nPos per-location
%               preferred-phase vectors (from phase_progression.mat).
%               NOT a phase coherence — a circular-concentration measure.
%   'coherence' PHASE COHERENCE.  R = coherence over ALL trials after each
%               trial is de-rotated by its own location's predicted lag. At
%               k = 0 this IS the Phase_coherence/ measure, locations pooled.
%               Needs ph_all_sess.mat and one SLURM job per channel.
% Order matters only for figure row order. Set to a single name to run one.
ESTIMATORS = {'phase','coherence'};

% Force re-submission of the per-channel SLURM jobs that build the per-location
% trial sums S(p,f) used by the 'coherence' estimator, even if a cache exists
% (set true after changing dv, nPerm, or functions/trial_position_sums_chan.m).
RECOMPUTE_TRIAL_SUMS = true;

if RF_VALID_ONLY
    rf_tag  = '_validRF';
    rf_note = 'RF: Valid_Gaussian only';
else
    rf_tag  = '';
    rf_note = 'RF: all centres incl. Extrapolated';
end
% Per-estimator naming used in filenames, panel titles and console output. The
% short label goes in panel titles (space is tight); the long one in sgtitles.
EST_INFO = struct( ...
    'phase',     struct('tag','_phase',     'short','PHASE ALIGNMENT', ...
        'long','phase estimator = PHASE ALIGNMENT: R = resultant of the per-location preferred-phase vectors (NOT a phase coherence)'), ...
    'coherence', struct('tag','_coherence', 'short','PHASE COHERENCE', ...
        'long','coherence estimator = PHASE COHERENCE: R = coherence over all trials, each de-rotated by its own location (the Phase_coherence/ measure at k=0)'));
for ei = 1:numel(ESTIMATORS)
    if ~isfield(EST_INFO, ESTIMATORS{ei}), error('Unknown ESTIMATOR: %s', ESTIMATORS{ei}); end
end
tag  = rf_tag;                         % shared suffix (grids figure + .mat)
note = rf_note;                        % shown in the figure sgtitles

out_dir = fullfile(base,'Plots','scanning','stimulus_loc_wave','cp10_till_100', dv);
res_dir = fullfile(base,'results_combined','scanning','stimulus_loc_wave','cp10_till_100', dv);
if ~exist(out_dir,'dir'), mkdir(out_dir); end
if ~exist(res_dir,'dir'), mkdir(res_dir); end

% Three combination modes, all in VISUAL degrees / deg-per-second. They differ
% ONLY in the distance each de-rotates by; null, thresholds and figures are
% shared. Writing the de-rotation as a table over (channel c, location p):
%
%   visual          : INCOHERENT — de-rotate by k*d_p. The same number for
%                     every channel at a given location. |resultant| per
%                     channel over locations, then averaged (per-channel phase
%                     offset discarded; robust but blind to per-electrode
%                     delays).
%   visual_coherent : COHERENT — de-rotate by k*(d_p + D_c), D_c = electrode
%                     RF eccentricity (deg). D_c has no p in it, so its row of
%                     the table is FLAT: a channel-constant phase, which a
%                     per-channel |.| would delete exactly — hence this mode
%                     MUST sum coherently across channels. Tests ONE wave
%                     sweeping outward from the fovea, independent of where
%                     the stimulus is. Assumes a common phase reference, so a
%                     null is ambiguous (offsets can only LOWER R).
%   visual_arrival  : INCOHERENT — de-rotate by k*(d_p + a(c,p)), where
%                     a(c,p) = visual-field separation between electrode c's
%                     RF centre and stimulus location p (deg). This is the
%                     ARRIVAL-TIME model: the stimulus lands on its
%                     retinotopic patch and the information spreads outward at
%                     the same speed v, so electrode c hears about location p
%                     after a(c,p)/v. Each channel's row DIPS at the location
%                     that falls on its own RF, so the term varies across
%                     locations within a channel and SURVIVES the per-channel
%                     |.| — the only mode that can express "the wave starts
%                     where the stimulus is" while staying robust to
%                     per-channel offsets. Note |.|_c projects out the row
%                     mean of a, so this tests the position-VARYING part of
%                     the arrival delay.
metrics = {'visual','visual_coherent','visual_arrival'};
Vsets   = {V_VISUAL, V_VISUAL, V_VISUAL};
Vunit   = {'deg/s','deg/s','deg/s'};
cmode   = {'incoherent','coherent','incoherent'};

G = struct();   % G.(metric)(ia): grid, threshold, sig, freq, speeds

%% ─── Per animal ──────────────────────────────────────────────────────
for ia = 1:numel(animals)
    animalName = animals{ia};
    pp = fullfile(base, ['results_' animalName], 'scanning', ...
        'phase_progression','cp10_till_100', dv, 'phase_progression.mat');
    if ~isfile(pp), warning('No data for %s — skipping.', animalName); continue; end
    S = load(pp, 'pref_phase','coh_mag','freq','positions');
    freq = S.freq(:); positions = S.positions(:);
    [nCh,nFreq,nPos] = size(S.pref_phase);
    fprintf('\n=== %s / %s : %d ch, %d freq, %d pos ===\n', animalName, upper(dv), nCh, nFreq, nPos);

    % coherence-significance mask [nCh x nFreq], identical to planar wave
    if strcmp(CH_FILTER,'significant')
        coh_root = fullfile(base, ['results_' animalName], 'phase_coherence', ...
            'complex','cp10_till_100', dv, 'all_loc_difflev');
        coh_sig = load_coh_sig_mask(coh_root, nCh, nFreq, COH_SIG_ALPHA);
    else
        coh_sig = true(nCh, nFreq);
    end

    f_use = find(freq >= FREQ_RANGE(1) & freq <= FREQ_RANGE(2));
    fHz   = freq(f_use);

    % ── Put stimulus and RF centres in ONE frame: fixation-centred degrees ──
    % Stimulus: trialinfo col 16 (== the 'positions' values) and col 17 are
    % ALREADY fixation-centred pixels, so only /ppd is needed.
    [x_px, y_px] = target_pix_per_pos(positions, base, animalName);
    ppd      = PIX_PER_DEG.(animalName);          % pixels per degree (this animal)
    stim_deg = [x_px(:), y_px(:)] / ppd;          % nPos x 2, (0,0) = fovea
    ecc      = hypot(stim_deg(:,1), stim_deg(:,2));   % eccentricity (deg)
    d_vis    = ecc - min(ecc);                    % linear visual distance (deg)
    fprintf('  ppd=%.2f px/deg | eccentricity %.2f-%.2f deg\n', ppd, min(ecc), max(ecc));

    % Electrode RF centres in the SAME frame (the summary table holds SCREEN
    % pixels, so the screen centre is subtracted inside elec_rf_deg before
    % /ppd). Used by the coherent mode only; referenced like d_p (subtract
    % min — a constant offset is a global phase rotation and cannot change R).
    [rf_deg, rf_valid] = elec_rf_deg(base, animalName, nCh, RF_DATE, ppd, SCREEN_XY);
    ecc_elec = hypot(rf_deg(:,1), rf_deg(:,2));   % eccentricity (deg), NaN if no centre
    if RF_VALID_ONLY
        ecc_elec(~rf_valid) = NaN;                % drop Extrapolated / No_Data centres
    end
    Dc = ecc_elec - min(ecc_elec, [], 'omitnan');
    fprintf('  electrode RF ecc [%s]: %d/%d channels usable (%d Valid_Gaussian) | %.2f-%.2f deg\n', ...
        rf_note, sum(isfinite(ecc_elec)), nCh, sum(rf_valid), ...
        min(ecc_elec,[],'omitnan'), max(ecc_elec,[],'omitnan'));

    % ARRIVAL-DISTANCE table a(c,p) for the visual_arrival mode [nCh x nPos]:
    % visual-field separation between electrode c's RF centre and stimulus
    % location p — a straight difference, legitimate only because both are now
    % in the SAME fixation-centred frame (see the COORDINATE FRAME note). Row c
    % dips at the location that lands on channel c's RF. No per-row referencing
    % is applied: a row-constant offset is deleted by the per-channel |.|
    % anyway, so it cannot affect R.
    a_arr = hypot(rf_deg(:,1) - stim_deg(:,1).', rf_deg(:,2) - stim_deg(:,2).');
    if RF_VALID_ONLY, a_arr(~rf_valid, :) = NaN; end
    fprintf('  arrival distance a(c,p): %d/%d channels usable | %.2f-%.2f deg\n', ...
        sum(all(isfinite(a_arr),2)), nCh, min(a_arr(:),[],'omitnan'), max(a_arr(:),[],'omitnan'));

    % 'coherence' estimator: fetch the per-location complex trial sums S(p,f)
    % (one SLURM job per channel, cached on disk). The 'phase' estimator never
    % touches these — it works entirely from phase_progression.mat.
    if any(strcmp(ESTIMATORS,'coherence'))
        [Ssum, Sperm, Wch] = get_trial_sums(base, animalName, dv, nCh, nPerm, RECOMPUTE_TRIAL_SUMS);
    end

  for ei = 1:numel(ESTIMATORS)
    est = ESTIMATORS{ei};
    fprintf('  --- estimator: %s (%s) ---\n', est, EST_INFO.(est).short);
    for mi = 1:numel(metrics)
        d = d_vis(:).'; Vs = Vsets{mi};
        % The ONLY place the two estimators diverge. Everything below this
        % switch — thresholds, standardisation, gain, pooling — is identical.
        %   'phase'     sums nPos per-location preferred-phase vectors
        %   'coherence' sums all trials (via the cached per-location sums S)
        switch [est '/' metrics{mi}]
            case 'phase/visual'
                [Robs, R0, Rnull_max, Gnull_max, Rnull, Gnull] = align_grid(S.pref_phase, S.coh_mag, coh_sig, ...
                    f_use, fHz, d, Vs, MIN_LOC, nPerm);
            case 'phase/visual_coherent'
                [Robs, R0, Rnull_max, Gnull_max, Rnull, Gnull] = align_grid_coherent(S.pref_phase, S.coh_mag, coh_sig, ...
                    f_use, fHz, d, Dc, Vs, MIN_LOC, nPerm);
            case 'phase/visual_arrival'
                [Robs, R0, Rnull_max, Gnull_max, Rnull, Gnull] = align_grid_arrival(S.pref_phase, S.coh_mag, coh_sig, ...
                    f_use, fHz, d, a_arr, Vs, MIN_LOC, nPerm);
            case 'coherence/visual'
                [Robs, R0, Rnull_max, Gnull_max, Rnull, Gnull] = coh_grid(Ssum, Sperm, Wch, coh_sig, ...
                    f_use, fHz, d, Vs, nPerm);
            case 'coherence/visual_coherent'
                [Robs, R0, Rnull_max, Gnull_max, Rnull, Gnull] = coh_grid_coherent(Ssum, Sperm, Wch, coh_sig, ...
                    f_use, fHz, d, Dc, Vs, nPerm);
            case 'coherence/visual_arrival'
                [Robs, R0, Rnull_max, Gnull_max, Rnull, Gnull] = coh_grid_arrival(Ssum, Sperm, Wch, coh_sig, ...
                    f_use, fHz, d, a_arr, Vs, nPerm);
            otherwise
                error('Unhandled ESTIMATOR/mode: %s', [est '/' metrics{mi}]);
        end
        thr = quantile(Rnull_max, 1-alpha);
        sig = Robs >= thr;

        % ── Standardise against this animal's OWN null, cell by cell ──────
        % The animals sit on very different scales (R0 ~0.80 vs ~0.93, and
        % thresholds differ 4-6x), so a raw cross-animal average would simply
        % be dominated by the animal with the larger numbers. Converting each
        % cell to a z-score against its own permutation distribution puts both
        % animals in comparable units, which is what the pooled test below
        % averages. The same permutations are used to centre/scale and to build
        % the pooled threshold — standard practice for a standardised
        % ("pseudo-t") max-stat permutation test.
        % 'omitnan' is essential: mean/std do NOT skip NaN by default, so a
        % single NaN permutation at one cell would make that cell's z NaN for
        % EVERY permutation, and the NaN then spreads through the cross-animal
        % sum below (NaN + finite = NaN).
        muR = mean(Rnull,3,'omitnan'); sdR = std(Rnull,0,3,'omitnan');
        muG = mean(Gnull,3,'omitnan'); sdG = std(Gnull,0,3,'omitnan');
        zR_obs  = (Robs  - muR) ./ max(sdR, eps);
        zR_null = (Rnull - muR) ./ max(sdR, eps);
        zG_obs  = ((Robs - R0(:)) - muG) ./ max(sdG, eps);
        zG_null = (Gnull - muG) ./ max(sdG, eps);
        fprintf('  %-16s: finite cells in standardised grid = %d/%d\n', ...
            metrics{mi}, sum(isfinite(zR_obs(:))), numel(zR_obs));

        % de-rotation gain vs the un-rotated ("actual") coherence + its
        % max-stat null: dR(f,v)=R(f,v)-R0(f); a wave shows a significant
        % positive gain at a (near-)constant best speed across frequencies.
        gain     = Robs - R0(:);                 % nF × nV (broadcast per row)
        thr_gain = quantile(Gnull_max, 1-alpha);
        sig_gain = gain >= thr_gain;
        [gpk, vbest] = max(Robs, [], 2);         % best speed index per freq
        vbest_speed  = Vs(vbest(:));             % best de-rotation speed per freq
        vbest_sig    = (gpk - R0(:)) >= thr_gain;% is that gain significant?

        m = metrics{mi};
        G.(est).(m)(ia).animal = animalName;
        G.(est).(m)(ia).estimator = est;
        G.(est).(m)(ia).R = Robs; G.(est).(m)(ia).thr = thr;
        G.(est).(m)(ia).sig = sig; G.(est).(m)(ia).fHz = fHz;
        G.(est).(m)(ia).speeds = Vs;
        G.(est).(m)(ia).R0 = R0(:);              % no-rotation baseline per freq
        G.(est).(m)(ia).gain = gain;             % de-rotation gain grid
        G.(est).(m)(ia).thr_gain = thr_gain;
        G.(est).(m)(ia).sig_gain = sig_gain;
        G.(est).(m)(ia).vbest = vbest_speed;     % best speed per freq
        G.(est).(m)(ia).vbest_sig = vbest_sig;
        G.(est).(m)(ia).zR_obs = zR_obs;         % standardised grids, for pooling
        G.(est).(m)(ia).zG_obs = zG_obs;
        G.(est).(m)(ia).zR_null = zR_null;       % nF x nV x nPerm (stripped before save)
        G.(est).(m)(ia).zG_null = zG_null;

        [rmax,ix] = max(Robs(:)); [fi,vi] = ind2sub(size(Robs), ix);
        fprintf('    %-16s: peak R=%.4f at f=%.1f Hz, v=%.1f %s | R0 range %.4f-%.4f | thr=%.4f | sig cells=%d\n', ...
            m, rmax, fHz(fi), Vs(vi), Vunit{mi}, min(R0,[],'omitnan'), max(R0,[],'omitnan'), thr, sum(sig(:)));
        [gmax,gix] = max(gain(:)); [gfi,gvi] = ind2sub(size(gain), gix);
        fprintf('    %-16s: peak gain dR=%.4f at f=%.1f Hz, v=%.1f %s | thr_gain=%.4f | sig-gain cells=%d | freqs w/ sig gain=%d/%d\n', ...
            m, gmax, fHz(gfi), Vs(gvi), Vunit{mi}, thr_gain, sum(sig_gain(:)), sum(vbest_sig), numel(vbest_sig));
    end
  end
end
% Animals that produced grids (same set for every estimator/mode).
valid = find(arrayfun(@(s) ~isempty(s.animal), G.(ESTIMATORS{1}).(metrics{1})));

%% ─── Combine animals: (a) REPLICATION, (b) POOLED standardised test ──
% Two different criteria, reported side by side. They answer different
% questions and are NOT interchangeable:
%
%   REPLICATION  repl = sig(hermes) AND sig(klecks), cell by cell. Each animal
%                was already max-stat corrected over its whole grid, so a
%                replicated cell is very strong evidence. It cannot be driven
%                by one animal — but it also cannot aggregate weak evidence,
%                and it has low power. PRIMARY criterion.
%
%   POOLED       average the STANDARDISED grids across animals and test that
%                average against its own max-stat null, built by pairing
%                permutation b of one animal with permutation b of the other
%                (legitimate because the animals are independent, so the
%                pairing is arbitrary and samples the product null):
%                    Zobs(f,v)   = mean_a z_a(f,v)
%                    Znull(f,v,b)= mean_a z_a_null(f,v,b)
%                    thr_pool    = (1-alpha) quantile of max over grid of Znull
%                This AGGREGATES evidence, so it has more power — but with only
%                two animals a "pooled significant" cell can be ~100% one
%                animal. SECONDARY criterion; always report which animal drives
%                it (the per-animal z at that cell is printed below).
%
% NOTE: with n=2 neither criterion supports a population-level inference —
% between-animal variance is not estimable. Replication is simply the more
% conservative descriptive rule.
% Done SEPARATELY FOR EACH ESTIMATOR — the two never share a null or a
% threshold, because their R values live on different scales.
C = struct();
for ei = 1:numel(ESTIMATORS)
  est = ESTIMATORS{ei};
  fprintf('\n--- combining animals | estimator: %s (%s) ---\n', est, EST_INFO.(est).short);
  for mi = 1:numel(metrics)
    m = metrics{mi};
    Rsum = 0; Gsum = 0;
    ZRsum = 0; ZGsum = 0; ZRnull = 0; ZGnull = 0;
    repl = true(size(G.(est).(m)(valid(1)).R)); repl_gain = repl;
    for ia = valid
        Rsum = Rsum + G.(est).(m)(ia).R;
        Gsum = Gsum + G.(est).(m)(ia).gain;
        repl = repl & G.(est).(m)(ia).sig;
        repl_gain = repl_gain & G.(est).(m)(ia).sig_gain;
        ZRsum  = ZRsum  + G.(est).(m)(ia).zR_obs;    ZGsum  = ZGsum  + G.(est).(m)(ia).zG_obs;
        ZRnull = ZRnull + G.(est).(m)(ia).zR_null;   ZGnull = ZGnull + G.(est).(m)(ia).zG_null;
    end
    nA = numel(valid);
    C.(est).(m).R = Rsum/nA; C.(est).(m).repl = repl;
    C.(est).(m).gain = Gsum/nA; C.(est).(m).repl_gain = repl_gain;
    C.(est).(m).fHz = G.(est).(m)(valid(1)).fHz; C.(est).(m).speeds = Vsets{mi};

    % pooled standardised grids + their max-stat thresholds.
    % max() omits NaN per column, so skipped frequencies are harmless; but if a
    % whole animal's grid is NaN (every frequency skipped for that mode) every
    % column is NaN, quantile returns NaN, and the panel/threshold are
    % meaningless. Guard that explicitly instead of letting NaN through.
    C.(est).(m).ZR = ZRsum/nA;  C.(est).(m).ZG = ZGsum/nA;
    maxR = max(reshape(ZRnull/nA, [], nPerm), [], 1).';   % nPerm × 1
    maxG = max(reshape(ZGnull/nA, [], nPerm), [], 1).';
    tR = quantile(maxR, 1-alpha);  tG = quantile(maxG, 1-alpha);
    if ~isfinite(tR) || ~isfinite(tG)
        warning(['%s/%s: pooled null is entirely NaN — at least one animal produced no finite ' ...
                 'cells for this mode (every frequency skipped? too few usable channels?). ' ...
                 'Pooled test disabled for this mode; replication is unaffected.'], est, m);
        if ~isfinite(tR), tR = Inf; end     % Inf => nothing significant, never NaN
        if ~isfinite(tG), tG = Inf; end
    end
    C.(est).(m).thr_pool_R = tR;                    C.(est).(m).thr_pool_G = tG;
    C.(est).(m).sig_pool_R = C.(est).(m).ZR >= tR;  C.(est).(m).sig_pool_G = C.(est).(m).ZG >= tG;
    C.(est).(m).pooled_ok  = isfinite(tR) && isfinite(tG);

    % zR and zG are algebraically identical (see the header note on why the gain
    % figures carry no pooled column). Assert it rather than trust it silently.
    dz = max(abs(C.(est).(m).ZR(:) - C.(est).(m).ZG(:)), [], 'omitnan');
    if isfinite(dz) && dz > 1e-8
        warning('%s/%s: pooled zR and zG differ by %.3g — the R0-invariance assumption is violated.', est, m, dz);
    end

    fprintf('%-16s: BOTH-animal (replication) R=%d gain=%d | POOLED z R=%d (thr %.2f) gain=%d (thr %.2f)\n', ...
        m, sum(repl(:)), sum(repl_gain(:)), sum(C.(est).(m).sig_pool_R(:)), C.(est).(m).thr_pool_R, ...
        sum(C.(est).(m).sig_pool_G(:)), C.(est).(m).thr_pool_G);

    % For any pooled-significant gain cell, show how lopsided it is: a pooled
    % hit carried by one animal is NOT a combined result.
    if any(C.(est).(m).sig_pool_G(:))
        [~, ix] = max(C.(est).(m).ZG(:) .* double(C.(est).(m).sig_pool_G(:)));
        [pf, pv] = ind2sub(size(C.(est).(m).ZG), ix);
        za = arrayfun(@(ia) G.(est).(m)(ia).zG_obs(pf,pv), valid);
        fprintf('%-16s:   peak pooled gain cell f=%.2f Hz v=%.1f: per-animal z = [%s] (%s)\n', ...
            m, C.(est).(m).fHz(pf), C.(est).(m).speeds(pv), ...
            strjoin(compose('%.2f', za), ', '), strjoin({G.(est).(m)(valid).animal}, ', '));
    end
  end
end

%% ─── Figure: freq × speed alignment heatmaps ─────────────────────────
% MATLAB's default TeX interpreter renders '_' as a subscript, which mangles
% the mode names ('visual_coherent' -> "visual" + subscript c) and rf_note
% ('Valid_Gaussian'). Escape underscores rather than switching the interpreter
% off, so the intentional TeX (\DeltaR, R_0, 2\pi) still renders. Titles are
% cell arrays = one line per element, which stops them overlapping.
esc = @(s) strrep(s, '_', '\_');

% Map a speed in physical units onto the index-based y axis (shared by both
% figures, so the best-speed line can be drawn on the R grid as well as on the
% gain grid).
speeds_to_idx = @(sp,sv) interp1(sv, 1:numel(sv), sp, 'linear', NaN);

%% ─── FIGURE 1: the DE-ROTATION R grids, both estimators ──────────────
% Rows = 3 modes x 2 estimators (estimator block named in every panel title,
% because R means a DIFFERENT THING in the two blocks and the values are NOT
% comparable between them). Cols = animals + mean/replication + pooled z.
% Colour limits are set PER ESTIMATOR BLOCK from that block's own data — a
% fixed [0 1] scale renders the phase-coherence block (R ~ 0.08) as flat blue.
nEst = numel(ESTIMATORS);
ncol = numel(valid)+2;                       % animals + mean + pooled z
nrow = nEst*numel(metrics);
f1 = figure('Visible','off','Position',[40 40 380*ncol 320*nrow]);
for ei = 1:nEst
  est = ESTIMATORS{ei};
  % shared colour range across all modes/animals of THIS estimator
  rmax_e = 0;
  for mi = 1:numel(metrics)
      rmax_e = max(rmax_e, max(C.(est).(metrics{mi}).R(:),[],'omitnan'));
      for ia = valid, rmax_e = max(rmax_e, max(G.(est).(metrics{mi})(ia).R(:),[],'omitnan')); end
  end
  if ~isfinite(rmax_e) || rmax_e <= 0, rmax_e = 1; end
  for mi = 1:numel(metrics)
    m   = metrics{mi};
    row = (ei-1)*numel(metrics) + mi;
    for k = 1:numel(valid)
        ia = valid(k);
        ax = subplot(nrow, ncol, (row-1)*ncol + k); hold(ax,'on');
        imagesc(ax, G.(est).(m)(ia).fHz, 1:numel(G.(est).(m)(ia).speeds), G.(est).(m)(ia).R.');
        set(ax,'YDir','normal'); axis(ax,'tight');
        contour(ax, G.(est).(m)(ia).fHz, 1:numel(G.(est).(m)(ia).speeds), double(G.(est).(m)(ia).sig.'), [0.5 0.5], 'w','LineWidth',1.2);
        yt = round(linspace(1,numel(G.(est).(m)(ia).speeds),5));
        set(ax,'YTick',yt,'YTickLabel',compose('%.0f',G.(est).(m)(ia).speeds(yt)));
        caxis(ax,[0 rmax_e]); colorbar(ax);
        xlabel(ax,'Frequency (Hz)'); ylabel(ax,sprintf('speed (%s)',Vunit{mi}));
        title(ax, {sprintf('%s — %s', G.(est).(m)(ia).animal, esc(m)), ...
                   sprintf('R  [%s]', EST_INFO.(est).short)}, 'FontSize',9);
    end
    % mean R, contour = REPLICATION (significant in both animals)
    ax = subplot(nrow, ncol, (row-1)*ncol + ncol-1); hold(ax,'on');
    imagesc(ax, C.(est).(m).fHz, 1:numel(C.(est).(m).speeds), C.(est).(m).R.');
    set(ax,'YDir','normal'); axis(ax,'tight');
    contour(ax, C.(est).(m).fHz, 1:numel(C.(est).(m).speeds), double(C.(est).(m).repl.'), [0.5 0.5], 'w','LineWidth',1.4);
    yt = round(linspace(1,numel(C.(est).(m).speeds),5));
    set(ax,'YTick',yt,'YTickLabel',compose('%.0f',C.(est).(m).speeds(yt)));
    caxis(ax,[0 rmax_e]); colorbar(ax);
    xlabel(ax,'Frequency (Hz)'); ylabel(ax,sprintf('speed (%s)',Vunit{mi}));
    title(ax, {sprintf('%s  [%s]', esc(m), EST_INFO.(est).short), 'mean R (replication outlined)'}, 'FontSize',9);

    % POOLED standardised z, contour = its own max-stat threshold. Units are z,
    % NOT R — auto-scaled, so do not compare its colours with the panels left of
    % it. Drawn only here: the pooled gain z is algebraically identical (header).
    ax = subplot(nrow, ncol, (row-1)*ncol + ncol); hold(ax,'on');
    imagesc(ax, C.(est).(m).fHz, 1:numel(C.(est).(m).speeds), C.(est).(m).ZR.');
    set(ax,'YDir','normal'); axis(ax,'tight');
    contour(ax, C.(est).(m).fHz, 1:numel(C.(est).(m).speeds), double(C.(est).(m).sig_pool_R.'), [0.5 0.5], 'w','LineWidth',1.4);
    set(ax,'YTick',yt,'YTickLabel',compose('%.0f',C.(est).(m).speeds(yt)));
    if ~any(isfinite(C.(est).(m).ZR(:))), caxis(ax,[0 1]); end   % all-NaN: sane range
    colorbar(ax);
    xlabel(ax,'Frequency (Hz)'); ylabel(ax,sprintf('speed (%s)',Vunit{mi}));
    if C.(est).(m).pooled_ok, ttl_p = sprintf('POOLED z (thr %.2f)', C.(est).(m).thr_pool_R);
    else,                     ttl_p = 'POOLED z — n/a (see warning)'; end
    title(ax, {sprintf('%s  [%s]', esc(m), EST_INFO.(est).short), ttl_p}, 'FontSize',9);
  end
end
% Plain ASCII in the sgtitle: MATLAB's TeX symbols (\Delta, subscripts) export to
% PDF with broken advance widths, leaving visible gaps mid-word.
sgtitle({['DE-ROTATION R: coherence/alignment after de-rotating locations by k*d, k = 2*pi*f/v   (' esc(note) ')'], ...
         'ROWS 1-3 = PHASE ALIGNMENT estimator (R = resultant of the per-location preferred-phase vectors — NOT a phase coherence)', ...
         'ROWS 4-6 = PHASE COHERENCE estimator (R = coherence over all trials, each de-rotated by its own location — the Phase\_coherence/ measure at k=0)', ...
         'R IS NOT COMPARABLE BETWEEN THE TWO BLOCKS: colour limits are set per block, thresholds are per block.', ...
         'white contour = significant (max-stat corrected over the whole grid)   |   col 3 = mean, contour = replication   |   col 4 = pooled standardised z (z units, auto-scaled)'}, ...
         'FontSize',9);
set(f1,'PaperPositionMode','auto'); pos=get(f1,'Position'); set(f1,'PaperUnits','points','PaperSize',pos(3:4));
saveas(f1, fullfile(out_dir, ['stimulus_loc_grids' tag '.pdf']));

%% ─── FIGURES 2..n: de-rotation GAIN, ONE FIGURE PER ESTIMATOR ────────
% dR(f,v) = R(f,v) - R0(f): the speeds at which de-rotating beats not rotating.
% A wave = significant positive gain (white contour) along a near-CONSTANT best
% speed (black line) across frequencies; a rising black line means constant
% wavenumber, i.e. a fixed phase offset rather than a propagating wave.
% NO pooled column here — the pooled gain z equals the pooled R z exactly, and
% is drawn once on the grids figure (see the header note).
ncol_g = numel(valid)+1;                     % animals + mean/replication
for ei = 1:nEst
  est = ESTIMATORS{ei};
  fg = figure('Visible','off','Position',[40 40 380*ncol_g 320*numel(metrics)]);
  for mi = 1:numel(metrics)
    m = metrics{mi};
    gmax_m = max(C.(est).(m).gain(:),[],'omitnan');
    for ia = valid, gmax_m = max(gmax_m, max(G.(est).(m)(ia).gain(:),[],'omitnan')); end
    if ~isfinite(gmax_m) || gmax_m <= 0, gmax_m = eps; end
    for k = 1:numel(valid)
        ia = valid(k);
        ax = subplot(numel(metrics), ncol_g, (mi-1)*ncol_g + k); hold(ax,'on');
        imagesc(ax, G.(est).(m)(ia).fHz, 1:numel(G.(est).(m)(ia).speeds), G.(est).(m)(ia).gain.');
        set(ax,'YDir','normal'); axis(ax,'tight');
        contour(ax, G.(est).(m)(ia).fHz, 1:numel(G.(est).(m)(ia).speeds), double(G.(est).(m)(ia).sig_gain.'), [0.5 0.5], 'w','LineWidth',1.2);
        plot(ax, G.(est).(m)(ia).fHz, speeds_to_idx(G.(est).(m)(ia).vbest, G.(est).(m)(ia).speeds), 'k-','LineWidth',1.3);
        yt = round(linspace(1,numel(G.(est).(m)(ia).speeds),5));
        set(ax,'YTick',yt,'YTickLabel',compose('%.0f',G.(est).(m)(ia).speeds(yt)));
        caxis(ax,[0 gmax_m]); colorbar(ax);
        xlabel(ax,'Frequency (Hz)'); ylabel(ax,sprintf('speed (%s)',Vunit{mi}));
        title(ax, {sprintf('%s — %s', G.(est).(m)(ia).animal, esc(m)), ...
                   sprintf('gain dR (R0 %.3f-%.3f)', min(G.(est).(m)(ia).R0,[],'omitnan'), max(G.(est).(m)(ia).R0,[],'omitnan'))}, 'FontSize',9);
    end
    % mean gain, contour = REPLICATION
    ax = subplot(numel(metrics), ncol_g, (mi-1)*ncol_g + ncol_g); hold(ax,'on');
    imagesc(ax, C.(est).(m).fHz, 1:numel(C.(est).(m).speeds), C.(est).(m).gain.');
    set(ax,'YDir','normal'); axis(ax,'tight');
    contour(ax, C.(est).(m).fHz, 1:numel(C.(est).(m).speeds), double(C.(est).(m).repl_gain.'), [0.5 0.5], 'w','LineWidth',1.4);
    yt = round(linspace(1,numel(C.(est).(m).speeds),5));
    set(ax,'YTick',yt,'YTickLabel',compose('%.0f',C.(est).(m).speeds(yt)));
    caxis(ax,[0 gmax_m]); colorbar(ax);
    xlabel(ax,'Frequency (Hz)'); ylabel(ax,sprintf('speed (%s)',Vunit{mi}));
    title(ax, {esc(m), 'mean gain (replication outlined)'}, 'FontSize',9);
  end
  sgtitle({sprintf('DE-ROTATION GAIN dR = R(f,v) - R0(f)   |   %s   (%s)', upper(EST_INFO.(est).short), esc(note)), ...
           esc(EST_INFO.(est).long), ...
           'R0 = the same statistic at k=0 (nothing de-rotated). dR asks: does assuming a wave buy any alignment/coherence over assuming none?', ...
           'white contour = significant   |   black line = best-fit speed per frequency   |   FLAT line = real wave, RISING = constant wavenumber (not a wave)', ...
           'no pooled column: the pooled gain z is identical to the pooled R z, drawn once on the grids figure'}, ...
           'FontSize',9);
  set(fg,'PaperPositionMode','auto'); pos=get(fg,'Position'); set(fg,'PaperUnits','points','PaperSize',pos(3:4));
  saveas(fg, fullfile(out_dir, ['stimulus_loc_gain' EST_INFO.(est).tag tag '.pdf']));
end

% Drop the per-permutation null grids before saving — they are only needed to
% build the pooled threshold above and would add ~100 MB to the .mat.
% The standardised OBSERVED grids (zR_obs/zG_obs) are kept: they are what tells
% you which animal drives any pooled-significant cell.
for ei = 1:numel(ESTIMATORS)
    for mi = 1:numel(metrics)
        for ia = valid
            G.(ESTIMATORS{ei}).(metrics{mi})(ia).zR_null = [];
            G.(ESTIMATORS{ei}).(metrics{mi})(ia).zG_null = [];
        end
    end
end

% BOTH estimators live in one file: results.G.(estimator).(mode)(animal) and
% results.C.(estimator).(mode). Only the RF toggle is in the filename.
results = struct('G',G,'C',C,'animals',{animals},'dv',dv, ...
    'FREQ_RANGE',FREQ_RANGE,'V_VISUAL',V_VISUAL, ...
    'PIX_PER_DEG',PIX_PER_DEG,'RF_DATE',RF_DATE,'SCREEN_XY',SCREEN_XY, ...
    'RF_VALID_ONLY',RF_VALID_ONLY,'ESTIMATORS',{ESTIMATORS},'EST_INFO',EST_INFO, ...
    'metrics',{metrics},'cmode',{cmode},'Vunit',{Vunit},'nPerm',nPerm,'alpha',alpha);
save(fullfile(res_dir, ['stimulus_loc_wave' tag '.mat']),'results','-v7.3');
fprintf('\nSaved under %s :\n  stimulus_loc_grids%s.pdf   (de-rotation R, both estimators)\n', out_dir, tag);
for ei = 1:numel(ESTIMATORS)
    fprintf('  stimulus_loc_gain%s%s.pdf   (gain dR, %s)\n', ...
        EST_INFO.(ESTIMATORS{ei}).tag, tag, EST_INFO.(ESTIMATORS{ei}).short);
end

%% =====================================================================
%% Helpers
%%
%%  align_grid* = ESTIMATOR 'phase'      PHASE ALIGNMENT
%%      Input : pref_phase + coh_mag, i.e. nPos PRE-AVERAGED per-location
%%              vectors (the trials are already gone).
%%      R     : weighted RESULTANT LENGTH over those locations.
%%              NOT a phase coherence — a circular-concentration measure.
%%              R0 runs high (~0.8-0.93) because a handful of angles are
%%              easily concentrated.
%%      Uses  : S(p,f)/n_p implicitly, weighted by |c_p| ~ 1/sqrt(n_p)
%%              -> fewer-trial locations get LARGER weight.
%%
%%  coh_grid*  = ESTIMATOR 'coherence'   PHASE COHERENCE
%%      Input : S(p,f) = raw complex sums over the TRIALS at each location.
%%      R     : a genuine PHASE COHERENCE over all trials, each de-rotated by
%%              its own location's predicted lag. At k=0 it equals
%%              Phase_coherence/functions/phase_coherence.m with locations
%%              pooled (normaliser /W instead of /nTrials, so R is in [0,1]).
%%      Uses  : S(p,f) raw -> every trial counts exactly once.
%%
%%  Everything else is shared: both return the same six outputs, both define
%%  R0 as their own statistic at k=0, both feed the identical gain / null /
%%  threshold / pooling code. R values are NOT comparable between the two.
%% =====================================================================
function [Robs, R0, Rnull_max, Gnull_max, Rnull, Gnull] = align_grid(pref, coh, coh_sig, f_use, fHz, d, Vs, MIN_LOC, nPerm)
% ESTIMATOR 'phase' (PHASE ALIGNMENT), mode 'visual'. INCOHERENT across channels.
% R(f,v)  = mean over coherence-sig channels of the DE-ROTATED resultant
%           across stimulus locations (de-rotate each location by 2pi f d/v).
% R0(f)   = the same resultant with NO rotation (k=0, i.e. v->inf): the
%           "actual" location coherence with no wave assumed. This is the
%           baseline the de-rotation gain  dR(f,v)=R(f,v)-R0(f)  is measured
%           against — dR isolates the wave-specific phase ramp and is
%           invariant to locations merely sharing a preferred phase.
% Returns observed grid + no-rotation baseline + per-perm grid maxima of
%   R  (absolute null)  and  of  dR=R-R0  (gain null). R0 does not use d, so
% it is unchanged by the location<->distance shuffle -> a clean gain null.
nF = numel(f_use); nV = numel(Vs); nPos = numel(d);
Robs      = nan(nF, nV);
R0        = nan(nF, 1);
Rnull_max = nan(nPerm, 1);
Gnull_max = nan(nPerm, 1);
Rnull     = nan(nF, nV, nPerm);   % full null grid, kept for cross-animal pooling
Gnull     = nan(nF, nV, nPerm);

% Pre-shuffle the location order once per permutation (shared across all
% channels and all (f,v) cells -> synchronised null + grid-wide max-stat).
perms = zeros(nPerm, nPos);
for b = 1:nPerm, perms(b,:) = randperm(nPos); end
null_max_running  = -inf(nPerm,1);
gnull_max_running = -inf(nPerm,1);

for fi = 1:nF
    f = f_use(fi); fh = fHz(fi);
    PHI = squeeze(pref(:, f, :));     % nCh × nPos
    W   = squeeze(coh(:, f, :));      % nCh × nPos
    M   = isfinite(PHI) & isfinite(W) & (W > 0);
    W(~M) = 0; PHI(~M) = 0;
    useC = coh_sig(:, f) & (sum(M,2) >= MIN_LOC);
    if ~any(useC)
        continue
    end
    A   = W .* exp(1i*PHI);           % nCh × nPos (weighted phase vectors)
    sw  = sum(W, 2);                  % nCh × 1
    A   = A(useC,:); sw = sw(useC);

    % no-rotation baseline (k=0): "actual" location coherence, no wave assumed
    R0(fi) = mean( abs(sum(A,2)) ./ max(sw,eps) );

    for vi = 1:nV
        k = 2*pi*fh / Vs(vi);                       % rad per unit distance
        % observed: de-rotate by θ_p = k·d_p
        Robs(fi,vi) = mean( abs(A * exp(-1i*k*d(:))) ./ max(sw,eps) );
        % (per channel: row of A dotted with the de-rotation phasors)
    end

    % null: same channels/weights, shuffled location<->distance mapping.
    % R0 is unaffected by the shuffle (k=0 uses no d), so subtracting it
    % gives a clean null on the de-rotation GAIN dR = R - R0.
    for b = 1:nPerm
        dp = d(perms(b,:));
        rb = -inf;
        for vi = 1:nV
            k = 2*pi*fh / Vs(vi);
            rcell = mean( abs(A * exp(-1i*k*dp(:))) ./ max(sw,eps) );
            Rnull(fi,vi,b) = rcell;
            Gnull(fi,vi,b) = rcell - R0(fi);
            if rcell > rb, rb = rcell; end
        end
        if rb > null_max_running(b),          null_max_running(b)  = rb;         end
        if rb - R0(fi) > gnull_max_running(b), gnull_max_running(b) = rb - R0(fi); end
    end
end
[Rnull_max, Gnull_max] = null_guard(null_max_running, gnull_max_running, nPerm);
end

function [Robs, R0, Rnull_max, Gnull_max, Rnull, Gnull] = align_grid_coherent(pref, coh, coh_sig, f_use, fHz, d, Dc, Vs, MIN_LOC, nPerm)
% ESTIMATOR 'phase' (PHASE ALIGNMENT), mode 'visual_coherent'.
% COHERENT combination across channels: test ONE cortical wave in which BOTH
% the stimulus-location eccentricity d_p AND the electrode RF eccentricity
% Dc_c contribute the same wave phase k*distance (k = 2*pi*f/v):
%   R(f,v) = | Σ_{c,p} w_cp e^{i(φ_cp − k d_p − k Dc_c)} | / Σ w_cp
%          = | (e^{-ik·Dc})ᵀ · A · e^{-ik·d} | / Σ w_cp .
% Unlike align_grid this does NOT take |·| per channel, so the per-electrode
% delay is now VISIBLE (but unmodelled per-channel offsets can only lower R —
% never fabricate it — so a positive band is trustworthy, a null ambiguous).
% Channels without a finite Dc are dropped. R0 = the k=0 coherent resultant
% (no wave). Null shuffles BOTH d→location and Dc→channel with a SINGLE
% shared shuffle per permutation (synchronised across the whole grid).
nF = numel(f_use); nV = numel(Vs); nPos = numel(d); nCh = numel(Dc);
Robs      = nan(nF, nV);
R0        = nan(nF, 1);
Rnull_max = nan(nPerm, 1);
Gnull_max = nan(nPerm, 1);
Rnull     = nan(nF, nV, nPerm);   % full null grid, kept for cross-animal pooling
Gnull     = nan(nF, nV, nPerm);
Dc = Dc(:);

% one shared shuffle per permutation, over locations AND over channels,
% reused across every (f,v) cell -> synchronised null + grid-wide max-stat.
% NaN-SAFE: channels with no finite Dc (RF_VALID_ONLY drops the Extrapolated
% centres) are EXCLUDED from the observed set by useC, but an unrestricted
% randperm would shuffle one of them INTO a used slot, making that permutation
% NaN. NaN > rb is false, so the running max would stay -inf and every cell
% would come out "significant" against thr = -inf. Permute only among the
% finite-geometry channels and leave the NaN slots where they are.
okC = isfinite(Dc(:)).'; iOK = find(okC); iNo = find(~okC);
permP = zeros(nPerm, nPos); permC = zeros(nPerm, nCh);
for b = 1:nPerm
    permP(b,:)   = randperm(nPos);
    permC(b,iOK) = iOK(randperm(numel(iOK)));
    permC(b,iNo) = iNo;
end
null_max_running  = -inf(nPerm,1);
gnull_max_running = -inf(nPerm,1);

for fi = 1:nF
    f = f_use(fi); fh = fHz(fi);
    PHI = squeeze(pref(:, f, :));     % nCh × nPos
    W   = squeeze(coh(:, f, :));      % nCh × nPos
    M   = isfinite(PHI) & isfinite(W) & (W > 0);
    W(~M) = 0; PHI(~M) = 0;
    useC = coh_sig(:, f) & (sum(M,2) >= MIN_LOC) & isfinite(Dc);
    if sum(useC) < 2, continue; end
    A    = W .* exp(1i*PHI);          % nCh × nPos (weighted phase vectors)
    A    = A(useC,:);                 % nUse × nPos (keeps per-channel phase)
    swtot = sum(W(useC,:), 'all');
    Dcu  = Dc(useC);

    % no-rotation baseline (k=0): fully coherent resultant, no wave assumed
    R0(fi) = abs(sum(A(:))) / max(swtot, eps);

    for vi = 1:nV
        k = 2*pi*fh / Vs(vi);                            % rad per unit distance
        % Z = row(1×nUse) · A(nUse×nPos) · col(nPos×1) = Σ_{c,p} w e^{i(φ−k(d+Dc))}
        Z = exp(-1i*k*Dcu(:)).' * A * exp(-1i*k*d(:));
        Robs(fi,vi) = abs(Z) / max(swtot, eps);
    end

    % null: de-correlate the wave by shuffling both assignments. Dc is
    % permuted over ALL channels then subset to the used set (shared across
    % frequencies); R0 (k=0) is invariant to the shuffle -> clean gain null.
    for b = 1:nPerm
        dp   = d(permP(b,:));
        Dcpb = Dc(permC(b,:));                           % nCh (col), then restrict
        Dcpb = Dcpb(useC);
        rb = -inf;
        for vi = 1:nV
            k = 2*pi*fh / Vs(vi);
            Z = exp(-1i*k*Dcpb(:)).' * A * exp(-1i*k*dp(:));
            rcell = abs(Z) / max(swtot, eps);
            Rnull(fi,vi,b) = rcell;
            Gnull(fi,vi,b) = rcell - R0(fi);
            if rcell > rb, rb = rcell; end
        end
        if rb > null_max_running(b),          null_max_running(b)  = rb;          end
        if rb - R0(fi) > gnull_max_running(b), gnull_max_running(b) = rb - R0(fi); end
    end
end
[Rnull_max, Gnull_max] = null_guard(null_max_running, gnull_max_running, nPerm);
end

function [Robs, R0, Rnull_max, Gnull_max, Rnull, Gnull] = align_grid_arrival(pref, coh, coh_sig, f_use, fHz, d, a, Vs, MIN_LOC, nPerm)
% ESTIMATOR 'phase' (PHASE ALIGNMENT), mode 'visual_arrival'.
% ARRIVAL-TIME model, INCOHERENT across channels.
%   R_c(f,v) = | Σ_p w_cp e^{i(φ_cp − k(d_p + a(c,p)))} | / Σ_p w_cp ,
%   R(f,v)   = mean over coherence-significant channels of R_c ,   k = 2πf/v.
%
% a [nCh x nPos] is the visual-field separation (deg) between electrode c's RF
% centre and stimulus location p. Unlike align_grid_coherent's D_c, a(c,p)
% VARIES across locations within a channel (it dips where the stimulus lands on
% that channel's RF), so it is NOT removed by the per-channel |.| — which is
% why this model can be fitted incoherently and therefore needs no
% common-reference assumption. Only the position-varying part of a is tested:
% |.|_c deletes the row mean of a exactly as it deletes any channel constant.
%
% Vs = v is shared by BOTH terms (the scan across locations and the spread from
% the driven patch), so no extra free parameter is introduced — the grid stays
% frequency x speed.
%
% Null: ONE synchronised shuffle per permutation, reused across every (f,v)
% cell, permuting the geometry two ways at once —
%   columns  location <-> (d_p, a(:,p))  : breaks "phase tracks the stimulus"
%   rows     channel  <-> a(c,:)         : breaks "the dip is at the RIGHT
%                                           electrode"
% R0 (k=0) uses neither d nor a, so it is shuffle-invariant -> clean gain null.
% (R0 here is numerically identical to align_grid's R0, as it must be.)
nF = numel(f_use); nV = numel(Vs); nPos = numel(d); nChAll = size(a,1);
Robs      = nan(nF, nV);
R0        = nan(nF, 1);
Rnull_max = nan(nPerm, 1);
Gnull_max = nan(nPerm, 1);
Rnull     = nan(nF, nV, nPerm);   % full null grid, kept for cross-animal pooling
Gnull     = nan(nF, nV, nPerm);

% NaN-SAFE: same issue as align_grid_coherent — channels whose arrival row
% a(c,:) is not finite (RF_VALID_ONLY) must not be shuffled into a used slot,
% or the whole permutation goes NaN and thr collapses to -inf.
okC = all(isfinite(a),2).'; iOK = find(okC); iNo = find(~okC);
permP = zeros(nPerm, nPos); permC = zeros(nPerm, nChAll);
for b = 1:nPerm
    permP(b,:)   = randperm(nPos);
    permC(b,iOK) = iOK(randperm(numel(iOK)));
    permC(b,iNo) = iNo;
end
null_max_running  = -inf(nPerm,1);
gnull_max_running = -inf(nPerm,1);

for fi = 1:nF
    f = f_use(fi); fh = fHz(fi);
    PHI = squeeze(pref(:, f, :));     % nCh × nPos
    W   = squeeze(coh(:, f, :));      % nCh × nPos
    M   = isfinite(PHI) & isfinite(W) & (W > 0);
    W(~M) = 0; PHI(~M) = 0;
    useC = coh_sig(:, f) & (sum(M,2) >= MIN_LOC) & all(isfinite(a), 2);
    if ~any(useC), continue; end
    A  = W .* exp(1i*PHI);            % nCh × nPos (weighted phase vectors)
    sw = sum(W, 2);
    A  = A(useC,:); sw = sw(useC);

    % no-rotation baseline (k=0): "actual" location coherence, no wave assumed
    R0(fi) = mean( abs(sum(A,2)) ./ max(sw,eps) );

    % Per speed, cache the de-rotation phasors on the FULL channel set. Row
    % permutation commutes with the elementwise exp, so the null can reuse
    % these by indexing instead of recomputing exp() 1000x per cell.
    Ea = cell(1,nV); Ed = cell(1,nV);
    for vi = 1:nV
        k = 2*pi*fh / Vs(vi);                 % rad per degree
        Ea{vi} = exp(-1i * k * a);            % nChAll × nPos
        Ed{vi} = exp(-1i * k * d(:).');       % 1 × nPos
        E = Ea{vi}(useC,:) .* Ed{vi};
        Robs(fi,vi) = mean( abs(sum(A .* E, 2)) ./ max(sw,eps) );
    end

    for b = 1:nPerm
        cc = permC(b,:); pp = permP(b,:);
        rb = -inf;
        for vi = 1:nV
            Ep    = Ea{vi}(cc,:);                        % rows reassigned
            Ep    = Ep(useC, pp) .* Ed{vi}(pp);          % columns relabelled
            rcell = mean( abs(sum(A .* Ep, 2)) ./ max(sw,eps) );
            Rnull(fi,vi,b) = rcell;
            Gnull(fi,vi,b) = rcell - R0(fi);
            if rcell > rb, rb = rcell; end
        end
        if rb > null_max_running(b),           null_max_running(b)  = rb;          end
        if rb - R0(fi) > gnull_max_running(b), gnull_max_running(b) = rb - R0(fi); end
    end
end
[Rnull_max, Gnull_max] = null_guard(null_max_running, gnull_max_running, nPerm);
end

% =====================================================================
% TRIAL-level estimators. Same three modes, same outputs, same nulls as the
% align_grid* family above — only R is defined differently:
%
%   'phase'     PHASE ALIGNMENT — R = resultant of the per-location PREFERRED
%               PHASES (coh_mag-weighted). Not a coherence.
%   'coherence' PHASE COHERENCE — R = |c|, c = (1/W)*SUM over TRIALS
%               y*exp(i(phi - k*dist)). A coherence in the Phase_coherence/ sense.
%
% Everything works off the per-location complex sums produced by
% functions/trial_position_sums_chan.m:
%       S(c,p,f) = SUM over trials at location p of  y_t*exp(i*phi_tf)
%       W(c)     = SUM over trials of |y_t|
% because the de-rotation depends on a trial only through its location:
%       c_c(f,v) = ( 1/W_c ) * SUM_p S(c,p,f)*exp(-i*k*d_p) .
% At k=0 this is exactly the phase coherence of the phase_coherence/ pipeline
% with all locations pooled, so R0 here is a coherence in the familiar sense.
% =====================================================================
function [Robs, R0, Rnull_max, Gnull_max, Rnull, Gnull] = ...
        coh_grid(Ssum, Sperm, Wch, coh_sig, f_use, fHz, d, Vs, nPerm)
% ESTIMATOR 'coherence' (PHASE COHERENCE), mode 'visual'. The trial-level twin
% of align_grid: same de-rotation, same modes, same null — but R here is a real
% phase coherence over trials, not a resultant of pre-averaged location vectors.
% R0 = |sum_p S(p,f)| / W is exactly the Phase_coherence/ number, locations
% pooled. INCOHERENT across channels (|.| per channel, then averaged).
nF = numel(f_use); nV = numel(Vs);
[Robs, R0, Rnull, Gnull] = deal(nan(nF,nV), nan(nF,1), nan(nF,nV,nPerm), nan(nF,nV,nPerm));
null_max_running = -inf(nPerm,1); gnull_max_running = -inf(nPerm,1);

for fi = 1:nF
    f = f_use(fi); fh = fHz(fi);
    useC = coh_sig(:,f) & isfinite(Wch) & (Wch > 0);
    if ~any(useC), continue; end
    Sf = squeeze(Ssum(:,:,f));                       % nCh × nPos
    Wu = Wch(useC);
    R0(fi) = mean( abs(sum(Sf(useC,:),2)) ./ Wu );   % k = 0 -> pooled coherence

    E = exp(-1i * (2*pi*fh./Vs(:).') .* d(:));       % nPos × nV
    Zo = Sf(useC,:) * E;                             % nUse × nV
    Robs(fi,:) = mean( abs(Zo) ./ Wu, 1 );

    % null: all permutations in one matrix product
    Sp = double(reshape(permute(Sperm(useC,:,f,:), [1 4 2 3]), [], size(Sperm,2)));
    Zn = reshape(abs(Sp * E), sum(useC), nPerm, nV);  % (nUse, nPerm, nV)
    Rn = squeeze(mean(Zn ./ Wu, 1));                  % nPerm × nV
    Rnull(fi,:,:) = Rn.';  Gnull(fi,:,:) = (Rn - R0(fi)).';
    rb = max(Rn, [], 2);                              % nPerm × 1 (max over speed)
    null_max_running  = max(null_max_running,  rb);
    gnull_max_running = max(gnull_max_running, rb - R0(fi));
end
[Rnull_max, Gnull_max] = null_guard(null_max_running, gnull_max_running, nPerm);
end

function [Robs, R0, Rnull_max, Gnull_max, Rnull, Gnull] = ...
        coh_grid_coherent(Ssum, Sperm, Wch, coh_sig, f_use, fHz, d, Dc, Vs, nPerm)
% ESTIMATOR 'coherence' (PHASE COHERENCE), mode 'visual_coherent'.
% COHERENT across channels — the trial-level twin of align_grid_coherent.
% Channels are summed complex (per-channel phase kept), each de-rotated by its
% own k*Dc, so the per-electrode delay is testable. Same shuffle logic: the
% location relabelling is already synchronised across channels inside Sperm,
% and Dc is additionally permuted over channels here.
nF = numel(f_use); nV = numel(Vs); nCh = numel(Dc); Dc = Dc(:);
[Robs, R0, Rnull, Gnull] = deal(nan(nF,nV), nan(nF,1), nan(nF,nV,nPerm), nan(nF,nV,nPerm));
null_max_running = -inf(nPerm,1); gnull_max_running = -inf(nPerm,1);

okC = isfinite(Dc).'; iOK = find(okC); iNo = find(~okC);
permC = zeros(nPerm, nCh);
for b = 1:nPerm
    permC(b,iOK) = iOK(randperm(numel(iOK)));   % NaN slots stay put (see align_grid_coherent)
    permC(b,iNo) = iNo;
end

for fi = 1:nF
    f = f_use(fi); fh = fHz(fi);
    useC = coh_sig(:,f) & isfinite(Wch) & (Wch > 0) & isfinite(Dc);
    if sum(useC) < 2, continue; end
    Sf = squeeze(Ssum(:,:,f)); Wtot = sum(Wch(useC));
    R0(fi) = abs(sum(sum(Sf(useC,:)))) / Wtot;

    for vi = 1:nV
        k  = 2*pi*fh / Vs(vi);
        ed = exp(-1i*k*d(:));  u = exp(-1i*k*Dc(useC));
        Robs(fi,vi) = abs( u.' * (Sf(useC,:) * ed) ) / Wtot;
        Zc = squeeze(sum(Sperm(useC,:,f,:) .* ed.', 2));      % nUse × nPerm
        Dp = Dc(permC(:,useC).');                              % nUse × nPerm
        Rn = abs(sum(double(Zc) .* exp(-1i*k*Dp), 1)).' / Wtot; % nPerm × 1
        Rnull(fi,vi,:) = Rn;  Gnull(fi,vi,:) = Rn - R0(fi);
    end
    rb = max(squeeze(Rnull(fi,:,:)), [], 1).';
    null_max_running  = max(null_max_running,  rb);
    gnull_max_running = max(gnull_max_running, rb - R0(fi));
end
[Rnull_max, Gnull_max] = null_guard(null_max_running, gnull_max_running, nPerm);
end

function [Robs, R0, Rnull_max, Gnull_max, Rnull, Gnull] = ...
        coh_grid_arrival(Ssum, Sperm, Wch, coh_sig, f_use, fHz, d, a, Vs, nPerm)
% ESTIMATOR 'coherence' (PHASE COHERENCE), mode 'visual_arrival'.
% ARRIVAL model, INCOHERENT — the trial-level twin of align_grid_arrival.
% De-rotates by k*(d_p + a(c,p)); a varies across locations within a channel,
% so it survives the per-channel |.| exactly as in the phase estimator.
nF = numel(f_use); nV = numel(Vs); nChAll = size(a,1);
[Robs, R0, Rnull, Gnull] = deal(nan(nF,nV), nan(nF,1), nan(nF,nV,nPerm), nan(nF,nV,nPerm));
null_max_running = -inf(nPerm,1); gnull_max_running = -inf(nPerm,1);

okC = all(isfinite(a),2).'; iOK = find(okC); iNo = find(~okC);
permC = zeros(nPerm, nChAll);
for b = 1:nPerm
    permC(b,iOK) = iOK(randperm(numel(iOK)));
    permC(b,iNo) = iNo;
end

for fi = 1:nF
    f = f_use(fi); fh = fHz(fi);
    useC = coh_sig(:,f) & isfinite(Wch) & (Wch > 0) & all(isfinite(a),2);
    if ~any(useC), continue; end
    Sf = squeeze(Ssum(:,:,f)); Wu = Wch(useC);
    R0(fi) = mean( abs(sum(Sf(useC,:),2)) ./ Wu );

    for vi = 1:nV
        k  = 2*pi*fh / Vs(vi);
        Ea = exp(-1i*k*a);  ed = exp(-1i*k*d(:)).';        % nChAll × nPos, 1 × nPos
        Robs(fi,vi) = mean( abs(sum(Sf(useC,:) .* Ea(useC,:) .* ed, 2)) ./ Wu );
        % null: rows of a reassigned across channels, locations already shuffled
        Ep = Ea(permC(:,useC).', :);                        % (nUse*nPerm) × nPos
        Sp = double(reshape(permute(Sperm(useC,:,f,:), [1 4 2 3]), [], size(Sperm,2)));
        Zn = abs(sum(Sp .* Ep .* ed, 2));
        Rn = mean(reshape(Zn, sum(useC), nPerm) ./ Wu, 1).';
        Rnull(fi,vi,:) = Rn;  Gnull(fi,vi,:) = Rn - R0(fi);
    end
    rb = max(squeeze(Rnull(fi,:,:)), [], 1).';
    null_max_running  = max(null_max_running,  rb);
    gnull_max_running = max(gnull_max_running, rb - R0(fi));
end
[Rnull_max, Gnull_max] = null_guard(null_max_running, gnull_max_running, nPerm);
end

function [Rnull_max, Gnull_max] = null_guard(rm, gm, nPerm)
% Shared by BOTH helper families (align_grid* and coh_grid*): never return
% -inf/NaN silently. A -inf threshold marks every cell significant.
if ~all(isfinite(rm))
    warning(['%d/%d permutations produced no finite statistic — the max-stat ' ...
             'threshold is INVALID. Usual cause: no usable channels at any frequency.'], ...
             sum(~isfinite(rm)), nPerm);
end
Rnull_max = rm; Gnull_max = gm;
end

function [Ssum, Sperm, Wch] = get_trial_sums(base, animalName, dv, nCh, nPerm, force)
% Build (via one SLURM job per channel) and aggregate the per-location complex
% sums S(c,p,f) the 'coherence' estimator runs on. Re-uses saved sums when present,
% so re-running the analysis does not re-submit the cluster jobs.
data_dir = fullfile(base, ['results_' animalName], 'multi_lin_reg', 'cp10_till_100');
sum_dir  = fullfile(base, ['results_' animalName], 'scanning', 'trial_position_sums', 'cp10_till_100', dv);
if ~exist(sum_dir,'dir'), mkdir(sum_dir); end

need = force;
if ~need
    for ch = 1:nCh
        if ~isfile(fullfile(sum_dir, num2str(ch), 'trial_position_sums.mat')), need = true; break; end
    end
end
if need
    fprintf('  submitting %d SLURM jobs for trial-level position sums...\n', nCh);
    cfg = cell(1, nCh);
    for ch = 1:nCh
        cfg{ch}.ichan = ch;  cfg{ch}.nPerm = nPerm;  cfg{ch}.dv = dv;
        cfg{ch}.infile = data_dir;  cfg{ch}.outfile = sum_dir;
        cfg{ch}.perm_seed_base = 2025;      % SHARED across channels (synchronised null)
    end
    slurmfun(@trial_position_sums_chan, cfg, ...
        'partition','8GB', 'stopOnError',false, 'useUserPath',true);
else
    fprintf('  re-using cached trial-level position sums in %s\n', sum_dir);
end

% Peek at the first available channel for the dimensions, then preallocate:
% Sperm is the big one (nCh x nPos x nFreq x nPerm, single complex).
first = '';
for ch = 1:nCh
    fch = fullfile(sum_dir, num2str(ch), 'trial_position_sums.mat');
    if isfile(fch), first = fch; break; end
end
if isempty(first), error('no trial-level position sums found under %s', sum_dir); end
D0 = load(first, 'S_obs');
[nPos, nFreq] = size(D0.S_obs);
Ssum  = complex(nan(nCh, nPos, nFreq));
Sperm = complex(nan(nCh, nPos, nFreq, nPerm, 'single'));
Wch   = nan(nCh,1);
for ch = 1:nCh
    fch = fullfile(sum_dir, num2str(ch), 'trial_position_sums.mat');
    if ~isfile(fch), warning('missing trial sums for channel %d', ch); continue; end
    D = load(fch);
    Ssum(ch,:,:)    = D.S_obs;
    Sperm(ch,:,:,:) = D.S_perm;
    Wch(ch)         = D.W;
end
fprintf('  trial sums: %d/%d channels loaded | Sperm %.0f MB\n', ...
    sum(isfinite(Wch)), nCh, numel(Sperm)*8/1e6);
end

function [xy_deg, valid_gauss] = elec_rf_deg(base, animal, nCh, rf_date, ppd, screen_xy)
% Per-channel electrode RF centre in FIXATION-CENTRED DEGREES:
%   xy_deg [nCh x 2] = (x, y) with (0,0) = fovea — the SAME frame the stimulus
%   positions are converted to, so the two may be differenced safely.
%
% Source: the Gaussian-fit centre written by RF_Mapping/mapping_lfp.m into
%   Plots/RF_Mapping/<animal>/loc_RF_map/gaussian_overlap/
%       <animal>_<date>_..._channel_target_summary.txt
% RF_Center_X/Y there are SCREEN pixels — chan_loc_mua.m adds centerX/centerY
% when writing them — so the screen centre is SUBTRACTED here before /ppd.
% The stimulus values (trialinfo col 16/17) are already fixation-centred and
% therefore need no such offset; mixing the two frames without this
% subtraction is off by (840, 525) px. Rows with no fit hold '-' -> NaN.
% Channel index is 1:1 with the phase-progression array (64 rows).
%
% valid_gauss(ch) is TRUE only where Status == 'Valid_Gaussian', i.e. a real
% 2-D Gaussian fit to that channel's RF map. Status 'Extrapolated' means the
% fit FAILED and the centre was filled in from the median of the channel's
% 8-connected array neighbours, or from a plane fit across the array
% (RF_Mapping/mapping_lfp.m:513) — inferred under a retinotopic-smoothness
% assumption, not measured. 'No_Data' rows are NaN and never valid.
xy_deg = nan(nCh,2); valid_gauss = false(nCh,1);
odir = fullfile(base,'Plots','RF_Mapping',animal,'loc_RF_map','gaussian_overlap');
rf_file = fullfile(odir, sprintf('%s_%s_rfmapping_bar_1_channel_target_summary.txt', animal, rf_date));
if ~isfile(rf_file)                                   % fall back to any bar_* for that date
    d = dir(fullfile(odir, sprintf('%s_%s_*channel_target_summary.txt', animal, rf_date)));
    if isempty(d), warning('RF summary for %s (%s) not found — no electrode centres.', animal, rf_date); return; end
    rf_file = fullfile(d(1).folder, d(1).name);
end
t = readtable(rf_file, 'Delimiter','\t');
x = t.RF_Center_X; y = t.RF_Center_Y;                 % may import as cell if '-' present
if iscell(x), x = str2double(x); end
if iscell(y), y = str2double(y); end
st = t.Status;                                        % char/cell/string/categorical
if ~iscell(st), st = cellstr(string(st)); end
cx = screen_xy(1)/2; cy = screen_xy(2)/2;             % fixation = screen centre (840, 525)
n = min(nCh, height(t));
xy_deg(1:n,1) = (x(1:n) - cx) / ppd;                  % NaN passes through
xy_deg(1:n,2) = (y(1:n) - cy) / ppd;
valid_gauss(1:n) = strcmp(strtrim(st(1:n)), 'Valid_Gaussian') & ...
                   isfinite(xy_deg(1:n,1)) & isfinite(xy_deg(1:n,2));
end

function [x_px, y_px] = target_pix_per_pos(positions, base, animal)
% Target location per stimulus position in screen pixels:
%   x = trialinfo col 16 (== the 'positions' values themselves),
%   y = trialinfo col 17 (looked up from ph_all_sess.mat; 0 if unavailable).
% (units confirmed in code/RF_Mapping/chan_loc_mua.m: x/y_target_pix =
%  trialinfo(:,16)/(:,17), fixation-centred pixels.)
x_px = positions(:);
y_px = zeros(numel(positions),1);
f = fullfile(base, ['results_' animal], 'multi_lin_reg','cp10_till_100','ph_all_sess.mat');
if ~isfile(f)
    warning('ph_all_sess.mat not found for %s — using x only (y=0).', animal); return
end
try
    Sd = load(f,'ph_comb'); ti = Sd.ph_comb.trialinfo;
    for i = 1:numel(positions)
        r = ti(:,16) == positions(i);
        if any(r), y_px(i) = median(ti(r,17),'omitnan'); end
    end
catch ME
    warning('Could not read trialinfo col 17 (%s) — using x only (y=0).', ME.message);
    y_px = zeros(numel(positions),1);
end
end

function sig = load_coh_sig_mask(coh_root, nCh, nFreq, alpha)
% Per-channel coherence-significance mask [nCh x nFreq] (max-stat across
% frequency), identical to traveling_planar_wave.m / origin script.
sig = false(nCh, nFreq);
for ch = 1:nCh
    cfile = fullfile(coh_root, num2str(ch), 'coherence.mat');
    pfile = fullfile(coh_root, num2str(ch), 'coh_perm_complex.mat');
    if ~isfile(cfile) || ~isfile(pfile), continue; end
    Cc = load(cfile, 'coh', 'coh_complex');
    Pp = load(pfile, 'coh_perm_complex');
    if isfield(Cc,'coh'),             obs = Cc.coh(:)';
    elseif isfield(Cc,'coh_complex'), obs = abs(Cc.coh_complex(:))';
    else, continue; end
    perm_mag = abs(Pp.coh_perm_complex);
    if any(isnan(obs)) || any(isnan(perm_mag(:))), continue; end
    thr = quantile(max(perm_mag,[],2), 1-alpha);
    n = min(numel(obs), nFreq);
    sig(ch,1:n) = obs(1:n) >= thr;
end
end
