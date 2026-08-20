% =====================================================================
% Traveling-wave test across STIMULUS LOCATIONS. De-rotate each location's
% phase by theta(p) = 2*pi*f*d(p)/v and sweep a (frequency x speed) grid for
% the parameters that make the locations align.
%
% d(p) = stimulus eccentricity in deg, v in deg/s. A genuine wave gives a band
% of high resultant R at a frequency-independent speed; a chance per-frequency
% alignment does not. Gain dR = R(f,v) - R0(f), with R0 the same statistic at
% k = 0 (no wave assumed), isolates the wave-specific ramp.
%
% THREE MODES -- what distance goes into k*d. Null, thresholds and figures are
% shared; only the de-rotation differs.
%   visual           R_c = |SUM_p w_p e^{i(phi_p - k d_p)}| / SUM_p w_p, then
%                    averaged over coherence-significant channels. The
%                    per-channel |.| cancels each channel's constant offset,
%                    and with it any per-electrode delay.
%   visual_coherent  de-rotate by k*(d_p + Dc_c), Dc_c = electrode RF
%                    eccentricity, and sum coherently across channels, so the
%                    per-electrode delay is visible. Tests one wave sweeping
%                    outward from the fovea, independent of the stimulus.
%                    Assumes a common phase reference; unmodelled offsets can
%                    only lower R, so a band is trustworthy but a null is
%                    ambiguous.
%   visual_arrival   de-rotate by k*(d_p + a(c,p)), a = visual-field
%                    separation between channel c's RF centre and location p:
%                    the stimulus lands on its retinotopic patch and spreads
%                    at the same speed v. a(c,p) varies across locations
%                    within a channel, so it survives the per-channel |.| and
%                    can be fitted incoherently -- the only mode expressing
%                    "the wave starts where the stimulus is" while staying
%                    robust to per-channel offsets. |.|_c projects out the row
%                    mean, so only the position-varying part is tested.
%
% TWO ESTIMATORS -- what is vector-summed. Both run in one pass and share the
% de-rotation, the modes, R0, dR, the null and the thresholds.
%   'phase'      resultant of the nPos per-location preferred-phase vectors
%                from phase_progression.m, weighted by coh_mag:
%                R_c = |SUM_p coh_mag e^{i(pref_phase - k d_p)}| / SUM_p coh_mag.
%                Not a coherence: the trials were averaged away upstream, and
%                coh_mag ~ 1/sqrt(n_p) gives locations with fewer trials more
%                weight. R0 ~ 0.80-0.93 here, leaving a ceiling of 0.07-0.20.
%   'coherence'  one trial-level coherence, every trial de-rotated by the
%                prediction for its own location:
%                c(f,v) = SUM_t y_t e^{i(phi_tf - k d_p(t))} / SUM_t |y_t|.
%                At k = 0 this is exactly phase_coherence.m with locations
%                pooled. R0 ~ 0.06-0.08, so there is headroom for a gain.
%   R is not comparable between the two -- judge each against its own
%   threshold. A lower R0 under 'coherence' is headroom, not a weaker result.
%
%   Implementation: the rotation depends on a trial only through its location,
%   so both estimators collapse onto per-location complex sums
%   S(p,f) = SUM over trials at p of y_t*e^{i phi_tf}, and the three modes are
%   re-weightings of S. functions/trial_position_sums_chan.m produces S for the
%   observed labels and all nPerm shuffles, one SLURM job per channel; the
%   cache keeps the name 'trial_position_sums' (only the estimator was renamed
%   'trial' -> 'coherence').
%
% Per animal throughout, never pooled: grid -> animals -> mean grid plus
% "significant in both" replication. Null: one synchronised shuffle of the
% location<->distance assignment (plus electrode<->Dc for coherent,
% electrode<->a(c,:) for arrival) shared across the grid, max over the grid ->
% max-statistic threshold. R0 is shuffle-invariant, so dR gets a clean null too.
%
% COORDINATE FRAME. Stimulus positions and RF centres are stored with different
% origins, so both are converted to fixation-centred degrees, (0,0) = fovea:
%   stimulus   trialinfo col-16/17 are already fixation-centred pixels
%              -> xy_deg = [col16 col17] / ppd
%   RF centre  *_channel_target_summary.txt holds SCREEN pixels
%              -> xy_deg = ([RF_Center_X RF_Center_Y] - [840 525]) / ppd
% ppd = PIX_PER_DEG, the per-animal rig calibration. Both then give
% eccentricity = hypot(x,y), so electrode-minus-stimulus is frame-safe.
% (Resolved from code/RF_Mapping/chan_loc_mua.m.)
%
% RF_VALID_ONLY keeps only channels whose RF centre is a real Gaussian fit;
% 'Extrapolated' centres were filled from array neighbours or a plane fit after
% a failed fit, i.e. inferred rather than measured. Affects visual_coherent and
% visual_arrival only. Output names carry a '_validRF' tag when it is on, so
% both variants coexist on disk.
%
% OUTPUT -- three figures, in
% Plots/scanning/stimulus_loc_wave/cp10_till_100/<dv>/
%   stimulus_loc_grids.pdf       R. 6 rows = 3 modes x 2 estimators (estimator
%                                named in every title, since R differs between
%                                them), 4 cols = animals, mean/replication,
%                                pooled z.
%   stimulus_loc_gain_<est>.pdf  dR for one estimator, 3 x 3.
%   results_combined/scanning/stimulus_loc_wave/cp10_till_100/<dv>/
%     stimulus_loc_wave.mat      both estimators: results.G.(est).(mode)(animal)
% The gain figures carry no pooled column: R0 depends on neither speed nor
% shuffle, so zG = (R - R0 - (muR - R0))/sdR = zR (verified on saved output,
% max|zR - zG| = 1.75e-13). It would duplicate the pooled R panel.
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
MIN_TRIALS_FRAC = 0.5;             % drop stimulus positions sampled with fewer than this
                                   % fraction of the median trials/position (0 = keep all)
FREQ_RANGE = [2 100];               % Hz swept (low-freq band the hypothesis lives in)
NV         = 30;                   % # speed steps
V_VISUAL   = logspace(log10(1),   log10(200), NV);   % deg/s
rng(2025);

% Experiment calibration
PIX_PER_DEG   = struct('hermes', 53.24, 'klecks', 50.56);   % pixels per degree (per-animal rig calibration)

% Electrode RF centres (COHERENT and ARRIVAL modes). Read from the per-channel
% summary table written by RF_Mapping/mapping_lfp.m -- the same file
% cortical_planar_wave_PGD.m uses. RF_Center_X/Y are screen pixels, so
% elec_rf_deg subtracts the screen centre before /ppd; see COORDINATE FRAME
% above. Both arrays are 8x8 and the summary's 64 rows match the
% phase-progression channel index 1:1, so no remap is needed.
RF_DATE   = '20170829';   % RF session both animals share
SCREEN_XY = [1680 1050];  % screen pixels; fixation/fovea at the centre

% Use ONLY channels whose RF centre is a real Gaussian fit (Status ==
% 'Valid_Gaussian')? 'Extrapolated' centres were interpolated from array
% neighbours after a failed fit, not measured (see header). Affects
% visual_coherent (via Dc) and visual_arrival (via a(c,p)); plain visual never
% uses RF centres and is unchanged either way.
RF_VALID_ONLY = false;

% Which estimators to run -- see TWO ESTIMATORS in the header. Both are
% computed in one pass; the grids figure shows them side by side and each gets
% its own gain figure. 'coherence' additionally needs ph_all_sess.mat and one
% SLURM job per channel. Order sets figure row order; a single name runs one.
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

% Three modes, all in visual degrees / deg-per-second, differing only in the
% distance de-rotated by (see the header for the full statement):
%   visual          k*d_p, |resultant| per channel then averaged -- robust to
%                   per-channel offsets, blind to per-electrode delays.
%   visual_coherent k*(d_p + D_c). D_c has no p in it, so a per-channel |.|
%                   would delete it exactly -- this mode must sum coherently.
%   visual_arrival  k*(d_p + a(c,p)). a dips at the location falling on the
%                   channel's own RF, so it varies within a channel and
%                   survives the per-channel |.|.
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
    [S, pos_keep] = drop_undersampled_positions(S, pp, MIN_TRIALS_FRAC);
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

    % ARRIVAL-DISTANCE table a(c,p) for visual_arrival [nCh x nPos]: visual-field
    % separation between channel c's RF centre and location p. A straight
    % difference, legitimate only because both are now in the same
    % fixation-centred frame. Row c dips at the location landing on its RF. No
    % per-row referencing: the per-channel |.| deletes a row constant anyway.
    a_arr = hypot(rf_deg(:,1) - stim_deg(:,1).', rf_deg(:,2) - stim_deg(:,2).');
    if RF_VALID_ONLY, a_arr(~rf_valid, :) = NaN; end
    fprintf('  arrival distance a(c,p): %d/%d channels usable | %.2f-%.2f deg\n', ...
        sum(all(isfinite(a_arr),2)), nCh, min(a_arr(:),[],'omitnan'), max(a_arr(:),[],'omitnan'));

    % 'coherence' estimator: fetch the per-location complex trial sums S(p,f)
    % (one SLURM job per channel, cached on disk). The 'phase' estimator never
    % touches these — it works entirely from phase_progression.mat.
    if any(strcmp(ESTIMATORS,'coherence'))
        [Ssum, Sperm, Wch] = get_trial_sums(base, animalName, dv, nCh, nPerm, RECOMPUTE_TRIAL_SUMS);
        % The cache is built over ALL positions, so drop the same ones here —
        % dim 2 has to stay aligned with d/d_vis/a_arr, which are now nPos-long.
        % Wch is a global sum(|y|) over all trials and is left as-is: it is a
        % per-channel scale that the R statistic normalises out.
        Ssum  = Ssum(:, pos_keep, :);
        Sperm = Sperm(:, pos_keep, :, :);
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
        % The animals sit on very different scales (R0 ~0.80 vs ~0.93,
        % thresholds 4-6x apart), so a raw cross-animal average would be
        % dominated by the larger numbers. The same permutations centre/scale
        % here and build the pooled threshold below -- a standardised
        % ("pseudo-t") max-stat test.
        % 'omitnan' is essential: mean/std do not skip NaN, so one NaN
        % permutation would make that cell NaN for every permutation and the
        % NaN would spread through the cross-animal sum.
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
% Two criteria, reported side by side, answering different questions.
%
%   REPLICATION  repl = sig(hermes) AND sig(klecks), cell by cell. Each animal
%                is already max-stat corrected over its grid, so a replicated
%                cell is strong evidence and cannot be driven by one animal.
%                It cannot aggregate weak evidence and has low power. PRIMARY.
%
%   POOLED       average the standardised grids and test against their own
%                max-stat null, pairing permutation b of one animal with
%                permutation b of the other (legitimate: the animals are
%                independent, so the pairing samples the product null):
%                    Zobs(f,v)    = mean_a z_a(f,v)
%                    Znull(f,v,b) = mean_a z_a_null(f,v,b)
%                    thr_pool     = (1-alpha) quantile of max over grid
%                More power, but with two animals a significant cell can be
%                almost entirely one of them. SECONDARY -- the per-animal z at
%                that cell is printed below.
%
% With n=2 neither supports a population-level inference; between-animal
% variance is not estimable. Run separately per estimator -- the two never
% share a null or a threshold, their R values being on different scales.
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
% Escape underscores rather than switching the TeX interpreter off, so mode
% names survive ('visual_coherent' would render as a subscript) while the
% intentional TeX (\DeltaR, R_0, 2\pi) still renders. Cell-array titles put one
% line per element, which stops them overlapping.
esc = @(s) strrep(s, '_', '\_');

% Map a speed in physical units onto the index-based y axis (shared by both
% figures, so the best-speed line can be drawn on the R grid as well as on the
% gain grid).
speeds_to_idx = @(sp,sv) interp1(sv, 1:numel(sv), sp, 'linear', NaN);

%% ─── FIGURE 1: the DE-ROTATION R grids, both estimators ──────────────
% Rows = 3 modes x 2 estimators (named in every title -- R means a different
% thing in the two blocks). Cols = animals + mean/replication + pooled z.
% Colour limits are per estimator block: a fixed [0 1] scale renders the
% phase-coherence block (R ~ 0.08) as flat blue.
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
% A wave = significant positive gain (white contour) along a near-constant best
% speed (black line); a rising black line means constant wavenumber, i.e. a
% fixed phase offset rather than propagation.
% No pooled column -- the pooled gain z equals the pooled R z (see header).
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
%%  align_grid*  ESTIMATOR 'phase'. Input: pref_phase + coh_mag, i.e. nPos
%%               pre-averaged per-location vectors. R = weighted resultant
%%               length over those locations, not a coherence; R0 runs high
%%               (~0.8-0.93) because a handful of angles concentrate easily.
%%               Weighted by |c_p| ~ 1/sqrt(n_p), so fewer-trial locations
%%               count more.
%%  coh_grid*    ESTIMATOR 'coherence'. Input: S(p,f), the raw complex sums
%%               over trials at each location. R = a genuine phase coherence
%%               over all trials, each de-rotated by its own location's lag;
%%               at k=0 it equals phase_coherence.m with locations pooled.
%%               Every trial counts once.
%%
%%  Otherwise identical: same six outputs, R0 defined as each one's own
%%  statistic at k=0, same gain / null / threshold / pooling code. R values are
%%  not comparable between the two.
%% =====================================================================
function [S, pos_keep] = drop_undersampled_positions(S, pp, min_frac)
% Drop stimulus positions sampled with far fewer trials than the rest.
% Relevant to the 'phase' estimator only: as the Helpers note above says, its
% weight |c_p| ~ 1/sqrt(n_p), so a position with few trials counts MORE, not
% less. pos_keep indexes back into the ORIGINAL position numbering and must be
% applied to anything else built over all positions (Ssum/Sperm).
pos_keep = 1:size(S.pref_phase,3);
if min_frac <= 0, return; end
w = load(pp, 'n_pos');
if ~isfield(w,'n_pos') || isempty(w.n_pos), return; end
trials = double(max(w.n_pos, [], 1));          % dead channels sit at 0 -> use max
keep   = trials >= min_frac * median(trials);
if all(keep), return; end
for pd = find(~keep)
    fprintf('  DROP position %d (coord %g): %d trials vs median %d — under-sampled\n', ...
        pd, S.positions(pd), trials(pd), round(median(trials)));
end
pos_keep     = find(keep);
S.pref_phase = S.pref_phase(:,:,keep);
S.coh_mag    = S.coh_mag(:,:,keep);
S.positions  = S.positions(keep);
end

function [Robs, R0, Rnull_max, Gnull_max, Rnull, Gnull] = align_grid(pref, coh, coh_sig, f_use, fHz, d, Vs, MIN_LOC, nPerm)
% ESTIMATOR 'phase', mode 'visual'. Incoherent across channels.
% R(f,v) = mean over coherence-sig channels of the de-rotated resultant across
%          stimulus locations.
% R0(f)  = the same resultant unrotated (k=0, v->inf): the baseline the gain
%          dR = R - R0 is measured against, so dR is invariant to locations
%          merely sharing a preferred phase.
% Returns the observed grid, the baseline, and per-perm grid maxima of R and of
% dR. R0 does not use d, so the location<->distance shuffle leaves it
% unchanged -> a clean gain null.
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
% ESTIMATOR 'phase', mode 'visual_coherent'. One cortical wave in which both
% the stimulus eccentricity d_p and the electrode RF eccentricity Dc_c
% contribute k*distance, k = 2*pi*f/v:
%   R(f,v) = | Σ_{c,p} w_cp e^{i(φ_cp − k d_p − k Dc_c)} | / Σ w_cp
%          = | (e^{-ik·Dc})ᵀ · A · e^{-ik·d} | / Σ w_cp .
% No per-channel |·|, so the per-electrode delay is visible -- but unmodelled
% offsets can only lower R, so a band is trustworthy and a null ambiguous.
% Channels without a finite Dc are dropped. R0 = the k=0 coherent resultant.
% Null shuffles d→location and Dc→channel together, one shared shuffle per
% permutation.
nF = numel(f_use); nV = numel(Vs); nPos = numel(d); nCh = numel(Dc);
Robs      = nan(nF, nV);
R0        = nan(nF, 1);
Rnull_max = nan(nPerm, 1);
Gnull_max = nan(nPerm, 1);
Rnull     = nan(nF, nV, nPerm);   % full null grid, kept for cross-animal pooling
Gnull     = nan(nF, nV, nPerm);
Dc = Dc(:);

% one shared shuffle per permutation, over locations AND channels, reused
% across every (f,v) cell -> synchronised null + grid-wide max-stat.
% NaN-safe: channels with no finite Dc are excluded by useC, but an
% unrestricted randperm would shuffle one into a used slot and make that
% permutation NaN. NaN > rb is false, so the running max would stay -inf and
% every cell would come out significant. Permute only among the
% finite-geometry channels.
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
% ESTIMATOR 'phase', mode 'visual_arrival'. Arrival-time model, incoherent
% across channels.
%   R_c(f,v) = | Σ_p w_cp e^{i(φ_cp − k(d_p + a(c,p)))} | / Σ_p w_cp ,
%   R(f,v)   = mean over coherence-significant channels of R_c ,   k = 2πf/v.
%
% a [nCh x nPos] is the visual-field separation (deg) between channel c's RF
% centre and location p. Unlike Dc_c it varies across locations within a
% channel (it dips where the stimulus lands on that channel's RF), so the
% per-channel |.| does not remove it and no common-reference assumption is
% needed. |.|_c deletes the row mean, so only the position-varying part is
% tested. One v serves both terms, so the grid stays frequency x speed.
%
% Null: one synchronised shuffle per permutation, permuting the geometry two
% ways -- columns location <-> (d_p, a(:,p)), breaking "phase tracks the
% stimulus"; rows channel <-> a(c,:), breaking "the dip is at the right
% electrode". R0 (k=0) uses neither d nor a, so it is shuffle-invariant and
% numerically identical to align_grid's R0.
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
% TRIAL-level estimators. Same three modes, outputs and nulls as the align_grid*
% family; only R differs -- here it is a genuine phase coherence over trials,
% each de-rotated by its own location's predicted lag.
%
% Everything works off the per-location complex sums from
% functions/trial_position_sums_chan.m,
%       S(c,p,f) = SUM over trials at location p of y_t*exp(i*phi_tf)
%       W(c)     = SUM over trials of |y_t|
% because the de-rotation depends on a trial only through its location:
%       c_c(f,v) = ( 1/W_c ) * SUM_p S(c,p,f)*exp(-i*k*d_p) .
% At k=0 this is the phase_coherence/ measure with all locations pooled.
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
% ESTIMATOR 'coherence', mode 'visual_coherent'. The trial-level twin of
% align_grid_coherent: channels summed complex, each de-rotated by its own
% k*Dc, so the per-electrode delay is testable. The location relabelling is
% already synchronised across channels inside Sperm; Dc is permuted here.
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
% frequency). Kept identical to the copies in cortical_planar_wave_PGD.m
% and cortical_planar_wave_derotation.m — change all three together.
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
