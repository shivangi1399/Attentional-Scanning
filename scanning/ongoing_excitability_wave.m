% =====================================================================
% ONGOING EXCITABILITY PROPAGATION -- does an EXCITABILITY signal move across
% the array in the ONGOING, pre-target LFP?
%
% WHY THIS EXISTS
%   Every wave analysis in this folder so far measures PHASE, and the whole
%   paper turns on the fact that phase alone cannot separate propagation from a
%   fixed spatial offset: at one frequency cos(2*pi*f*t - k*x) is a wave at
%   v = 2*pi*f/k for any k. cortical_planar_wave_derotation.m settled that
%   across frequency and found a diagonal ridge -- constant wavenumber, i.e. a
%   frozen offset, not propagation.
%
%   erp_latency_wave.m is the one existing script that reads TIME instead of
%   phase, but it reads the time of the EVOKED response, and its own header
%   says so: a null there does not refute a wave in ongoing activity. Nothing
%   in the folder yet asks whether ONGOING excitability propagates. That is the
%   substance of the definition Fries (2023) actually uses -- "waves of enhanced
%   neuronal excitability" sweeping a representation -- so it is the gap this
%   script fills.
%
% THE KEY IDEA -- WHY AN ENVELOPE AND NOT A PHASE
%   This script cross-correlates an AMPLITUDE ENVELOPE (or MUA), not a phase.
%   An envelope has no carrier, so its peak lag is a TIME in milliseconds, with
%   no k-versus-v ambiguity and nothing to rescale with frequency. The entire
%   fixed-offset confound that killed the phase analyses cannot arise here:
%       real propagation  ->  tau grows linearly with distance, one speed
%       fixed offset      ->  no tau at all; a static phase map has zero lag
%   That is the point of running it. A positive result here would be the first
%   evidence in this dataset for definition-3 propagation (something moves);
%   a null closes the last route by which the negative phase results could be
%   an artefact of measuring phase.
%
% ESTIMATOR -- leave-one-out mean reference, the same logic as erp_latency_wave.m
%   Per trial, each channel's envelope is z-scored and cross-correlated against
%   the mean envelope of all OTHER channels (leave-one-out, so a channel does
%   not correlate with itself through the reference). Cross-spectra are summed
%   over trials and inverse-transformed ONCE at the end, so the correlogram is
%   averaged before the peak is picked -- picking a peak per trial and averaging
%   the peaks would average a noise floor instead. The peak is refined by
%   parabolic interpolation, giving sub-sample lag resolution, which matters:
%   0.1-0.8 m/s across 2.8 mm of array is only 3.5-28 ms end to end.
%
%   tau_c is then referenced to the array mean (mean 0 by construction, so only
%   the SPREAD is meaningful) and a plane is fitted in lag space,
%       tau_c ~ s . r_c ,   s = slowness vector (ms/mm),  speed = 1/|s|
%   The fit's R2 is the test statistic; |s| gives the speed and angle(s) the
%   propagation direction, which is then compared with the PGD direction.
%
%   DO_PAIRWISE additionally fits all 2016 electrode pairs (tau_cd against
%   (r_d - r_c)), which is more robust but O(nChan^2) in the cross-spectrum
%   accumulation. Off by default; use it with MAX_SESS set small.
%
% TWO WINDOWS, and the comparison is the point
%   PRE_WIN   ongoing activity before the target -- the question being asked.
%   POST_WIN  the evoked transient -- a POSITIVE CONTROL. Muller et al. (2014)
%             and Zanos et al. (2015) both report propagation in exactly this
%             epoch, so if the method cannot find a lag gradient in POST_WIN it
%             is not sensitive enough to interpret a null in PRE_WIN. Read
%             POST first; PRE is only interpretable if POST is positive.
%
% NULLS -- three, because a smooth lag map is easy to produce by accident
%   (1) POSITION shuffle. Shuffle which array coordinate each channel sits at
%       and refit the plane. Same null as cortical_planar_wave_PGD.m and the
%       de-rotation scripts. Kills "any spatially smooth lag field looks planar".
%   (2) TRIAL shuffle. Correlate channel c on trial t against the reference
%       built from trial t' /= t. Anything stimulus-locked, any common drift and
%       any fixed filter delay survives this; only genuine within-trial
%       co-fluctuation does not. This is the strong control.
%   (3) SIGN-FLIP of the lag axis is not needed -- the plane fit is signed and
%       the position shuffle already covers it.
%   Both nulls are max-stat over the same quantity that is reported.
%
% WHAT A POSITIVE LOOKS LIKE
%   R2 above both nulls, AND a speed inside SPEED_OK (0.1-0.8 m/s, the
%   unmyelinated horizontal-fibre range of Muller et al. 2018 -- the band that
%   licenses calling a wave first-order/mesoscopic), AND replication in both
%   animals. A speed far outside that band is not a mesoscopic wave even if R2
%   is significant; it is reported, flagged, and must be discussed as such.
%
% WHAT THIS IS NOT
%   Not a phase analysis and not comparable to PGD numbers. Not a test of the
%   scanning hypothesis: it asks whether anything propagates across CORTEX, not
%   whether it sweeps the stimulus representation. The retinotopic version of
%   the question is the 'retinotopic' axis mode below, which projects the same
%   lags onto RF eccentricity instead of array distance -- that is the scanning
%   readout, and it is the one that matters for RAS.
%
% DATA  <base>/results_<animal>/<session>/clean_lfp.mat -> `clean_data`
%       <base>/results_<animal>/<session>/clean_mua.mat -> `clean_mua`
%   FieldTrip raw structures: .trial {1 x nTrials} each [nChan x nTime], .time
%   (0 = stimulus event), .label, .trialinfo, .fsample. V4 channels map onto
%   canonical slots 1..64 by label number -- see functions/v4_channel_slots.m.
%
% OUTPUT
%   Plots/scanning/ongoing_excitability/cp10_till_100/<SIGNAL>/ongoing_excitability_wave.pdf
%   results_combined/scanning/ongoing_excitability/cp10_till_100/<SIGNAL>/ongoing_excitability_wave.mat
% =====================================================================

clearvars; close all; clc
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/scanning/functions

%% --- Settings -------------------------------------------------------
animals    = {'hermes','klecks'};
base       = '/mnt/hpc/projects/MWSampling/4Shivangi';
grid_rows  = 8; grid_cols = 8;
SPACING_MM = 0.4;                  % electrode pitch (mm)
nChTot     = grid_rows*grid_cols;
rng(2025);

% Which excitability signal. Run one per invocation; the output folder is
% named after it so the variants do not overwrite each other.
%   'gamma_env'  LFP band-pass GAMMA_HZ -> Hilbert -> |.|   (the standard proxy)
%   'bb_env'     LFP band-pass BB_HZ    -> Hilbert -> |.|   (broadband/high-gamma)
%   'mua'        clean_mua trial data, rectified is not needed (already MUA),
%                smoothed with SMOOTH_MS -- the most direct excitability measure
SIGNAL     = 'gamma_env';
GAMMA_HZ   = [40 90];
BB_HZ      = [60 150];
SMOOTH_MS  = 5;                    % envelope/MUA smoothing before correlation

% Windows relative to the event at t = 0 in clean_data.time.
% PRE closes well before the earliest critical time (54 ms in monkey H) so no
% evoked activity leaks in; the analysis windows elsewhere in the paper close
% 100 ms after onset, and POST_WIN matches that.
PRE_WIN    = [-0.60 -0.05];        % s -- the ongoing question
POST_WIN   = [ 0.02  0.12];        % s -- the positive control
WINDOWS    = {'pre','post'};

MAX_LAG_MS = 25;                   % lag search half-width
MIN_CH     = 40;                   % channels a trial must have to be used
MIN_TRIALS = 200;                  % per animal, to attempt an estimate
MAX_SESS   = Inf;                  % set small (e.g. 3) for a quick run
DO_PAIRWISE = false;               % all-pairs fit as well (slow -- see header)

nPerm      = 1000;                 % position shuffles
nPermTrial = 200;                  % trial shuffles (each is a full re-estimate)
alpha      = 0.05;

% Plausible mesoscopic propagation band, Muller et al. (2018): unmyelinated
% long-range horizontal fibres of the superficial layers.
SPEED_OK   = [10 80];              % cm/s

AXIS_MODES = {'cortical','retinotopic'};
PGD_BAND   = [13.3 23.3];          % Hz -- the replicated beta band, for the
                                   % direction comparison only
% Retinotopic axis calibration -- same constants as erp_latency_wave.m
PIX_PER_DEG   = struct('hermes', 53.24, 'klecks', 50.56);
RF_DATE       = '20170829';
SCREEN_XY     = [1680 1050];
RF_VALID_ONLY = true;

dv_tag  = 'cp10_till_100';         % kept for path consistency with the folder
out_dir = fullfile(base,'Plots','scanning','ongoing_excitability', dv_tag, SIGNAL);
res_dir = fullfile(base,'results_combined','scanning','ongoing_excitability', dv_tag, SIGNAL);
if ~exist(out_dir,'dir'), mkdir(out_dir); end
if ~exist(res_dir,'dir'), mkdir(res_dir); end

% Array coordinates -- SAME index<->position map as the PGD, de-rotation and
% ERP-latency scripts. Do not change one without the others.
ch_col = ceil((1:nChTot)' / grid_rows);
ch_row = grid_rows - mod((1:nChTot)' - 1, grid_rows);
XY     = [ch_col, ch_row] * SPACING_MM;         % nCh x 2, mm
XY     = XY - mean(XY,1);                       % centre, so the fit intercept
                                                % is the array-mean lag

% =====================================================================
% ASSUMPTIONS TO CHECK BEFORE TRUSTING THE OUTPUT
% =====================================================================
% (1) clean_data.time is zeroed on the same STIMULUS event the phase pipeline
%     used. A different zero shifts both windows together; PRE_WIN would then
%     no longer be guaranteed target-free, which WOULD matter here. The script
%     prints the first session's time axis -- eyeball it.
% (2) Trials whose time axis does not cover a window are skipped, not padded.
%     The count is printed per window; a large drop means the epoching is
%     shorter than assumed and PRE_WIN must be shortened.
% (3) Filtering is zero-phase (filtfilt) throughout. A causal filter would
%     impose a group delay, which is the very quantity being measured. Any
%     per-channel delay from the acquisition hardware is common to all
%     channels and cancels in the array-mean reference; a per-channel hardware
%     delay would not, and cannot be detected from these data.
% (4) MUA and LFP come from the same electrodes but different files; both are
%     mapped through v4_channel_slots so slot i is the same electrode in both.
% =====================================================================

fprintf('ONGOING EXCITABILITY PROPAGATION | signal = %s\n', SIGNAL);
fprintf('  PRE  %.3f..%.3f s (the question)\n', PRE_WIN(1), PRE_WIN(2));
fprintf('  POST %.3f..%.3f s (positive control -- read this first)\n', POST_WIN(1), POST_WIN(2));

A = struct();

%% --- Per animal -----------------------------------------------------
for ia = 1:numel(animals)
    animalName = animals{ia};
    sess = list_sessions(base, animalName);
    if isempty(sess)
        warning('No sessions found for %s -- skipping.', animalName); continue
    end
    nS = min(numel(sess), MAX_SESS);
    fprintf('\n##### %s : %d sessions #####\n', animalName, nS);

    A(ia).animal = animalName;

    for iw = 1:numel(WINDOWS)
        win_name = WINDOWS{iw};
        switch win_name
            case 'pre',  WIN = PRE_WIN;
            case 'post', WIN = POST_WIN;
        end

        % ---- accumulate cross-spectra over every trial of every session ----
        ACC = [];                       % lazily sized once fs is known
        for is = 1:nS
            [X, fs, nUsed, nSkip, tinfo] = load_envelopes(base, animalName, sess{is}, ...
                SIGNAL, GAMMA_HZ, BB_HZ, SMOOTH_MS, WIN, nChTot, MIN_CH);
            if isempty(X), continue, end
            if isempty(ACC), ACC = init_acc(size(X,2), nChTot, DO_PAIRWISE); ACC.fs = fs; end
            ACC = accumulate_xspec(ACC, X, DO_PAIRWISE);
            fprintf('  . %-45s %5d trials used, %4d skipped\n', ...
                shortname(sess{is}), nUsed, nSkip);
        end
        nGot = 0; if ~isempty(ACC), nGot = ACC.nTrials; end
        if nGot < MIN_TRIALS
            warning('  %s/%s: only %d trials -- window skipped.', animalName, win_name, nGot);
            continue
        end

        % ---- observed lags and plane fit ---------------------------------
        maxlagS = round(MAX_LAG_MS/1000 * ACC.fs);
        tau     = lags_from_acc(ACC, maxlagS);              % s, nCh x 1
        W       = ACC.nPerCh(:);                            % weight = trials seen
        [sl, R2, dirRad, speed_cms] = fit_lag_plane(tau, XY, W);

        % ---- null 1: position shuffle ------------------------------------
        R2n = nan(nPerm,1);
        for p = 1:nPerm
            pm = randperm(nChTot);
            [~, R2n(p)] = fit_lag_plane(tau, XY(pm,:), W(pm));
        end
        p_pos = (1 + sum(R2n >= R2)) / (1 + nPerm);

        % ---- null 2: trial shuffle ---------------------------------------
        % Re-estimate tau with the reference taken from a DIFFERENT trial. The
        % accumulator stores what is needed to do this cheaply; see the helper.
        R2t = nan(nPermTrial,1);
        if isfield(ACC,'canShuffleTrials') && ACC.canShuffleTrials
            for p = 1:nPermTrial
                tau_s = lags_from_acc(ACC, maxlagS, true, p);
                [~, R2t(p)] = fit_lag_plane(tau_s, XY, W);
            end
            p_trl = (1 + sum(R2t >= R2)) / (1 + nPermTrial);
        else
            p_trl = NaN;
        end

        % ---- retinotopic projection: the SCANNING readout ----------------
        RET = struct('slope',NaN,'p',NaN,'speed',NaN,'n',0);
        ppd = PIX_PER_DEG.(animalName);
        [rf_deg, rf_valid] = elec_rf_deg(base, animalName, nChTot, RF_DATE, ppd, SCREEN_XY);
        ecc = sqrt(sum(rf_deg.^2, 2));                       % deg from fovea
        useR = isfinite(ecc) & isfinite(tau);
        if RF_VALID_ONLY, useR = useR & rf_valid(:); end
        if sum(useR) >= 10
            cf = polyfit(ecc(useR), tau(useR), 1);            % s per deg
            rr = corr(ecc(useR), tau(useR));
            % null: same position shuffle, applied to the eccentricity labels
            rn = nan(nPerm,1);
            iu = find(useR);
            for p = 1:nPerm
                pm = iu(randperm(numel(iu)));
                rn(p) = corr(ecc(pm), tau(useR));
            end
            RET.slope = cf(1); RET.r = rr; RET.n = sum(useR);
            RET.p     = (1 + sum(abs(rn) >= abs(rr))) / (1 + nPerm);
            RET.speed = 1/cf(1);                              % deg/s
        end

        % ---- direction agreement with the phase gradient ------------------
        th_pgd  = pgd_direction(base, animalName, 'lfp', PGD_BAND);
        d_ang   = NaN;
        if isfinite(th_pgd), d_ang = rad2deg(angdiff_(dirRad, th_pgd)); end

        % ---- optional all-pairs fit --------------------------------------
        PW = struct('speed',NaN,'R2',NaN,'dir',NaN);
        if DO_PAIRWISE && isfield(ACC,'CS')
            [PW.speed, PW.dir, PW.R2] = fit_lag_pairs(ACC, XY, maxlagS);
        end

        % ---- report -------------------------------------------------------
        okSpeed = abs(speed_cms) >= SPEED_OK(1) && abs(speed_cms) <= SPEED_OK(2);
        fprintf('\n  --- %s / %s (%d trials) ---\n', animalName, upper(win_name), ACC.nTrials);
        fprintf('    lag spread across channels   %.2f ms (sd), range %.2f ms\n', ...
            1e3*std(tau,'omitnan'), 1e3*(max(tau)-min(tau)));
        fprintf('    lag-plane R2                 %.3f   (position-shuffle p = %.4f)\n', R2, p_pos);
        if isfinite(p_trl)
            fprintf('                                        (trial-shuffle    p = %.4f)\n', p_trl);
        end
        fprintf('    speed                        %.1f cm/s  [%s]\n', speed_cms, ...
            ternary(okSpeed,'inside the 10-80 cm/s mesoscopic band','OUTSIDE the mesoscopic band'));
        fprintf('    direction                    %.0f deg (array coords)', rad2deg(dirRad));
        if isfinite(d_ang), fprintf('   | %+.0f deg from the PGD axis', d_ang); end
        fprintf('\n');
        if RET.n > 0
            fprintf('    retinotopic (SCANNING) slope %+.3f ms/deg  r = %+.2f  p = %.4f  (n = %d ch)\n', ...
                1e3*RET.slope, RET.r, RET.p, RET.n);
        end
        if DO_PAIRWISE
            fprintf('    all-pairs fit                %.1f cm/s, R2 %.3f\n', PW.speed, PW.R2);
        end
        verdict = wave_verdict(R2, p_pos, p_trl, okSpeed, alpha);
        fprintf('    --> %s\n', verdict);

        A(ia).(win_name) = struct('tau',tau,'W',W,'R2',R2,'p_pos',p_pos,'p_trl',p_trl, ...
            'slowness',sl,'speed_cms',speed_cms,'dir',dirRad,'dir_vs_pgd',d_ang, ...
            'R2_null_pos',R2n,'R2_null_trl',R2t,'nTrials',ACC.nTrials, ...
            'RET',RET,'PW',PW,'speed_ok',okSpeed,'verdict',verdict,'fs',ACC.fs);
    end
end

%% --- Replication across animals -------------------------------------
fprintf('\n===== REPLICATION =====\n');
for iw = 1:numel(WINDOWS)
    w = WINDOWS{iw}; ok = true; spd = [];
    for ia = 1:numel(A)
        if ~isfield(A(ia),w) || isempty(A(ia).(w)), ok = false; continue, end
        s = A(ia).(w);
        ok = ok && s.p_pos < alpha && (isnan(s.p_trl) || s.p_trl < alpha) && s.speed_ok;
        spd(end+1) = s.speed_cms; %#ok<SAGROW>
    end
    if numel(spd) == numel(animals) && ok
        fprintf('  %-4s : propagation replicated in both animals (%.0f and %.0f cm/s)\n', ...
            upper(w), spd(1), spd(2));
    else
        fprintf('  %-4s : NOT replicated -- no ongoing propagation claim is supported\n', upper(w));
    end
end

%% --- Figure ---------------------------------------------------------
nA = numel(A);
fg = new_fig(1500, 420*max(nA,1));
for ia = 1:nA
    for iw = 1:numel(WINDOWS)
        w = WINDOWS{iw};
        if ~isfield(A(ia),w) || isempty(A(ia).(w)), continue, end
        S = A(ia).(w);

        % (a) lag map on the array + fitted gradient arrow
        ax = subplot(nA*2, 3, (ia-1)*6 + (iw-1)*3 + 1);
        M = nan(grid_rows, grid_cols);
        for c = 1:nChTot, M(ch_row(c), ch_col(c)) = 1e3*S.tau(c); end
        imagesc(ax, M); axis(ax,'image'); colormap(ax, parula); cb = colorbar(ax);
        cb.Label.String = 'lag (ms)';
        hold(ax,'on');
        q = 2.2*[cos(S.dir) sin(S.dir)];
        quiver(ax, (grid_cols+1)/2, (grid_rows+1)/2, q(1), -q(2), 0, 'k','LineWidth',2,'MaxHeadSize',2);
        title(ax, {sprintf('%s -- %s', esc(A(ia).animal), upper(w)), ...
                   sprintf('lag map, %.1f cm/s', S.speed_cms)}, 'FontSize',9);

        % (b) lag vs distance along the fitted axis
        ax = subplot(nA*2, 3, (ia-1)*6 + (iw-1)*3 + 2); hold(ax,'on');
        d = XY(:,1)*cos(S.dir) + XY(:,2)*sin(S.dir);
        plot(ax, d, 1e3*S.tau, 'o','MarkerFaceColor',[.2 .4 .8],'MarkerEdgeColor','none');
        xx = linspace(min(d),max(d),10);
        plot(ax, xx, 1e3*(xx*norm(S.slowness)), 'k-','LineWidth',1.4);
        xlabel(ax,'distance along fitted axis (mm)'); ylabel(ax,'lag (ms)');
        title(ax, sprintf('R^2 = %.3f, p_{pos} = %.4f', S.R2, S.p_pos),'FontSize',9);
        grid(ax,'on');

        % (c) observed R2 against both nulls
        ax = subplot(nA*2, 3, (ia-1)*6 + (iw-1)*3 + 3); hold(ax,'on');
        histogram(ax, S.R2_null_pos, 30, 'FaceColor',[.7 .7 .7],'EdgeColor','none', ...
            'Normalization','pdf','DisplayName','position shuffle');
        if ~all(isnan(S.R2_null_trl))
            histogram(ax, S.R2_null_trl, 20, 'FaceColor',[.9 .6 .2],'EdgeColor','none', ...
                'FaceAlpha',.6,'Normalization','pdf','DisplayName','trial shuffle');
        end
        xline(ax, S.R2, 'r-','LineWidth',2,'DisplayName','observed');
        xlabel(ax,'lag-plane R^2'); legend(ax,'Location','best','FontSize',7); grid(ax,'on');
        title(ax, S.verdict, 'FontSize', 8);
    end
end
sgtitle({sprintf('ONGOING EXCITABILITY PROPAGATION (%s): does an excitability signal MOVE across V4?', SIGNAL), ...
         'An envelope lag is a TIME, so the fixed-offset confound of the phase analyses cannot arise here.', ...
         'Read POST (positive control) first: if the method finds no lag gradient there, a null in PRE means nothing.'}, ...
         'FontSize', 9);
save_pdf(fg, fullfile(out_dir,'ongoing_excitability_wave.pdf'));

results = struct('A',A,'animals',{animals},'SIGNAL',SIGNAL,'PRE_WIN',PRE_WIN, ...
    'POST_WIN',POST_WIN,'GAMMA_HZ',GAMMA_HZ,'BB_HZ',BB_HZ,'SMOOTH_MS',SMOOTH_MS, ...
    'MAX_LAG_MS',MAX_LAG_MS,'SPACING_MM',SPACING_MM,'nPerm',nPerm, ...
    'nPermTrial',nPermTrial,'alpha',alpha,'SPEED_OK',SPEED_OK,'PGD_BAND',PGD_BAND);
save(fullfile(res_dir,'ongoing_excitability_wave.mat'),'results','-v7.3');
fprintf('\nSaved under %s\n', out_dir);

%% ===================== local functions ==============================

function s = list_sessions(base, animal)
% Session folders are <base>/results_<animal>/<animal>_<date>_attentional-sampling_*
% Same contract as the local of the same name in erp_latency_wave.m: FULL PATHS.
d = dir(fullfile(base, ['results_' animal], [animal '_*attentional-sampling*']));
d = d([d.isdir]);
s = arrayfun(@(x) fullfile(x.folder, x.name), d, 'uni', 0);
end

function n = shortname(p)
[~, n] = fileparts(p);
end

function [X, fs, nUsed, nSkip, tinfo] = load_envelopes(base, animal, sessDir, ...
        SIGNAL, GAMMA_HZ, BB_HZ, SMOOTH_MS, WIN, nChTot, MIN_CH) %#ok<INUSL>
% -> X : nTrials cell-free 3-D array [nCh x nSamp x nTrials], z-scored per
%        channel within trial, NaN rows where the channel is absent.
X = []; fs = NaN; nUsed = 0; nSkip = 0; tinfo = [];

switch SIGNAL
    case 'mua'
        fp = fullfile(sessDir, 'clean_mua.mat');
        if ~isfile(fp), return, end
        S = load(fp, 'clean_mua'); D = S.clean_mua;
    otherwise
        fp = fullfile(sessDir, 'clean_lfp.mat');
        if ~isfile(fp), return, end
        S = load(fp, 'clean_data'); D = S.clean_data;
end
fs = D.fsample;
[rows, slots, ~] = v4_channel_slots(D.label, animal, nChTot);
if isempty(rows), return, end
tinfo = D.trialinfo;

% Filter design once per session. Zero-phase: a causal filter's group delay
% IS the quantity being measured (see assumption 3 in the header).
switch SIGNAL
    case 'gamma_env', bnd = GAMMA_HZ;
    case 'bb_env',    bnd = BB_HZ;
    case 'mua',       bnd = [];
    otherwise, error('Unknown SIGNAL: %s', SIGNAL);
end
if ~isempty(bnd)
    [bb, aa] = butter(3, bnd/(fs/2), 'bandpass');
end
nsm = max(1, round(SMOOTH_MS/1000*fs));
ker = ones(1, nsm)/nsm;

tvec  = WIN(1) : 1/fs : WIN(2);
nSamp = numel(tvec);
buf   = nan(nChTot, nSamp, numel(D.trial));

for it = 1:numel(D.trial)
    tt = D.time{it};
    if tt(1) > WIN(1) || tt(end) < WIN(2), nSkip = nSkip + 1; continue, end
    raw = D.trial{it}(rows,:);
    % filter on the FULL trial, then cut -- filtering the cut window would put
    % edge transients inside the analysis window
    if ~isempty(bnd)
        for r = 1:size(raw,1)
            if all(isfinite(raw(r,:))), raw(r,:) = filtfilt(bb, aa, raw(r,:)); end
        end
        env = abs(hilbert(raw.')).';
    else
        env = raw;                     % MUA is already an activity measure
    end
    env = conv2(env, ker, 'same');
    seg = interp1(tt, env.', tvec, 'linear').';       % nRows x nSamp
    seg = detrend(seg.').';                           % remove per-trial ramp
    sd  = std(seg, 0, 2);
    seg = (seg - mean(seg,2)) ./ max(sd, eps);        % z-score within trial
    good = all(isfinite(seg),2) & sd > 0;
    if sum(good) < MIN_CH, nSkip = nSkip + 1; continue, end
    nUsed = nUsed + 1;
    buf(slots(good), :, nUsed) = seg(good,:);
end
X = buf(:,:,1:nUsed);
end

function ACC = init_acc(nSamp, nChTot, doPair)
nfft = 2^nextpow2(2*nSamp);
ACC = struct('nfft',nfft,'nSamp',nSamp,'nChTot',nChTot, ...
    'CSref', zeros(nChTot, nfft), ...         % channel x reference cross-spectrum
    'nPerCh', zeros(nChTot,1), 'nTrials', 0, ...
    'canShuffleTrials', true, 'Xstore', {{}}, 'fs', NaN);
if doPair
    ACC.CS = zeros(nChTot, nChTot, nfft);
end
end

function ACC = accumulate_xspec(ACC, X, doPair)
% X : nCh x nSamp x nTrials. Accumulates, over trials, the cross-spectrum of
% every channel against the LEAVE-ONE-OUT mean of the other channels, so a
% channel never correlates with itself through the reference. Summing spectra
% and inverse-transforming once at the end averages the CORRELOGRAM before the
% peak is picked -- picking per trial and averaging peaks averages noise.
nfft = ACC.nfft;
for it = 1:size(X,3)
    x = X(:,:,it);
    good = all(isfinite(x),2);
    if sum(good) < 2, continue, end
    g  = find(good);
    Xf = fft(x(g,:), nfft, 2);                    % nGood x nfft
    tot = sum(Xf, 1);
    ref = (tot - Xf) / max(numel(g)-1, 1);        % leave-one-out reference
    ACC.CSref(g,:)  = ACC.CSref(g,:) + Xf .* conj(ref);
    ACC.nPerCh(g)   = ACC.nPerCh(g) + 1;
    ACC.nTrials     = ACC.nTrials + 1;
    if doPair
        Xp = complex(zeros(ACC.nChTot, 1, nfft));
        Xp(g,1,:) = reshape(Xf, numel(g), 1, nfft);
        ACC.CS = ACC.CS + pagemtimes(Xp, 'none', Xp, 'ctranspose');
    end
    % Keep a bounded reservoir of trials so the trial-shuffle null can be
    % re-estimated without re-reading the data. Bounded, not complete: the
    % null only needs enough trials to break the within-trial pairing.
    if numel(ACC.Xstore) < 400
        ACC.Xstore{end+1} = single(x);
    end
end
end

function tau = lags_from_acc(ACC, maxlagS, shuffleTrials, seed)
% Peak lag per channel, in SECONDS, referenced to the array mean.
if nargin < 3, shuffleTrials = false; end
if shuffleTrials
    CS = shuffled_csref(ACC, seed);
else
    CS = ACC.CSref;
end
xc  = real(ifft(CS, [], 2));                       % nCh x nfft, circular
xc  = fftshift(xc, 2);
mid = ACC.nfft/2 + 1;
sel = (mid-maxlagS) : (mid+maxlagS);
xc  = xc(:, sel);
nCh = size(xc,1);
tau = nan(nCh,1);
for c = 1:nCh
    v = xc(c,:);
    if ~any(isfinite(v)) || all(v == 0), continue, end
    [~, im] = max(v);
    if im > 1 && im < numel(v)
        % parabolic refinement -- 0.1-0.8 m/s across 2.8 mm is only 3.5-28 ms
        % end to end, so sub-sample resolution is not a luxury here
        d = (v(im-1) - v(im+1)) / (2*(v(im-1) - 2*v(im) + v(im+1)));
    else
        d = 0;
    end
    if ~isfinite(d), d = 0; end
    tau(c) = ((im - 1 - maxlagS) + d) / ACC.fs;
end
tau = tau - mean(tau, 'omitnan');                  % only the SPREAD is meaningful
end

function CS = shuffled_csref(ACC, seed)
% Trial-shuffle null: correlate each channel against a reference built from a
% DIFFERENT trial. Anything stimulus-locked, any common drift and any fixed
% filter delay survives this; genuine within-trial co-fluctuation does not.
rng(9000 + seed);
n = numel(ACC.Xstore);
CS = zeros(ACC.nChTot, ACC.nfft);
if n < 2, return, end
ord = randperm(n);
for it = 1:n
    x = double(ACC.Xstore{it});
    y = double(ACC.Xstore{ord(it)});
    if it == ord(it), continue, end
    g = find(all(isfinite(x),2) & all(isfinite(y),2));
    if numel(g) < 2, continue, end
    Xf = fft(x(g,:), ACC.nfft, 2);
    Yf = fft(y(g,:), ACC.nfft, 2);
    tot = sum(Yf,1);
    ref = (tot - Yf) / max(numel(g)-1, 1);
    CS(g,:) = CS(g,:) + Xf .* conj(ref);
end
end

function [sl, R2, dirRad, speed_cms] = fit_lag_plane(tau, XY, W)
% Least-squares plane in LAG space: tau_c ~ s . r_c. s is the slowness vector
% in s/mm; speed = 1/|s|, direction = angle(s).
use = isfinite(tau) & all(isfinite(XY),2) & W > 0;
sl = [0;0]; R2 = 0; dirRad = NaN; speed_cms = NaN;
if sum(use) < 6, return, end
Xd = [XY(use,:), ones(sum(use),1)];
w  = W(use) / max(sum(W(use)), eps);
Aw = Xd .* sqrt(w);  bw = tau(use) .* sqrt(w);
c  = Aw \ bw;
sl = c(1:2);
pred = Xd * c;
ss_res = sum(w .* (tau(use) - pred).^2);
ss_tot = sum(w .* (tau(use) - sum(w.*tau(use))).^2);
R2 = 1 - ss_res / max(ss_tot, eps);
nsl = norm(sl);
if nsl > 0
    dirRad    = atan2(sl(2), sl(1));
    speed_cms = (1/nsl) / 10;          % (mm/s) -> cm/s
end
end

function [speed_cms, dirRad, R2] = fit_lag_pairs(ACC, XY, maxlagS)
% All-pairs version: tau_cd fitted against (r_d - r_c). More robust than the
% mean reference but needs ACC.CS, which is O(nChan^2 * nfft) to accumulate.
speed_cms = NaN; dirRad = NaN; R2 = NaN;
if ~isfield(ACC,'CS'), return, end
n = ACC.nChTot; mid = ACC.nfft/2 + 1;
dd = []; tt = [];
for c = 1:n-1
    Xc = real(ifft(squeeze(ACC.CS(c, c+1:n, :)), [], 2));
    Xc = fftshift(Xc, 2);
    Xc = Xc(:, (mid-maxlagS):(mid+maxlagS));
    for j = 1:size(Xc,1)
        d2 = c + j;
        v = Xc(j,:);
        if all(v == 0) || ~any(isfinite(v)), continue, end
        [~, im] = max(v);
        tt(end+1,1) = (im - 1 - maxlagS)/ACC.fs; %#ok<AGROW>
        dd(end+1,:) = XY(d2,:) - XY(c,:);        %#ok<AGROW>
    end
end
if numel(tt) < 20, return, end
s = [dd, ones(numel(tt),1)] \ tt;
pred = [dd, ones(numel(tt),1)] * s;
R2 = 1 - sum((tt-pred).^2) / max(sum((tt-mean(tt)).^2), eps);
if norm(s(1:2)) > 0
    dirRad = atan2(s(2), s(1));
    speed_cms = (1/norm(s(1:2)))/10;
end
end

function v = wave_verdict(R2, p_pos, p_trl, okSpeed, alpha) %#ok<INUSL>
if p_pos >= alpha
    v = 'NO propagation: lag field is no more planar than a position shuffle';
elseif isfinite(p_trl) && p_trl >= alpha
    v = 'NOT within-trial: survives the trial shuffle, so it is common/stimulus-locked structure';
elseif ~okSpeed
    v = 'planar lag gradient, but speed OUTSIDE the mesoscopic band -- not a first-order wave';
else
    v = 'PROPAGATION: planar lag gradient at a mesoscopic speed';
end
end

function th = pgd_direction(base, animal, dv, band)
% Circular-mean propagation direction over `band`, from the saved PGD results.
% Same reader as the local of this name in erp_latency_wave.m.
th = NaN;
fp = fullfile(base,'results_combined','scanning','planar_wave_existence', ...
    'cp10_till_100', dv, 'planar_wave_existence.mat');
if ~isfile(fp), return; end
S = load(fp,'results'); A = S.results.A;
ia = find(strcmp({A.animal}, animal), 1);
if isempty(ia), return; end
sel = A(ia).freq >= band(1) & A(ia).freq <= band(2) & A(ia).sig(:);
if ~any(sel), sel = A(ia).freq >= band(1) & A(ia).freq <= band(2); end
th = angle(mean(exp(1i*A(ia).DIR(sel)), 'omitnan'));
end

function d = angdiff_(a, b)
d = angle(exp(1i*(a-b)));
end

function s = esc(s), s = strrep(s, '_', '\_'); end

function out = ternary(c,a,b), if c, out=a; else, out=b; end, end

function f = new_fig(w, h)
ss = get(0,'ScreenSize');
w = min(w, ss(3)-80); h = min(h, ss(4)-120);
f = figure('Units','pixels','Position',[40 40 w h],'Color','w');
end

function save_pdf(fig, fname)
drawnow;
set(fig,'Color','w','InvertHardcopy','off');
set(fig,'Units','inches'); p = get(fig,'Position');
set(fig,'PaperUnits','inches','PaperSize',[p(3) p(4)], ...
        'PaperPosition',[0 0 p(3) p(4)],'PaperPositionMode','manual');
set(fig,'Units','pixels');
try
    exportgraphics(fig, fname, 'ContentType','vector','BackgroundColor','white');
catch
    print(fig, fname, '-dpdf', '-painters', '-r300');
end
end
