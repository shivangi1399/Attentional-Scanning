% =====================================================================
% ERP LATENCY TEST — does the phase gradient correspond to a real TIME delay?
%
% WHY THIS SCRIPT EXISTS
% ----------------------
% Every test in this folder so far measures PHASE. Phase alone cannot separate
% a propagating signal from a fixed spatial phase offset, because at a SINGLE
% frequency the two are the same object:
%
%     signal(x,t) = cos(2*pi*f*t - k*x)     is a wave at v = 2*pi*f/k
%                                            for ANY k, real or artefactual.
%
% The separation is a cross-frequency one:
%
%     real propagation  -> fixed TIME delay tau = d/v
%                       -> phase lag = 2*pi*f*tau  GROWS with f  ->  k ∝ f
%     fixed offset      -> fixed PHASE lag, same at every f      ->  k const
%
% This script measures tau DIRECTLY, in the time domain, from the
% stimulus-locked ERP. No phase wrapping, no frequency ambiguity, no k-vs-v
% confound. It then asks the decisive question:
%
%     does the measured tau PREDICT the observed phase gradient?
%         predicted d(phi)/d(distance) at frequency f  =  2*pi*f * (dtau/ddistance)
%     MATCH    -> the phase tilt IS conduction delay -> real traveling wave
%     NO MATCH -> phase tilt exists with no corresponding time delay
%                 -> not propagation
%
% TWO AXES, TWO DIFFERENT QUESTIONS (do not conflate them)
% -------------------------------------------------------
%   AXIS_MODE = 'cortical'  latency vs position along the array (mm), projected
%                           on the propagation direction theta that
%                           cortical_planar_wave_PGD.m already fitted.
%                           -> tests CORTICAL propagation (the §2/§4 question)
%   AXIS_MODE = 'retinotopic'  latency vs RF eccentricity (deg), foveal ->
%                           peripheral, using elec_rf_deg (same frame as
%                           stimulus_loc_traveling_wave.m).
%                           -> tests the SCANNING/retinotopic question (§5)
% Both are run by default.
%
% THE LEAD THIS IS AIMED AT
% -------------------------
% The slope test on planar_wave_existence.mat (see README) found ONE band that
% behaves like a real wave: hermes 13.3-25.5 Hz, v ~ 28 cm/s, slope of
% log(v) vs log(f) = -0.10, k rising in proportion. It does not replicate in
% klecks and de-rotation did not confirm it. If that band is real, THIS test
% should return a cortical latency slope near 1/(28 cm/s) = 3.6 ms/mm along the
% same direction. If it returns ~0 ms/mm, the beta band is a fixed offset too.
%
% CAVEAT
% -----------------------
% ERP latency measures propagation of the EVOKED response. The phase analyses
% measure ONGOING oscillation. These need not be the same process. A MATCH is
% strong evidence for propagation; a NULL does not by itself refute a wave in
% ongoing activity. Report it that way.
%
% DATA
% ----
%   <base>/results_<animal>/<session>/clean_lfp.mat  -> variable `clean_data`,
%       a FieldTrip raw structure:
%           .trial      {1 x nTrials} each [nChan x nTime]
%           .time       {1 x nTrials} each [1 x nTime], 0 = stimulus event
%           .label      {nChan x 1}
%           .trialinfo  [nTrials x nCols]  col 16/17 = target x/y (fix-centred px)
%           .fsample    Hz
%   V4 channels are mapped onto the canonical 1..64 slots BY LABEL NUMBER
%   ('V4-n' -> slot n for hermes, n-64 for klecks), so channel i here is
%   channel i in phase_progression.mat. The offsets are verified against that
%   file's empty channels — see v4_channel_slots.
%
%   NOT by sorted position: sessions are missing different subsets of the 64
%   V4 channels (hermes 60/64, klecks 62/64, and WHICH ones varies), so a
%   positional sort would put a different electrode in row i in different
%   sessions. Absent channels are left NaN for that session and simply
%   contribute no trials, exactly as Phase_combine_sessions.m does when it
%   builds the 64-channel phase array. See v4_channel_slots.
%
% Output:
%   Plots/scanning/erp_latency/cp10_till_100/<dv>/erp_latency.pdf
%   results_combined/scanning/erp_latency/cp10_till_100/<dv>/erp_latency.mat
% =====================================================================

clearvars; close all; clc
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/scanning/functions

%% ─── Settings ────────────────────────────────────────────────────────
animals    = {'hermes','klecks'};
dv         = 'lfp';
base       = '/mnt/hpc/projects/MWSampling/4Shivangi';
grid_rows  = 8; grid_cols = 8;
SPACING_MM = 0.4;                  % electrode pitch (mm)
nPerm      = 1000;
alpha      = 0.05;
rng(2025);

% ERP window, relative to the event at t = 0 in clean_data.time
ERP_WIN    = [-0.05 0.30];         % s
BASE_WIN   = [-0.05 0.00];         % s, for baseline subtraction
LP_HZ      = 100;                  % low-pass the ERP before latency estimation

% Latency estimator:
%   'xcorr'    lag of peak cross-correlation against the array-mean ERP.
%              Most robust — uses the whole waveform, not one landmark.
%   'halfpeak' time at which |ERP| first reaches 50% of its peak
%   'peak'     time of the ERP extremum
LAT_METHOD = 'xcorr';
MAX_LAG_MS = 30;                   % search window for the lag (ms)
MIN_TRIALS = 20;                   % per channel, to accept an ERP

AXIS_MODES = {'cortical','retinotopic'};

% Propagation direction per animal, in ARRAY coordinates (radians), taken from
% cortical_planar_wave_PGD.m. Left empty = read from planar_wave_existence.mat
% for the band in PGD_BAND below. Set explicitly to override.
PGD_BAND   = [13.3 25.5];          % Hz — the hermes wave-like band (see header)
THETA_FIX  = struct('hermes', [], 'klecks', []);

% Retinotopic axis calibration (only for AXIS_MODE = 'retinotopic')
PIX_PER_DEG = struct('hermes', 53.24, 'klecks', 50.56);
RF_DATE     = '20170829';
SCREEN_XY   = [1680 1050];
RF_VALID_ONLY = true;              % Valid_Gaussian RF centres only

out_dir = fullfile(base,'Plots','scanning','erp_latency','cp10_till_100', dv);
res_dir = fullfile(base,'results_combined','scanning','erp_latency','cp10_till_100', dv);
if ~exist(out_dir,'dir'), mkdir(out_dir); end
if ~exist(res_dir,'dir'), mkdir(res_dir); end

% Array coordinates — SAME index<->position map as the PGD and de-rotation
% scripts. Do not change one without the others.
nChTot = grid_rows*grid_cols;
ch_col = ceil((1:nChTot)' / grid_rows);
ch_row = grid_rows - mod((1:nChTot)' - 1, grid_rows);
XY     = [ch_col, ch_row] * SPACING_MM;         % nCh x 2, mm

% =====================================================================
% ASSUMPTIONS THAT MUST BE CHECKED BEFORE TRUSTING ANY OUTPUT
% =====================================================================
% (1) clean_data.time is zeroed on the STIMULUS event that phase_progression.m
%     also used. If it is zeroed on trial start or on the response instead,
%     every latency here is offset by a constant — harmless for the SLOPE
%     (which is all that matters) but the absolute numbers are meaningless.
%     -> the script prints the time axis of the first session; eyeball it.
% (2) Channel slot matches phase_progression.mat 1:1. VERIFIED, not assumed:
%     the slots absent from every session here reproduce that file's
%     zero-trial channels exactly, in both animals (see v4_channel_slots).
%     Expected coverage, for comparison with what the script prints:
%
%        hermes   31 sessions   57-61 of 64 channels each
%                 never present: 48, 64      (zero trials in phase_progression)
%                 sometimes missing: 22, 37, 43, 49, 55, 56, 63
%        klecks   25 sessions   61-62 of 64 channels each
%                 never present: 45          (zero trials in phase_progression)
%                 slot 64 exists in ONE session only (20170906) -> 68 trials
%
%     A DIFFERENT set here means the label convention has changed upstream;
%     stop and re-derive it rather than reinterpreting the numbers.
%     -> sanity check: the ERP should look like an ERP on most channels, and
%        the channel-mean ERP should have a clear evoked deflection.
% (3) Sessions are pooled by concatenating trials. If sessions differ in
%     electrode drift this adds noise but no bias in the slope.
% =====================================================================

L = struct();

%% ─── Per animal: build ERPs and estimate per-channel latency ─────────
for ia = 1:numel(animals)
    animalName = animals{ia};
    sess = list_sessions(base, animalName);
    if isempty(sess)
        warning('No sessions found for %s — skipping.', animalName); continue
    end
    fprintf('\n=== %s : %d sessions ===\n', animalName, numel(sess));

    ERPsum = []; ERPn = zeros(nChTot,1); tvec = [];
    nSessUsed = 0; slotSeen = false(nChTot,1); covLog = [];
    for is = 1:numel(sess)
        fp = fullfile(sess{is}, 'clean_lfp.mat');
        if ~isfile(fp), continue; end
        S = load(fp, 'clean_data');  D = S.clean_data;

        % Sessions are missing different subsets of the 64 V4 channels, so
        % map each present channel onto its CANONICAL slot rather than
        % assuming the session holds all 64 in order. See v4_channel_slots.
        [rows, slots, nBad] = v4_channel_slots(D.label, animalName, nChTot);
        if isempty(rows)
            warning('  %s: no V4 channels — session skipped.', sess{is});
            continue
        end
        if nBad > 0
            warning('  %s: %d V4 label(s) outside the canonical 1..%d range — ignored.', ...
                sess{is}, nBad, nChTot);
        end
        nSessUsed = nSessUsed + 1;
        slotSeen(slots) = true;
        covLog(end+1,:) = [is numel(slots)]; %#ok<AGROW>

        if isempty(tvec)
            fs   = D.fsample;
            tvec = ERP_WIN(1) : 1/fs : ERP_WIN(2);
            ERPsum = zeros(nChTot, numel(tvec));
            fprintf('  fs=%g Hz | session time axis: %.3f .. %.3f s (0 = event)\n', ...
                fs, D.time{1}(1), D.time{1}(end));
        end
        for it = 1:numel(D.trial)
            tt = D.time{it};
            if tt(1) > ERP_WIN(1) || tt(end) < ERP_WIN(2), continue; end
            % interpolate only the channels this session has, then scatter
            % them into their canonical slots (absent slots stay NaN)
            seg = nan(nChTot, numel(tvec));
            seg(slots,:) = interp1(tt, D.trial{it}(rows,:).', tvec, 'linear').';
            bb  = tvec >= BASE_WIN(1) & tvec <= BASE_WIN(2);
            seg = seg - mean(seg(:,bb), 2, 'omitnan');                     % baseline
            good = all(isfinite(seg), 2);
            ERPsum(good,:) = ERPsum(good,:) + seg(good,:);
            ERPn(good)     = ERPn(good) + 1;
        end
    end
    if isempty(ERPsum), warning('  no usable trials for %s', animalName); continue; end

    % Channel coverage — READ THIS before trusting the map. A slot that is
    % absent from every session can never produce a latency, and a slot with
    % few sessions is noisier than its neighbours.
    fprintf('  sessions used: %d/%d | V4 channels per session: %d–%d of %d\n', ...
        nSessUsed, numel(sess), min(covLog(:,2)), max(covLog(:,2)), nChTot);
    if ~all(slotSeen)
        fprintf('  canonical slots never present in ANY session: %s\n', ...
            mat2str(find(~slotSeen).'));
    end

    ERP = ERPsum ./ max(ERPn,1);            % nCh x nT, channel-wise ERP
    ERP(ERPn < MIN_TRIALS, :) = NaN;
    ERP = lowpass_erp(ERP, 1/mean(diff(tvec)), LP_HZ);
    fprintf('  channels with >=%d trials: %d/%d\n', MIN_TRIALS, sum(ERPn>=MIN_TRIALS), nChTot);

    % ── per-channel latency (s), relative to the array-mean ERP ────────
    lat = erp_latency(ERP, tvec, LAT_METHOD, MAX_LAG_MS/1000);
    fprintf('  latency spread across channels: %.2f ms (sd), range %.2f ms\n', ...
        1e3*std(lat,'omitnan'), 1e3*(max(lat)-min(lat)));

    L(ia).animal = animalName; L(ia).ERP = ERP; L(ia).t = tvec;
    L(ia).lat = lat; L(ia).nTrials = ERPn;

    %% ── Regress latency on each distance axis ──────────────────────────
    for am = 1:numel(AXIS_MODES)
        mode = AXIS_MODES{am};
        switch mode
            case 'cortical'
                th = THETA_FIX.(animalName);
                if isempty(th), th = pgd_direction(base, animalName, dv, PGD_BAND); end
                if isnan(th)
                    warning('  no PGD direction for %s — cortical axis skipped.', animalName);
                    continue
                end
                dist = XY(:,1)*cos(th) + XY(:,2)*sin(th);      % mm along propagation
                unit = 'mm';  unit_v = 'cm/s';
            case 'retinotopic'
                ppd = PIX_PER_DEG.(animalName);
                [rf_deg, rf_valid] = elec_rf_deg(base, animalName, nChTot, RF_DATE, ppd, SCREEN_XY);
                dist = hypot(rf_deg(:,1), rf_deg(:,2));         % RF eccentricity, deg
                if RF_VALID_ONLY, dist(~rf_valid) = NaN; end
                unit = 'deg'; unit_v = 'deg/s';
        end

        use = isfinite(lat) & isfinite(dist);
        if sum(use) < 8
            warning('  %s/%s: only %d usable channels — skipped.', animalName, mode, sum(use));
            continue
        end
        % slope in s per unit distance; speed = 1/slope
        cf    = polyfit(dist(use), lat(use), 1);         % [slope intercept]
        slope = cf(1);                                   % s per unit
        speed = 1/slope;                                 % unit per s
        if strcmp(mode,'cortical'), speed = speed/10; end   % mm/s -> cm/s

        % Null: shuffle latency across channels (breaks the spatial relation,
        % preserves the latency distribution), take |slope| as the statistic.
        nullslope = nan(nPerm,1);
        idx = find(use);
        for ib = 1:nPerm
            lp = lat; lp(idx) = lat(idx(randperm(numel(idx))));
            cfp = polyfit(dist(use), lp(use), 1);
            nullslope(ib) = cfp(1);
        end
        pval = mean(abs(nullslope) >= abs(slope));
        r    = corr(dist(use), lat(use));

        fprintf('  [%-11s] slope = %+.3f ms/%s  -> speed %+.1f %s | r = %+.2f | p = %.3f  %s\n', ...
            mode, 1e3*slope, unit, speed, unit_v, r, pval, ternary(pval<alpha,'*',''));

        L(ia).(mode) = struct('dist',dist,'slope',slope,'coef',cf,'speed',speed, ...
            'r',r,'p',pval,'nullslope',nullslope,'unit',unit,'unit_v',unit_v,'use',use);
    end
end
valid = find(arrayfun(@(s) ~isempty(s.animal), L));

%% ─── THE DECISIVE COMPARISON: does tau predict the phase gradient? ───
% A measured latency gradient dtau/dd (s per mm) predicts a phase gradient
%       k_pred(f) = 2*pi*f * dtau/dd        rad/mm
% that GROWS LINEARLY WITH FREQUENCY. Compare against k_corr(f) that
% cortical_planar_wave_PGD.m measured. If the observed k is flat while k_pred
% rises, the phase tilt is not a time delay.
fprintf('\n================ PHASE-GRADIENT PREDICTION ================\n');
for ia = valid
    if ~isfield(L(ia),'cortical') || isempty(L(ia).cortical), continue; end
    [fHz, k_obs, sig] = pgd_wavenumber(base, L(ia).animal, dv);
    if isempty(fHz), continue; end
    k_pred = 2*pi*fHz * L(ia).cortical.slope;          % rad/mm, from ERP latency
    L(ia).k_obs = k_obs; L(ia).k_pred = k_pred; L(ia).fHz = fHz; L(ia).pgd_sig = sig;
    ok = isfinite(k_obs) & isfinite(k_pred) & sig;
    if sum(ok) >= 3
        rr = corr(k_pred(ok), k_obs(ok));
        fprintf('%-8s: over PGD-significant freqs, corr(k_pred, k_obs) = %+.2f  (n=%d)\n', ...
            L(ia).animal, rr, sum(ok));
        fprintf('          k_obs  %.3f-%.3f rad/mm   |  k_pred %.3f-%.3f rad/mm\n', ...
            min(k_obs(ok)), max(k_obs(ok)), min(k_pred(ok)), max(k_pred(ok)));
        fprintf('          -> high corr AND similar magnitude = the tilt IS a time delay\n');
        fprintf('          -> flat k_obs while k_pred rises   = fixed phase offset\n');
    end
end

%% ─── Figure ──────────────────────────────────────────────────────────
esc = @(s) strrep(s,'_','\_');
nA  = numel(valid);
f1 = figure('Visible','off','Position',[40 40 420*3 320*max(nA,1)]);
for k = 1:nA
    ia = valid(k);
    % col 1: the ERPs themselves (sanity check — these must look like ERPs)
    ax = subplot(nA,3,(k-1)*3+1); hold(ax,'on');
    plot(ax, L(ia).t*1e3, L(ia).ERP.', 'Color',[.7 .7 .7 .5]);
    plot(ax, L(ia).t*1e3, mean(L(ia).ERP,1,'omitnan'), 'k','LineWidth',1.5);
    xline(ax,0,'r:'); xlabel(ax,'time (ms)'); ylabel(ax,'ERP (z)');
    title(ax,{esc(L(ia).animal),'channel ERPs (black = mean)'},'FontSize',9);

    % col 2/3: latency vs each distance axis
    for am = 1:numel(AXIS_MODES)
        mode = AXIS_MODES{am};
        if ~isfield(L(ia),mode) || isempty(L(ia).(mode)), continue; end
        Q = L(ia).(mode);
        ax = subplot(nA,3,(k-1)*3+1+am); hold(ax,'on');
        plot(ax, Q.dist(Q.use), 1e3*L(ia).lat(Q.use), 'o','MarkerFaceColor',[.2 .4 .8], ...
            'MarkerEdgeColor','none');
        xx = linspace(min(Q.dist(Q.use)), max(Q.dist(Q.use)), 10);
        plot(ax, xx, 1e3*polyval(Q.coef, xx), 'k-','LineWidth',1.4);
        xlabel(ax,sprintf('distance (%s)',Q.unit)); ylabel(ax,'ERP latency (ms)');
        title(ax,{sprintf('%s — %s', esc(L(ia).animal), esc(mode)), ...
            sprintf('%+.2f ms/%s | v=%+.1f %s | p=%.3f', 1e3*Q.slope, Q.unit, Q.speed, Q.unit_v, Q.p)}, ...
            'FontSize',9);
        grid(ax,'on');
    end
end
sgtitle({'ERP LATENCY TEST: does the phase gradient correspond to a real TIME delay?', ...
         'cortical axis = propagation along the array (tests cortical wave)  |  retinotopic axis = RF eccentricity (tests the scanning hypothesis)', ...
         'a real wave gives a NON-ZERO latency slope whose 1/slope matches the phase-derived speed'}, ...
         'FontSize',9);
set(f1,'PaperPositionMode','auto'); pos=get(f1,'Position');
set(f1,'PaperUnits','points','PaperSize',pos(3:4));
saveas(f1, fullfile(out_dir,'erp_latency.pdf'));

results = struct('L',L,'animals',{animals},'dv',dv,'ERP_WIN',ERP_WIN, ...
    'LAT_METHOD',LAT_METHOD,'SPACING_MM',SPACING_MM,'nPerm',nPerm,'alpha',alpha, ...
    'PGD_BAND',PGD_BAND,'RF_VALID_ONLY',RF_VALID_ONLY);
save(fullfile(res_dir,'erp_latency.mat'),'results','-v7.3');
fprintf('\nSaved under %s\n', out_dir);

%% =====================================================================
%% Helpers
%% =====================================================================
function s = list_sessions(base, animal)
% Session folders are <base>/results_<animal>/<animal>_<date>_attentional-sampling_*
d = dir(fullfile(base, ['results_' animal], [animal '_*attentional-sampling*']));
d = d([d.isdir]);
s = arrayfun(@(x) fullfile(x.folder, x.name), d, 'uni', 0);
end

function [rows, slots, nBad] = v4_channel_slots(labels, animal, nChTot)
% Map this session's V4 channels onto the CANONICAL 1..nChTot slots.
%
% WHY NOT A PLAIN SORT. Sessions do not all contain all 64 V4 channels —
% cleaning drops a few, and WHICH ones differs per session (hermes: 60/64,
% missing e.g. 43 in one session and 37 in another; klecks: 62/64). Sorting
% the present labels and using position would therefore put a different
% physical electrode in row i in different sessions, and silently average
% channel 43 onto channel 44 when pooling. Map by LABEL NUMBER instead.
%
% This mirrors how the phase pipeline builds its 64-channel array
% (Phase_analysis/masters_code/Phase_combine_sessions.m): chan_orig = 1:64,
% each present channel written to its own slot, absent channels left NaN —
% never compacted. Slot i here is channel i in phase_progression.mat.
%
%   hermes   label V4-n     -> slot n         (labels run   1..64,  blank ->  64)
%   klecks   label V4-n     -> slot n - 64    (labels run  65..128, blank -> 128)
%
% THE KLECKS OFFSET IS -64, NOT -63. clean_lfp.mat already holds the
% INCREMENTED numbering (65..128), not mapping_lfp.m's pre-increment 64..127.
% Verified against phase_progression.mat, which is what these slots must line
% up with — three independent facts, all exact:
%
%   hermes  channels with zero trials              = {48, 64}
%           slots absent from every session here   = {48, 64}      match
%   klecks  V4-109 is absent from all 25 sessions -> slot 45
%           the ONLY klecks channel with zero trials = 45          match
%   klecks  V4-128 appears in one session (20170906) -> slot 64
%           channel 64 has only 68 trials, vs ~1400 elsewhere      match
%
% With -63 every klecks channel would sit one slot off, silently rotating the
% array by one electrode. Re-run that check if the upstream pipeline changes.
%
%   rows   row indices into D.label / D.trial for the channels we keep
%   slots  canonical slot each of those rows maps to (same length as rows)
%   nBad   V4 labels whose number fell outside 1..nChTot (should be 0)
labels = labels(:);
isV4 = startsWith(strtrim(labels), 'V4-');
V4i  = find(isV4);
nums = str2double(erase(strtrim(labels(isV4)), 'V4-'));   % ' 1' -> 1, '' -> NaN
switch lower(animal)
    case 'hermes'
        nums(isnan(nums)) = 64;
        slots = nums;
    case 'klecks'
        nums(isnan(nums)) = 128;
        slots = nums - 64;
    otherwise
        error('Unknown animal: %s', animal);
end
keep  = slots >= 1 & slots <= nChTot & isfinite(slots);
nBad  = sum(~keep);
rows  = V4i(keep);
slots = slots(keep);
if numel(unique(slots)) ~= numel(slots)
    error('Duplicate canonical slots in %s — label convention is wrong.', animal);
end
end

function E = lowpass_erp(E, fs, fcut)
% Zero-phase low-pass. Zero-phase matters: a causal filter would introduce a
% group delay, which is exactly the quantity being measured.
if isempty(fcut) || fcut >= fs/2, return; end
n = 4; [b,a] = butter(n, fcut/(fs/2), 'low');
for c = 1:size(E,1)
    if all(isfinite(E(c,:))), E(c,:) = filtfilt(b,a,E(c,:)); end
end
end

function lat = erp_latency(E, t, method, maxlag)
% Per-channel latency in SECONDS, referenced to the array-mean ERP (so the
% mean is 0 by construction and only the SPREAD across channels is meaningful
% — which is all the regression uses).
nCh = size(E,1); lat = nan(nCh,1);
ref = mean(E, 1, 'omitnan');
dt  = mean(diff(t));
switch method
    case 'xcorr'
        mx = round(maxlag/dt);
        for c = 1:nCh
            if ~all(isfinite(E(c,:))), continue; end
            [xc, lags] = xcorr(E(c,:) - mean(E(c,:)), ref - mean(ref), mx, 'coeff');
            [~, im] = max(xc);
            lat(c) = lags(im) * dt;
        end
    case 'halfpeak'
        for c = 1:nCh
            if ~all(isfinite(E(c,:))), continue; end
            [pk, ip] = max(abs(E(c,:)));
            j = find(abs(E(c,1:ip)) >= 0.5*pk, 1, 'first');
            if ~isempty(j), lat(c) = t(j); end
        end
        lat = lat - mean(lat,'omitnan');
    case 'peak'
        for c = 1:nCh
            if ~all(isfinite(E(c,:))), continue; end
            [~, ip] = max(abs(E(c,:))); lat(c) = t(ip);
        end
        lat = lat - mean(lat,'omitnan');
    otherwise
        error('Unknown LAT_METHOD: %s', method);
end
end

function th = pgd_direction(base, animal, dv, band)
% Circular-mean propagation direction over `band`, from the saved PGD results.
th = NaN;
fp = fullfile(base,'results_combined','scanning','planar_wave_existence', ...
    'cp10_till_100', dv, 'planar_wave_existence.mat');
if ~isfile(fp), warning('planar_wave_existence.mat not found'); return; end
S = load(fp,'results'); A = S.results.A;
ia = find(strcmp({A.animal}, animal), 1);
if isempty(ia), return; end
sel = A(ia).freq >= band(1) & A(ia).freq <= band(2) & A(ia).sig(:);
if ~any(sel), sel = A(ia).freq >= band(1) & A(ia).freq <= band(2); end
th = angle(mean(exp(1i*A(ia).DIR(sel)), 'omitnan'));
fprintf('  PGD direction over %.1f-%.1f Hz = %.0f deg (array coords)\n', ...
    band(1), band(2), rad2deg(th));
end

function [fHz, k_obs, sig] = pgd_wavenumber(base, animal, dv)
% Per-frequency measured wavenumber (rad/mm) + the cluster-significant mask.
fHz = []; k_obs = []; sig = [];
fp = fullfile(base,'results_combined','scanning','planar_wave_existence', ...
    'cp10_till_100', dv, 'planar_wave_existence.mat');
if ~isfile(fp), return; end
S = load(fp,'results'); A = S.results.A;
ia = find(strcmp({A.animal}, animal), 1);
if isempty(ia), return; end
fHz = A(ia).freq(:); k_obs = A(ia).k_corr(:); sig = logical(A(ia).sig(:));
end

function out = ternary(c,a,b), if c, out=a; else, out=b; end, end
