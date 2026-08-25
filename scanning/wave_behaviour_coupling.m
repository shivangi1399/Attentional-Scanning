% =====================================================================
% SINGLE-TRIAL WAVES AND DETECTION -- do the POSITION and TIMING of a wave
% before the target predict whether the target is seen?
%
% WHY THIS EXISTS
%   Every wave analysis in this folder so far asks whether a wave EXISTS. None
%   asks whether it MATTERS. Davis et al. (2020, Nature) is the template and
%   the strongest form the claim can take: in marmoset MT the timing and
%   position of a spontaneous wave before target onset predicted the magnitude
%   of the evoked response and whether the target was detected, while spatially
%   DISORGANISED fluctuations of the same amplitude did not. That contrast --
%   organised versus merely large -- is the whole result, and it is why the
%   control in section B4 below is not optional.
%
%   This is also the only analysis in the folder that could still return a
%   POSITIVE for rhythmic attentional scanning. cortical_planar_wave_derotation.m
%   and stimulus_loc_traveling_wave.m rejected a PROPAGATING phase ramp in the
%   trial average. A wave can nevertheless be behaviourally relevant on single
%   trials while cancelling in that average -- which is exactly what
%   checks/check_single_trial_planes.m was built to detect, and it found the
%   trial-to-trial direction concentration to be low, i.e. the average was
%   discarding a great deal. This script uses that single-trial structure
%   instead of averaging it away, and asks whether what is discarded is coupled
%   to behaviour.
%
% WHAT IS MEASURED, PER TRIAL AND PER TIME POINT
%   The band-passed analytic signal across the array is reduced to five numbers:
%     PGD(t)     phase-gradient directionality -- how planar the phase map is.
%                Same definition as cortical_planar_wave_PGD.m, computed here on
%                a single trial at a single time rather than on pref_phase.
%     dir(t)     propagation direction of the best-fitting plane.
%     k(t)       its wavenumber (rad/mm).
%     phi_tgt(t) the wave's LOCAL PHASE AT THE TARGET'S CORTICAL PATCH -- where
%                the wavefront sits relative to the tissue that will receive the
%                target. This is the "position" half of "position and timing"
%                and is the quantity a scanning account actually predicts.
%     amp(t)     array-mean band amplitude -- the CONFOUND, carried everywhere.
%
%   Sign convention for (k, dir) is the one used by
%   checks/check_single_trial_planes.m and cortical_planar_wave_derotation.m:
%   the plane is fitted as exp(-i k.r) against the observed phasors, so the
%   maximising k points along the direction of propagation. Do not mix
%   conventions between scripts.
%
% THE FOUR TESTS
%   B1  WAVE STRENGTH.   Is PGD at target onset higher on hits than misses?
%   B2  WAVE DIRECTION.  Is the wave travelling TOWARD the target's cortical
%                        patch more often on hits? Statistic cos(dir - dir_tgt).
%   B3  WAVE POSITION.   Does hit probability depend on the wave's local phase
%                        at the target patch? Measured with the same complex
%                        resultant the rest of the paper uses, so the number is
%                        comparable to the phase-coherence results rather than
%                        being a new kind of statistic.
%   B4  THE DAVIS CONTROL. Repeat B1 with array-mean AMPLITUDE in place of PGD.
%                        If amplitude predicts detection just as well there is
%                        no wave claim to make, only a power claim. B1 and B3
%                        are therefore ALSO reported with amplitude partialled
%                        out of the outcome -- the same amplitude control the
%                        regression analyses already apply.
%
% NULLS
%   Hit/miss labels are permuted WITHIN SESSION, so session composition, trial
%   counts and any session-level difference in wave statistics are preserved
%   exactly. Note deliberately that trials are NOT stratified by difficulty:
%   the permutation already preserves N and the design, and stratifying would
%   change the estimator rather than the null. Difficulty (trialinfo col 18) is
%   carried per trial so it can be used as a covariate if wanted.
%
% THE PGD FLOOR
%   A plane fitted to ~60 noisy single-trial phases returns a non-zero PGD by
%   construction. That floor is estimated honestly, by recomputing PGD on a
%   subsample of trials with the electrode-to-coordinate assignment shuffled --
%   the same null as cortical_planar_wave_PGD.m. It does NOT depend on the
%   outcome, so it cannot manufacture a hit/miss difference; it can only
%   compress one, which makes B1 conservative. The floor is what the event
%   threshold sits above, and it is plotted so the margin is visible.
%
% BY-PRODUCT PANEL (DO_EVENT_RATE)
%   Wave EVENTS -- PGD above that shuffle threshold for at least MIN_DUR_MS --
%   are counted in the pre-target window and their inter-event intervals
%   histogrammed. RAS requires a scan to recur about once per theta cycle, so
%   the theta range is marked. This is DESCRIPTIVE only: it falls out of the
%   event detection the script needs anyway, it has no test attached, and it
%   must be reported as a description rather than as evidence of theta-rhythmic
%   scanning.
%
% DATA  <base>/results_<animal>/<session>/clean_lfp.mat -> `clean_data`
%   FieldTrip raw: .trial {1 x nTrials} [nChan x nTime], .time (0 = stimulus
%   event), .label, .trialinfo, .fsample.
%   trialinfo columns used: 16/17 target x/y (fixation-centred px), 18
%   difficulty, 20 outcome (1 = hit, 5 = miss).
%   V4 channels map onto canonical slots by label -- functions/v4_channel_slots.m.
%
% OUTPUT
%   Plots/scanning/wave_behaviour/cp10_till_100/<band>/wave_behaviour_coupling.pdf
%   results_combined/scanning/wave_behaviour/cp10_till_100/<band>/wave_behaviour_coupling.mat
% =====================================================================

clearvars; close all; clc
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/scanning/functions

%% --- Settings -------------------------------------------------------
animals    = {'hermes','klecks'};
base       = '/mnt/hpc/projects/MWSampling/4Shivangi';
grid_rows  = 8; grid_cols = 8;
SPACING_MM = 0.4;
nChTot     = grid_rows*grid_cols;
rng(2025);

% Bands. Defaults are the two ranges in which the planar gradient was
% significant in BOTH animals (cortical_planar_wave_PGD.m): theta 4.4-6.1 Hz
% and beta 13.3-23.3 Hz. Run one per invocation; output folders are named after
% BAND_NAME so the variants do not overwrite each other.
BAND_NAME  = 'theta';
BANDS      = struct('theta', [4.4 6.1], 'beta', [13.3 23.3]);
BAND       = BANDS.(BAND_NAME);

% Analysis window, relative to the event at t = 0. Everything closes before the
% earliest critical time (54 ms in monkey H), so the wave is measured on
% genuinely PRE-target activity and cannot be contaminated by the evoked
% response whose consequences are being predicted.
WIN            = [-0.60 -0.02];    % s, window PGD is tracked over
T_ANALYSIS     = -0.02;            % s, the single time B1-B4 are evaluated at
TIME_STEP_MS   = 10;               % PGD tracked on this grid, not per sample

MIN_CH     = 40;                   % channels a trial must have to be used
MIN_TRIALS = 200;                  % per animal
MAX_SESS   = Inf;                  % set small (e.g. 3) for a quick run
nPerm      = 1000;
alpha      = 0.05;

% Plane-fit grid, same construction as checks/check_single_trial_planes.m
K_MAX      = 2.5;                  % rad/mm
N_K        = 16;
N_THETA    = 24;

% Target-patch definition: electrodes are weighted by how close their RF centre
% is to the target, and the weighted centroid of their ARRAY positions is the
% cortical patch the target lands on.
SIGMA_DEG   = 1.5;                 % RF-to-target distance constant (deg)
MIN_PATCH_W = 1e-3;                % reject a trial if no electrode is near

% PGD floor / event detection
N_FLOOR_TRIALS = 300;              % trials used for the shuffled-position floor
DO_EVENT_RATE  = true;
MIN_DUR_MS     = 30;
EVENT_PCTILE   = 95;               % of the shuffled-position floor

% Retinotopy calibration -- same constants as erp_latency_wave.m
PIX_PER_DEG   = struct('hermes', 53.24, 'klecks', 50.56);
RF_DATE       = '20170829';
SCREEN_XY     = [1680 1050];
RF_VALID_ONLY = true;

dv_tag  = 'cp10_till_100';
out_dir = fullfile(base,'Plots','scanning','wave_behaviour', dv_tag, BAND_NAME);
res_dir = fullfile(base,'results_combined','scanning','wave_behaviour', dv_tag, BAND_NAME);
if ~exist(out_dir,'dir'), mkdir(out_dir); end
if ~exist(res_dir,'dir'), mkdir(res_dir); end

% Array coordinates -- SAME index<->position map as every other script here.
ch_col = ceil((1:nChTot)' / grid_rows);
ch_row = grid_rows - mod((1:nChTot)' - 1, grid_rows);
XY     = [ch_col, ch_row] * SPACING_MM;
XY     = XY - mean(XY,1);
LIN    = sub2ind([grid_rows grid_cols], ch_row, ch_col);   % channel -> grid cell

% =====================================================================
% ASSUMPTIONS TO CHECK BEFORE TRUSTING THE OUTPUT
% =====================================================================
% (1) t = 0 in clean_data.time is target onset, the same event the phase
%     pipeline used. If it is not, WIN is not pre-target and B1-B4 measure the
%     evoked response instead. The script prints the first session's time axis
%     -- eyeball it before reading any p-value.
% (2) Outcome coding 1 = hit, 5 = miss in trialinfo col 20, as used throughout
%     Correlation_analysis/. Trials with any other code are dropped and the
%     count printed; a large count means the coding changed.
% (3) The target patch depends on RF centres, 18/64 and 21/64 of which were
%     filled from neighbours upstream. RF_VALID_ONLY = true uses only
%     Valid_Gaussian fits; repeat with it false to confirm the answer does not
%     hinge on the filled channels.
% (4) Filtering is zero-phase and is applied to the FULL trial before cutting,
%     so no edge transient and no group delay enters the analysis window.
% =====================================================================

fprintf('SINGLE-TRIAL WAVES AND DETECTION | band = %s (%.1f-%.1f Hz)\n', ...
    BAND_NAME, BAND(1), BAND(2));
fprintf('  wave tracked over %.3f..%.3f s, tests evaluated at %.3f s\n', ...
    WIN(1), WIN(2), T_ANALYSIS);

[E, kmag, kdir] = build_k_grid(XY, K_MAX, N_K, N_THETA);
tgrid = WIN(1) : TIME_STEP_MS/1000 : WIN(2);
[~, iT] = min(abs(tgrid - T_ANALYSIS));

Acell = cell(1, numel(animals));

%% --- Per animal -----------------------------------------------------
for ia = 1:numel(animals)
    animalName = animals{ia};
    sess = list_sessions(base, animalName);
    if isempty(sess)
        warning('No sessions found for %s -- skipping.', animalName); continue
    end
    nS = min(numel(sess), MAX_SESS);
    fprintf('\n##### %s : %d sessions #####\n', animalName, nS);

    ppd = PIX_PER_DEG.(animalName);
    [rf_deg, rf_valid] = elec_rf_deg(base, animalName, nChTot, RF_DATE, ppd, SCREEN_XY);
    if RF_VALID_ONLY, rf_deg(~rf_valid, :) = NaN; end

    PGD = []; DIRW = []; PHI_T = []; AMP = []; HIT = []; SESS = [];
    DIRT = []; DL = []; FLOOR = [];
    nDrop = 0; printedAxis = false;

    for is = 1:nS
        [Z, tinfo, fs, nUsed, nSkip, t0str] = load_analytic(base, animalName, sess{is}, ...
            BAND, WIN, tgrid, nChTot, MIN_CH);
        if isempty(Z), continue, end
        if ~printedAxis
            fprintf('  fs = %g Hz | session time axis: %s (0 = event)\n', fs, t0str);
            printedAxis = true;
        end
        nDrop = nDrop + nSkip;

        % ---- per-trial wave descriptors + the shuffled-position floor ----
        nFl = min(N_FLOOR_TRIALS, size(Z,3));
        [pgd_t, dir_t, k_t, phi0_t, amp_t, pgd_fl] = wave_descriptors(Z, E, kmag, kdir, ...
            LIN, grid_rows, grid_cols, nFl);

        % ---- target patch on the array, per trial ------------------------
        tgt_deg = tinfo(:,[16 17]) / ppd;               % fixation-centred deg
        [patch_xy, patch_ok] = target_patch(tgt_deg, rf_deg, XY, SIGMA_DEG, MIN_PATCH_W);
        dir_tgt = atan2(patch_xy(:,2), patch_xy(:,1));  % array centre -> patch

        % local wave phase AT the patch: phi = phi0 + k . r_patch
        nTr = size(pgd_t,1);
        phi_at_tgt = nan(nTr,1);
        kv = [k_t(:,iT).*cos(dir_t(:,iT)), k_t(:,iT).*sin(dir_t(:,iT))];
        okp = patch_ok & isfinite(phi0_t(:,iT)) & all(isfinite(kv),2);
        phi_at_tgt(okp) = angle(exp(1i*(phi0_t(okp,iT) + sum(kv(okp,:).*patch_xy(okp,:),2))));

        % ---- outcome ------------------------------------------------------
        oc  = tinfo(:,20);
        hit = nan(size(oc));
        hit(oc == 1) = 1; hit(oc == 5) = 0;
        nDrop = nDrop + sum(~isfinite(hit));
        keep = isfinite(hit) & okp & isfinite(pgd_t(:,iT));

        PGD   = [PGD;   pgd_t(keep,:)];            %#ok<AGROW>
        DIRW  = [DIRW;  dir_t(keep,iT)];           %#ok<AGROW>
        PHI_T = [PHI_T; phi_at_tgt(keep)];         %#ok<AGROW>
        AMP   = [AMP;   amp_t(keep,:)];            %#ok<AGROW>
        DIRT  = [DIRT;  dir_tgt(keep)];            %#ok<AGROW>
        HIT   = [HIT;   hit(keep)];                %#ok<AGROW>
        DL    = [DL;    tinfo(keep,18)];           %#ok<AGROW>
        SESS  = [SESS;  repmat(is, sum(keep), 1)]; %#ok<AGROW>
        FLOOR = [FLOOR; pgd_fl(:)];                %#ok<AGROW>
        fprintf('  . %-45s %5d trials used, %4d skipped\n', ...
            shortname(sess{is}), sum(keep), nSkip);
    end

    if numel(HIT) < MIN_TRIALS
        warning('  %s: only %d usable trials -- skipped.', animalName, numel(HIT)); continue
    end
    fprintf('  pooled: %d trials (%d hits, %d misses), %d dropped\n', ...
        numel(HIT), sum(HIT==1), sum(HIT==0), nDrop);

    % =================================================================
    % B1  wave STRENGTH at target onset
    % =================================================================
    s_pgd = PGD(:,iT);
    [d_pgd, p_pgd] = perm_stat(s_pgd, HIT, SESS, nPerm, 'diff');

    % =================================================================
    % B2  wave DIRECTION relative to the target's cortical patch
    % =================================================================
    align = cos(DIRW - DIRT);
    [d_ali, p_ali] = perm_stat(align, HIT, SESS, nPerm, 'diff');

    % =================================================================
    % B3  wave POSITION: local wave phase at the patch vs hit probability
    % =================================================================
    [R_phi, p_phi] = phase_outcome_coherence(PHI_T, HIT, SESS, nPerm);

    % =================================================================
    % B4  THE DAVIS CONTROL, then B1/B3 again with amplitude removed
    % =================================================================
    s_amp = AMP(:,iT);
    [d_amp, p_amp] = perm_stat(s_amp, HIT, SESS, nPerm, 'diff');

    hit_r = residualise(HIT, s_amp, SESS);          % outcome, amplitude out
    [d_pgd_r, p_pgd_r] = perm_stat(s_pgd, hit_r, SESS, nPerm, 'corr');
    [R_phi_r, p_phi_r] = phase_outcome_coherence(PHI_T, hit_r, SESS, nPerm);

    % =================================================================
    % by-product: wave events and their recurrence interval
    % =================================================================
    EV = struct('rate',NaN,'iei',[],'thr',NaN);
    if DO_EVENT_RATE && ~isempty(FLOOR)
        EV.thr = prctile(FLOOR, EVENT_PCTILE);
        [EV.rate, EV.iei] = event_stats(PGD, tgrid, EV.thr, MIN_DUR_MS);
        fprintf('    wave events: %.2f /trial-second above PGD %.3f (shuffled floor p%.0f),\n', ...
            EV.rate, EV.thr, EVENT_PCTILE);
        fprintf('                 median inter-event interval %.0f ms\n', ...
            1e3*median(EV.iei,'omitnan'));
    end

    % ---- report -------------------------------------------------------
    fprintf('\n  --- %s : single-trial waves vs detection (%s) ---\n', animalName, BAND_NAME);
    fprintf('    B1 strength   PGD  hit-miss %+0.4f   p = %.4f\n', d_pgd, p_pgd);
    fprintf('    B2 direction  cos  hit-miss %+0.4f   p = %.4f  (toward target patch)\n', d_ali, p_ali);
    fprintf('    B3 position   R    %.4f              p = %.4f  (wave phase at patch)\n', R_phi, p_phi);
    fprintf('    B4 CONTROL    amp  hit-miss %+0.4f   p = %.4f  <-- must be weaker than B1/B3\n', d_amp, p_amp);
    fprintf('       B1 with amplitude partialled out  r = %+0.4f  p = %.4f\n', d_pgd_r, p_pgd_r);
    fprintf('       B3 with amplitude partialled out  R = %.4f    p = %.4f\n', R_phi_r, p_phi_r);
    verdict = behaviour_verdict(p_pgd_r, p_ali, p_phi_r, p_amp, alpha);
    fprintf('    --> %s\n', verdict);

    Acell{ia} = struct('animal',animalName,'nTrials',numel(HIT),'nHit',sum(HIT==1), ...
        'tgrid',tgrid,'iT',iT,'BAND',BAND,'BAND_NAME',BAND_NAME, ...
        'd_pgd',d_pgd,'p_pgd',p_pgd,'d_ali',d_ali,'p_ali',p_ali, ...
        'R_phi',R_phi,'p_phi',p_phi,'d_amp',d_amp,'p_amp',p_amp, ...
        'd_pgd_r',d_pgd_r,'p_pgd_r',p_pgd_r,'R_phi_r',R_phi_r,'p_phi_r',p_phi_r, ...
        'PGD',PGD,'AMP',AMP,'HIT',HIT,'PHI_T',PHI_T,'align',align,'DL',DL, ...
        'SESS',SESS,'FLOOR',FLOOR,'EV',EV,'verdict',verdict);
end

A = [Acell{~cellfun(@isempty, Acell)}];

%% --- Replication across animals -------------------------------------
fprintf('\n===== REPLICATION =====\n');
if numel(A) == numel(animals)
    rep = all([A.p_pgd_r] < alpha) || all([A.p_phi_r] < alpha) || all([A.p_ali] < alpha);
    ctl = all([A.p_amp] >= alpha | abs([A.d_amp]) < abs([A.d_pgd]));
    if rep && ctl
        fprintf('  Wave-behaviour coupling replicated in both animals, and it is not amplitude.\n');
    elseif rep
        fprintf('  Coupling replicated, BUT amplitude explains as much -- report it as a power effect.\n');
    else
        fprintf('  NOT replicated: no single-trial wave-behaviour coupling in these data.\n');
    end
else
    fprintf('  Only %d of %d animals produced a result -- no replication statement.\n', ...
        numel(A), numel(animals));
end

%% --- Figure ---------------------------------------------------------
nA = numel(A);
fg = new_fig(1500, 420*max(nA,1));
for j = 1:nA
    S = A(j);

    % (a) PGD time course, hits vs misses, against the shuffled floor
    ax = subplot(nA, 4, (j-1)*4 + 1); hold(ax,'on');
    plot(ax, 1e3*S.tgrid, mean(S.PGD(S.HIT==1,:),1,'omitnan'), 'LineWidth',1.6,'DisplayName','hit');
    plot(ax, 1e3*S.tgrid, mean(S.PGD(S.HIT==0,:),1,'omitnan'), 'LineWidth',1.6,'DisplayName','miss');
    if ~isempty(S.FLOOR)
        yline(ax, median(S.FLOOR,'omitnan'), 'k:','LineWidth',1.2,'DisplayName','shuffled floor');
    end
    xline(ax, 1e3*S.tgrid(S.iT), 'k--','HandleVisibility','off');
    xlabel(ax,'time from target (ms)'); ylabel(ax,'PGD');
    title(ax, {sprintf('%s -- B1 wave strength', esc(S.animal)), ...
               sprintf('hit-miss %+.4f, p = %.4f', S.d_pgd, S.p_pgd)},'FontSize',9);
    legend(ax,'Location','best','FontSize',7); grid(ax,'on');

    % (b) direction alignment to the target patch
    ax = subplot(nA, 4, (j-1)*4 + 2); hold(ax,'on');
    histogram(ax, S.align(S.HIT==1), 20,'Normalization','pdf','FaceAlpha',.55, ...
        'EdgeColor','none','DisplayName','hit');
    histogram(ax, S.align(S.HIT==0), 20,'Normalization','pdf','FaceAlpha',.55, ...
        'EdgeColor','none','DisplayName','miss');
    xlabel(ax,'cos(wave dir - target dir)'); ylabel(ax,'pdf');
    title(ax, {'B2 direction toward target patch', ...
               sprintf('hit-miss %+.4f, p = %.4f', S.d_ali, S.p_ali)},'FontSize',9);
    legend(ax,'Location','best','FontSize',7); grid(ax,'on');

    % (c) hit rate vs wave phase at the patch
    ax = subplot(nA, 4, (j-1)*4 + 3);
    nb = 12; edges = linspace(-pi,pi,nb+1); ctr = edges(1:end-1)+diff(edges)/2;
    hr = nan(nb,1);
    for b = 1:nb
        sel = S.PHI_T >= edges(b) & S.PHI_T < edges(b+1);
        if any(sel), hr(b) = mean(S.HIT(sel)); end
    end
    bar(ax, ctr, hr, 1, 'FaceColor',[.45 .35 .65],'EdgeColor','none');
    yline(ax, mean(S.HIT), 'k--');
    xlabel(ax,'wave phase at target patch (rad)'); ylabel(ax,'hit rate');
    title(ax, {'B3 wave position', ...
               sprintf('R = %.4f, p = %.4f (amp out: p = %.4f)', S.R_phi, S.p_phi, S.p_phi_r)},'FontSize',9);
    xlim(ax,[-pi pi]); grid(ax,'on');

    % (d) the Davis control + the by-product event-interval histogram
    ax = subplot(nA, 4, (j-1)*4 + 4); hold(ax,'on');
    if ~isempty(S.EV.iei)
        histogram(ax, 1e3*S.EV.iei, 30, 'FaceColor',[.6 .6 .6],'EdgeColor','none', ...
            'DisplayName','inter-event interval');
        xline(ax, 1000/8, 'r-','LineWidth',1.5,'DisplayName','8 Hz');
        xline(ax, 1000/4, 'r--','LineWidth',1.5,'DisplayName','4 Hz');
        xlabel(ax,'inter-event interval (ms)'); ylabel(ax,'count');
        legend(ax,'Location','best','FontSize',7);
    end
    title(ax, {sprintf('B4 control: amp hit-miss %+.4f (p = %.4f)', S.d_amp, S.p_amp), ...
               'by-product: wave-event recurrence, theta marked'},'FontSize',8);
    grid(ax,'on');
end
sgtitle({sprintf('SINGLE-TRIAL WAVES AND DETECTION (%s band): do wave POSITION and TIMING predict seeing the target?', BAND_NAME), ...
         'B4 is the load-bearing control: if array-mean amplitude predicts detection as well as PGD does, this is a power effect, not a wave effect.', ...
         'The inter-event histogram is DESCRIPTIVE -- no test is attached and it is not evidence of theta-rhythmic scanning.'}, ...
         'FontSize', 9);
save_pdf(fg, fullfile(out_dir,'wave_behaviour_coupling.pdf'));

results = struct('A',A,'animals',{animals},'BAND',BAND,'BAND_NAME',BAND_NAME, ...
    'WIN',WIN,'T_ANALYSIS',T_ANALYSIS,'TIME_STEP_MS',TIME_STEP_MS,'K_MAX',K_MAX, ...
    'N_K',N_K,'N_THETA',N_THETA,'SIGMA_DEG',SIGMA_DEG,'nPerm',nPerm,'alpha',alpha, ...
    'RF_VALID_ONLY',RF_VALID_ONLY,'MIN_DUR_MS',MIN_DUR_MS,'EVENT_PCTILE',EVENT_PCTILE, ...
    'N_FLOOR_TRIALS',N_FLOOR_TRIALS);
save(fullfile(res_dir,'wave_behaviour_coupling.mat'),'results','-v7.3');
fprintf('\nSaved under %s\n', out_dir);

%% ===================== local functions ==============================

function s = list_sessions(base, animal)
d = dir(fullfile(base, ['results_' animal], [animal '_*attentional-sampling*']));
d = d([d.isdir]);
s = arrayfun(@(x) fullfile(x.folder, x.name), d, 'uni', 0);
end

function n = shortname(p), [~, n] = fileparts(p); end

function [Z, tinfo, fs, nUsed, nSkip, t0str] = load_analytic(base, animal, sessDir, ...
        BAND, WIN, tgrid, nChTot, MIN_CH) %#ok<INUSL>
% -> Z : nCh x nTime x nTrials complex analytic signal on tgrid, NaN where a
%        channel is absent. Filtering is applied to the FULL trial and only
%        then cut, so edge transients stay outside the analysis window, and it
%        is zero-phase so no group delay enters a phase measurement.
Z = []; tinfo = []; fs = NaN; nUsed = 0; nSkip = 0; t0str = '';
fp = fullfile(sessDir, 'clean_lfp.mat');
if ~isfile(fp), return, end
S = load(fp, 'clean_data'); D = S.clean_data;
fs = D.fsample;
t0str = sprintf('%.3f .. %.3f s', D.time{1}(1), D.time{1}(end));
[rows, slots, ~] = v4_channel_slots(D.label, animal, nChTot);
if isempty(rows), return, end

[bb, aa] = butter(3, BAND/(fs/2), 'bandpass');
nT  = numel(tgrid);
buf = complex(nan(nChTot, nT, numel(D.trial)));
ti  = nan(numel(D.trial), size(D.trialinfo,2));

for it = 1:numel(D.trial)
    tt = D.time{it};
    if tt(1) > WIN(1) || tt(end) < WIN(2), nSkip = nSkip + 1; continue, end
    raw = D.trial{it}(rows,:);
    ok  = all(isfinite(raw), 2);
    if sum(ok) < MIN_CH, nSkip = nSkip + 1; continue, end
    fl  = filtfilt(bb, aa, raw(ok,:).');           % nTime x nGood
    an  = hilbert(fl);
    seg = interp1(tt, an, tgrid, 'linear').';      % nGood x nT, complex
    nUsed = nUsed + 1;
    buf(slots(ok), :, nUsed) = seg;
    ti(nUsed,:) = D.trialinfo(it,:);
end
Z = buf(:,:,1:nUsed);
tinfo = ti(1:nUsed,:);
end

function [E, kmag, kdir] = build_k_grid(XY, K_MAX, N_K, N_THETA)
% De-rotation design over (wavenumber x direction) plus the k = 0 cell. Same
% construction as checks/check_single_trial_planes.m, so a fitted direction
% means the same thing in both scripts.
ks = linspace(K_MAX/N_K, K_MAX, N_K);
th = linspace(0, 2*pi, N_THETA+1); th(end) = [];
kx = 0; ky = 0; kmag = 0; kdir = 0;
for i = 1:numel(ks)
    kx   = [kx,   ks(i)*cos(th)]; %#ok<AGROW>
    ky   = [ky,   ks(i)*sin(th)]; %#ok<AGROW>
    kmag = [kmag, repmat(ks(i),1,N_THETA)]; %#ok<AGROW>
    kdir = [kdir, th]; %#ok<AGROW>
end
E = exp(-1i*(kx(:)*XY(:,1).' + ky(:)*XY(:,2).'));      % grid x chan
kmag = kmag(:); kdir = kdir(:);
end

function [pgd, dirw, kw, phi0, amp, pgd_floor] = wave_descriptors(Z, E, kmag, kdir, ...
        LIN, gr, gc, nFloor)
% Per trial and per time point: PGD, best-plane direction and wavenumber, the
% plane's phase offset, and the array-mean amplitude. pgd_floor is the same PGD
% recomputed on the first nFloor trials with the electrode-to-coordinate
% assignment shuffled -- the chance level for this array and this much phase
% scatter, and the threshold the event detection sits above.
[nCh, nT, nTr] = size(Z);
pgd = nan(nTr, nT); dirw = nan(nTr, nT); kw = nan(nTr, nT);
phi0 = nan(nTr, nT); amp = nan(nTr, nT); pgd_floor = nan(nFloor, nT);
for it = 1:nTr
    zt = Z(:,:,it);                                    % nCh x nT
    good = all(isfinite(zt), 2);
    if sum(good) < 20, continue, end
    amp(it,:) = mean(abs(zt(good,:)), 1);

    % --- PGD from local finite differences, vectorised over time -------
    PH = nan(gr*gc, nT);
    PH(LIN(good), :) = angle(zt(good,:));
    pgd(it,:) = pgd_from_maps(reshape(PH, gr, gc, nT));

    % --- best plane over the (k, theta) grid ---------------------------
    u  = zt(good,:) ./ max(abs(zt(good,:)), eps);      % unit phasors
    Eg = E(:, good);
    Rg = abs(Eg * u) / sum(good);                      % grid x nT
    [~, b] = max(Rg, [], 1);
    dirw(it,:) = kdir(b);
    kw(it,:)   = kmag(b);
    % each time point against ITS OWN best plane -- the diagonal of Eg(b,:)*u,
    % taken directly so the off-diagonal nT x nT product is never formed
    phi0(it,:) = angle(sum(Eg(b,:).' .* u, 1));

    % --- the floor, on the first nFloor trials -------------------------
    if it <= nFloor
        pm = randperm(nCh);
        PHs = nan(gr*gc, nT);
        gs = pm(good);
        PHs(LIN(gs), :) = angle(zt(good,:));
        pgd_floor(it,:) = pgd_from_maps(reshape(PHs, gr, gc, nT));
    end
end
end

function v = pgd_from_maps(PH)
% PGD = |mean(grad phi)| / mean(|grad phi|) per time slice of a gr x gc x nT
% phase array. Wrapped differences, no 2-D unwrap -- identical definition to
% cortical_planar_wave_PGD.m.
dx = angle(exp(1i*(PH(:,2:end,:) - PH(:,1:end-1,:))));
dy = angle(exp(1i*(PH(2:end,:,:) - PH(1:end-1,:,:))));
gx = dx(1:end-1,:,:); gy = dy(:,1:end-1,:);
g  = gx + 1i*gy;
nT = size(PH,3);
v  = nan(1, nT);
for j = 1:nT
    gj = g(:,:,j); gj = gj(isfinite(gj));
    if numel(gj) < 10, continue, end
    v(j) = abs(mean(gj)) / max(mean(abs(gj)), eps);
end
end

function [patch_xy, ok] = target_patch(tgt_deg, rf_deg, XY, sigma, minw)
% Weighted centroid, in ARRAY coordinates, of the electrodes whose receptive
% fields sit near the target: the cortical patch the target lands on.
nTr = size(tgt_deg,1);
patch_xy = nan(nTr,2); ok = false(nTr,1);
for it = 1:nTr
    if ~all(isfinite(tgt_deg(it,:))), continue, end
    d2 = sum((rf_deg - tgt_deg(it,:)).^2, 2);
    w  = exp(-d2 / (2*sigma^2));
    w(~isfinite(w)) = 0;
    sw = sum(w);
    if sw < minw, continue, end
    patch_xy(it,:) = (w.' * XY) / sw;
    ok(it) = true;
end
end

function [d, p] = perm_stat(x, y, sess, nPerm, kind)
% Hit-minus-miss difference ('diff') or correlation ('corr') between x and the
% outcome y, with y permuted WITHIN SESSION so session composition and N are
% preserved exactly. Trials are deliberately NOT stratified by difficulty --
% the permutation already preserves the design.
use = isfinite(x) & isfinite(y);
x = x(use); y = y(use); sess = sess(use);
d = stat_(x, y, kind);
us  = unique(sess);
idx = arrayfun(@(k) find(sess == k), us, 'uni', 0);
dn  = nan(nPerm,1);
for pp = 1:nPerm
    ys = y;
    for k = 1:numel(idx)
        ii = idx{k};
        ys(ii) = y(ii(randperm(numel(ii))));
    end
    dn(pp) = stat_(x, ys, kind);
end
p = (1 + sum(abs(dn) >= abs(d))) / (1 + nPerm);
end

function s = stat_(x, y, kind)
switch kind
    case 'diff', s = mean(x(y==1)) - mean(x(y==0));
    case 'corr', s = corr(x, y);
    otherwise, error('Unknown statistic: %s', kind);
end
end

function [R, p] = phase_outcome_coherence(phi, y, sess, nPerm)
% The same complex resultant the rest of the paper uses: weight each trial's
% phasor by its mean-removed outcome and take the normalised magnitude, so a
% dependence of outcome on phase produces a large R and no dependence cancels.
% Works unchanged for a binary outcome and for a residualised continuous one.
use = isfinite(phi) & isfinite(y);
phi = phi(use); y = y(use); sess = sess(use);
w = y - mean(y);
R = abs(sum(w .* exp(1i*phi))) / max(sum(abs(w)), eps);
us  = unique(sess);
idx = arrayfun(@(k) find(sess == k), us, 'uni', 0);
Rn  = nan(nPerm,1);
for pp = 1:nPerm
    ys = y;
    for k = 1:numel(idx)
        ii = idx{k};
        ys(ii) = y(ii(randperm(numel(ii))));
    end
    ws = ys - mean(ys);
    Rn(pp) = abs(sum(ws .* exp(1i*phi))) / max(sum(abs(ws)), eps);
end
p = (1 + sum(Rn >= R)) / (1 + nPerm);
end

function r = residualise(y, x, sess)
% Remove the linear contribution of amplitude from the outcome, within session.
r = nan(size(y));
us = unique(sess);
for k = 1:numel(us)
    m = sess == us(k) & isfinite(x) & isfinite(y);
    if sum(m) < 10, continue, end
    Xd = [x(m), ones(sum(m),1)];
    r(m) = y(m) - Xd * (Xd \ y(m));
end
end

function [rate, iei] = event_stats(PGD, tgrid, thr, minDurMs)
% Wave events = PGD above the shuffled-position threshold for >= minDurMs.
dt   = mean(diff(tgrid));
minN = max(1, round(minDurMs/1000/dt));
nEv  = 0; iei = [];
for it = 1:size(PGD,1)
    v = PGD(it,:) > thr;
    v(~isfinite(PGD(it,:))) = false;
    st = find(diff([false v]) == 1);
    en = find(diff([v false]) == -1);
    st = st((en - st + 1) >= minN);
    nEv = nEv + numel(st);
    if numel(st) > 1, iei = [iei, diff(tgrid(st))]; end %#ok<AGROW>
end
rate = nEv / (size(PGD,1) * (tgrid(end) - tgrid(1)));
iei  = iei(:);
end

function v = behaviour_verdict(p_pgd_r, p_ali, p_phi_r, p_amp, alpha)
anyWave = (p_pgd_r < alpha) || (p_ali < alpha) || (p_phi_r < alpha);
if ~anyWave
    v = 'NO coupling: single-trial wave strength, direction and position do not predict detection';
elseif p_amp < alpha && p_pgd_r >= alpha && p_phi_r >= alpha
    v = 'AMPLITUDE effect, not a wave effect -- the control explains it';
else
    v = 'COUPLING: wave organisation predicts detection beyond array-mean amplitude';
end
end

function s = esc(s), s = strrep(s, '_', '\_'); end

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
