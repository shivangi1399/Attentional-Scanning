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
% Two distance metrics (the sketch's "cortical" vs "visual"):
%   VISUAL   : d(p) = eccentricity in degrees (linear in stimulus pos);
%              speed in deg/s.
%   CORTICAL : d(p) = cortical distance via a cortical-magnification map
%              M(E) = 1/(E+E2)  ->  x(E) ∝ ln(E+E2); speed in cm/s. The
%              log compression is what makes it differ from VISUAL.
%
% Pipeline, PER ANIMAL (never pooled), then combined:
%   per channel:  R_c(f,v) = | Σ_p w_p e^{i(φ_p - θ_p)} | / Σ_p w_p
%                 (w_p = coherence magnitude; resultant is invariant to the
%                 channel's own phase offset, so channels can be averaged).
%   channels  ->  R_grid(f,v) = mean over coherence-significant channels.
%   animals   ->  mean grid + "significant in both" replication.
% Significance: permute the location<->distance assignment with a SINGLE
%   shuffle shared across channels (synchronised null, as in
%   phase_progression_chan.m), recompute the whole grid, take its max ->
%   max-statistic threshold that corrects across the entire (f,v) grid.
%
% Units (resolved from code/RF_Mapping/chan_loc_mua.m):
%   The stimulus 'positions' are trialinfo col-16 = x target location in
%   fixation-centred SCREEN PIXELS (col-17 = y, loaded here for the 2D
%   eccentricity). The fovea/fixation is the origin (0,0), so a location's
%   eccentricity = hypot(x,y) pixels, and degrees = pixels / PIX_PER_DEG,
%   where PIX_PER_DEG is the PER-ANIMAL pixels-per-degree (rig calibration;
%   the raw .RF sessInfo did not store ppd). VISUAL distance d(p) =
%   eccentricity (deg); CORTICAL distance via CMF M(E)=1/(E+E2) ->
%   x(E)=MM_PER_LOGDEG·ln(E+E2) (cm).
%   The absolute scale only relabels the speed axis; it does NOT change the
%   diagonal-vs-horizontal (wave-or-not) conclusion.
%
% Output:
%   Plots/scanning/phase_alignment_wave/cp10_till_100/<dv>/
%       phase_alignment_grids.pdf
%   results_combined/scanning/phase_alignment_wave/cp10_till_100/<dv>/
%       phase_alignment.mat
% =====================================================================

clearvars; close all; clc

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
NV         = 30;                   % # speed steps per metric
V_VISUAL   = logspace(log10(1),   log10(200), NV);   % deg/s
V_CORTICAL = logspace(log10(5),   log10(200), NV);   % cm/s
rng(2025);

% Experiment calibration
PIX_PER_DEG   = struct('hermes', 53.24, 'klecks', 50.56);   % pixels per degree (per-animal rig calibration)
% Cortical-magnification model:  M(E) = 1/(E+E2)  [mm/deg],
%   integrated to cortical position  x(E) = MM_PER_LOGDEG * ln(E+E2)  [mm].
%   Inverse-linear CMF form: Cowey & Rolls (1974) Exp Brain Res 21:447-454;
%     Rovamo & Virsu (1979) Exp Brain Res 37:495-510.
%   Complex-log / ln(E+E2) mapping: Schwartz (1980) Vision Res 20:645-669.
%   Foundational monkey V1 CMF: Daniel & Whitteridge (1961) J Physiol 159:203-221.
E2_DEG        = 1.0;   % CMF E2 (deg): eccentricity at which M halves.
                       %   Macaque V1 E2 ~= 0.7-1.0 deg: Van Essen, Newsome &
                       %   Maunsell (1984) Vision Res 24:429-448. E2 concept:
                       %   Levi, Klein & Aitsebaomo (1985) Vision Res 25:963-977.
MM_PER_LOGDEG = 8.0;   % cortical mm per natural-log-degree (V4-ish scale).
                       %   Macaque V4 visuotopy/magnification: Gattass, Sousa &
                       %   Gross (1988) J Neurosci 8:1831-1845; Motter (2009)
                       %   J Neurosci 29:5749-5757. (Only rescales the speed
                       %   axis; does not affect the wave-vs-no-wave pattern.)

out_dir = fullfile(base,'Plots','scanning','phase_alignment_wave','cp10_till_100', dv);
res_dir = fullfile(base,'results_combined','scanning','phase_alignment_wave','cp10_till_100', dv);
if ~exist(out_dir,'dir'), mkdir(out_dir); end
if ~exist(res_dir,'dir'), mkdir(res_dir); end

% ppd is known per animal, so positions -> degrees and both metrics run.
metrics = {'visual','cortical'}; Vsets = {V_VISUAL, V_CORTICAL}; Vunit = {'deg/s','cm/s'};

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

    % stimulus location per position in screen pixels: x = col 16 (== the
    % 'positions' values), y = col 17 (from trialinfo). Fovea = (0,0), so
    % eccentricity = hypot(x,y); convert to degrees with the animal's ppd.
    [x_px, y_px] = target_pix_per_pos(positions, base, animalName);
    ppd = PIX_PER_DEG.(animalName);           % pixels per degree (this animal)
    ecc = hypot(x_px, y_px) / ppd;            % eccentricity in degrees
    d_vis = ecc - min(ecc);                   % linear visual distance (deg)
    x_mm  = MM_PER_LOGDEG * log(ecc + E2_DEG);   % cortical position (mm)
    d_cort = (x_mm - min(x_mm)) / 10;         % cortical distance (cm)
    Dsets = {d_vis, d_cort};
    fprintf('  ppd=%.2f px/deg | eccentricity %.2f-%.2f deg\n', ppd, min(ecc), max(ecc));

    for mi = 1:numel(metrics)
        d = Dsets{mi}(:).'; Vs = Vsets{mi};
        [Robs, Rnull_max] = align_grid(S.pref_phase, S.coh_mag, coh_sig, ...
            f_use, fHz, d, Vs, MIN_LOC, nPerm);
        thr = quantile(Rnull_max, 1-alpha);
        sig = Robs >= thr;

        G.(metrics{mi})(ia).animal = animalName;
        G.(metrics{mi})(ia).R = Robs; G.(metrics{mi})(ia).thr = thr;
        G.(metrics{mi})(ia).sig = sig; G.(metrics{mi})(ia).fHz = fHz;
        G.(metrics{mi})(ia).speeds = Vs;

        [rmax,ix] = max(Robs(:)); [fi,vi] = ind2sub(size(Robs), ix);
        fprintf('  %-8s: peak R=%.3f at f=%.1f Hz, v=%.1f %s | thr=%.3f | sig cells=%d\n', ...
            metrics{mi}, rmax, fHz(fi), Vs(vi), Vunit{mi}, thr, sum(sig(:)));
    end
end
valid = find(arrayfun(@(s) ~isempty(s.animal), G.(metrics{1})));

%% ─── Combine animals (mean grid + replication) ───────────────────────
C = struct();
for mi = 1:numel(metrics)
    m = metrics{mi};
    Rsum = 0; repl = true(size(G.(m)(valid(1)).R));
    for ia = valid
        Rsum = Rsum + G.(m)(ia).R;
        repl = repl & G.(m)(ia).sig;
    end
    C.(m).R = Rsum/numel(valid); C.(m).repl = repl;
    C.(m).fHz = G.(m)(valid(1)).fHz; C.(m).speeds = Vsets{mi};
    fprintf('%-8s: cells significant in BOTH animals = %d\n', m, sum(repl(:)));
end

%% ─── Figure: freq × speed alignment heatmaps ─────────────────────────
ncol = numel(valid)+1;
f1 = figure('Visible','off','Position',[40 40 360*ncol 300*numel(metrics)]);
for mi = 1:numel(metrics)
    m = metrics{mi};
    for k = 1:numel(valid)
        ia = valid(k);
        ax = subplot(numel(metrics), ncol, (mi-1)*ncol + k); hold(ax,'on');
        imagesc(ax, G.(m)(ia).fHz, 1:numel(G.(m)(ia).speeds), G.(m)(ia).R.');
        set(ax,'YDir','normal'); axis(ax,'tight');
        contour(ax, G.(m)(ia).fHz, 1:numel(G.(m)(ia).speeds), double(G.(m)(ia).sig.'), [0.5 0.5], 'w','LineWidth',1.2);
        yt = round(linspace(1,numel(G.(m)(ia).speeds),5));
        set(ax,'YTick',yt,'YTickLabel',compose('%.0f',G.(m)(ia).speeds(yt)));
        caxis(ax,[0 1]); colorbar(ax);
        xlabel(ax,'Frequency (Hz)'); ylabel(ax,sprintf('%s speed (%s)',m,Vunit{mi}));
        title(ax,sprintf('%s — %s (sig outlined)', G.(m)(ia).animal, m),'FontSize',9);
    end
    ax = subplot(numel(metrics), ncol, (mi-1)*ncol + ncol); hold(ax,'on');
    imagesc(ax, C.(m).fHz, 1:numel(C.(m).speeds), C.(m).R.');
    set(ax,'YDir','normal'); axis(ax,'tight');
    contour(ax, C.(m).fHz, 1:numel(C.(m).speeds), double(C.(m).repl.'), [0.5 0.5], 'w','LineWidth',1.4);
    yt = round(linspace(1,numel(C.(m).speeds),5));
    set(ax,'YTick',yt,'YTickLabel',compose('%.0f',C.(m).speeds(yt)));
    caxis(ax,[0 1]); colorbar(ax);
    xlabel(ax,'Frequency (Hz)'); ylabel(ax,sprintf('%s speed (%s)',m,Vunit{mi}));
    title(ax,sprintf('mean (repl. outlined) — %s', m),'FontSize',9);
end
sgtitle('Phase-alignment traveling-wave grid: R after de-rotating locations by 2\pi f d/v','FontSize',11);
set(f1,'PaperPositionMode','auto'); pos=get(f1,'Position'); set(f1,'PaperUnits','points','PaperSize',pos(3:4));
saveas(f1, fullfile(out_dir,'phase_alignment_grids.pdf'));

results = struct('G',G,'C',C,'animals',{animals},'dv',dv, ...
    'FREQ_RANGE',FREQ_RANGE,'V_VISUAL',V_VISUAL,'V_CORTICAL',V_CORTICAL, ...
    'PIX_PER_DEG',PIX_PER_DEG,'E2_DEG',E2_DEG,'MM_PER_LOGDEG',MM_PER_LOGDEG, ...
    'metrics',{metrics},'Vunit',{Vunit},'nPerm',nPerm,'alpha',alpha);
save(fullfile(res_dir,'phase_alignment.mat'),'results','-v7.3');
fprintf('\nSaved figure + results under %s\n', out_dir);

%% =====================================================================
%% Helpers
%% =====================================================================
function [Robs, Rnull_max] = align_grid(pref, coh, coh_sig, f_use, fHz, d, Vs, MIN_LOC, nPerm)
% R(f,v) = mean over coherence-sig channels of the de-rotated resultant
% across stimulus locations. Returns observed grid + per-perm grid maxima
% (synchronised location shuffle shared across channels).
nF = numel(f_use); nV = numel(Vs); nPos = numel(d);
Robs      = nan(nF, nV);
Rnull_max = nan(nPerm, 1);

% Pre-shuffle the location order once per permutation (shared across all
% channels and all (f,v) cells -> synchronised null + grid-wide max-stat).
perms = zeros(nPerm, nPos);
for b = 1:nPerm, perms(b,:) = randperm(nPos); end
null_max_running = -inf(nPerm,1);

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

    for vi = 1:nV
        k = 2*pi*fh / Vs(vi);                       % rad per unit distance
        % observed: de-rotate by θ_p = k·d_p
        Robs(fi,vi) = mean( abs(A * exp(-1i*k*d(:))) ./ max(sw,eps) );
        % (per channel: row of A dotted with the de-rotation phasors)
    end

    % null: same channels/weights, shuffled location<->distance mapping
    for b = 1:nPerm
        dp = d(perms(b,:));
        rb = -inf;
        for vi = 1:nV
            k = 2*pi*fh / Vs(vi);
            rcell = mean( abs(A * exp(-1i*k*dp(:))) ./ max(sw,eps) );
            if rcell > rb, rb = rcell; end
        end
        if rb > null_max_running(b), null_max_running(b) = rb; end
    end
end
Rnull_max = null_max_running;
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
