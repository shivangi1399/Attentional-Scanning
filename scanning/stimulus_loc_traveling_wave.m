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
% eccentricity (deg). Two ways to combine channels are computed:
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
% Also computed for each mode: the de-rotation GAIN dR = R(f,v) − R0(f),
% where R0 is the un-rotated (k=0) coherence — the "actual" coherence with no
% wave assumed; dR isolates the wave-specific phase ramp.
%
% Pipeline, PER ANIMAL (never pooled): grid -> animals -> mean grid +
%   "significant in both" replication.
% Significance: permute the location<->distance (and, for the coherent mode,
%   electrode<->Dc) assignment with a SINGLE synchronised shuffle shared
%   across the grid, recompute, take its max -> max-statistic threshold that
%   corrects across the entire (f,v) grid. R0 is shuffle-invariant, giving a
%   clean max-stat null on the gain dR too.
%
% Units (resolved from code/RF_Mapping/chan_loc_mua.m):
%   The stimulus 'positions' are trialinfo col-16 = x target location in
%   fixation-centred SCREEN PIXELS (col-17 = y, loaded here for the 2D
%   eccentricity). The fovea/fixation is the origin (0,0), so a location's
%   eccentricity = hypot(x,y) pixels, and degrees = pixels / PIX_PER_DEG,
%   where PIX_PER_DEG is the PER-ANIMAL pixels-per-degree (rig calibration;
%   the raw .RF sessInfo did not store ppd). VISUAL distance d(p) =
%   eccentricity (deg).
%
% Output (rows = the two modes, cols = animals + combined):
%   Plots/scanning/phase_alignment_wave/cp10_till_100/<dv>/
%       phase_alignment_grids.pdf   (absolute R)
%       phase_alignment_gain.pdf    (de-rotation gain dR + best-speed line)
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
NV         = 30;                   % # speed steps
V_VISUAL   = logspace(log10(1),   log10(200), NV);   % deg/s
rng(2025);

% Experiment calibration
PIX_PER_DEG   = struct('hermes', 53.24, 'klecks', 50.56);   % pixels per degree (per-animal rig calibration)

% Electrode RF centres (for the COHERENT mode). Each channel's RF-centre
% eccentricity in DEGREES comes from the Gaussian-fit centre already written
% by RF_Mapping/mapping_lfp.m into the per-channel summary table (the SAME
% file traveling_planar_wave.m reads). Its RF_Center_X/Y are in screen pixels
% with fixation at the screen centre, so eccentricity (deg) =
% hypot(X-cx, Y-cy)/ppd — same fovea-origin convention as the stimulus d_p.
% Both animals are 8x8 = 64-channel arrays and the summary has 64 rows 1:1
% with the phase-progression channel index, so no channel remap is needed.
RF_DATE   = '20170829';   % RF session both animals share
SCREEN_XY = [1680 1050];  % screen pixels; fixation/fovea at the centre

out_dir = fullfile(base,'Plots','scanning','phase_alignment_wave','cp10_till_100', dv);
res_dir = fullfile(base,'results_combined','scanning','phase_alignment_wave','cp10_till_100', dv);
if ~exist(out_dir,'dir'), mkdir(out_dir); end
if ~exist(res_dir,'dir'), mkdir(res_dir); end

% Two combination modes, both in VISUAL degrees / deg-per-second:
%   visual          : INCOHERENT — |resultant| per channel over stimulus
%                     locations, then averaged (per-channel phase offset
%                     discarded; robust but blind to per-electrode delays).
%   visual_coherent : COHERENT — sum complex phase vectors across ALL
%                     channel x location pairs after de-rotating each by
%                     k*(d_p + D_c), where D_c = electrode RF eccentricity
%                     (deg). Tests ONE cortical wave in visual space; a
%                     per-electrode delay is now visible (assumes phases
%                     share a common reference — a null is thus ambiguous).
metrics = {'visual','visual_coherent'};
Vsets   = {V_VISUAL, V_VISUAL};
Vunit   = {'deg/s','deg/s'};
cmode   = {'incoherent','coherent'};

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
    fprintf('  ppd=%.2f px/deg | eccentricity %.2f-%.2f deg\n', ppd, min(ecc), max(ecc));

    % electrode RF-centre eccentricity (deg) for the coherent mode; referenced
    % like d_p (subtract min). NaN where a channel has no valid RF centre.
    ecc_elec = elec_ecc_deg(base, animalName, nCh, RF_DATE, ppd, SCREEN_XY);
    Dc = ecc_elec - min(ecc_elec, [], 'omitnan');
    fprintf('  electrode RF ecc: %d/%d channels with a centre | %.2f-%.2f deg\n', ...
        sum(isfinite(ecc_elec)), nCh, min(ecc_elec,[],'omitnan'), max(ecc_elec,[],'omitnan'));

    for mi = 1:numel(metrics)
        d = d_vis(:).'; Vs = Vsets{mi};
        if strcmp(cmode{mi},'coherent')
            [Robs, R0, Rnull_max, Gnull_max] = align_grid_coherent(S.pref_phase, S.coh_mag, coh_sig, ...
                f_use, fHz, d, Dc, Vs, MIN_LOC, nPerm);
        else
            [Robs, R0, Rnull_max, Gnull_max] = align_grid(S.pref_phase, S.coh_mag, coh_sig, ...
                f_use, fHz, d, Vs, MIN_LOC, nPerm);
        end
        thr = quantile(Rnull_max, 1-alpha);
        sig = Robs >= thr;

        % de-rotation gain vs the un-rotated ("actual") coherence + its
        % max-stat null: dR(f,v)=R(f,v)-R0(f); a wave shows a significant
        % positive gain at a (near-)constant best speed across frequencies.
        gain     = Robs - R0(:);                 % nF × nV (broadcast per row)
        thr_gain = quantile(Gnull_max, 1-alpha);
        sig_gain = gain >= thr_gain;
        [gpk, vbest] = max(Robs, [], 2);         % best speed index per freq
        vbest_speed  = Vs(vbest(:));             % best de-rotation speed per freq
        vbest_sig    = (gpk - R0(:)) >= thr_gain;% is that gain significant?

        G.(metrics{mi})(ia).animal = animalName;
        G.(metrics{mi})(ia).R = Robs; G.(metrics{mi})(ia).thr = thr;
        G.(metrics{mi})(ia).sig = sig; G.(metrics{mi})(ia).fHz = fHz;
        G.(metrics{mi})(ia).speeds = Vs;
        G.(metrics{mi})(ia).R0 = R0(:);              % no-rotation baseline per freq
        G.(metrics{mi})(ia).gain = gain;             % de-rotation gain grid
        G.(metrics{mi})(ia).thr_gain = thr_gain;
        G.(metrics{mi})(ia).sig_gain = sig_gain;
        G.(metrics{mi})(ia).vbest = vbest_speed;     % best speed per freq
        G.(metrics{mi})(ia).vbest_sig = vbest_sig;

        [rmax,ix] = max(Robs(:)); [fi,vi] = ind2sub(size(Robs), ix);
        fprintf('  %-8s: peak R=%.3f at f=%.1f Hz, v=%.1f %s | thr=%.3f | sig cells=%d\n', ...
            metrics{mi}, rmax, fHz(fi), Vs(vi), Vunit{mi}, thr, sum(sig(:)));
        [gmax,gix] = max(gain(:)); [gfi,gvi] = ind2sub(size(gain), gix);
        fprintf('  %-8s: peak gain dR=%.3f at f=%.1f Hz, v=%.1f %s | thr_gain=%.3f | sig-gain cells=%d | freqs w/ sig gain=%d/%d\n', ...
            metrics{mi}, gmax, fHz(gfi), Vs(gvi), Vunit{mi}, thr_gain, sum(sig_gain(:)), sum(vbest_sig), numel(vbest_sig));
    end
end
valid = find(arrayfun(@(s) ~isempty(s.animal), G.(metrics{1})));

%% ─── Combine animals (mean grid + replication) ───────────────────────
C = struct();
for mi = 1:numel(metrics)
    m = metrics{mi};
    Rsum = 0; Gsum = 0;
    repl = true(size(G.(m)(valid(1)).R)); repl_gain = repl;
    for ia = valid
        Rsum = Rsum + G.(m)(ia).R;
        Gsum = Gsum + G.(m)(ia).gain;
        repl = repl & G.(m)(ia).sig;
        repl_gain = repl_gain & G.(m)(ia).sig_gain;
    end
    C.(m).R = Rsum/numel(valid); C.(m).repl = repl;
    C.(m).gain = Gsum/numel(valid); C.(m).repl_gain = repl_gain;
    C.(m).fHz = G.(m)(valid(1)).fHz; C.(m).speeds = Vsets{mi};
    fprintf('%-8s: cells significant in BOTH animals = %d | gain-significant in BOTH = %d\n', ...
        m, sum(repl(:)), sum(repl_gain(:)));
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

%% ─── Figure: de-rotation GAIN over the un-rotated ("actual") coherence ─
% dR(f,v)=R(f,v)-R0(f): the speeds at which de-rotating beats no rotation.
% A wave = significant positive gain (contour) along a near-constant best
% speed (black line) across frequencies.
f2 = figure('Visible','off','Position',[40 40 360*ncol 300*numel(metrics)]);
for mi = 1:numel(metrics)
    m = metrics{mi};
    speeds_to_idx = @(sp,sv) interp1(sv, 1:numel(sv), sp, 'linear', NaN);
    gmax_m = max(C.(m).gain(:));
    for ia = valid, gmax_m = max(gmax_m, max(G.(m)(ia).gain(:))); end
    for k = 1:numel(valid)
        ia = valid(k);
        ax = subplot(numel(metrics), ncol, (mi-1)*ncol + k); hold(ax,'on');
        imagesc(ax, G.(m)(ia).fHz, 1:numel(G.(m)(ia).speeds), G.(m)(ia).gain.');
        set(ax,'YDir','normal'); axis(ax,'tight');
        contour(ax, G.(m)(ia).fHz, 1:numel(G.(m)(ia).speeds), double(G.(m)(ia).sig_gain.'), [0.5 0.5], 'w','LineWidth',1.2);
        plot(ax, G.(m)(ia).fHz, speeds_to_idx(G.(m)(ia).vbest, G.(m)(ia).speeds), 'k-','LineWidth',1.3);
        yt = round(linspace(1,numel(G.(m)(ia).speeds),5));
        set(ax,'YTick',yt,'YTickLabel',compose('%.0f',G.(m)(ia).speeds(yt)));
        caxis(ax,[0 max(gmax_m,eps)]); colorbar(ax);
        xlabel(ax,'Frequency (Hz)'); ylabel(ax,sprintf('%s speed (%s)',m,Vunit{mi}));
        title(ax,sprintf('%s — %s gain \\DeltaR (sig outlined, best-v line)', G.(m)(ia).animal, m),'FontSize',9);
    end
    ax = subplot(numel(metrics), ncol, (mi-1)*ncol + ncol); hold(ax,'on');
    imagesc(ax, C.(m).fHz, 1:numel(C.(m).speeds), C.(m).gain.');
    set(ax,'YDir','normal'); axis(ax,'tight');
    contour(ax, C.(m).fHz, 1:numel(C.(m).speeds), double(C.(m).repl_gain.'), [0.5 0.5], 'w','LineWidth',1.4);
    yt = round(linspace(1,numel(C.(m).speeds),5));
    set(ax,'YTick',yt,'YTickLabel',compose('%.0f',C.(m).speeds(yt)));
    caxis(ax,[0 max(gmax_m,eps)]); colorbar(ax);
    xlabel(ax,'Frequency (Hz)'); ylabel(ax,sprintf('%s speed (%s)',m,Vunit{mi}));
    title(ax,sprintf('mean gain (repl. outlined) — %s', m),'FontSize',9);
end
sgtitle('De-rotation gain \DeltaR = R(f,v) - R_0(f): does de-rotating beat the un-rotated coherence?','FontSize',11);
set(f2,'PaperPositionMode','auto'); pos=get(f2,'Position'); set(f2,'PaperUnits','points','PaperSize',pos(3:4));
saveas(f2, fullfile(out_dir,'phase_alignment_gain.pdf'));

results = struct('G',G,'C',C,'animals',{animals},'dv',dv, ...
    'FREQ_RANGE',FREQ_RANGE,'V_VISUAL',V_VISUAL, ...
    'PIX_PER_DEG',PIX_PER_DEG,'RF_DATE',RF_DATE,'SCREEN_XY',SCREEN_XY, ...
    'metrics',{metrics},'cmode',{cmode},'Vunit',{Vunit},'nPerm',nPerm,'alpha',alpha);
save(fullfile(res_dir,'phase_alignment.mat'),'results','-v7.3');
fprintf('\nSaved figures + results under %s\n', out_dir);

%% =====================================================================
%% Helpers
%% =====================================================================
function [Robs, R0, Rnull_max, Gnull_max] = align_grid(pref, coh, coh_sig, f_use, fHz, d, Vs, MIN_LOC, nPerm)
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
            if rcell > rb, rb = rcell; end
        end
        if rb > null_max_running(b),          null_max_running(b)  = rb;         end
        if rb - R0(fi) > gnull_max_running(b), gnull_max_running(b) = rb - R0(fi); end
    end
end
Rnull_max = null_max_running;
Gnull_max = gnull_max_running;
end

function [Robs, R0, Rnull_max, Gnull_max] = align_grid_coherent(pref, coh, coh_sig, f_use, fHz, d, Dc, Vs, MIN_LOC, nPerm)
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
Dc = Dc(:);

% one shared shuffle per permutation, over locations AND over channels,
% reused across every (f,v) cell -> synchronised null + grid-wide max-stat.
permP = zeros(nPerm, nPos); permC = zeros(nPerm, nCh);
for b = 1:nPerm, permP(b,:) = randperm(nPos); permC(b,:) = randperm(nCh); end
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
            if rcell > rb, rb = rcell; end
        end
        if rb > null_max_running(b),          null_max_running(b)  = rb;          end
        if rb - R0(fi) > gnull_max_running(b), gnull_max_running(b) = rb - R0(fi); end
    end
end
Rnull_max = null_max_running;
Gnull_max = gnull_max_running;
end

function ecc = elec_ecc_deg(base, animal, nCh, rf_date, ppd, screen_xy)
% Per-channel electrode RF-centre eccentricity (degrees) from the Gaussian-fit
% centre already written by RF_Mapping/mapping_lfp.m into the per-channel
% summary table (Plots/RF_Mapping/<animal>/loc_RF_map/gaussian_overlap/
% <animal>_<date>_..._channel_target_summary.txt) — the SAME file
% traveling_planar_wave.m reads. RF_Center_X/Y are screen pixels with fixation
% at the screen centre, so eccentricity = hypot(X-cx, Y-cy)/ppd (deg). Rows
% with no fit hold '-' -> NaN. Channel index is 1:1 with the phase-progression
% array (64 rows). NaN where a channel has no centre.
ecc = nan(nCh,1);
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
cx = screen_xy(1)/2; cy = screen_xy(2)/2;             % fixation = screen centre
n = min(nCh, height(t));
ecc(1:n) = hypot(x(1:n) - cx, y(1:n) - cy) / ppd;     % NaN passes through
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
