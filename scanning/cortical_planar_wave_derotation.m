% =====================================================================
% Cortical planar traveling-wave detection by PHASE-VECTOR ALIGNMENT
% across ELECTRODES, swept over frequency × speed × direction
% ("rotate-to-overlap" plane fit) — the cortical-space twin of
% stimulus_loc_traveling_wave.m.
%
% Idea : each electrode c sits at a physical position r_c = (x_c,y_c) mm on
% the 8x8 array and has a preferred-phase vector e^{iφ_c}. If a PLANAR
% traveling wave of speed v and direction θ crosses the array, the phase at
% each electrode is delayed relative to a reference by
%       Δt_c = d_c(θ)/v ,   d_c(θ) = x_c cosθ + y_c sinθ   (mm along the
%                                                            wave axis)
% i.e. rotated by  k·d_c(θ),  k = 2π f / v  (rad/mm). De-rotating each
% electrode's phase vector by k·d_c(θ) makes them OVERLAP (resultant R → 1)
% — but only at the wave's true (f, v, θ). So we sweep a grid of
% (frequency × speed × direction) and read off where the de-rotated vectors
% align. A GENUINE wave shows high R along a (near-)CONSTANT speed ACROSS
% frequencies (speed is frequency-independent, i.e. k ∝ f) — the key
% discriminator from a chance per-frequency alignment.
%
% Structural difference from the stim-loc test: there the "distance" d_p was
% a given scalar (stimulus eccentricity along the radial axis) so only speed
% was swept. Here the wave direction on the 2-D array is unknown, so we also
% sweep DIRECTION θ (equivalently, fit the 2-D wavevector k). R is reported
% maxed over θ; the best θ is recorded.
%
% COHERENT across electrodes: we sum the complex phase vectors keeping each
% electrode's own phase (no per-channel |·|), because cortical electrodes
% share a common reference (all measured against the same stimulus/clock), so
% the per-electrode delay IS the wave and is directly testable:
%       R(f,v,θ) = | Σ_c w_c e^{i(φ_c − k d_c(θ))} | / Σ_c w_c ,
% with w_c = coherence magnitude. Positions are collapsed into one map per
% frequency (magnitude-weighted, as in the collapsed PGD test).
%
% TWO readouts (both as in stimulus_loc_traveling_wave.m):
%   (1) ABSOLUTE R with a max-statistic permutation threshold that corrects
%       over the whole (f,v,θ) grid.
%   (2) INCREASE IN COHERENCE (de-rotation GAIN) dR = R(f,v,θ) − R0(f),
%       where R0 is the un-rotated (k=0) coherence — the "actual" electrode
%       coherence with NO wave assumed (pure synchrony). dR isolates the
%       wave-specific phase ramp and is invariant to electrodes merely
%       sharing a common phase. Its own max-stat null is clean because R0 is
%       unaffected by the electrode<->position shuffle (k=0 uses no d).
%
% Significance: permute the electrode<->position assignment (shuffle which
% array coordinate each (phase,weight) sits at) with a single synchronised
% shuffle shared across the whole grid, recompute, take its max -> max-stat
% threshold correcting across the entire (f,v,θ) grid, for both R and dR.
%
% Pipeline, PER ANIMAL (never pool channels across animals): grid -> animals
%   -> mean grid + "significant in BOTH" replication (for R and for the gain).
%
% Units: physical array coordinates in mm (pitch SPACING_MM); speed in cm/s
% (converted to mm/s internally); k in rad/mm. SPEED_OK flags the
% physiologically plausible cortical band.
%
% Output:
%   Plots/scanning/planar_wave_derotation/cp10_till_100/<dv>/
%       derotation_R_grids.pdf         (absolute R)
%       derotation_gain_grids.pdf      (increase in coherence dR + best-speed line)
%       derotation_direction_speed.pdf (best-fit direction & speed vs frequency)
%   results_combined/scanning/planar_wave_derotation/cp10_till_100/<dv>/
%       planar_wave_derotation.mat
% =====================================================================

% =====================================================================
% TWO ESTIMATORS (setting: ESTIMATOR = 'phase' | 'trial')
% =====================================================================
% Both build the same object — one complex number per (electrode, frequency),
% whose phase is that electrode's preferred phase and whose magnitude says how
% reliable it is — and then run the identical de-rotation, null, gain and
% pooling machinery. They differ only in HOW that number is formed:
%
%   'phase' (default, unchanged): collapse the per-location coherences from
%       phase_progression.mat with a magnitude-weighted circular sum,
%           z_c(f) = SUM_p coh_mag(c,p,f) * exp(i*pref_phase(c,p,f)) .
%       Each term is a per-location MEAN over trials, so a location with few
%       trials contributes a term whose magnitude is inflated (|c_p| ~ 1/sqrt(n_p)).
%
%   'trial': one phase coherence over ALL trials, in the phase_coherence/ sense,
%           z_c(f) = ( 1/W_c ) * SUM over all trials  y_t * exp(i*phi_tf) ,
%           W_c = SUM_t |y_t| ,
%       so |z_c| is that channel's overall coherence and every trial counts once.
%       Computed from the cached per-location sums S(c,p,f) that
%       functions/trial_position_sums_chan.m writes (one SLURM job per channel),
%       via z_c(f) = SUM_p S(c,p,f) / W_c.
%
%   NOTE this script only needs the OBSERVED sums: its null shuffles the
%   electrode <-> array-coordinate assignment, not trial labels, so the
%   permutation part of the cached worker output is not used here (~0.3 MB per
%   animal rather than the ~286 MB the stimulus-location script needs).
%   The cache is shared with stimulus_loc_traveling_wave.m.
%
% =====================================================================
% COMBINING ANIMALS: replication AND pooled standardised z
% =====================================================================
%   REPLICATION (primary): sig in BOTH animals, cell by cell. Cannot be driven
%     by one animal; cannot aggregate weak evidence; low power.
%   POOLED (secondary): z-score each animal's grid against its OWN permutation
%     null, average the z's, and test that against a paired-permutation
%     max-stat null. Aggregates evidence, so more power — but with two animals
%     a pooled hit can be ~100% one animal, so the per-animal z at the peak
%     cell is printed alongside.
%   Standardising first is essential: the animals sit on very different scales
%   (R0 ~0.95 vs ~0.99 here), so a raw average would just track the larger one.
%   ONE pooled statistic is reported, not two: because gain = R - R0(f) and R0
%   is constant across speed and across permutations within a frequency, the
%   constant cancels in both the numerator and the null, so the standardised
%   gain grid is IDENTICAL to the standardised R grid. Verified numerically.
%   Hence the pooled column appears once, on the R figure.
% =====================================================================

clearvars; close all; clc

%% ─── Dependencies (only the 'trial' estimator needs slurmfun) ────────
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/scanning/functions

%% ─── Settings ────────────────────────────────────────────────────────
animals    = {'hermes','klecks'};
dv         = 'lfp';
base       = '/mnt/hpc/projects/MWSampling/4Shivangi';
grid_rows  = 8; grid_cols = 8;
SPACING_MM = 0.4;                 % electrode pitch (mm)
nPerm      = 1000;
alpha      = 0.05;
CH_FILTER  = 'significant';       % coherence-sig channels only (as in planar/stim-loc)
COH_SIG_ALPHA = 0.05;
MIN_CH     = 8;                   % min reliable channels to attempt a plane fit
FREQ_RANGE = [2 100];             % Hz swept
NV         = 30;                  % # speed steps
V_CORTICAL = logspace(log10(1), log10(300), NV);   % cortical wave speed (cm/s)
NTHETA     = 24;                  % # propagation-direction steps (15° each)
THETA      = linspace(0, 2*pi, NTHETA+1); THETA(end) = [];
SPEED_OK   = [5 100];             % plausible cortical wave speed (cm/s), for flagging
rng(2025);

% Which estimator forms the per-electrode complex map (see header).
% 'phase' = collapse the per-location coherences (original, unchanged).
% 'trial' = one phase coherence over all trials; needs the cached SLURM sums.
ESTIMATOR = 'phase';
RECOMPUTE_TRIAL_SUMS = false;     % force re-submission of the per-channel jobs
switch ESTIMATOR
    case 'phase', tag = '';        est_note = 'estimator: phase (per-location preferred phases)';
    case 'trial', tag = '_trial';  est_note = 'estimator: trial (phase coherence over all trials)';
    otherwise,    error('Unknown ESTIMATOR: %s', ESTIMATOR);
end

out_dir = fullfile(base,'Plots','scanning','planar_wave_derotation','cp10_till_100', dv);
res_dir = fullfile(base,'results_combined','scanning','planar_wave_derotation','cp10_till_100', dv);
if ~exist(out_dir,'dir'), mkdir(out_dir); end
if ~exist(res_dir,'dir'), mkdir(res_dir); end

% electrode grid coordinates (mm); same index<->position map as the PGD script
ch_col = ceil((1:(grid_rows*grid_cols))' / grid_rows);
ch_row = grid_rows - mod((1:(grid_rows*grid_cols))' - 1, grid_rows);
XY = [ch_col, ch_row] * SPACING_MM;      % [nCh x 2] (x,y) in mm

G = struct();   % G(ia): grids, thresholds, sig, gain, best speed/dir

%% ─── Per animal ──────────────────────────────────────────────────────
for ia = 1:numel(animals)
    animalName = animals{ia};
    pp = fullfile(base, ['results_' animalName], 'scanning', ...
        'phase_progression','cp10_till_100', dv, 'phase_progression.mat');
    if ~isfile(pp), warning('No data for %s — skipping.', animalName); continue; end
    S = load(pp, 'pref_phase','coh_mag','freq','positions');
    freq = S.freq(:); [nCh,nFreq,nPos] = size(S.pref_phase);
    fprintf('\n=== %s / %s : %d ch, %d freq, %d pos ===\n', animalName, upper(dv), nCh, nFreq, nPos);

    % coherence-significance mask [nCh x nFreq], identical to planar/stim-loc
    if strcmp(CH_FILTER,'significant')
        coh_root = fullfile(base, ['results_' animalName], 'phase_coherence', ...
            'complex','cp10_till_100', dv, 'all_loc_difflev');
        coh_sig = load_coh_sig_mask(coh_root, nCh, nFreq, COH_SIG_ALPHA);
        fprintf('  coherence-significant channels: %d–%d across freqs (median %d)\n', ...
            min(sum(coh_sig,1)), max(sum(coh_sig,1)), round(median(sum(coh_sig,1))));
    else
        coh_sig = true(nCh, nFreq);
    end

    f_use = find(freq >= FREQ_RANGE(1) & freq <= FREQ_RANGE(2));
    fHz   = freq(f_use);
    Vs_mm = V_CORTICAL * 10;                 % cm/s -> mm/s (k = 2*pi*f / v_mm)

    % ── Per-electrode complex map Zmap(ch,freq): the only thing the two
    %    estimators do differently (see header) ─────────────────────────
    switch ESTIMATOR
        case 'phase'
            Zmap = complex(nan(nCh, nFreq));
            for f = 1:nFreq
                PHI = squeeze(S.pref_phase(:,f,:));  W = squeeze(S.coh_mag(:,f,:));
                M = isfinite(PHI) & isfinite(W) & (W > 0);  W(~M) = 0; PHI(~M) = 0;
                Zmap(:,f) = sum(W .* exp(1i*PHI), 2);      % magnitude-weighted circular sum
            end
        case 'trial'
            % z_c(f) = SUM_p S(c,p,f) / W_c = coherence over ALL trials.
            % Only the observed sums are needed; the null here is spatial.
            [Ssum, Wch] = get_trial_sums_obs(base, animalName, dv, nCh, nPerm, RECOMPUTE_TRIAL_SUMS);
            Zmap = squeeze(sum(Ssum, 2)) ./ Wch;           % nCh × nFreq
            fprintf('  trial estimator: |z| range %.3f–%.3f (per-channel coherence)\n', ...
                min(abs(Zmap(:)),[],'omitnan'), max(abs(Zmap(:)),[],'omitnan'));
    end

    [Robs, R0, Rnull_max, Gnull_max, DIRbest, Rnull, Gnull] = derotate_grid( ...
        Zmap, coh_sig, f_use, fHz, XY(1:nCh,:), Vs_mm, THETA, MIN_CH, nPerm);

    % (1) absolute R: max-stat threshold over the whole (f,v,θ) grid
    thr = quantile(Rnull_max, 1-alpha);
    sig = Robs >= thr;

    % (2) increase in coherence (de-rotation gain) dR = R - R0, with its
    %     own max-stat null (R0 is shuffle-invariant -> clean gain null)
    gain     = Robs - R0(:);                  % nF × nV (broadcast per row)
    thr_gain = quantile(Gnull_max, 1-alpha);
    sig_gain = gain >= thr_gain;

    % best speed per frequency + is that gain significant?
    [gpk, vbest] = max(Robs, [], 2);
    vbest_speed  = V_CORTICAL(vbest(:));      % cm/s
    vbest_dir    = arrayfun(@(fi) DIRbest(fi, vbest(fi)), (1:numel(fHz))');
    vbest_sig    = (gpk - R0(:)) >= thr_gain;

    G(ia).animal = animalName;
    G(ia).R = Robs; G(ia).thr = thr; G(ia).sig = sig;
    G(ia).R0 = R0(:); G(ia).gain = gain; G(ia).thr_gain = thr_gain; G(ia).sig_gain = sig_gain;
    G(ia).fHz = fHz; G(ia).speeds = V_CORTICAL; G(ia).DIRbest = DIRbest;
    G(ia).vbest = vbest_speed; G(ia).vbest_dir = vbest_dir; G(ia).vbest_sig = vbest_sig;

    % ── Standardise against this animal's OWN null, cell by cell ─────────
    % 'omitnan' matters: mean/std do NOT skip NaN, so one NaN permutation would
    % make that cell NaN for every permutation and the NaN would then spread
    % through the cross-animal sum (NaN + finite = NaN).
    mu = mean(Rnull,3,'omitnan');  sd = std(Rnull,0,3,'omitnan');
    G(ia).z_obs  = (Robs  - mu) ./ max(sd, eps);
    G(ia).z_null = (Rnull - mu) ./ max(sd, eps);
    fprintf('  finite cells in standardised grid = %d/%d\n', ...
        sum(isfinite(G(ia).z_obs(:))), numel(G(ia).z_obs));

    [rmax,ix] = max(Robs(:)); [fi,vi] = ind2sub(size(Robs), ix);
    fprintf('  peak R=%.3f at f=%.1f Hz, v=%.1f cm/s, dir=%.0f° | thr=%.3f | sig cells=%d\n', ...
        rmax, fHz(fi), V_CORTICAL(vi), rad2deg(DIRbest(fi,vi)), thr, sum(sig(:)));
    [gmax,gix] = max(gain(:)); [gfi,gvi] = ind2sub(size(gain), gix);
    okspeed = V_CORTICAL(gvi) >= SPEED_OK(1) && V_CORTICAL(gvi) <= SPEED_OK(2);
    fprintf('  peak gain dR=%.3f at f=%.1f Hz, v=%.1f cm/s %s | thr_gain=%.3f | sig-gain cells=%d | freqs w/ sig gain=%d/%d\n', ...
        gmax, fHz(gfi), V_CORTICAL(gvi), ternary(okspeed,'(plausible)','(NON-PHYSICAL)'), ...
        thr_gain, sum(sig_gain(:)), sum(vbest_sig), numel(vbest_sig));
end
valid = find(arrayfun(@(s) ~isempty(s.animal), G));

%% ─── Combine animals (mean grid + replication) ───────────────────────
Rsum = 0; Gsum = 0; Zsum = 0; Znull = 0;
repl = true(size(G(valid(1)).R)); repl_gain = repl;
for ia = valid
    Rsum = Rsum + G(ia).R;
    Gsum = Gsum + G(ia).gain;
    repl = repl & G(ia).sig;
    repl_gain = repl_gain & G(ia).sig_gain;
    Zsum  = Zsum  + G(ia).z_obs;
    Znull = Znull + G(ia).z_null;
end
nA = numel(valid);
C.R = Rsum/nA; C.repl = repl;
C.gain = Gsum/nA; C.repl_gain = repl_gain;
C.fHz = G(valid(1)).fHz; C.speeds = V_CORTICAL;

% POOLED standardised test. Permutation b of one animal is paired with
% permutation b of the other — legitimate because the animals are independent,
% so the pairing is arbitrary and samples the product null. max() omits NaN per
% column, so skipped frequencies are harmless; an all-NaN null (no usable
% channels anywhere) is caught explicitly rather than yielding thr = NaN.
C.Z  = Zsum/nA;
maxZ = max(reshape(Znull/nA, [], nPerm), [], 1).';
tZ   = quantile(maxZ, 1-alpha);
if ~isfinite(tZ)
    warning(['pooled null is entirely NaN — at least one animal produced no finite ' ...
             'cells (too few usable channels?). Pooled test disabled; replication unaffected.']);
    tZ = Inf;
end
C.thr_pool = tZ;  C.sig_pool = C.Z >= tZ;  C.pooled_ok = isfinite(tZ);

fprintf('\n================ COMBINED ================\n');
fprintf('cells significant in BOTH animals (absolute R)      = %d\n', sum(repl(:)));
fprintf('cells gain-significant in BOTH animals (increase-R)  = %d\n', sum(repl_gain(:)));
fprintf('cells significant in the POOLED standardised test    = %d (thr z = %.2f)\n', ...
    sum(C.sig_pool(:)), C.thr_pool);
if any(C.sig_pool(:))
    [~, ix] = max(C.Z(:) .* double(C.sig_pool(:)));
    [pf, pv] = ind2sub(size(C.Z), ix);
    za = arrayfun(@(ia) G(ia).z_obs(pf,pv), valid);
    fprintf('  peak pooled cell f=%.2f Hz v=%.1f cm/s: per-animal z = [%s] (%s)\n', ...
        C.fHz(pf), C.speeds(pv), strjoin(compose('%.2f', za), ', '), strjoin({G(valid).animal}, ', '));
end

%% ─── Figure 1: freq × speed ABSOLUTE R heatmaps ──────────────────────
% R maxed over direction; significant cells outlined. A wave = high R along a
% (near-)constant speed across frequencies.
esc  = @(s) strrep(s, '_', '\_');   % TeX renders '_' as a subscript otherwise
ncol = numel(valid)+2;              % animals + mean/replication + pooled z
f1 = figure('Visible','off','Position',[40 40 380*ncol 340]);
for k = 1:numel(valid)
    ia = valid(k);
    ax = subplot(1, ncol, k); hold(ax,'on');
    imagesc(ax, G(ia).fHz, 1:numel(G(ia).speeds), G(ia).R.');
    set(ax,'YDir','normal'); axis(ax,'tight');
    contour(ax, G(ia).fHz, 1:numel(G(ia).speeds), double(G(ia).sig.'), [0.5 0.5], 'w','LineWidth',1.2);
    speed_yticks(ax, G(ia).speeds, SPEED_OK);
    caxis(ax,[0 1]); colorbar(ax);
    xlabel(ax,'Frequency (Hz)'); ylabel(ax,'wave speed (cm/s)');
    title(ax,sprintf('%s — R (sig outlined)', G(ia).animal),'FontSize',9);
end
ax = subplot(1, ncol, ncol-1); hold(ax,'on');
imagesc(ax, C.fHz, 1:numel(C.speeds), C.R.');
set(ax,'YDir','normal'); axis(ax,'tight');
contour(ax, C.fHz, 1:numel(C.speeds), double(C.repl.'), [0.5 0.5], 'w','LineWidth',1.4);
speed_yticks(ax, C.speeds, SPEED_OK);
caxis(ax,[0 1]); colorbar(ax);
xlabel(ax,'Frequency (Hz)'); ylabel(ax,'wave speed (cm/s)');
title(ax, {'mean', 'contour = replication'},'FontSize',9);

% POOLED standardised z. Units are z, NOT R — auto-scaled, so its colours are
% not comparable with the panels to the left. One pooled panel only: the
% standardised gain grid is identical to this one (see header).
ax = subplot(1, ncol, ncol); hold(ax,'on');
imagesc(ax, C.fHz, 1:numel(C.speeds), C.Z.');
set(ax,'YDir','normal'); axis(ax,'tight');
contour(ax, C.fHz, 1:numel(C.speeds), double(C.sig_pool.'), [0.5 0.5], 'w','LineWidth',1.4);
speed_yticks(ax, C.speeds, SPEED_OK);
if ~any(isfinite(C.Z(:))), caxis(ax,[0 1]); end
colorbar(ax);
xlabel(ax,'Frequency (Hz)'); ylabel(ax,'wave speed (cm/s)');
if C.pooled_ok, ttl_p = sprintf('POOLED z (thr %.2f)', C.thr_pool);
else,           ttl_p = 'POOLED z - n/a (see warning)'; end
title(ax, {'pooled standardised', ttl_p},'FontSize',9);

sgtitle({['Cortical planar-wave de-rotation grid: R after de-rotating electrodes by k*d(theta)   (' esc(est_note) ')'], ...
         'white contour = significant   |   col 3 = mean, contour = replication   |   col 4 = pooled standardised z'}, ...
         'FontSize',9);
set(f1,'PaperPositionMode','auto'); pos=get(f1,'Position'); set(f1,'PaperUnits','points','PaperSize',pos(3:4));
saveas(f1, fullfile(out_dir, ['derotation_R_grids' tag '.pdf']));

%% ─── Figure 2: INCREASE IN COHERENCE (de-rotation gain) dR ───────────
% dR(f,v)=R(f,v)-R0(f): speeds at which de-rotating BEATS no rotation. A wave
% = significant positive gain (contour) along a near-constant best speed
% (black line) across frequencies.
gmax_m = max(C.gain(:));
for ia = valid, gmax_m = max(gmax_m, max(G(ia).gain(:))); end
ncol_g = numel(valid)+1;    % animals + mean; the pooled panel lives on figure 1 only
f2 = figure('Visible','off','Position',[40 40 380*ncol_g 340]);
for k = 1:numel(valid)
    ia = valid(k);
    ax = subplot(1, ncol_g, k); hold(ax,'on');
    imagesc(ax, G(ia).fHz, 1:numel(G(ia).speeds), G(ia).gain.');
    set(ax,'YDir','normal'); axis(ax,'tight');
    contour(ax, G(ia).fHz, 1:numel(G(ia).speeds), double(G(ia).sig_gain.'), [0.5 0.5], 'w','LineWidth',1.2);
    plot(ax, G(ia).fHz, interp1(G(ia).speeds, 1:numel(G(ia).speeds), G(ia).vbest, 'linear', NaN), 'k-','LineWidth',1.3);
    speed_yticks(ax, G(ia).speeds, SPEED_OK);
    caxis(ax,[0 max(gmax_m,eps)]); colorbar(ax);
    xlabel(ax,'Frequency (Hz)'); ylabel(ax,'wave speed (cm/s)');
    title(ax, {G(ia).animal, 'gain dR'},'FontSize',9);
end
ax = subplot(1, ncol_g, ncol_g); hold(ax,'on');
imagesc(ax, C.fHz, 1:numel(C.speeds), C.gain.');
set(ax,'YDir','normal'); axis(ax,'tight');
contour(ax, C.fHz, 1:numel(C.speeds), double(C.repl_gain.'), [0.5 0.5], 'w','LineWidth',1.4);
speed_yticks(ax, C.speeds, SPEED_OK);
caxis(ax,[0 max(gmax_m,eps)]); colorbar(ax);
xlabel(ax,'Frequency (Hz)'); ylabel(ax,'wave speed (cm/s)');
title(ax, {'mean gain', 'contour = replication'},'FontSize',9);
sgtitle({['Increase in coherence dR = R(f,v) - R0(f): does de-rotating beat the un-rotated coherence?   (' esc(est_note) ')'], ...
         'white contour = significant   |   black line = best-fit speed   |   the pooled standardised test is on the R-grid figure (it is identical for R and dR)'}, ...
         'FontSize',9);
set(f2,'PaperPositionMode','auto'); pos=get(f2,'Position'); set(f2,'PaperUnits','points','PaperSize',pos(3:4));
saveas(f2, fullfile(out_dir, ['derotation_gain_grids' tag '.pdf']));

%% ─── Figure 3: best-fit propagation DIRECTION vs frequency ────────────
% At each frequency the best-fit plane wave has a direction (taken at the
% peak-R speed). Grey = all frequencies; filled = frequencies whose gain is
% significant (a real wave). A genuine wave = a CONSISTENT direction across
% the significant band (points cluster); scattered directions => no coherent
% wave. Bottom row: cross-animal direction agreement per frequency.
cols = lines(numel(animals));
f3 = figure('Visible','off','Position',[40 40 380*numel(valid) 560]);
for k = 1:numel(valid)
    ia = valid(k);
    fr = G(ia).fHz; dg = mod(rad2deg(G(ia).vbest_dir(:)), 360); sg = G(ia).vbest_sig(:);
    % top: direction vs frequency (deg), significant highlighted
    ax = subplot(2, numel(valid), k); hold(ax,'on');
    plot(ax, fr, dg, 'o', 'Color',[.75 .75 .75], 'MarkerSize',4, 'DisplayName','all freqs');
    if any(sg)
        plot(ax, fr(sg), dg(sg), 'o', 'Color',cols(ia,:), 'MarkerFaceColor',cols(ia,:), ...
            'MarkerSize',6, 'DisplayName','sig gain');
        Rdir = abs(mean(exp(1i*G(ia).vbest_dir(sg))));    % direction consistency 0..1
        ttl = sprintf('%s — prop. direction (sig: R_{dir}=%.2f)', G(ia).animal, Rdir);
    else
        ttl = sprintf('%s — prop. direction (no sig band)', G(ia).animal);
    end
    ylim(ax,[0 360]); set(ax,'YTick',0:90:360);
    xlabel(ax,'Frequency (Hz)'); ylabel(ax,'best-fit direction (deg)');
    title(ax, ttl,'FontSize',9); grid(ax,'on'); legend(ax,'Location','best','FontSize',7);
    % bottom: best speed vs frequency (constant speed => real wave), sig marked
    ax2 = subplot(2, numel(valid), numel(valid)+k); hold(ax2,'on');
    plot(ax2, fr, G(ia).vbest(:), 'o', 'Color',[.75 .75 .75], 'MarkerSize',4);
    if any(sg), plot(ax2, fr(sg), G(ia).vbest(sg), 'o', 'Color',cols(ia,:),'MarkerFaceColor',cols(ia,:),'MarkerSize',6); end
    yline(ax2, SPEED_OK(1), 'k:'); yline(ax2, SPEED_OK(2), 'k:');
    set(ax2,'YScale','log'); ylim(ax2,[V_CORTICAL(1) V_CORTICAL(end)]);
    xlabel(ax2,'Frequency (Hz)'); ylabel(ax2,'best-fit speed (cm/s)');
    title(ax2,'best speed (flat across freq = real wave; dotted = plausible band)','FontSize',8);
    grid(ax2,'on');
end
sgtitle('Best-fit planar-wave direction & speed vs frequency (filled = significant increase in coherence)','FontSize',11);
set(f3,'PaperPositionMode','auto'); pos=get(f3,'Position'); set(f3,'PaperUnits','points','PaperSize',pos(3:4));
saveas(f3, fullfile(out_dir, ['derotation_direction_speed' tag '.pdf']));

results = struct('G',G,'C',C,'animals',{animals},'dv',dv, ...
    'FREQ_RANGE',FREQ_RANGE,'V_CORTICAL',V_CORTICAL,'THETA',THETA, ...
    'SPACING_MM',SPACING_MM,'SPEED_OK',SPEED_OK,'nPerm',nPerm,'alpha',alpha,'ESTIMATOR',ESTIMATOR);
% drop the per-permutation null grids before saving (~100 MB otherwise)
for ia = valid, G(ia).z_null = []; end
results.G = G;
save(fullfile(res_dir, ['planar_wave_derotation' tag '.mat']),'results','-v7.3');
fprintf('\nSaved figures + results under %s\n', out_dir);

%% =====================================================================
%% Helpers
%% =====================================================================
function [Robs, R0, Rnull_max, Gnull_max, DIRbest, Rnull, Gnull] = derotate_grid(Zmap, coh_sig, f_use, fHz, XY, Vs_mm, Th, MIN_CH, nPerm)
% Coherent cortical planar-wave plane fit by de-rotation, over the
% (frequency × speed × direction) grid.
%   R(f,v,θ) = | Σ_c w_c e^{i(φ_c − k d_c(θ)) } | / Σ_c w_c ,  k = 2πf/v_mm,
%              d_c(θ) = x_c cosθ + y_c sinθ  (mm).
% Positions are collapsed per frequency (magnitude-weighted circular sum).
% Returns R maxed over direction (Robs), the k=0 baseline R0, per-perm grid
% maxima of R (absolute null) and of dR=R−R0 (gain null), and the best
% direction per (f,v). Null shuffles the electrode<->position assignment
% (permutes each (phase,weight) to a random array coordinate); R0 is
% shuffle-invariant (k=0 uses no coordinate) -> clean gain null.
nF = numel(f_use); nV = numel(Vs_mm);
Robs = nan(nF, nV); R0 = nan(nF, 1); DIRbest = nan(nF, nV);
Rnull_max = -inf(nPerm,1); Gnull_max = -inf(nPerm,1);
Rnull = nan(nF, nV, nPerm); Gnull = nan(nF, nV, nPerm);   % kept for cross-animal pooling
ct = [cos(Th(:)), sin(Th(:))].';          % 2 × nTheta

for fi = 1:nF
    f = f_use(fi); fh = fHz(fi);
    z    = Zmap(:, f);                    % nCh × 1 per-electrode complex map
    useC = coh_sig(:, f) & isfinite(z) & (abs(z) > 0);
    if sum(useC) < MIN_CH, continue; end

    A  = z(useC);                          % nUse × 1 (weighted phase vectors)
    sw = sum(abs(A));                      % Σ w_c
    P  = XY(useC, :);                       % nUse × 2 (mm)
    D  = P * ct;                            % nUse × nTheta (d along each dir)

    % no-rotation baseline (k=0): "actual" electrode coherence (synchrony)
    R0(fi) = abs(sum(A)) / max(sw, eps);

    % observed grid + per-speed de-rotation phasor cache (reused for the null)
    Ecache = cell(1, nV);
    for vi = 1:nV
        k = 2*pi*fh / Vs_mm(vi);            % rad/mm
        E = exp(-1i * k * D);               % nUse × nTheta
        Ecache{vi} = E;
        Zt = abs(A.' * E) / max(sw, eps);   % 1 × nTheta
        [rmx, ti] = max(Zt);
        Robs(fi,vi) = rmx; DIRbest(fi,vi) = Th(ti);
    end

    % null: permute (phase,weight) across electrode positions; sw invariant.
    for b = 1:nPerm
        Ap = A(randperm(numel(A)));
        rb = -inf;
        for vi = 1:nV
            m = max(abs(Ap.' * Ecache{vi})) / max(sw, eps);
            Rnull(fi,vi,b) = m;                 % full grid kept for pooling
            Gnull(fi,vi,b) = m - R0(fi);
            if m > rb, rb = m; end
        end
        if rb > Rnull_max(b),          Rnull_max(b) = rb;          end
        if rb - R0(fi) > Gnull_max(b), Gnull_max(b) = rb - R0(fi); end
    end
end
if ~all(isfinite(Rnull_max))
    warning(['%d/%d permutations produced no finite statistic — the max-stat ' ...
             'threshold is INVALID (thr=-inf makes every cell "significant").'], ...
             sum(~isfinite(Rnull_max)), nPerm);
end
end

function [Ssum, Wch] = get_trial_sums_obs(base, animalName, dv, nCh, nPerm, force)
% Observed per-location complex sums S(c,p,f) and per-channel normalisers W(c)
% for the 'trial' estimator. One SLURM job per channel via
% functions/trial_position_sums_chan.m, cached on disk and SHARED with
% stimulus_loc_traveling_wave.m. Only S_obs is read here: this script's null
% shuffles electrode <-> array coordinate, not trial labels, so the (large)
% permutation part of each file is never loaded.
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
        cfg{ch}.perm_seed_base = 2025;
    end
    slurmfun(@trial_position_sums_chan, cfg, ...
        'partition','8GB', 'stopOnError',false, 'useUserPath',true);
else
    fprintf('  re-using cached trial-level position sums in %s\n', sum_dir);
end

first = '';
for ch = 1:nCh
    fch = fullfile(sum_dir, num2str(ch), 'trial_position_sums.mat');
    if isfile(fch), first = fch; break; end
end
if isempty(first), error('no trial-level position sums found under %s', sum_dir); end
D0 = load(first, 'S_obs');
[nPos, nFreq] = size(D0.S_obs);
Ssum = complex(nan(nCh, nPos, nFreq));  Wch = nan(nCh,1);
for ch = 1:nCh
    fch = fullfile(sum_dir, num2str(ch), 'trial_position_sums.mat');
    if ~isfile(fch), warning('missing trial sums for channel %d', ch); continue; end
    D = load(fch, 'S_obs', 'W');       % NOT S_perm — not needed here
    Ssum(ch,:,:) = D.S_obs;  Wch(ch) = D.W;
end
fprintf('  trial sums: %d/%d channels loaded (observed only, %.1f MB)\n', ...
    sum(isfinite(Wch)), nCh, numel(Ssum)*16/1e6);
end

function speed_yticks(ax, speeds, SPEED_OK)
% label the (index-based) speed axis with actual cm/s + mark the plausible band
yt = round(linspace(1, numel(speeds), 5));
set(ax, 'YTick', yt, 'YTickLabel', compose('%.0f', speeds(yt)));
lo = interp1(speeds, 1:numel(speeds), SPEED_OK(1), 'linear', NaN);
hi = interp1(speeds, 1:numel(speeds), SPEED_OK(2), 'linear', NaN);
if isfinite(lo), yline(ax, lo, 'w:', 'LineWidth', 0.8); end
if isfinite(hi), yline(ax, hi, 'w:', 'LineWidth', 0.8); end
end

function [sig, freq] = load_coh_sig_mask(coh_root, nCh, nFreq, alpha)
% Per-channel coherence-significance mask [nCh x nFreq] from the
% phase_coherence/complex pipeline. Per channel, the max-statistic
% (freq-corrected) threshold is the (1-alpha) quantile of the per-perm max
% over frequency of |coh_perm_complex|; the channel is significant at freq f
% where |coh_complex(f)| exceeds that single threshold. (Identical to the
% mask used by cortical_planar_wave_PGD.m and stimulus_loc_traveling_wave.m.)
sig = false(nCh, nFreq); freq = [];
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

function out = ternary(c,a,b), if c, out=a; else, out=b; end, end
