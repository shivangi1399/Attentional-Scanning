% =====================================================================
% Cortical planar traveling-wave detection by phase-vector alignment across
% ELECTRODES, swept over frequency × speed × direction -- the cortical-space
% twin of stimulus_loc_traveling_wave.m.
%
% Electrode c sits at r_c = (x_c,y_c) mm on the 8x8 array with preferred-phase
% vector e^{iφ_c}. A planar wave of speed v and direction θ delays it by
%   Δt_c = d_c(θ)/v ,  d_c(θ) = x_c cosθ + y_c sinθ   (mm along the wave axis)
% i.e. rotates it by k·d_c(θ), k = 2πf/v (rad/mm). De-rotating makes the
% vectors overlap (R -> 1), but only at the true (f,v,θ), so we sweep the grid
% and read off where they align. A genuine wave gives high R along a
% near-constant speed across frequencies (k ∝ f); a chance per-frequency
% alignment does not.
%
% Unlike the stim-loc test, the wave direction on the 2-D array is unknown, so
% θ is swept too (equivalently, the 2-D wavevector is fitted). R is reported
% maxed over θ, with the best θ recorded.
%
% Coherent across electrodes: the complex vectors are summed keeping each
% electrode's own phase, because cortical electrodes share a common reference,
% so the per-electrode delay IS the wave and is directly testable:
%   R(f,v,θ) = | Σ_c w_c e^{i(φ_c − k d_c(θ))} | / Σ_c w_c ,  w_c = |Zmap|.
% Positions are collapsed into one map per frequency (magnitude-weighted, as in
% the collapsed PGD test).
%
% TWO READOUTS
%   (1) absolute R, max-stat corrected over the whole (f,v,θ) grid;
%   (2) de-rotation gain dR = R(f,v,θ) − R0(f), R0 being the unrotated (k=0)
%       coherence -- pure synchrony, no wave assumed. dR isolates the
%       wave-specific ramp and is invariant to electrodes merely sharing a
%       phase. Its null is clean because k=0 uses no d, so the
%       electrode<->position shuffle leaves R0 unchanged.
%
% Null: shuffle which array coordinate each (phase,weight) sits at, one
% synchronised shuffle shared across the grid, max over the grid -> max-stat
% threshold for both R and dR. Per animal throughout, never pooled across
% animals: grid -> animals -> mean grid + "significant in both" replication.
%
% Units: array coordinates in mm (pitch SPACING_MM), speed in cm/s (mm/s
% internally), k in rad/mm. SPEED_OK flags the plausible cortical band.
%
% TWO ESTIMATORS -- both run in one pass. Both build Zmap(c,f), one complex
% number per (electrode, frequency) whose angle is that electrode's preferred
% phase and whose magnitude is its weight; the de-rotation, null, gain,
% thresholds and pooling are then identical. They differ only in how Zmap is
% formed, and so in what |Zmap| means.
%   'phase'      magnitude-weighted circular sum of the per-location coherences
%                from phase_progression.mat,
%                z_c(f) = SUM_p coh_mag(c,p,f) e^{i pref_phase(c,p,f)}.
%                |z_c| is a resultant length over locations, not a coherence;
%                each term is a per-location mean, so a location with fewer
%                trials contributes an inflated magnitude, |c_p| ~ 1/sqrt(n_p).
%   'coherence'  one phase coherence over all trials,
%                z_c(f) = SUM_t y_t e^{i phi_tf} / SUM_t |y_t|, so |z_c| is
%                that channel's coherence in [0,1] and every trial counts once.
%                Built from the cached per-location sums as SUM_p S(c,p,f)/W_c
%                -- summing over p just undoes the grouping, since this script
%                does not use the stimulus locations at all.
%   R is not comparable between the two; judge each against its own threshold.
%
%   Switching estimator here does NOT lower R0, unlike in the stimulus-location
%   script. The final statistic is always a resultant ACROSS ELECTRODES,
%     R = | SUM_c w_c e^{i(phi_c - k d_c)} | / SUM_c w_c ,  w_c = |Zmap(c,f)| ,
%   whose normalisation divides |Zmap| straight back out, so the estimator
%   changes the weights and the angles but never the scale of R. Measured R0
%   medians: 0.947 -> 0.953 (hermes), 0.9916 -> 0.9916 (klecks). The
%   near-saturated baseline, and the small gain ceiling that follows from it,
%   is a property of this analysis and caveats both estimators.
%
%   Only the observed sums are needed: the null shuffles electrode <-> array
%   coordinate, not trial labels, so the cached worker's permutation part is
%   never loaded (~0.3 MB per animal instead of ~286 MB). The cache is shared
%   with stimulus_loc_traveling_wave.m and keeps the name 'trial_position_sums'.
%
% COMBINING ANIMALS
%   REPLICATION (primary)  significant in both animals, cell by cell. Cannot be
%     driven by one animal, cannot aggregate weak evidence, low power.
%   POOLED (secondary)     z-score each animal's grid against its own null,
%     average, and test against a paired-permutation max-stat null. More power,
%     but with two animals a hit can be almost entirely one of them, so the
%     per-animal z at the peak cell is printed alongside. Standardising first
%     is essential -- R0 ~0.95 vs ~0.99, so a raw average would track the
%     larger animal.
%   One pooled statistic, not two: R0 is constant across speed and across
%   permutations within a frequency, so it cancels in both numerator and null
%   and the standardised gain grid is identical to the standardised R grid
%   (verified numerically). The pooled column appears once, on the R figure.
%
% OUTPUT -- one shared R-grid figure plus a gain and a direction/speed figure
% per estimator, in Plots/scanning/planar_wave_derotation/cp10_till_100/<dv>/
%   derotation_R_grids.pdf        row 1 = 'phase', row 2 = 'coherence' (named
%                                 in every title; colour limits and thresholds
%                                 per row). Cols = animals + mean/replication
%                                 + pooled z.
%   derotation_gain_grids_<est>.pdf        dR, animals + mean/replication
%   derotation_direction_speed_<est>.pdf   best-fit direction & speed vs freq
%   results_combined/scanning/planar_wave_derotation/cp10_till_100/<dv>/
%     planar_wave_derotation.mat  results.G.(est)(animal), results.C.(est)
% =====================================================================

clearvars; close all; clc

%% ─── Dependencies (only the 'coherence' estimator needs slurmfun) ────
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

% Which estimators form the per-electrode complex map Zmap(c,f) -- see header.
% Both run in one pass: the R-grid figure shows them side by side and each gets
% its own gain figure. 'coherence' needs the cached SLURM sums (observed part).
ESTIMATORS = {'phase','coherence'};
RECOMPUTE_TRIAL_SUMS = true;     % force re-submission of the per-channel jobs

% Per-estimator naming used in filenames, panel titles and console output.
EST_INFO = struct( ...
    'phase',     struct('tag','_phase',     'short','PHASE ALIGNMENT', ...
        'long','phase estimator = PHASE ALIGNMENT: |Zmap| = resultant of the per-location preferred phases (NOT a coherence)'), ...
    'coherence', struct('tag','_coherence', 'short','PHASE COHERENCE', ...
        'long','coherence estimator = PHASE COHERENCE: |Zmap| = each channel''s coherence over all trials (the Phase_coherence/ measure)'));
for ei = 1:numel(ESTIMATORS)
    if ~isfield(EST_INFO, ESTIMATORS{ei}), error('Unknown ESTIMATOR: %s', ESTIMATORS{ei}); end
end
tag = '';   % no RF toggle in this script; the estimator goes in the gain names

out_dir = fullfile(base,'Plots','scanning','planar_wave_derotation','cp10_till_100', dv);
res_dir = fullfile(base,'results_combined','scanning','planar_wave_derotation','cp10_till_100', dv);
if ~exist(out_dir,'dir'), mkdir(out_dir); end
if ~exist(res_dir,'dir'), mkdir(res_dir); end

% electrode grid coordinates (mm); same index<->position map as the PGD script
ch_col = ceil((1:(grid_rows*grid_cols))' / grid_rows);
ch_row = grid_rows - mod((1:(grid_rows*grid_cols))' - 1, grid_rows);
XY = [ch_col, ch_row] * SPACING_MM;      % [nCh x 2] (x,y) in mm

G = struct();   % G.(estimator)(ia): grids, thresholds, sig, gain, best speed/dir

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

    % 'coherence' estimator: per-channel totals from the cached trial sums.
    % Loaded ONCE per animal and reused across estimators (cheap: observed only).
    if any(strcmp(ESTIMATORS,'coherence'))
        [Ssum, Wch] = get_trial_sums_obs(base, animalName, dv, nCh, nPerm, RECOMPUTE_TRIAL_SUMS);
    end

  for ei = 1:numel(ESTIMATORS)
    est = ESTIMATORS{ei};
    fprintf('  --- estimator: %s (%s) ---\n', est, EST_INFO.(est).short);

    % ── Per-electrode complex map Zmap(ch,freq): THE ONLY THING THE TWO
    %    ESTIMATORS DO DIFFERENTLY. Everything after this switch — the
    %    de-rotation, null, gain, thresholds, pooling — is identical.
    %      'phase'     |Zmap| = resultant length over locations (NOT a coherence)
    %      'coherence' |Zmap| = that channel's phase coherence over all trials
    switch est
        case 'phase'
            Zmap = complex(nan(nCh, nFreq));
            for f = 1:nFreq
                PHI = squeeze(S.pref_phase(:,f,:));  W = squeeze(S.coh_mag(:,f,:));
                M = isfinite(PHI) & isfinite(W) & (W > 0);  W(~M) = 0; PHI(~M) = 0;
                Zmap(:,f) = sum(W .* exp(1i*PHI), 2);      % magnitude-weighted circular sum
            end
            fprintf('    PHASE ALIGNMENT: |Zmap| range %.3f-%.3f (resultant length, not a coherence)\n', ...
                min(abs(Zmap(:)),[],'omitnan'), max(abs(Zmap(:)),[],'omitnan'));
        case 'coherence'
            % z_c(f) = SUM_p S(c,p,f) / W_c = coherence over ALL trials.
            % Summing over p just undoes the per-location grouping; this script
            % does not use the stimulus locations. Only the observed sums are
            % needed — the null here is spatial (electrode <-> coordinate).
            Zmap = squeeze(sum(Ssum, 2)) ./ Wch;           % nCh × nFreq
            fprintf('    PHASE COHERENCE: |Zmap| range %.3f-%.3f (per-channel coherence, in [0,1])\n', ...
                min(abs(Zmap(:)),[],'omitnan'), max(abs(Zmap(:)),[],'omitnan'));
        otherwise
            error('Unhandled ESTIMATOR: %s', est);
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

    G.(est)(ia).animal = animalName;  G.(est)(ia).estimator = est;
    G.(est)(ia).R = Robs; G.(est)(ia).thr = thr; G.(est)(ia).sig = sig;
    G.(est)(ia).R0 = R0(:); G.(est)(ia).gain = gain;
    G.(est)(ia).thr_gain = thr_gain; G.(est)(ia).sig_gain = sig_gain;
    G.(est)(ia).fHz = fHz; G.(est)(ia).speeds = V_CORTICAL; G.(est)(ia).DIRbest = DIRbest;
    G.(est)(ia).vbest = vbest_speed; G.(est)(ia).vbest_dir = vbest_dir; G.(est)(ia).vbest_sig = vbest_sig;

    % ── Standardise against this animal's OWN null, cell by cell ─────────
    % 'omitnan' matters: mean/std do NOT skip NaN, so one NaN permutation would
    % make that cell NaN for every permutation and the NaN would then spread
    % through the cross-animal sum (NaN + finite = NaN).
    mu = mean(Rnull,3,'omitnan');  sd = std(Rnull,0,3,'omitnan');
    G.(est)(ia).z_obs  = (Robs  - mu) ./ max(sd, eps);
    G.(est)(ia).z_null = (Rnull - mu) ./ max(sd, eps);
    fprintf('    finite cells in standardised grid = %d/%d\n', ...
        sum(isfinite(G.(est)(ia).z_obs(:))), numel(G.(est)(ia).z_obs));

    [rmax,ix] = max(Robs(:)); [fi,vi] = ind2sub(size(Robs), ix);
    fprintf('    peak R=%.4f at f=%.1f Hz, v=%.1f cm/s, dir=%.0f° | R0 %.4f-%.4f | thr=%.4f | sig cells=%d\n', ...
        rmax, fHz(fi), V_CORTICAL(vi), rad2deg(DIRbest(fi,vi)), ...
        min(R0,[],'omitnan'), max(R0,[],'omitnan'), thr, sum(sig(:)));
    [gmax,gix] = max(gain(:)); [gfi,gvi] = ind2sub(size(gain), gix);
    okspeed = V_CORTICAL(gvi) >= SPEED_OK(1) && V_CORTICAL(gvi) <= SPEED_OK(2);
    fprintf('    peak gain dR=%.4f at f=%.1f Hz, v=%.1f cm/s %s | thr_gain=%.4f | sig-gain cells=%d | freqs w/ sig gain=%d/%d\n', ...
        gmax, fHz(gfi), V_CORTICAL(gvi), ternary(okspeed,'(plausible)','(NON-PHYSICAL)'), ...
        thr_gain, sum(sig_gain(:)), sum(vbest_sig), numel(vbest_sig));
  end
end
valid = find(arrayfun(@(s) ~isempty(s.animal), G.(ESTIMATORS{1})));

%% ─── Combine animals (mean grid + replication + pooled z) ────────────
% Done SEPARATELY FOR EACH ESTIMATOR — the two never share a null or a
% threshold, because their R values live on different scales.
C = struct();
for ei = 1:numel(ESTIMATORS)
  est = ESTIMATORS{ei};
  Rsum = 0; Gsum = 0; Zsum = 0; Znull = 0;
  repl = true(size(G.(est)(valid(1)).R)); repl_gain = repl;
  for ia = valid
      Rsum = Rsum + G.(est)(ia).R;
      Gsum = Gsum + G.(est)(ia).gain;
      repl = repl & G.(est)(ia).sig;
      repl_gain = repl_gain & G.(est)(ia).sig_gain;
      Zsum  = Zsum  + G.(est)(ia).z_obs;
      Znull = Znull + G.(est)(ia).z_null;
  end
  nA = numel(valid);
  C.(est).R = Rsum/nA; C.(est).repl = repl;
  C.(est).gain = Gsum/nA; C.(est).repl_gain = repl_gain;
  C.(est).fHz = G.(est)(valid(1)).fHz; C.(est).speeds = V_CORTICAL;

  % POOLED standardised test. Permutation b of one animal is paired with
  % permutation b of the other — legitimate because the animals are independent,
  % so the pairing is arbitrary and samples the product null. max() omits NaN per
  % column, so skipped frequencies are harmless; an all-NaN null (no usable
  % channels anywhere) is caught explicitly rather than yielding thr = NaN.
  C.(est).Z  = Zsum/nA;
  maxZ = max(reshape(Znull/nA, [], nPerm), [], 1).';
  tZ   = quantile(maxZ, 1-alpha);
  if ~isfinite(tZ)
      warning(['%s: pooled null is entirely NaN — at least one animal produced no finite ' ...
               'cells (too few usable channels?). Pooled test disabled; replication unaffected.'], est);
      tZ = Inf;
  end
  C.(est).thr_pool = tZ;  C.(est).sig_pool = C.(est).Z >= tZ;
  C.(est).pooled_ok = isfinite(tZ);

  fprintf('\n========= COMBINED | estimator: %s (%s) =========\n', est, EST_INFO.(est).short);
  fprintf('cells significant in BOTH animals (absolute R)       = %d\n', sum(repl(:)));
  fprintf('cells gain-significant in BOTH animals (increase-R)   = %d\n', sum(repl_gain(:)));
  fprintf('cells significant in the POOLED standardised test     = %d (thr z = %.2f)\n', ...
      sum(C.(est).sig_pool(:)), C.(est).thr_pool);
  if any(C.(est).sig_pool(:))
      [~, ix] = max(C.(est).Z(:) .* double(C.(est).sig_pool(:)));
      [pf, pv] = ind2sub(size(C.(est).Z), ix);
      za = arrayfun(@(ia) G.(est)(ia).z_obs(pf,pv), valid);
      fprintf('  peak pooled cell f=%.2f Hz v=%.1f cm/s: per-animal z = [%s] (%s)\n', ...
          C.(est).fHz(pf), C.(est).speeds(pv), strjoin(compose('%.2f', za), ', '), ...
          strjoin({G.(est)(valid).animal}, ', '));
  end
end

%% ─── FIGURE 1: the DE-ROTATION R grids, both estimators ──────────────
% R maxed over direction; significant cells outlined. A wave = high R along a
% near-constant speed across frequencies. Rows = estimators; R is not
% comparable between them, so colour limits are per row (a fixed [0 1] scale
% renders the phase-coherence row as flat blue).
esc  = @(s) strrep(s, '_', '\_');   % TeX renders '_' as a subscript otherwise
nEst = numel(ESTIMATORS);
ncol = numel(valid)+2;              % animals + mean/replication + pooled z
f1 = figure('Visible','off','Position',[40 40 380*ncol 340*nEst]);
for ei = 1:nEst
  est = ESTIMATORS{ei};
  rmax_e = max(C.(est).R(:),[],'omitnan');
  for ia = valid, rmax_e = max(rmax_e, max(G.(est)(ia).R(:),[],'omitnan')); end
  if ~isfinite(rmax_e) || rmax_e <= 0, rmax_e = 1; end
  for k = 1:numel(valid)
      ia = valid(k);
      ax = subplot(nEst, ncol, (ei-1)*ncol + k); hold(ax,'on');
      imagesc(ax, G.(est)(ia).fHz, 1:numel(G.(est)(ia).speeds), G.(est)(ia).R.');
      set(ax,'YDir','normal'); axis(ax,'tight');
      contour(ax, G.(est)(ia).fHz, 1:numel(G.(est)(ia).speeds), double(G.(est)(ia).sig.'), [0.5 0.5], 'w','LineWidth',1.2);
      speed_yticks(ax, G.(est)(ia).speeds, SPEED_OK);
      caxis(ax,[0 rmax_e]); colorbar(ax);
      xlabel(ax,'Frequency (Hz)'); ylabel(ax,'wave speed (cm/s)');
      title(ax, {sprintf('%s — R (sig outlined)', G.(est)(ia).animal), ...
                 sprintf('[%s]', EST_INFO.(est).short)},'FontSize',9);
  end
  ax = subplot(nEst, ncol, (ei-1)*ncol + ncol-1); hold(ax,'on');
  imagesc(ax, C.(est).fHz, 1:numel(C.(est).speeds), C.(est).R.');
  set(ax,'YDir','normal'); axis(ax,'tight');
  contour(ax, C.(est).fHz, 1:numel(C.(est).speeds), double(C.(est).repl.'), [0.5 0.5], 'w','LineWidth',1.4);
  speed_yticks(ax, C.(est).speeds, SPEED_OK);
  caxis(ax,[0 rmax_e]); colorbar(ax);
  xlabel(ax,'Frequency (Hz)'); ylabel(ax,'wave speed (cm/s)');
  title(ax, {sprintf('mean R  [%s]', EST_INFO.(est).short), 'contour = replication'},'FontSize',9);

  % POOLED standardised z. Units are z, NOT R — auto-scaled, so its colours are
  % not comparable with the panels to the left. One pooled panel per estimator:
  % the standardised gain grid is identical to this one (see header).
  ax = subplot(nEst, ncol, (ei-1)*ncol + ncol); hold(ax,'on');
  imagesc(ax, C.(est).fHz, 1:numel(C.(est).speeds), C.(est).Z.');
  set(ax,'YDir','normal'); axis(ax,'tight');
  contour(ax, C.(est).fHz, 1:numel(C.(est).speeds), double(C.(est).sig_pool.'), [0.5 0.5], 'w','LineWidth',1.4);
  speed_yticks(ax, C.(est).speeds, SPEED_OK);
  if ~any(isfinite(C.(est).Z(:))), caxis(ax,[0 1]); end
  colorbar(ax);
  xlabel(ax,'Frequency (Hz)'); ylabel(ax,'wave speed (cm/s)');
  if C.(est).pooled_ok, ttl_p = sprintf('POOLED z (thr %.2f)', C.(est).thr_pool);
  else,                 ttl_p = 'POOLED z - n/a (see warning)'; end
  title(ax, {sprintf('pooled standardised  [%s]', EST_INFO.(est).short), ttl_p},'FontSize',9);
end
sgtitle({'DE-ROTATION R: coherence/alignment after de-rotating electrodes by k*d(theta), k = 2*pi*f/v', ...
         'ROW 1 = PHASE ALIGNMENT estimator (|Zmap| = resultant of the per-location preferred phases — NOT a coherence)', ...
         'ROW 2 = PHASE COHERENCE estimator (|Zmap| = each channel''s coherence over all trials — the Phase\_coherence/ measure)', ...
         'R IS NOT COMPARABLE BETWEEN ROWS: colour limits and thresholds are per row.', ...
         'white contour = significant   |   col 3 = mean, contour = replication   |   col 4 = pooled standardised z (z units, auto-scaled)'}, ...
         'FontSize',9);
set(f1,'PaperPositionMode','auto'); pos=get(f1,'Position'); set(f1,'PaperUnits','points','PaperSize',pos(3:4));
saveas(f1, fullfile(out_dir, ['derotation_R_grids' tag '.pdf']));

%% ─── FIGURES 2..n: INCREASE IN COHERENCE (gain) dR, ONE PER ESTIMATOR ─
% dR(f,v) = R(f,v) - R0(f): speeds at which de-rotating beats no rotation. A
% wave = significant positive gain (contour) along a near-constant best speed
% (black line); a rising line means constant wavenumber, i.e. a fixed phase
% offset rather than propagation.
% No pooled column -- the pooled gain z equals the pooled R z (see header).
ncol_g = numel(valid)+1;    % animals + mean
for ei = 1:nEst
  est = ESTIMATORS{ei};
  gmax_m = max(C.(est).gain(:),[],'omitnan');
  for ia = valid, gmax_m = max(gmax_m, max(G.(est)(ia).gain(:),[],'omitnan')); end
  if ~isfinite(gmax_m) || gmax_m <= 0, gmax_m = eps; end
  fg = figure('Visible','off','Position',[40 40 380*ncol_g 360]);
  for k = 1:numel(valid)
      ia = valid(k);
      ax = subplot(1, ncol_g, k); hold(ax,'on');
      imagesc(ax, G.(est)(ia).fHz, 1:numel(G.(est)(ia).speeds), G.(est)(ia).gain.');
      set(ax,'YDir','normal'); axis(ax,'tight');
      contour(ax, G.(est)(ia).fHz, 1:numel(G.(est)(ia).speeds), double(G.(est)(ia).sig_gain.'), [0.5 0.5], 'w','LineWidth',1.2);
      plot(ax, G.(est)(ia).fHz, interp1(G.(est)(ia).speeds, 1:numel(G.(est)(ia).speeds), G.(est)(ia).vbest, 'linear', NaN), 'k-','LineWidth',1.3);
      speed_yticks(ax, G.(est)(ia).speeds, SPEED_OK);
      caxis(ax,[0 gmax_m]); colorbar(ax);
      xlabel(ax,'Frequency (Hz)'); ylabel(ax,'wave speed (cm/s)');
      title(ax, {G.(est)(ia).animal, sprintf('gain dR (R0 %.3f-%.3f)', ...
          min(G.(est)(ia).R0,[],'omitnan'), max(G.(est)(ia).R0,[],'omitnan'))},'FontSize',9);
  end
  ax = subplot(1, ncol_g, ncol_g); hold(ax,'on');
  imagesc(ax, C.(est).fHz, 1:numel(C.(est).speeds), C.(est).gain.');
  set(ax,'YDir','normal'); axis(ax,'tight');
  contour(ax, C.(est).fHz, 1:numel(C.(est).speeds), double(C.(est).repl_gain.'), [0.5 0.5], 'w','LineWidth',1.4);
  speed_yticks(ax, C.(est).speeds, SPEED_OK);
  caxis(ax,[0 gmax_m]); colorbar(ax);
  xlabel(ax,'Frequency (Hz)'); ylabel(ax,'wave speed (cm/s)');
  title(ax, {'mean gain', 'contour = replication'},'FontSize',9);
  sgtitle({sprintf('DE-ROTATION GAIN dR = R(f,v) - R0(f)   |   %s', upper(EST_INFO.(est).short)), ...
           esc(EST_INFO.(est).long), ...
           'R0 = the same statistic at k=0 (nothing de-rotated). dR asks: does assuming a wave buy any alignment/coherence over assuming none?', ...
           'white contour = significant   |   black line = best-fit speed   |   FLAT line = real wave, RISING = constant wavenumber (not a wave)', ...
           'no pooled column: the pooled gain z is identical to the pooled R z, drawn once on the grids figure'}, ...
           'FontSize',9);
  set(fg,'PaperPositionMode','auto'); pos=get(fg,'Position'); set(fg,'PaperUnits','points','PaperSize',pos(3:4));
  saveas(fg, fullfile(out_dir, ['derotation_gain_grids' EST_INFO.(est).tag tag '.pdf']));
end

%% ─── Figure 3: best-fit propagation DIRECTION vs frequency ────────────
% Direction of the best-fit plane wave at each frequency, taken at the peak-R
% speed. Grey = all frequencies, filled = frequencies with significant gain. A
% genuine wave = a consistent direction across the significant band; scattered
% directions mean no coherent wave. Bottom row: cross-animal agreement per
% frequency. One figure per estimator, since direction and speed are read off
% that estimator's own grid.
cols = lines(numel(animals));
for ei = 1:nEst
  est = ESTIMATORS{ei};
  f3 = figure('Visible','off','Position',[40 40 380*numel(valid) 580]);
  for k = 1:numel(valid)
    ia = valid(k);
    fr = G.(est)(ia).fHz; dg = mod(rad2deg(G.(est)(ia).vbest_dir(:)), 360); sg = G.(est)(ia).vbest_sig(:);
    % top: direction vs frequency (deg), significant highlighted
    ax = subplot(2, numel(valid), k); hold(ax,'on');
    plot(ax, fr, dg, 'o', 'Color',[.75 .75 .75], 'MarkerSize',4, 'DisplayName','all freqs');
    if any(sg)
        plot(ax, fr(sg), dg(sg), 'o', 'Color',cols(ia,:), 'MarkerFaceColor',cols(ia,:), ...
            'MarkerSize',6, 'DisplayName','sig gain');
        Rdir = abs(mean(exp(1i*G.(est)(ia).vbest_dir(sg))));    % direction consistency 0..1
        ttl = sprintf('%s — prop. direction (sig: R_{dir}=%.2f)', G.(est)(ia).animal, Rdir);
    else
        ttl = sprintf('%s — prop. direction (no sig band)', G.(est)(ia).animal);
    end
    ylim(ax,[0 360]); set(ax,'YTick',0:90:360);
    xlabel(ax,'Frequency (Hz)'); ylabel(ax,'best-fit direction (deg)');
    title(ax, ttl,'FontSize',9); grid(ax,'on'); legend(ax,'Location','best','FontSize',7);
    % bottom: best speed vs frequency (constant speed => real wave), sig marked
    ax2 = subplot(2, numel(valid), numel(valid)+k); hold(ax2,'on');
    plot(ax2, fr, G.(est)(ia).vbest(:), 'o', 'Color',[.75 .75 .75], 'MarkerSize',4);
    if any(sg), plot(ax2, fr(sg), G.(est)(ia).vbest(sg), 'o', 'Color',cols(ia,:),'MarkerFaceColor',cols(ia,:),'MarkerSize',6); end
    yline(ax2, SPEED_OK(1), 'k:'); yline(ax2, SPEED_OK(2), 'k:');
    set(ax2,'YScale','log'); ylim(ax2,[V_CORTICAL(1) V_CORTICAL(end)]);
    xlabel(ax2,'Frequency (Hz)'); ylabel(ax2,'best-fit speed (cm/s)');
    title(ax2,'best speed (flat across freq = real wave; dotted = plausible band)','FontSize',8);
    grid(ax2,'on');
  end
  sgtitle({sprintf('Best-fit planar-wave direction & speed vs frequency   |   %s', upper(EST_INFO.(est).short)), ...
           'filled = significant increase in coherence   |   FLAT speed across frequency = real wave, RISING = constant wavenumber'}, ...
           'FontSize',10);
  set(f3,'PaperPositionMode','auto'); pos=get(f3,'Position'); set(f3,'PaperUnits','points','PaperSize',pos(3:4));
  saveas(f3, fullfile(out_dir, ['derotation_direction_speed' EST_INFO.(est).tag tag '.pdf']));
end

% drop the per-permutation null grids before saving (~100 MB otherwise)
for ei = 1:nEst
    for ia = valid, G.(ESTIMATORS{ei})(ia).z_null = []; end
end

% BOTH estimators live in one file: results.G.(estimator)(animal) and
% results.C.(estimator).
results = struct('G',G,'C',C,'animals',{animals},'dv',dv, ...
    'FREQ_RANGE',FREQ_RANGE,'V_CORTICAL',V_CORTICAL,'THETA',THETA, ...
    'SPACING_MM',SPACING_MM,'SPEED_OK',SPEED_OK,'nPerm',nPerm,'alpha',alpha, ...
    'ESTIMATORS',{ESTIMATORS},'EST_INFO',EST_INFO);
save(fullfile(res_dir, ['planar_wave_derotation' tag '.mat']),'results','-v7.3');
fprintf('\nSaved under %s :\n  derotation_R_grids%s.pdf   (de-rotation R, both estimators)\n', out_dir, tag);
for ei = 1:nEst
    fprintf('  derotation_gain_grids%s%s.pdf / derotation_direction_speed%s%s.pdf   (%s)\n', ...
        EST_INFO.(ESTIMATORS{ei}).tag, tag, EST_INFO.(ESTIMATORS{ei}).tag, tag, ...
        EST_INFO.(ESTIMATORS{ei}).short);
end

%% =====================================================================
%% Helpers
%% =====================================================================
function [Robs, R0, Rnull_max, Gnull_max, DIRbest, Rnull, Gnull] = derotate_grid(Zmap, coh_sig, f_use, fHz, XY, Vs_mm, Th, MIN_CH, nPerm)
% Coherent planar-wave plane fit by de-rotation over the (frequency × speed ×
% direction) grid.
%   R(f,v,θ) = | Σ_c w_c e^{i(φ_c − k d_c(θ)) } | / Σ_c w_c ,  k = 2πf/v_mm,
%              d_c(θ) = x_c cosθ + y_c sinθ  (mm).
% Positions are collapsed per frequency (magnitude-weighted circular sum).
% Returns R maxed over direction, the k=0 baseline R0, per-perm grid maxima of
% R and of dR, and the best direction per (f,v). The null permutes each
% (phase,weight) to a random array coordinate; R0 is shuffle-invariant, since
% k=0 uses no coordinate.
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
% for the 'coherence' estimator. One SLURM job per channel via
% functions/trial_position_sums_chan.m, cached and shared with
% stimulus_loc_traveling_wave.m. Only S_obs is read: this script's null shuffles
% electrode <-> array coordinate, so the large permutation part is never loaded.
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
% phase_coherence/complex pipeline. Per channel the max-stat threshold is the
% (1-alpha) quantile of the per-perm max over frequency of |coh_perm_complex|;
% the channel is significant where |coh_complex(f)| exceeds it. Identical to the
% mask used by cortical_planar_wave_PGD.m and stimulus_loc_traveling_wave.m.
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
