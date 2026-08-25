% =====================================================================
% Traveling wave & spatial phase analysis
%
% Tests whether the position-dependent preferred phase across channels is
% consistent with a wave propagating across cortex, and characterises the
% spatial structure of the phase pattern.
%
% SUPERSEDED: thresholds here come from the H2-H1 test, whereas coherence
% significance is the better criterion. Prefer cortical_planar_wave_PGD.m and
% cortical_wave_type_classification.m.
%
% Runs per (animal x DV). CH_FILTER: 'all' = all 64 channels, 'sig_only' =
% only channels significant in the H2-H1 test.
%
% NINE-PAGE PDF per (animal, DV)
%
% Spatial structure on the array
%   Page 1   Phase heatmap on the 8x8 grid, one cell per channel (HSV, -pi to
%            pi); white squares mark driven channels (RF overlaps stimulus). A
%            smooth colour gradient = spatially organised phase.
%   Page 4   The same heatmap for three representative positions (peripheral,
%            middle, foveal), with the Page-3 gradient arrow overlaid.
%   Page 5   Moran's I on cos(phase), queen-contiguity (8-connected) weights:
%            positive = neighbouring electrodes have similar phases, near 0 =
%            random. p from 1000 shuffles of phase across electrode locations.
%            Tests spatial organisation without assuming a wave model.
%
% Where and how does phase change with position?
%   Page 2A  Circular phase difference (deg) between two positions on the grid
%            -- peripheral vs foveal, peripheral vs middle, middle vs foveal.
%            Shows directly which channels change phase.
%   Page 2B  Circular variance of phase across positions, per channel (0 =
%            identical at every position, 1 = uniformly spread), on the grid;
%            plus a histogram for sig vs non-sig channels. Higher variance in
%            the sig channels validates H2-H1.
%   Page 2C  Phase trajectory per significant channel: x = position
%            (peripheral to foveal), y = unwrapped preferred phase aligned to 0
%            at position 1. Flat = no change, diverging = positions disagree,
%            similar trajectories in neighbours = a spatially coherent effect.
%   Page 6   Phase vs retinotopic distance (sig channels only): distance from
%            RF centre to the estimated stimulus location, i.e. does phase
%            follow the visual field map?
%
% Phase gradient characterisation
%   Page 3   2-D plane fit per position, phase = a*row + b*col + c on the
%            2D-unwrapped phases. R2 per position, gradient magnitude
%            sqrt(a^2+b^2) in deg/electrode, and the gradient direction (a,b)
%            as arrows, showing how it rotates across positions.
%
% Output:
%   Plots/scanning/phase_progression/cp10_till_100/<dv>/
%       <animal>_traveling_wave_<freq>[_sigonly].pdf
% =====================================================================
clearvars; close all; clc

%% ═══════════════════════════════════════════════════════════════════
%% USER SETTINGS — edit these to control what runs
%% ═══════════════════════════════════════════════════════════════════

% Which animals to run: {'klecks'}, {'hermes'}, or {'klecks','hermes'}
animals    = {'klecks','hermes'};

% Which DVs to run: {'lfp'}, {'mua'}, or {'lfp','mua'}
dv_list    = {'lfp','mua'};

% Frequency mode:
%   'single'  — use one frequency (set by single_freq below)
%   'band'    — average phase across a frequency band (set by band_range)
%               and mark channels significant at ANY freq in the band
FREQ_MODE   = 'single';        % <-- 'single' or 'band'
single_freq = 10;             % Hz, used when FREQ_MODE = 'single'
band_range  = [8 12];         % [lo hi] Hz, used when FREQ_MODE = 'band'

% Channel filter:
%   'all'      — use all 64 channels (original behavior)
%   'sig_only' — keep only channels significant in H2-H1 test;
%                non-significant channels are blanked (NaN) in all plots
CH_FILTER   = 'sig_only';    % <-- 'all' or 'sig_only'

%% ═══════════════════════════════════════════════════════════════════

dv_key_map = struct('lfp','LFP_ERP_ampl_all','mua','MUA_ERP_ampl_all');
nCh        = 64;
grid_rows  = 8;
grid_cols  = 8;
base       = '/mnt/hpc/projects/MWSampling/4Shivangi';
alpha      = 0.05;

for dv_idx = 1:numel(dv_list)
dv = dv_list{dv_idx};
dv_key = dv_key_map.(dv);

for a_idx = 1:numel(animals)
animalName = animals{a_idx};
fprintf('\n=== %s / %s / mode=%s ===\n', animalName, upper(dv), FREQ_MODE);

%% Load phase progression data
pp_file = fullfile(base, ['results_' animalName], 'scanning', ...
    'phase_progression','cp10_till_100', dv, 'phase_progression.mat');
if ~isfile(pp_file)
    warning('No phase progression data for %s/%s — skipping.', animalName, dv);
    continue
end
S = load(pp_file, 'pref_phase','coh_mag','freq','positions');
freq      = S.freq;
positions = S.positions;
nPos      = numel(positions);
nFreq     = numel(freq);
nCh       = size(S.pref_phase, 1);

% Rebuild grid mapping for actual nCh
ch_col = ceil((1:nCh)' / grid_rows);
ch_row = grid_rows - mod((1:nCh)' - 1, grid_rows);

%% Frequency selection + phase extraction
switch FREQ_MODE
    case 'single'
        [~, fi] = min(abs(freq - single_freq));
        f_hz = freq(fi);
        freq_label = sprintf('%.0fHz', f_hz);
        fprintf('Single frequency: index %d = %.1f Hz\n', fi, f_hz);

        phi = squeeze(S.pref_phase(:, fi, :));    % nCh × nPos
        mag = squeeze(S.coh_mag(:, fi, :));

        % H2-H1 significance at this single frequency
        sig_at_f = get_h2h1_sig_single(base, animalName, dv_key, nCh, nFreq, fi, alpha);

    case 'band'
        band_idx = find(freq >= band_range(1) & freq <= band_range(2));
        freq_label = sprintf('%g-%gHz', band_range(1), band_range(2));
        fprintf('Band mode: %d frequencies in %g–%g Hz\n', numel(band_idx), band_range(1), band_range(2));

        % Average phase across band using circular mean of complex vectors
        % per (channel, position)
        z_band = mean(S.coh_mag(:, band_idx, :) .* ...
                      exp(1i * S.pref_phase(:, band_idx, :)), 2);  % nCh × 1 × nPos
        phi = angle(reshape(z_band, [], nPos));   % nCh × nPos
        mag = abs(reshape(z_band, [], nPos));

        % Channel is significant if H2-H1 passes at ANY freq in band
        sig_at_f = get_h2h1_sig_band(base, animalName, dv_key, nCh, nFreq, band_idx, alpha);

    otherwise
        error('Unknown FREQ_MODE: %s', FREQ_MODE);
end

fprintf('  %d channels significant (H2-H1)\n', sum(sig_at_f));

% Keep both versions for pages that need all-ch vs sig-only comparison
phi_all = phi;
phi_sig = phi;
phi_sig(~sig_at_f, :) = NaN;

%% Apply channel filter
if strcmp(CH_FILTER, 'sig_only')
    phi(~sig_at_f, :) = NaN;
    mag(~sig_at_f, :) = NaN;
    ch_label = 'sig only';
    ch_suffix = '_sigonly';
    fprintf('  CH_FILTER=sig_only: blanked %d non-significant channels\n', sum(~sig_at_f));
else
    ch_label = 'all ch';
    ch_suffix = '';
end

%% Load RF data — identify driven channels per position
rf_sessions = struct('klecks','klecks_20170829_rfmapping_bar_1_channel_target_summary.txt', ...
                     'hermes','hermes_20170829_rfmapping_bar_1_channel_target_summary.txt');
rf_file = fullfile(base, 'Plots','RF_Mapping', animalName, 'loc_RF_map', ...
    'gaussian_overlap', rf_sessions.(animalName));
rf_tab = readtable(rf_file, 'Delimiter','\t');

rf_x    = rf_tab.RF_Center_X;
rf_y    = rf_tab.RF_Center_Y;
rf_fwhm = rf_tab.FWHM_radius;

% Parse Locations_Inside to get which positions each channel covers
% Locations are numbered 1..nPos (matching the position ordering)
ch_covers_pos = false(nCh, nPos);
for ch = 1:min(nCh, height(rf_tab))
    loc_str = rf_tab.Locations_Inside{ch};
    if strcmp(loc_str, 'none') || isempty(loc_str), continue; end
    locs = str2double(strsplit(loc_str, ','));
    locs = locs(~isnan(locs) & locs >= 1 & locs <= nPos);
    ch_covers_pos(ch, locs) = true;
end

% Channels with ANY stimulus overlap (Group A / driven)
driven_chs = find(any(ch_covers_pos, 2));
fprintf('Driven channels (RF overlaps stimuli): ');
fprintf('%d ', driven_chs); fprintf('\n');

%% Output setup
plot_dir = fullfile(base, 'Plots','scanning','phase_progression', ...
    'cp10_till_100', dv);
if ~exist(plot_dir,'dir'), mkdir(plot_dir); end
out_pdf = fullfile(plot_dir, sprintf('%s_traveling_wave_%s%s.pdf', animalName, freq_label, ch_suffix));
out_ps  = strrep(out_pdf, '.pdf', '.ps');
if isfile(out_pdf), delete(out_pdf); end
if isfile(out_ps),  delete(out_ps);  end

%% ═══════════════════════════════════════════════════════════════════
%% Page 1: Phase heatmap on the 8×8 array grid, per position
%% ═══════════════════════════════════════════════════════════════════
cols_pg1 = min(nPos, 5);
rows_pg1 = ceil(nPos / cols_pg1);

fig1 = figure('Visible','off','Units','centimeters', ...
    'Position',[1 1 max(24, 5*cols_pg1) max(16, 5*rows_pg1+2)]);
set(fig1,'PaperUnits','centimeters','PaperSize',fig1.Position(3:4), ...
    'PaperPosition',[0 0 fig1.Position(3:4)]);

for p = 1:nPos
    ax = subplot(rows_pg1, cols_pg1, p);

    % Build 8×8 phase image
    phase_grid = nan(grid_rows, grid_cols);
    for ch = 1:nCh
        phase_grid(ch_row(ch), ch_col(ch)) = phi(ch, p);
    end

    h_img = imagesc(ax, 1:grid_cols, 1:grid_rows, phase_grid);
    set(h_img, 'AlphaData', ~isnan(phase_grid));
    set(ax, 'YDir','normal', 'Color',[0.9 0.9 0.9]);
    colormap(ax, hsv);
    caxis(ax, [-pi pi]);
    set(ax, 'XTick', 1:grid_cols, 'YTick', 1:grid_rows, 'FontSize', 6);
    xlabel(ax, 'Col', 'FontSize', 6);
    ylabel(ax, 'Row', 'FontSize', 6);
    axis(ax, 'square');

    % Mark driven channels for this position with black squares
    for ch = 1:nCh
        if ch_covers_pos(ch, p)
            hold(ax, 'on');
            plot(ax, ch_col(ch), ch_row(ch), 'ks', 'MarkerSize', 10, ...
                'LineWidth', 1.5, 'MarkerFaceColor', 'none');
        end
    end

    title(ax, sprintf('pos %g', positions(p)), 'FontSize', 8);
end

cb = colorbar('Position',[0.93 0.15 0.015 0.7]);
cb.Ticks = [-pi -pi/2 0 pi/2 pi];
cb.TickLabels = {'-\pi','-\pi/2','0','\pi/2','\pi'};
cb.FontSize = 7;

sgtitle(sprintf('%s — %s — %s [%s] :: preferred phase on array grid\n(black squares = RF overlaps this position)', ...
    animalName, upper(dv), freq_label, ch_label), 'FontSize', 10, 'FontWeight','bold');

print(fig1, '-dpsc', '-append', '-painters', out_ps);
close(fig1);
fprintf('Page 1 done\n');

%% ═══════════════════════════════════════════════════════════════════
%% Page 2A: Phase change map — peripheral vs foveal on the 8x8 grid
%% ═══════════════════════════════════════════════════════════════════
fig2a = figure('Visible','off','Units','centimeters', ...
    'Position',[1 1 30 10]);
set(fig2a,'PaperUnits','centimeters','PaperSize',fig2a.Position(3:4), ...
    'PaperPosition',[0 0 fig2a.Position(3:4)]);

% Pick position pairs: first vs last, first vs middle, middle vs last
pair_idx = [1 nPos; 1 ceil(nPos/2); ceil(nPos/2) nPos];
diff_grids = cell(3,1);
for k = 1:3
    p1 = pair_idx(k,1); p2 = pair_idx(k,2);
    dg = nan(grid_rows, grid_cols);
    for ch = 1:nCh
        ph1 = phi(ch, p1); ph2 = phi(ch, p2);
        if ~isnan(ph1) && ~isnan(ph2)
            dg(ch_row(ch), ch_col(ch)) = rad2deg(angle(exp(1i*(ph2 - ph1))));
        end
    end
    diff_grids{k} = dg;
end
global_mx = max(cellfun(@(g) max(abs(g(:)),[],'omitnan'), diff_grids));
if ~isfinite(global_mx) || global_mx == 0, global_mx = 180; end

for k = 1:3
    p1 = pair_idx(k,1); p2 = pair_idx(k,2);
    ax = subplot(1,3,k); hold(ax,'on');
    h_img = imagesc(ax, 1:grid_cols, 1:grid_rows, diff_grids{k});
    set(h_img, 'AlphaData', ~isnan(diff_grids{k}));
    set(ax, 'YDir','normal', 'Color',[0.9 0.9 0.9]);
    colormap(ax, parula);
    caxis(ax, [-global_mx global_mx]);
    colorbar(ax, 'FontSize', 6);
    set(ax, 'XTick', 1:grid_cols, 'YTick', 1:grid_rows, 'FontSize', 6);
    xlabel(ax, 'Col', 'FontSize', 7); ylabel(ax, 'Row', 'FontSize', 7);
    axis(ax, 'square');
    title(ax, sprintf('pos %g -> %g', positions(p1), positions(p2)), ...
        'FontSize', 9, 'FontWeight','bold');
end

sgtitle(sprintf('%s -- %s -- %s [%s] :: phase change between positions (deg)\n(bright = large phase shift; gray = no data)', ...
    animalName, upper(dv), freq_label, ch_label), 'FontSize', 10, 'FontWeight','bold');

print(fig2a, '-dpsc', '-append', '-painters', out_ps);
close(fig2a);
fprintf('Page 2A done (phase change map)\n');

%% ═══════════════════════════════════════════════════════════════════
%% Page 2B: Circular variance map — which channels change phase?
%% ═══════════════════════════════════════════════════════════════════
fig2b = figure('Visible','off','Units','centimeters', ...
    'Position',[1 1 24 10]);
set(fig2b,'PaperUnits','centimeters','PaperSize',fig2b.Position(3:4), ...
    'PaperPosition',[0 0 fig2b.Position(3:4)]);

% Left: circular variance on the 8x8 grid
ax1 = subplot(1,2,1); hold(ax1,'on');
cvar_grid = nan(grid_rows, grid_cols);
cvar_vec  = nan(nCh, 1);
for ch = 1:nCh
    ph_ch = phi(ch, :);
    valid_ph = ph_ch(~isnan(ph_ch));
    if numel(valid_ph) >= 2
        R = abs(mean(exp(1i * valid_ph)));
        cvar_vec(ch) = 1 - R;
        cvar_grid(ch_row(ch), ch_col(ch)) = cvar_vec(ch);
    end
end

h_img = imagesc(ax1, 1:grid_cols, 1:grid_rows, cvar_grid);
set(h_img, 'AlphaData', ~isnan(cvar_grid));
set(ax1, 'YDir','normal', 'Color',[0.9 0.9 0.9]);
colormap(ax1, hot);
caxis(ax1, [0 1]);
cb = colorbar(ax1, 'FontSize', 6);
cb.Label.String = 'Circular variance';
set(ax1, 'XTick', 1:grid_cols, 'YTick', 1:grid_rows, 'FontSize', 6);
xlabel(ax1, 'Col', 'FontSize', 7); ylabel(ax1, 'Row', 'FontSize', 7);
axis(ax1, 'square');
title(ax1, 'Circular variance of phase across positions', 'FontSize', 9, 'FontWeight','bold');

% Right: histogram of circular variance, sig vs non-sig
ax2 = subplot(1,2,2); hold(ax2,'on');
cv_sig = cvar_vec(sig_at_f);
cv_nonsig = cvar_vec(~sig_at_f);
cv_sig = cv_sig(~isnan(cv_sig));
cv_nonsig = cv_nonsig(~isnan(cv_nonsig));
edges = linspace(0, 1, 15);
if ~isempty(cv_nonsig)
    histogram(ax2, cv_nonsig, edges, 'FaceColor', [0.6 0.6 0.6], ...
        'EdgeColor','none', 'FaceAlpha', 0.6, 'DisplayName', 'non-sig ch');
end
if ~isempty(cv_sig)
    histogram(ax2, cv_sig, edges, 'FaceColor', [0.85 0.25 0.20], ...
        'EdgeColor','none', 'FaceAlpha', 0.8, 'DisplayName', 'sig ch (H2-H1)');
end
xlabel(ax2, 'Circular variance', 'FontSize', 8);
ylabel(ax2, 'Channel count', 'FontSize', 8);
title(ax2, 'Sig vs non-sig channels', 'FontSize', 9, 'FontWeight','bold');
legend(ax2, 'Location','best', 'FontSize', 7);
grid(ax2, 'on'); box(ax2, 'on');

sgtitle(sprintf('%s -- %s -- %s [%s] :: phase variability across positions\n(high circ. variance = phase depends on position = H2 channel)', ...
    animalName, upper(dv), freq_label, ch_label), 'FontSize', 10, 'FontWeight','bold');

print(fig2b, '-dpsc', '-append', '-painters', out_ps);
close(fig2b);
fprintf('Page 2B done (circular variance map)\n');

%% ═══════════════════════════════════════════════════════════════════
%% Page 2C: Phase trajectory per significant channel across positions
%% ═══════════════════════════════════════════════════════════════════
sig_ch_2c = find(sig_at_f);
nSig_2c   = numel(sig_ch_2c);

fig2c = figure('Visible','off','Units','centimeters','Position',[1 1 22 14]);
set(fig2c,'PaperUnits','centimeters','PaperSize',fig2c.Position(3:4), ...
    'PaperPosition',[0 0 fig2c.Position(3:4)]);

ax = axes('Parent', fig2c); hold(ax, 'on');

if nSig_2c >= 1
    ch_colors = lines(nSig_2c);
    x = 1:nPos;

    for k = 1:nSig_2c
        ch = sig_ch_2c(k);
        ph_ch = phi_sig(ch, :);
        if all(isnan(ph_ch)), continue; end
        valid_idx = ~isnan(ph_ch);
        ph_uw = unwrap(ph_ch(valid_idx));
        ph_aligned = ph_uw - ph_uw(1);
        plot(ax, x(valid_idx), rad2deg(ph_aligned), '-o', ...
            'Color', ch_colors(k,:), 'LineWidth', 1.2, ...
            'MarkerSize', 4, 'MarkerFaceColor', ch_colors(k,:), ...
            'DisplayName', sprintf('ch%d', ch));
    end

    yline(ax, 0, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5, ...
        'HandleVisibility','off');

    set(ax, 'XTick', x, 'XTickLabel', ...
        arrayfun(@(p) sprintf('%g', positions(p)), x, 'UniformOutput', false));
    set(ax, 'XTickLabelRotation', 45);
    xlabel(ax, 'Stimulus position (peripheral -> foveal)', 'FontSize', 8);
    ylabel(ax, 'Phase, unwrapped & aligned to pos 1 (deg)', 'FontSize', 8);
    legend(ax, 'Location', 'bestoutside', 'FontSize', 6);
    grid(ax, 'on'); box(ax, 'on');
else
    text(ax, 0.5, 0.5, 'No significant channels', ...
        'Units','normalized','HorizontalAlignment','center','FontSize',11);
    set(ax, 'XTick',[],'YTick',[]); box(ax,'on');
end

sgtitle(sprintf('%s -- %s -- %s [sig only] :: phase trajectory per channel\n(unwrapped, aligned to 0 at position 1; fan = position-dependent phase)', ...
    animalName, upper(dv), freq_label), 'FontSize', 10, 'FontWeight','bold');

print(fig2c, '-dpsc', '-append', '-painters', out_ps);
close(fig2c);
fprintf('Page 2C done (phase trajectories)\n');

%% ═══════════════════════════════════════════════════════════════════
%% Page 3: 2D plane fit — phase(row,col) = a·row + b·col + c
%% ═══════════════════════════════════════════════════════════════════
R2_plane   = nan(nPos, 1);
R2_pval    = nan(nPos, 1);
grad_row   = nan(nPos, 1);   % a  (rad per row)
grad_col   = nan(nPos, 1);   % b  (rad per col)
wave_dir   = nan(nPos, 1);   % atan2(a, b) in degrees
wave_speed = nan(nPos, 1);   % sqrt(a² + b²) in rad per electrode
nPerm_plane = 1000;

for p = 1:nPos
    ph_p = phi(:, p);
    valid = ~isnan(ph_p);
    nv = sum(valid);
    if nv < 4, continue; end

    X_design = [ch_row(valid), ch_col(valid), ones(nv,1)];

    % Unwrap and fit plane
    ph_uw = unwrap_2d(ph_p, ch_row, ch_col, grid_rows, grid_cols);
    b_lin = X_design \ ph_uw(valid);
    grad_row(p) = b_lin(1);
    grad_col(p) = b_lin(2);
    wave_dir(p)   = atan2d(b_lin(1), b_lin(2));
    wave_speed(p) = sqrt(b_lin(1)^2 + b_lin(2)^2);

    ph_hat_lin = X_design * b_lin;
    SS_res = sum((ph_uw(valid) - ph_hat_lin).^2);
    SS_tot = sum((ph_uw(valid) - mean(ph_uw(valid))).^2);
    if SS_tot > 0
        R2_plane(p) = 1 - SS_res/SS_tot;
    else
        R2_plane(p) = 0;
    end

    % Permutation test: shuffle phases across valid channel locations
    R2_null = nan(nPerm_plane, 1);
    ph_vals = ph_uw(valid);
    for b = 1:nPerm_plane
        ph_shuf = ph_vals(randperm(nv));
        b_shuf  = X_design \ ph_shuf;
        ph_hat_shuf = X_design * b_shuf;
        SSr_s = sum((ph_shuf - ph_hat_shuf).^2);
        SSt_s = sum((ph_shuf - mean(ph_shuf)).^2);
        if SSt_s > 0
            R2_null(b) = 1 - SSr_s/SSt_s;
        else
            R2_null(b) = 0;
        end
    end
    R2_pval(p) = mean(R2_null >= R2_plane(p), 'omitnan');
end

fig3 = figure('Visible','off','Units','centimeters', ...
    'Position',[1 1 28 12]);
set(fig3,'PaperUnits','centimeters','PaperSize',fig3.Position(3:4), ...
    'PaperPosition',[0 0 fig3.Position(3:4)]);

% Left: R² per position (colored by significance)
ax1 = subplot(1,3,1); hold(ax1, 'on');
for p = 1:nPos
    if isnan(R2_plane(p)), continue; end
    if R2_pval(p) < alpha
        bc = [0.85 0.25 0.20];
    else
        bc = [0.3 0.5 0.75];
    end
    bar(ax1, p, R2_plane(p), 'FaceColor', bc, 'EdgeColor','none', ...
        'HandleVisibility','off');
    if R2_pval(p) < alpha
        text(ax1, p, R2_plane(p)+0.02, '*', 'HorizontalAlignment','center', ...
            'FontSize', 10, 'FontWeight','bold', 'Color', [0.85 0.25 0.20]);
    end
end
bar(ax1, nan, nan, 'FaceColor', [0.85 0.25 0.20], 'EdgeColor','none', ...
    'DisplayName', sprintf('p < %.2f', alpha));
bar(ax1, nan, nan, 'FaceColor', [0.3 0.5 0.75], 'EdgeColor','none', ...
    'DisplayName', 'n.s.');
legend(ax1, 'Location','best', 'FontSize', 6);
set(ax1, 'XTick', 1:nPos, 'XTickLabel', ...
    arrayfun(@(p) sprintf('%g',p), positions, 'UniformOutput',false));
set(ax1, 'XTickLabelRotation', 45, 'FontSize', 7);
xlabel(ax1, 'Stimulus position', 'FontSize', 8);
ylabel(ax1, 'R^2 (plane fit)', 'FontSize', 8);
title(ax1, sprintf('Plane fit quality (%d perms)', nPerm_plane), 'FontSize', 9, 'FontWeight','bold');
ylim(ax1, [0 max(0.5, max(R2_plane)*1.2)]);
grid(ax1, 'on'); box(ax1, 'on');

% Middle: gradient magnitude per position
ax2 = subplot(1,3,2); hold(ax2, 'on');
bar(ax2, 1:nPos, rad2deg(wave_speed), 'FaceColor', [0.5 0.7 0.3], 'EdgeColor','none');
set(ax2, 'XTick', 1:nPos, 'XTickLabel', ...
    arrayfun(@(p) sprintf('%g',p), positions, 'UniformOutput',false));
set(ax2, 'XTickLabelRotation', 45, 'FontSize', 7);
xlabel(ax2, 'Stimulus position', 'FontSize', 8);
ylabel(ax2, 'Gradient magnitude (deg/electrode)', 'FontSize', 8);
title(ax2, 'Wave strength', 'FontSize', 9, 'FontWeight','bold');
grid(ax2, 'on'); box(ax2, 'on');

% Right: wave direction arrows (quiver on a position axis)
ax3 = subplot(1,3,3); hold(ax3, 'on');
pos_cmap = [linspace(0.15, 0.85, nPos)', ...
            linspace(0.30, 0.15, nPos)', ...
            linspace(0.80, 0.20, nPos)'];
for p = 1:nPos
    if isnan(grad_row(p)), continue; end
    scale = 2;
    quiver(ax3, 0, 0, grad_col(p)*scale, grad_row(p)*scale, 0, ...
        'Color', pos_cmap(p,:), 'LineWidth', 2, 'MaxHeadSize', 0.5);
    text(ax3, grad_col(p)*scale*1.15, grad_row(p)*scale*1.15, ...
        sprintf('%g', positions(p)), 'FontSize', 6, ...
        'Color', pos_cmap(p,:), 'HorizontalAlignment','center');
end
axis(ax3, 'equal');
mx = max(abs([grad_row; grad_col])) * scale * 1.4;
if ~isfinite(mx) || mx == 0, mx = 1; end
xlim(ax3, [-mx mx]); ylim(ax3, [-mx mx]);
xlabel(ax3, 'Col gradient (rad/el)', 'FontSize', 8);
ylabel(ax3, 'Row gradient (rad/el)', 'FontSize', 8);
title(ax3, 'Wave direction per position', 'FontSize', 9, 'FontWeight','bold');
plot(ax3, [-mx mx], [0 0], 'k-', 'LineWidth', 0.3);
plot(ax3, [0 0], [-mx mx], 'k-', 'LineWidth', 0.3);
grid(ax3, 'on'); box(ax3, 'on');
set(ax3, 'FontSize', 7);

sgtitle(sprintf('%s — %s — %s [%s] :: 2D plane fit to phase gradient', ...
    animalName, upper(dv), freq_label, ch_label), 'FontSize', 10, 'FontWeight','bold');

print(fig3, '-dpsc', '-append', '-painters', out_ps);
close(fig3);
fprintf('Page 3 done\n');

%% ═══════════════════════════════════════════════════════════════════
%% Page 4: Summary — phase heatmap for 3 representative positions
%% with gradient arrows overlaid on the array
%% ═══════════════════════════════════════════════════════════════════
rep_pos = [1, ceil(nPos/2), nPos];   % peripheral, middle, foveal

fig4 = figure('Visible','off','Units','centimeters', ...
    'Position',[1 1 30 10]);
set(fig4,'PaperUnits','centimeters','PaperSize',fig4.Position(3:4), ...
    'PaperPosition',[0 0 fig4.Position(3:4)]);

for k = 1:3
    p = rep_pos(k);
    ax = subplot(1,3,k); hold(ax,'on');

    phase_grid = nan(grid_rows, grid_cols);
    for ch = 1:nCh
        phase_grid(ch_row(ch), ch_col(ch)) = phi(ch, p);
    end

    h_img = imagesc(ax, 1:grid_cols, 1:grid_rows, phase_grid);
    set(h_img, 'AlphaData', ~isnan(phase_grid));
    set(ax, 'YDir','normal', 'Color',[0.9 0.9 0.9]);
    colormap(ax, hsv); caxis(ax, [-pi pi]);

    % Overlay gradient arrow at array center
    if ~isnan(grad_row(p))
        arr_scale = 3;
        quiver(ax, grid_cols/2+0.5, grid_rows/2+0.5, ...
            grad_col(p)*arr_scale, grad_row(p)*arr_scale, 0, ...
            'k', 'LineWidth', 2.5, 'MaxHeadSize', 0.8);
    end

    % Mark driven channels
    for ch = 1:nCh
        if ch_covers_pos(ch, p)
            plot(ax, ch_col(ch), ch_row(ch), 'ks', 'MarkerSize', 12, ...
                'LineWidth', 1.5, 'MarkerFaceColor', 'none');
        end
    end

    set(ax, 'XTick', 1:grid_cols, 'YTick', 1:grid_rows, 'FontSize', 6);
    axis(ax, 'square');
    xlabel(ax, 'Col', 'FontSize', 7); ylabel(ax, 'Row', 'FontSize', 7);
    title(ax, sprintf('pos %g  (R²=%.2f)', positions(p), R2_plane(p)), ...
        'FontSize', 9, 'FontWeight','bold');
end

cb = colorbar('Position',[0.93 0.15 0.015 0.7]);
cb.Ticks = [-pi -pi/2 0 pi/2 pi];
cb.TickLabels = {'-\pi','-\pi/2','0','\pi/2','\pi'};
cb.FontSize = 7;

sgtitle(sprintf('%s — %s — %s [%s] :: phase gradient summary\n(black arrow = fitted wave direction; black squares = driven channels)', ...
    animalName, upper(dv), freq_label, ch_label), 'FontSize', 10, 'FontWeight','bold');

print(fig4, '-dpsc', '-append', '-painters', out_ps);
close(fig4);
fprintf('Page 4 done\n');

%% ═══════════════════════════════════════════════════════════════════
%% Page 5: Moran's I — spatial autocorrelation of phase per position
%% ═══════════════════════════════════════════════════════════════════
nPerm_moran = 1000;
I_obs  = nan(nPos, 1);
I_pval = nan(nPos, 1);

% Build spatial weights matrix (inverse-distance on the grid, queen contiguity)
% for the channels that have valid (non-NaN) phase.
% Recompute per position because sig_only filtering can leave different
% valid sets, but in practice the valid set is the same across positions
% (sig_at_f doesn't depend on position). Still, keep it general.

for p = 1:nPos
    [I_obs(p), I_pval(p)] = moran_i_phase(phi(:,p), ch_row, ch_col, nPerm_moran);
end

fig5 = figure('Visible','off','Units','centimeters', ...
    'Position',[1 1 28 12]);
set(fig5,'PaperUnits','centimeters','PaperSize',fig5.Position(3:4), ...
    'PaperPosition',[0 0 fig5.Position(3:4)]);

% Left panel: Moran's I per position
ax1 = subplot(1,2,1); hold(ax1,'on');
bar_colors = repmat([0.3 0.5 0.75], nPos, 1);
bar_colors(I_pval < alpha, :) = repmat([0.85 0.25 0.20], sum(I_pval < alpha), 1);
for p = 1:nPos
    bar(ax1, p, I_obs(p), 'FaceColor', bar_colors(p,:), 'EdgeColor','none');
end
set(ax1, 'XTick', 1:nPos, 'XTickLabel', ...
    arrayfun(@(p) sprintf('%g', p), positions, 'UniformOutput', false));
set(ax1, 'XTickLabelRotation', 45, 'FontSize', 7);
xlabel(ax1, 'Stimulus position', 'FontSize', 8);
ylabel(ax1, "Moran's I", 'FontSize', 8);
title(ax1, "Moran's I per position", 'FontSize', 9, 'FontWeight','bold');
yline(ax1, 0, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
grid(ax1, 'on'); box(ax1, 'on');
legend_entries = gobjects(2,1);
legend_entries(1) = bar(ax1, nan, nan, 'FaceColor', [0.85 0.25 0.20], 'EdgeColor','none');
legend_entries(2) = bar(ax1, nan, nan, 'FaceColor', [0.3 0.5 0.75], 'EdgeColor','none');
legend(ax1, legend_entries, {sprintf('p < %.2f', alpha), 'n.s.'}, ...
    'Location','best', 'FontSize', 7);

% Right panel: p-value per position (log scale)
ax2 = subplot(1,2,2); hold(ax2,'on');
for p = 1:nPos
    bar(ax2, p, -log10(max(I_pval(p), 1/nPerm_moran)), ...
        'FaceColor', bar_colors(p,:), 'EdgeColor','none');
end
yline(ax2, -log10(alpha), '--', 'Color', [0.85 0.15 0.15], 'LineWidth', 1);
set(ax2, 'XTick', 1:nPos, 'XTickLabel', ...
    arrayfun(@(p) sprintf('%g', p), positions, 'UniformOutput', false));
set(ax2, 'XTickLabelRotation', 45, 'FontSize', 7);
xlabel(ax2, 'Stimulus position', 'FontSize', 8);
ylabel(ax2, '-log_{10}(p)', 'FontSize', 8);
title(ax2, 'Permutation p-value', 'FontSize', 9, 'FontWeight','bold');
grid(ax2, 'on'); box(ax2, 'on');

sgtitle(sprintf("%s — %s — %s [%s] :: Moran's I (spatial autocorrelation of phase)\n(red = significant at p < %.2f, %d permutations)", ...
    animalName, upper(dv), freq_label, ch_label, alpha, nPerm_moran), ...
    'FontSize', 10, 'FontWeight','bold');

print(fig5, '-dpsc', '-append', '-painters', out_ps);
close(fig5);
fprintf('Page 5 done (Moran''s I)\n');

%% ═══════════════════════════════════════════════════════════════════
%% Page 6: Phase vs. retinotopic distance (Test 3, sig channels only)
%% ═══════════════════════════════════════════════════════════════════
n_rf = min(nCh, height(rf_tab));

fig6 = figure('Visible','off','Units','centimeters', ...
    'Position',[1 1 max(24, 5*cols_pg1) max(16, 5*rows_pg1+2)]);
set(fig6,'PaperUnits','centimeters','PaperSize',fig6.Position(3:4), ...
    'PaperPosition',[0 0 fig6.Position(3:4)]);

for p = 1:nPos
    ax = subplot(rows_pg1, cols_pg1, p); hold(ax, 'on');

    drv = find(ch_covers_pos(:, p));
    drv_rf = drv(drv <= n_rf);
    if isempty(drv_rf)
        text(ax, 0.5, 0.5, 'no driven ch', 'Units','normalized', ...
            'HorizontalAlignment','center', 'FontSize', 8);
        title(ax, sprintf('pos %g', positions(p)), 'FontSize', 8);
        continue
    end

    stim_x = mean(rf_x(drv_rf));
    stim_y = mean(rf_y(drv_rf));

    rf_dist = nan(nCh, 1);
    rf_dist(1:n_rf) = sqrt((rf_x(1:n_rf) - stim_x).^2 + ...
                            (rf_y(1:n_rf) - stim_y).^2);

    drv_sig = drv(sig_at_f(drv));
    if isempty(drv_sig)
        phi_drv_mean = angle(mean(exp(1i * phi_all(drv, p))));
    else
        phi_drv_mean = angle(mean(exp(1i * phi_sig(drv_sig, p))));
    end
    phi_rel = angle(exp(1i * (phi_sig(:,p) - phi_drv_mean)));

    valid = ~isnan(phi_rel) & ~isnan(rf_dist);

    for ch = 1:nCh
        if ~valid(ch), continue; end
        if ch_covers_pos(ch, p)
            col = [0.85 0.15 0.15]; mk = 8;
        else
            col = [0.2 0.4 0.75]; mk = 5;
        end
        scatter(ax, rf_dist(ch), rad2deg(phi_rel(ch)), mk^2, col, 'filled', ...
            'MarkerFaceAlpha', 0.6);
    end

    x = rf_dist(valid); y = phi_rel(valid);
    if numel(x) > 2
        p_fit = polyfit(x, y, 1);
        x_line = linspace(0, max(x), 50);
        plot(ax, x_line, rad2deg(polyval(p_fit, x_line)), 'k--', 'LineWidth', 1);
        y_hat = polyval(p_fit, x);
        SS_res = sum((y - y_hat).^2);
        SS_tot = sum((y - mean(y)).^2);
        R2 = max(0, 1 - SS_res / max(SS_tot, eps));
        text(ax, 0.95, 0.95, sprintf('R^2=%.2f', R2), ...
            'Units','normalized', 'HorizontalAlignment','right', ...
            'VerticalAlignment','top', 'FontSize', 6, 'BackgroundColor','w');
    end

    xlabel(ax, 'RF dist to stim (dva)', 'FontSize', 6);
    ylabel(ax, 'Phase rel. driven (deg)', 'FontSize', 6);
    title(ax, sprintf('pos %g', positions(p)), 'FontSize', 8);
    set(ax, 'FontSize', 6); grid(ax, 'on'); box(ax, 'on');
end

sgtitle(sprintf('%s -- %s -- %s [sig only] :: phase vs retinotopic distance\n(red = driven, blue = non-driven sig; dashed = linear fit)', ...
    animalName, upper(dv), freq_label), 'FontSize', 10, 'FontWeight','bold');

print(fig6, '-dpsc', '-append', '-painters', out_ps);
close(fig6);
fprintf('Page 6 done (phase vs retinotopic dist)\n');

%% Convert PS → PDF
[s, msg] = system(sprintf('ps2pdf "%s" "%s"', out_ps, out_pdf));
if s == 0 && isfile(out_pdf)
    delete(out_ps);
    fprintf('Saved %s\n', out_pdf);
else
    warning('ps2pdf failed (%d): %s', s, strtrim(msg));
end

end % animal loop
end % dv loop

fprintf('\nDone.\n');


%% =====================================================================
%% Helper: H2-H1 significance at a single frequency
%% =====================================================================
function sig = get_h2h1_sig_single(base, animalName, dv_key, nCh, nFreq, fi, alpha)
sig = false(nCh, 1);
H1f = fullfile(base,['results_' animalName],'multi_lin_reg','complex','cp10_till_100','multi_regression_channelwise_R2.mat');
H2f = fullfile(base,['results_' animalName],'multi_lin_reg','abs_per_pos','cp10_till_100','multi_regression_channelwise_R2_abs_per_pos.mat');
if ~isfile(H1f)||~isfile(H2f), return; end
S1 = load(H1f,'reg_results'); S2 = load(H2f,'reg_results');
if ~isfield(S1.reg_results,dv_key), return; end
obs_d = S2.reg_results.(dv_key).R2_phase - S1.reg_results.(dv_key).R2_phase;
h1r = fullfile(base,['results_' animalName],'multi_lin_reg','complex','cp10_till_100','perm_R',dv_key);
h2r = fullfile(base,['results_' animalName],'multi_lin_reg','abs_per_pos','cp10_till_100','perm_R_pos',dv_key);
for ch = 1:nCh
    try
        P1 = load(fullfile(h1r,num2str(ch),'per_channel_null.mat'),'null_R2_phase');
        P2 = load(fullfile(h2r,num2str(ch),'per_channel_null.mat'),'null_R2_phase');
        p1 = P1.null_R2_phase; p2 = P2.null_R2_phase;
        np = min(size(p1,1),size(p2,1));
        if np<1||size(p1,2)~=nFreq||size(p2,2)~=nFreq, continue; end
        dn = p2(1:np,:)-p1(1:np,:);
        thr = quantile(max(dn,[],2), 1-alpha);
        if isfinite(thr), sig(ch) = obs_d(ch,fi) >= thr; end
    catch
    end
end
end

%% =====================================================================
%% Helper: H2-H1 significance — union across a frequency band
%% =====================================================================
function sig = get_h2h1_sig_band(base, animalName, dv_key, nCh, nFreq, band_idx, alpha)
sig = false(nCh, 1);
H1f = fullfile(base,['results_' animalName],'multi_lin_reg','complex','cp10_till_100','multi_regression_channelwise_R2.mat');
H2f = fullfile(base,['results_' animalName],'multi_lin_reg','abs_per_pos','cp10_till_100','multi_regression_channelwise_R2_abs_per_pos.mat');
if ~isfile(H1f)||~isfile(H2f), return; end
S1 = load(H1f,'reg_results'); S2 = load(H2f,'reg_results');
if ~isfield(S1.reg_results,dv_key), return; end
obs_d = S2.reg_results.(dv_key).R2_phase - S1.reg_results.(dv_key).R2_phase;
h1r = fullfile(base,['results_' animalName],'multi_lin_reg','complex','cp10_till_100','perm_R',dv_key);
h2r = fullfile(base,['results_' animalName],'multi_lin_reg','abs_per_pos','cp10_till_100','perm_R_pos',dv_key);
for ch = 1:nCh
    try
        P1 = load(fullfile(h1r,num2str(ch),'per_channel_null.mat'),'null_R2_phase');
        P2 = load(fullfile(h2r,num2str(ch),'per_channel_null.mat'),'null_R2_phase');
        p1 = P1.null_R2_phase; p2 = P2.null_R2_phase;
        np = min(size(p1,1),size(p2,1));
        if np<1||size(p1,2)~=nFreq||size(p2,2)~=nFreq, continue; end
        dn = p2(1:np,:)-p1(1:np,:);
        thr = quantile(max(dn,[],2), 1-alpha);
        if isfinite(thr)
            sig(ch) = any(obs_d(ch, band_idx) >= thr);
        end
    catch
    end
end
end

%% =====================================================================
%% Helper: Moran's I for circular (phase) data on the electrode grid
%% =====================================================================
function [I, pval] = moran_i_phase(ph, ch_row, ch_col, nPerm)
% Computes Moran's I on cos(phase) using queen-contiguity weights
% (8-connected neighbors on the grid). Channels with NaN phase are
% excluded. Significance via random permutation of phases across the
% valid channel locations.

valid = find(~isnan(ph));
n = numel(valid);
if n < 4
    I = NaN; pval = NaN; return
end

r = ch_row(valid);
c = ch_col(valid);
x = cos(ph(valid));

% Queen-contiguity weight matrix: w_ij = 1 if Chebyshev distance <= 1
W = zeros(n);
for i = 1:n
    for j = i+1:n
        if max(abs(r(i)-r(j)), abs(c(i)-c(j))) <= 1
            W(i,j) = 1;
            W(j,i) = 1;
        end
    end
end
S0 = sum(W(:));
if S0 == 0
    I = NaN; pval = NaN; return
end

xbar = mean(x);
xc   = x - xbar;
I    = (n / S0) * (xc' * W * xc) / (xc' * xc);

% Permutation test
I_null = nan(nPerm, 1);
for b = 1:nPerm
    xp   = x(randperm(n));
    xpc  = xp - mean(xp);
    ss   = xpc' * xpc;
    if ss == 0, continue; end
    I_null(b) = (n / S0) * (xpc' * W * xpc) / ss;
end
pval = mean(I_null >= I, 'omitnan');
end

%% =====================================================================
%% Local helper: 2D unwrap
%% =====================================================================
function ph_uw = unwrap_2d(ph, ch_row, ch_col, nR, nC)
% Unwrap phases on a 2D grid by flood-filling from the center,
% adjusting each neighbor to be within π of its already-unwrapped
% neighbor. Falls back to raw phases for isolated channels.
ph_uw = ph;
grid = nan(nR, nC);
for k = 1:numel(ph)
    grid(ch_row(k), ch_col(k)) = ph(k);
end
visited = false(nR, nC);
uw_grid = grid;

% Start from center of array
sr = round(nR/2); sc = round(nC/2);
if isnan(grid(sr,sc))
    % Find nearest non-NaN to center
    for r = 1:nR; for c = 1:nC
        if ~isnan(grid(r,c)), sr = r; sc = c; break; end
    end; if ~isnan(grid(sr,sc)), break; end; end
end

queue = [sr, sc];
visited(sr, sc) = true;
while ~isempty(queue)
    r = queue(1,1); c = queue(1,2);
    queue(1,:) = [];
    nb = [r-1,c; r+1,c; r,c-1; r,c+1];
    for n = 1:4
        nr = nb(n,1); nc = nb(n,2);
        if nr < 1 || nr > nR || nc < 1 || nc > nC, continue; end
        if visited(nr,nc) || isnan(grid(nr,nc)), continue; end
        visited(nr,nc) = true;
        diff = uw_grid(nr,nc) - uw_grid(r,c);
        uw_grid(nr,nc) = uw_grid(nr,nc) - 2*pi*round(diff/(2*pi));
        queue = [queue; nr, nc]; %#ok<AGROW>
    end
end

for k = 1:numel(ph)
    if visited(ch_row(k), ch_col(k))
        ph_uw(k) = uw_grid(ch_row(k), ch_col(k));
    end
end
end
