% =====================================================================
% Assemble paper figure: position-dependent phase in V4 LFP
%
% Produces a single multi-page PDF with 5 panels (A–E), each on its
% own page, ready for arrangement in Illustrator / Inkscape.
%
% Panel A — H2 > H1: monkey-average paired differences (Reg R² LFP)
%           + per-hypothesis overlay showing H2 separates from H1 at
%           5–25 Hz.
%
% Panel B — Per-channel H2−H1 significance: heatmap (channel × freq)
%           for klecks + n_sig_channels bar plot.
%
% Panel C — Phase arrows at 10 Hz: 3 significant + 2 non-significant
%           channels from klecks, showing position-dependent spread
%           at significant channels.
%
% Panel D — Traveling wave: phase heatmaps on the 8×8 grid for 3
%           positions (peripheral, middle, foveal) with gradient
%           arrows, plus R² and direction summary.
%
% Panel E — RF overlap map: 8×8 grid showing which channels are
%           significant (H2−H1), which have RF overlap, and which are
%           significant WITHOUT RF overlap.
%
% Output:
%   Plots/sampling_compare/hypotheses/paper_figure_panels.pdf
% =====================================================================
clearvars; close all; clc

base = '/mnt/hpc/projects/MWSampling/4Shivangi';
alpha = 0.05;
nCh = 64;
grid_rows = 8;
grid_cols = 8;
animalName = 'klecks';
dv_key = 'LFP_ERP_ampl_all';
show_freq = 10;

out_dir = fullfile(base, 'Plots','sampling_compare','hypotheses');
out_pdf = fullfile(out_dir, 'paper_figure_panels.pdf');
out_ps  = strrep(out_pdf, '.pdf', '.ps');
if isfile(out_pdf), delete(out_pdf); end
if isfile(out_ps),  delete(out_ps);  end

%% Channel grid mapping (column-wise bottom-to-top)
ch_col = ceil((1:nCh)' / grid_rows);
ch_row = grid_rows - mod((1:nCh)' - 1, grid_rows);

%% Load shared data
% Phase progression
pp = load(fullfile(base, ['results_' animalName], 'scanning', ...
    'phase_progression','cp10_till_100','lfp','phase_progression.mat'), ...
    'pref_phase','coh_mag','freq','positions');
freq = pp.freq;
positions = pp.positions;
nPos = numel(positions);
nFreq = numel(freq);
[~, fi] = min(abs(freq - show_freq));
f_hz = freq(fi);

% H2−H1 per-channel significance
H1_obs = load(fullfile(base, ['results_' animalName], ...
    'multi_lin_reg','complex','cp10_till_100', ...
    'multi_regression_channelwise_R2.mat'), 'reg_results');
H2_obs = load(fullfile(base, ['results_' animalName], ...
    'multi_lin_reg','abs_per_pos','cp10_till_100', ...
    'multi_regression_channelwise_R2_abs_per_pos.mat'), 'reg_results');
obs_H1 = H1_obs.reg_results.(dv_key).R2_phase;
obs_H2 = H2_obs.reg_results.(dv_key).R2_phase;
obs_diff = obs_H2 - obs_H1;

h1_root = fullfile(base, ['results_' animalName], ...
    'multi_lin_reg','complex','cp10_till_100','perm_R', dv_key);
h2_root = fullfile(base, ['results_' animalName], ...
    'multi_lin_reg','abs_per_pos','cp10_till_100','perm_R_pos', dv_key);

sig_mask = false(nCh, nFreq);
pc_thr = nan(nCh, 1);
for ch = 1:nCh
    h1f = fullfile(h1_root, num2str(ch), 'per_channel_null.mat');
    h2f = fullfile(h2_root, num2str(ch), 'per_channel_null.mat');
    if ~isfile(h1f) || ~isfile(h2f), continue; end
    P1 = load(h1f, 'null_R2_phase');
    P2 = load(h2f, 'null_R2_phase');
    if ~isfield(P1,'null_R2_phase') || ~isfield(P2,'null_R2_phase'), continue; end
    p1 = P1.null_R2_phase; p2 = P2.null_R2_phase;
    n_pair = min(size(p1,1), size(p2,1));
    if n_pair < 1 || size(p1,2) ~= nFreq || size(p2,2) ~= nFreq, continue; end
    diff_perm = p2(1:n_pair,:) - p1(1:n_pair,:);
    thr = quantile(max(diff_perm, [], 2), 1 - alpha);
    pc_thr(ch) = thr;
    if isfinite(thr)
        sig_mask(ch,:) = obs_diff(ch,:) >= thr;
    end
end
sig_at_f = sig_mask(:, fi);
sig_chs = find(sig_at_f);

% RF data
rf_file = fullfile(base, 'Plots','RF_Mapping','klecks','loc_RF_map', ...
    'gaussian_overlap', ...
    'klecks_20170802_rfmapping_bar_1_channel_target_summary.txt');
rf_tab = readtable(rf_file, 'Delimiter','\t');
ch_has_rf = false(nCh, 1);
for ch = 1:min(nCh, height(rf_tab))
    loc_str = rf_tab.Locations_Inside{ch};
    if ~strcmp(loc_str, 'none') && ~isempty(loc_str)
        ch_has_rf(ch) = true;
    end
end

% Position colormap
pos_cmap = [linspace(0.15, 0.85, nPos)', ...
            linspace(0.30, 0.15, nPos)', ...
            linspace(0.80, 0.20, nPos)'];

%% ═══════════════════════════════════════════════════════════════════
%% Panel A — H2 > H1 monkey-average
%% ═══════════════════════════════════════════════════════════════════
% Load monkey-average paired differences
ma_file = fullfile(base, 'results_combined','multi_lin_reg','abs_per_pos', ...
    'cp10_till_100', dv_key, 'monkey_avg_results.mat');
ma = load(ma_file);

% H1 monkey-average
h1_ma_file = fullfile(base, 'results_combined','multi_lin_reg','complex', ...
    'cp10_till_100', dv_key, 'monkey_avg_results.mat');
if isfile(h1_ma_file)
    h1_ma = load(h1_ma_file);
end

% Channel-avg paired diff for Reg R²
h2_chan_avg = ma.monkey_avg_obs.phase;   % 1 × nFreq
h1_chan_avg = h1_ma.monkey_avg_obs.phase;
diff_obs = h2_chan_avg - h1_chan_avg;

% Null for paired diff
h2_perm = ma.perm_monkey_avg.phase;     % nPerm × nFreq
h1_perm = h1_ma.perm_monkey_avg.phase;
diff_null = h2_perm - h1_perm;
thr_diff = quantile(max(diff_null, [], 2), 1-alpha);
null_mean = mean(diff_null, 1);
sig_freqs = diff_obs >= thr_diff;

fig = figure('Visible','off','Units','centimeters','Position',[1 1 18 8]);
set(fig,'PaperUnits','centimeters','PaperSize',fig.Position(3:4), ...
    'PaperPosition',[0 0 fig.Position(3:4)]);

% Left: H1 vs H2 overlay
ax1 = subplot(1,2,1); hold(ax1,'on');
plot(ax1, freq, h1_chan_avg, 'Color',[0.2 0.4 0.8], 'LineWidth',1.8, ...
    'DisplayName','H1 (common phase)');
plot(ax1, freq, h2_chan_avg, 'Color',[0.2 0.7 0.3], 'LineWidth',1.8, ...
    'DisplayName','H2 (per-position phase)');
xlabel(ax1,'Frequency (Hz)','FontSize',9);
ylabel(ax1,'Partial R^2_{phase} (monkey avg)','FontSize',9);
legend(ax1,'Location','northeast','FontSize',7);
title(ax1,'A1: H1 vs H2','FontSize',10,'FontWeight','bold');
set(ax1,'FontSize',8); grid(ax1,'on'); box(ax1,'on');

% Right: paired difference H2−H1
ax2 = subplot(1,2,2); hold(ax2,'on');
col_diff = [0.2 0.65 0.4];
plot(ax2, freq, diff_obs, 'Color',col_diff, 'LineWidth',1.8, ...
    'DisplayName','H2 − H1');
plot(ax2, freq, null_mean, '--','Color',col_diff, 'LineWidth',0.8, ...
    'DisplayName','null mean');
if isfinite(thr_diff)
    yline(ax2, thr_diff, ':','Color',col_diff, 'LineWidth',1, ...
        'DisplayName','95% threshold');
end
if any(sig_freqs)
    yl = ylim(ax2); span = yl(2)-yl(1);
    scatter(ax2, freq(sig_freqs), repmat(yl(2)-0.05*span,1,sum(sig_freqs)), ...
        12, col_diff, 'filled', 'HandleVisibility','off');
end
xlabel(ax2,'Frequency (Hz)','FontSize',9);
ylabel(ax2,'\Delta R^2_{phase}','FontSize',9);
legend(ax2,'Location','northeast','FontSize',7);
title(ax2,'A2: Paired difference','FontSize',10,'FontWeight','bold');
set(ax2,'FontSize',8); grid(ax2,'on'); box(ax2,'on');

sgtitle('Panel A — Positions disagree on preferred phase (monkey avg, Reg R^2 LFP)', ...
    'FontSize',11,'FontWeight','bold');
print(fig, '-dpsc', '-append', '-painters', out_ps); close(fig);
fprintf('Panel A done\n');

%% ═══════════════════════════════════════════════════════════════════
%% Panel B — Per-channel heatmap + n sig channels
%% ═══════════════════════════════════════════════════════════════════
fig = figure('Visible','off','Units','centimeters','Position',[1 1 20 8]);
set(fig,'PaperUnits','centimeters','PaperSize',fig.Position(3:4), ...
    'PaperPosition',[0 0 fig.Position(3:4)]);

% Left: heatmap
ax1 = subplot(1,2,1);
img = imagesc(ax1, freq, 1:nCh, obs_diff);
set(ax1,'YDir','normal','Color',[1 1 1]);
colormap(ax1, parula);
v = obs_diff(~isnan(obs_diff));
if ~isempty(v)
    mx = max(abs(prctile(v,[2 98])));
    if mx > 0, caxis(ax1, [-mx mx]); end
end
A = 0.20 * ones(size(obs_diff));
A(sig_mask) = 1.0;
A(isnan(obs_diff)) = 0;
set(img, 'AlphaData', A);
cb = colorbar(ax1); cb.FontSize = 6;
cb.Label.String = '\Delta R^2_{phase}';
xlabel(ax1,'Frequency (Hz)','FontSize',9);
ylabel(ax1,'Channel','FontSize',9);
title(ax1,'B1: Per-channel H2−H1 (faded = n.s.)','FontSize',10,'FontWeight','bold');
set(ax1,'FontSize',8);

% Right: n sig channels vs freq
ax2 = subplot(1,2,2);
n_sig_per_freq = sum(sig_mask, 1);
bar(ax2, freq, n_sig_per_freq, 'FaceColor',[0.2 0.65 0.4], 'EdgeColor','none');
xlabel(ax2,'Frequency (Hz)','FontSize',9);
ylabel(ax2,'# sig channels','FontSize',9);
title(ax2,'B2: Significant channels vs freq','FontSize',10,'FontWeight','bold');
set(ax2,'FontSize',8); grid(ax2,'on'); box(ax2,'on');

sgtitle(sprintf('Panel B — klecks per-channel H2−H1 (Reg R^2 LFP)'), ...
    'FontSize',11,'FontWeight','bold');
print(fig, '-dpsc', '-append', '-painters', out_ps); close(fig);
fprintf('Panel B done\n');

%% ═══════════════════════════════════════════════════════════════════
%% Panel C — Phase arrows: 3 sig + 2 non-sig channels at 10 Hz
%% ═══════════════════════════════════════════════════════════════════
show_chs_sig = [17 42 51];
show_chs_ns  = [21 28];
show_chs = [show_chs_sig, show_chs_ns];
n_show = numel(show_chs);

phi = squeeze(pp.pref_phase(:, fi, :));
mag = squeeze(pp.coh_mag(:, fi, :));

fig = figure('Visible','off','Units','centimeters','Position',[1 1 5*n_show 5.5]);
set(fig,'PaperUnits','centimeters','PaperSize',fig.Position(3:4), ...
    'PaperPosition',[0 0 fig.Position(3:4)]);

rmax = max(mag(show_chs,:), [], 'all', 'omitnan');
if ~isfinite(rmax) || rmax == 0, rmax = 1; end

for k = 1:n_show
    ch = show_chs(k);
    ax = subplot(1, n_show, k); hold(ax,'on');

    is_sig = sig_at_f(ch);

    % Reference circle
    theta_c = linspace(0, 2*pi, 100);
    plot(ax, rmax*cos(theta_c), rmax*sin(theta_c), ...
        'Color',[0.85 0.85 0.85], 'LineWidth',0.3);
    plot(ax, [-rmax rmax],[0 0], 'Color',[0.9 0.9 0.9], 'LineWidth',0.2);
    plot(ax, [0 0],[-rmax rmax], 'Color',[0.9 0.9 0.9], 'LineWidth',0.2);

    if is_sig
        arrow_alpha = 1.0; lw = 1.5;
    else
        arrow_alpha = 0.3; lw = 0.8;
    end

    for p = 1:nPos
        if isnan(phi(ch,p)) || isnan(mag(ch,p)), continue; end
        x_end = mag(ch,p) * cos(phi(ch,p));
        y_end = mag(ch,p) * sin(phi(ch,p));
        col_p = pos_cmap(p,:);
        if ~is_sig, col_p = 0.5 + 0.5*col_p; end
        quiver(ax, 0, 0, x_end, y_end, 0, ...
            'Color',[col_p, arrow_alpha], 'LineWidth',lw, 'MaxHeadSize',0.6);
    end

    lim = rmax * 1.15;
    set(ax,'XLim',[-lim lim],'YLim',[-lim lim]);
    axis(ax,'square');
    set(ax,'XTick',[],'YTick',[]);

    if is_sig
        title(ax, sprintf('ch%d *', ch), 'FontSize',9, 'FontWeight','bold', ...
            'Color',[0.1 0.5 0.2]);
        set(ax,'Box','on','LineWidth',1.5, ...
            'XColor',[0.1 0.5 0.2],'YColor',[0.1 0.5 0.2]);
    else
        title(ax, sprintf('ch%d', ch), 'FontSize',9, 'Color',[0.5 0.5 0.5]);
        set(ax,'Box','on','LineWidth',0.3, ...
            'XColor',[0.8 0.8 0.8],'YColor',[0.8 0.8 0.8]);
    end
end

sgtitle(sprintf('Panel C — Phase vectors at %.0f Hz (blue=peripheral, red=foveal)', f_hz), ...
    'FontSize',11,'FontWeight','bold');
print(fig, '-dpsc', '-append', '-painters', out_ps); close(fig);
fprintf('Panel C done\n');

%% ═══════════════════════════════════════════════════════════════════
%% Panel D — Traveling wave: phase heatmaps + gradient summary
%% ═══════════════════════════════════════════════════════════════════
% Load RF overlap per position for driven channel markers
rf_tab2 = readtable(rf_file, 'Delimiter','\t');
ch_covers_pos = false(nCh, nPos);
for ch = 1:min(nCh, height(rf_tab2))
    loc_str = rf_tab2.Locations_Inside{ch};
    if strcmp(loc_str,'none') || isempty(loc_str), continue; end
    locs = str2double(strsplit(loc_str,','));
    locs = locs(~isnan(locs) & locs >= 1 & locs <= nPos);
    ch_covers_pos(ch, locs) = true;
end

rep_pos = [1, ceil(nPos/2), nPos];
phi_all = squeeze(pp.pref_phase(:, fi, :));

fig = figure('Visible','off','Units','centimeters','Position',[1 1 24 10]);
set(fig,'PaperUnits','centimeters','PaperSize',fig.Position(3:4), ...
    'PaperPosition',[0 0 fig.Position(3:4)]);

% Compute plane fit for each position (reuse traveling wave logic)
R2_plane = nan(nPos,1);
grad_row = nan(nPos,1);
grad_col = nan(nPos,1);
for p = 1:nPos
    ph_p = phi_all(:, p);
    valid = ~isnan(ph_p);
    if sum(valid) < 4, continue; end
    X_d = [ch_row(valid), ch_col(valid), ones(sum(valid),1)];
    ph_uw = unwrap_2d_local(ph_p, ch_row, ch_col, grid_rows, grid_cols);
    b_lin = X_d \ ph_uw(valid);
    grad_row(p) = b_lin(1);
    grad_col(p) = b_lin(2);
    ph_hat = X_d * b_lin;
    SS_res = sum((ph_uw(valid) - ph_hat).^2);
    SS_tot = sum((ph_uw(valid) - mean(ph_uw(valid))).^2);
    if SS_tot > 0, R2_plane(p) = 1 - SS_res/SS_tot; end
end

for k = 1:3
    p = rep_pos(k);
    ax = subplot(1,4,k); hold(ax,'on');
    phase_grid = nan(grid_rows, grid_cols);
    for ch = 1:nCh
        phase_grid(ch_row(ch), ch_col(ch)) = phi_all(ch, p);
    end
    imagesc(ax, 1:grid_cols, 1:grid_rows, phase_grid);
    set(ax,'YDir','normal','Color',[0.9 0.9 0.9]);
    colormap(ax, hsv); caxis(ax, [-pi pi]);

    % Gradient arrow
    if ~isnan(grad_row(p))
        quiver(ax, grid_cols/2+0.5, grid_rows/2+0.5, ...
            grad_col(p)*3, grad_row(p)*3, 0, ...
            'k','LineWidth',2.5,'MaxHeadSize',0.8);
    end

    % Driven channels
    for ch = 1:nCh
        if ch_covers_pos(ch, p)
            plot(ax, ch_col(ch), ch_row(ch), 'ws', 'MarkerSize',10, 'LineWidth',1.5);
        end
    end

    set(ax,'XTick',1:grid_cols,'YTick',1:grid_rows,'FontSize',6);
    axis(ax,'square');
    xlabel(ax,'Col','FontSize',7); ylabel(ax,'Row','FontSize',7);
    r2_str = 'n/a';
    if ~isnan(R2_plane(p)), r2_str = sprintf('%.2f', R2_plane(p)); end
    title(ax, sprintf('pos %g (R^2=%s)', positions(p), r2_str), ...
        'FontSize',8,'FontWeight','bold');
end

% 4th subplot: direction arrows
ax4 = subplot(1,4,4); hold(ax4,'on');
for p = 1:nPos
    if isnan(grad_row(p)), continue; end
    quiver(ax4, 0, 0, grad_col(p)*2, grad_row(p)*2, 0, ...
        'Color',pos_cmap(p,:), 'LineWidth',2, 'MaxHeadSize',0.5);
    text(ax4, grad_col(p)*2*1.15, grad_row(p)*2*1.15, ...
        sprintf('%g',positions(p)), 'FontSize',6, 'Color',pos_cmap(p,:), ...
        'HorizontalAlignment','center');
end
axis(ax4,'equal');
mx = max(abs([grad_row; grad_col])) * 2 * 1.4;
if ~isfinite(mx) || mx == 0, mx = 1; end
xlim(ax4,[-mx mx]); ylim(ax4,[-mx mx]);
plot(ax4,[-mx mx],[0 0],'k-','LineWidth',0.3);
plot(ax4,[0 0],[-mx mx],'k-','LineWidth',0.3);
xlabel(ax4,'Col gradient','FontSize',7);
ylabel(ax4,'Row gradient','FontSize',7);
title(ax4,'Wave direction','FontSize',8,'FontWeight','bold');
grid(ax4,'on'); box(ax4,'on'); set(ax4,'FontSize',6);

sgtitle(sprintf('Panel D — Phase gradient on array at %.0f Hz (black arrow = fitted direction)', f_hz), ...
    'FontSize',11,'FontWeight','bold');
print(fig, '-dpsc', '-append', '-painters', out_ps); close(fig);
fprintf('Panel D done\n');

%% ═══════════════════════════════════════════════════════════════════
%% Panel E — RF overlap map on 8×8 grid
%% ═══════════════════════════════════════════════════════════════════
fig = figure('Visible','off','Units','centimeters','Position',[1 1 18 12]);
set(fig,'PaperUnits','centimeters','PaperSize',fig.Position(3:4), ...
    'PaperPosition',[0 0 fig.Position(3:4)]);
ax = axes('Parent',fig,'Position',[0.08 0.08 0.55 0.82]); hold(ax,'on');

% Categories
cat_sig_rf    = sig_at_f & ch_has_rf;     % sig + RF overlap
cat_sig_norf  = sig_at_f & ~ch_has_rf;    % sig + NO RF overlap
cat_nosig_rf  = ~sig_at_f & ch_has_rf;    % not sig + RF overlap
cat_nosig     = ~sig_at_f & ~ch_has_rf;   % not sig + no RF

col_sig_rf   = [0.85 0.35 0.15];   % orange-red
col_sig_norf = [0.15 0.55 0.75];   % teal-blue
col_rf_only  = [0.75 0.75 0.55];   % muted yellow
col_none     = [0.85 0.85 0.85];   % gray

for ch = 1:nCh
    if cat_sig_rf(ch)
        c = col_sig_rf; lw = 2;
    elseif cat_sig_norf(ch)
        c = col_sig_norf; lw = 2;
    elseif cat_nosig_rf(ch)
        c = col_rf_only; lw = 1;
    else
        c = col_none; lw = 0.5;
    end
    scatter(ax, ch_col(ch), ch_row(ch), 250, c, 'filled', ...
        'MarkerEdgeColor','k', 'LineWidth',lw);
    text(ax, ch_col(ch), ch_row(ch), sprintf('%d',ch), ...
        'FontSize',5, 'HorizontalAlignment','center', ...
        'VerticalAlignment','middle', 'FontWeight','bold');
end

% Legend as manual text+markers on the right side of the figure
leg_x = 0.70; leg_y = 0.78; leg_dy = 0.08;
leg_cols  = {col_sig_rf, col_sig_norf, col_rf_only, col_none};
leg_names = {'Sig + RF overlap', 'Sig + no RF overlap', ...
             'Not sig + RF overlap', 'Not significant'};
for k = 1:4
    annotation(fig, 'ellipse', [leg_x, leg_y-(k-1)*leg_dy, 0.02, 0.025], ...
        'FaceColor', leg_cols{k}, 'Color', 'k', 'LineWidth', 0.8);
    annotation(fig, 'textbox', [leg_x+0.03, leg_y-(k-1)*leg_dy-0.005, 0.25, 0.03], ...
        'String', leg_names{k}, 'EdgeColor', 'none', 'FontSize', 7, ...
        'VerticalAlignment', 'middle');
end

set(ax,'XTick',1:grid_cols,'YTick',1:grid_rows,'FontSize',8);
xlim(ax,[0.5 grid_cols+0.5]); ylim(ax,[0.5 grid_rows+0.5]);
axis(ax,'square');
xlabel(ax,'Column','FontSize',9); ylabel(ax,'Row','FontSize',9);
title(ax, sprintf('Panel E — RF overlap vs H2−H1 significance at %.0f Hz', f_hz), ...
    'FontSize',10,'FontWeight','bold');
grid(ax,'on'); box(ax,'on');

print(fig, '-dpsc', '-append', '-painters', out_ps); close(fig);
fprintf('Panel E done\n');

%% Convert PS → PDF
[s, msg] = system(sprintf('ps2pdf "%s" "%s"', out_ps, out_pdf));
if s == 0 && isfile(out_pdf)
    delete(out_ps);
    fprintf('Saved %s\n', out_pdf);
else
    warning('ps2pdf failed (%d): %s', s, strtrim(msg));
end

%% =====================================================================
%% Local helper: 2D unwrap
%% =====================================================================
function ph_uw = unwrap_2d_local(ph, ch_row, ch_col, nR, nC)
ph_uw = ph;
grid = nan(nR, nC);
for k = 1:numel(ph)
    grid(ch_row(k), ch_col(k)) = ph(k);
end
visited = false(nR, nC);
uw_grid = grid;
sr = round(nR/2); sc = round(nC/2);
if isnan(grid(sr,sc))
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
        d = uw_grid(nr,nc) - uw_grid(r,c);
        uw_grid(nr,nc) = uw_grid(nr,nc) - 2*pi*round(d/(2*pi));
        queue = [queue; nr, nc]; %#ok<AGROW>
    end
end
for k = 1:numel(ph)
    if visited(ch_row(k), ch_col(k))
        ph_uw(k) = uw_grid(ch_row(k), ch_col(k));
    end
end
end
