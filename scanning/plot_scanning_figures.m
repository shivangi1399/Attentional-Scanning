% =====================================================================
% Plot all scanning / phase-progression figures.
%
% For each (animal × DV), produces ONE multi-page PDF at
%   Plots/scanning/phase_progression/cp10_till_100/<dv>/<animal>_scanning_figures.pdf
% with the following pages:
%
%   1. Main figure — preferred-phase trajectories across stimulus
%      position for every significant (channel × freq). Phases are
%      unwrapped across positions and aligned to 0 at position 1.
%      Thin gray = individual (ch × freq) trajectories. Thick blue =
%      circular mean trajectory across the significant set. Shaded
%      band = 95% bootstrap CI of the circular mean (also aligned).
%
%   2. Systematicity per channel — R_obs(ch, f) vs frequency, with the
%      channel's max-stat 95% threshold (dashed) and significant-freq
%      shading.
%
%   3. Preferred-phase heatmaps per channel (all freqs) — y = freq,
%      x = position, color = preferred phase (HSV, [-π, π]).
%
%   4. Preferred-phase heatmaps per channel (significant freqs only) —
%      non-significant rows blanked to white.
%
%   5. Polar histograms of step phases per channel (all freqs) — angle
%      = Δφ direction (peripheral→foveal), radius = count. Red radial
%      line = mean resultant vector; title shows R and mean step in
%      degrees.
%
%   6. Polar histograms of step phases per channel (significant freqs
%      only).
%
%   7. (LFP DV only) Band-restricted trajectory using only the
%      low-frequency band specified by LFP_LOW_BAND below. Strips
%      higher-freq noise from the main trajectory plot.
%
% Excluded channels: `excluded_channels` lists channels whose
% permutation null collapsed to all-zeros (e.g. ch37 in hermes — NaN
% input). Those channels are forced sig_mask = false so they don't
% contaminate pages 1, 4, 6, 7 or the cross-DV summary.
%
% After the per-(animal, DV) loop, a summary figure is written to
%   Plots/scanning/phase_progression/cp10_till_100/summary_sig_count_vs_freq.pdf
% showing the # of FWER-corrected significant channels at each
% frequency, one curve per (animal, DV). This is the fastest way to
% spot which spectral bands carry the effect.
%
% Per-channel thresholds (max-stat across freqs, MC-corrected within
% the channel) are also persisted alongside the data at:
%   results_<animal>/scanning/phase_progression/cp10_till_100/<dv>/
%       channel_thresholds.mat
% =====================================================================
clearvars; close all; clc

%% Settings
dv_types          = {'mua','lfp','RT','hit_miss'};
animals           = {'hermes','klecks'};
alpha             = 0.05;
nBins             = 24;             % 15° bins for polar histograms
excluded_channels = 37;             % null collapsed to 0 (NaN input)
LFP_LOW_BAND      = [0 10];         % Hz, for page-7 band-restricted plot

% Accumulator for the cross-DV summary figure built after the main loop.
summary = cell(numel(animals), numel(dv_types));

for a = 1:numel(animals)
    animalName = animals{a};
    for d = 1:numel(dv_types)
        dv = dv_types{d};

        in_file = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ...
            ['results_' animalName],'scanning','phase_progression','cp10_till_100', ...
            dv,'phase_progression.mat');
        if ~isfile(in_file)
            warning('No data for %s / %s — run phase_progression_per_chan.m first.', ...
                animalName, dv);
            continue
        end

        S = load(in_file, 'pref_phase','step_phase','R_obs','R_null','freq','positions');
        [nCh, ~, nPos] = size(S.pref_phase);
        if ~isfield(S,'positions') || isempty(S.positions)
            S.positions = (1:nPos)';
        end

        % Per-channel thresholds.
        %   thresh_chan_max  — max-stat across freqs within each channel
        %                      → FWER-corrected within-channel threshold.
        %   thresh_chan_perf — per (ch, freq) 95th pct of the null,
        %                      UNCORRECTED across frequencies. Kept on
        %                      disk for diagnostics only.
        % sig_mask uses the FWER-corrected threshold so every downstream
        % page reports the same standard as the dashed line on page 2.
        chan_max_null    = squeeze(max(S.R_null, [], 2));         % nCh × nPerm
        thresh_chan_max  = quantile(chan_max_null, 1-alpha, 2);   % nCh × 1
        thresh_chan_perf = quantile(S.R_null, 1-alpha, 3);        % nCh × nFreq
        sig_mask         = bsxfun(@ge, S.R_obs, thresh_chan_max); % nCh × nFreq

        % Hard-exclude artefact channels (null collapsed to 0).
        % Also auto-exclude any other channel whose max-stat threshold
        % is 0 or NaN — same diagnosis (empty/NaN data).
        bad_ch              = excluded_channels(:)';
        bad_ch              = unique([bad_ch, ...
            find(thresh_chan_max <= 0 | isnan(thresh_chan_max))']);
        sig_mask(bad_ch, :) = false;

        % Blank excluded channels in all per-channel arrays so the
        % grid pages (2, 3, 5) render them as empty panels rather than
        % spurious data.
        S.R_obs(bad_ch, :)         = NaN;
        S.pref_phase(bad_ch, :, :) = NaN;
        S.step_phase(bad_ch, :, :) = NaN;
        thresh_chan_max(bad_ch)    = NaN;

        save(fullfile(fileparts(in_file),'channel_thresholds.mat'), ...
            'thresh_chan_max','thresh_chan_perf','bad_ch','-v7.3');

        % Stash for the cross-DV summary plot.
        summary{a, d} = struct( ...
            'animal',  animalName, ...
            'dv',      dv, ...
            'freq',    S.freq(:)', ...
            'n_sig_f', sum(sig_mask, 1));

        plot_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots', ...
            'scanning','phase_progression','cp10_till_100', dv);
        if ~exist(plot_dir,'dir'), mkdir(plot_dir); end
        out_pdf = fullfile(plot_dir, sprintf('%s_scanning_figures.pdf', animalName));
        out_ps  = strrep(out_pdf, '.pdf', '.ps');
        if isfile(out_pdf), delete(out_pdf); end   % start fresh
        if isfile(out_ps),  delete(out_ps);  end

        tag = sprintf('%s — %s', animalName, upper(dv));

        % R2020 has exportgraphics but NOT its 'Append' option (R2022a+).

        % Page 1 — MAIN figure: phase-vs-position trajectories for all
        % significant (channel × frequency). This is the central panel.
        f = fig_phase_trajectory(S.pref_phase, sig_mask, S.positions, ...
            sprintf('%s :: preferred-phase trajectories (sig ch×freq, p<%.2f)', tag, alpha));
        print(f, '-dpsc', '-append', '-painters', out_ps);
        close(f);

        % Page 2 — systematicity per channel
        f = fig_systematicity_grid(S.R_obs, S.freq, thresh_chan_max, ...
            sprintf('%s :: R per channel (dashed = max-stat 95%%)', tag));
        print(f, '-dpsc', '-append', '-painters', out_ps);
        close(f);

        % Page 3 — preferred phase per channel, all freqs
        f = fig_pref_phase_grid(S.pref_phase, S.R_obs, [], S.freq, nPos, ...
            sprintf('%s :: preferred phase per channel (all freqs)', tag));
        print(f, '-dpsc', '-append', '-painters', out_ps);
        close(f);

        % Page 4 — preferred phase per channel, sig freqs only
        f = fig_pref_phase_grid(S.pref_phase, S.R_obs, sig_mask, S.freq, nPos, ...
            sprintf('%s :: preferred phase per channel (sig freqs, p<%.2f)', tag, alpha));
        print(f, '-dpsc', '-append', '-painters', out_ps);
        close(f);

        % Page 5 — step histograms per channel, all freqs
        f = fig_step_hist_grid(S.step_phase, [], nBins, ...
            sprintf('%s :: Δφ per channel (all freqs)', tag));
        print(f, '-dpsc', '-append', '-painters', out_ps);
        close(f);

        % Page 6 — step histograms per channel, sig freqs only
        f = fig_step_hist_grid(S.step_phase, sig_mask, nBins, ...
            sprintf('%s :: Δφ per channel (sig freqs, p<%.2f)', tag, alpha));
        print(f, '-dpsc', '-append', '-painters', out_ps);
        close(f);

        % Page 7 (LFP only) — trajectory restricted to the low-frequency
        % band. LFP is the DV where step-phase polar plots show a tight
        % near-zero direction across channels, so confining the
        % trajectory to that band should sharpen the population mean.
        if strcmp(dv, 'lfp')
            band_mask                 = S.freq(:)' >= LFP_LOW_BAND(1) ...
                                      & S.freq(:)' <= LFP_LOW_BAND(2);
            sig_mask_band             = sig_mask;
            sig_mask_band(:, ~band_mask) = false;
            f = fig_phase_trajectory(S.pref_phase, sig_mask_band, S.positions, ...
                sprintf('%s :: LFP %g–%g Hz trajectory (sig ch×freq, p<%.2f)', ...
                    tag, LFP_LOW_BAND(1), LFP_LOW_BAND(2), alpha));
            print(f, '-dpsc', '-append', '-painters', out_ps);
            close(f);
        end

        % Convert the multi-page PS to a single PDF and clean up.
        [conv_status, conv_msg] = system(sprintf('ps2pdf "%s" "%s"', out_ps, out_pdf));
        if conv_status == 0 && isfile(out_pdf)
            delete(out_ps);
            fprintf('Saved %s\n', out_pdf);
        else
            warning('ps2pdf failed for %s (status %d): %s\nLeaving PS at %s', ...
                out_pdf, conv_status, strtrim(conv_msg), out_ps);
        end
    end
end

%% Cross-DV summary: # sig channels vs frequency, one curve per (animal, DV)
summary_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots', ...
    'scanning','phase_progression','cp10_till_100');
if ~exist(summary_dir,'dir'), mkdir(summary_dir); end
sum_pdf = fullfile(summary_dir, 'summary_sig_count_vs_freq.pdf');
sum_ps  = strrep(sum_pdf, '.pdf', '.ps');
if isfile(sum_pdf), delete(sum_pdf); end
if isfile(sum_ps),  delete(sum_ps);  end

f = fig_sig_count_vs_freq(summary, animals, dv_types);
print(f, '-dpsc', '-painters', sum_ps);
close(f);
[conv_status, conv_msg] = system(sprintf('ps2pdf "%s" "%s"', sum_ps, sum_pdf));
if conv_status == 0 && isfile(sum_pdf)
    delete(sum_ps);
    fprintf('Saved %s\n', sum_pdf);
else
    warning('ps2pdf failed for summary (status %d): %s\nLeaving PS at %s', ...
        conv_status, strtrim(conv_msg), sum_ps);
end

fprintf('\nDone.\n');

%% =====================================================================
%% Local helpers — figure builders
%% =====================================================================
function fig = fig_sig_count_vs_freq(summary, animals, dv_types)
% Cross-DV summary: # FWER-significant channels at each frequency, one
% curve per (animal, DV). Color = DV, line style = animal.

nA = numel(animals);
nD = numel(dv_types);
fig = figure('Visible','off','Units','centimeters','Position',[1 1 18 12]);
set(fig,'PaperUnits','centimeters', ...
    'PaperSize',fig.Position(3:4), ...
    'PaperPosition',[0 0 fig.Position(3:4)]);
ax = axes('Parent', fig); hold(ax,'on');

dv_colors  = lines(nD);
animal_lty = {'-', '--', ':', '-.'};
plotted    = false;

for d = 1:nD
    for a = 1:nA
        S = summary{a, d};
        if isempty(S) || ~isfield(S,'freq') || isempty(S.freq), continue; end
        plot(ax, S.freq, S.n_sig_f, ...
            'Color', dv_colors(d,:), ...
            'LineStyle', animal_lty{min(a, numel(animal_lty))}, ...
            'LineWidth', 1.4, ...
            'DisplayName', sprintf('%s — %s', S.animal, upper(S.dv)));
        plotted = true;
    end
end

if ~plotted
    text(0.5, 0.5, 'No summary data available', ...
        'Units','normalized','HorizontalAlignment','center', ...
        'Parent', ax, 'FontSize', 11);
    set(ax,'XTick',[],'YTick',[]); box(ax,'on');
    title(ax, 'Sig channel count vs frequency','FontWeight','bold');
    return
end

xlabel(ax, 'Frequency (Hz)');
ylabel(ax, '# channels significant (FWER-corrected)');
title(ax,  '# sig channels vs frequency — all (animal, DV)', ...
    'FontWeight','bold');
grid(ax,'on'); box(ax,'on');
legend(ax,'Location','bestoutside','FontSize',8);
end



function fig = fig_phase_trajectory(pref_phase, sig_mask, positions, sgtitle_str)
% Main figure: preferred-phase trajectories vs stimulus position.
%
% pref_phase : nCh × nFreq × nPos, radians in (-π, π]
% sig_mask   : nCh × nFreq logical, true where (ch, f) is significant
% positions  : nPos × 1 numeric, plotted on x-axis
%
% For every significant (ch, f) pair the per-position phases are
% unwrapped and shifted so position 1 = 0 — this lets the *shape* of
% the trajectory (advance/retreat) be read off without each curve's
% absolute offset cluttering the panel.
% The mean trajectory is the circular mean across the significant set
% at each position, also unwrapped and aligned to 0 at position 1.
% The shaded band is a 95% bootstrap CI on that aligned mean.

[nCh, nFreq, nPos] = size(pref_phase);
fig = figure('Visible','off','Units','centimeters','Position',[1 1 18 13]);
set(fig,'PaperUnits','centimeters', ...
    'PaperSize',fig.Position(3:4), ...
    'PaperPosition',[0 0 fig.Position(3:4)]);
ax  = axes('Parent', fig); hold(ax,'on');

if isempty(sig_mask) || ~any(sig_mask(:))
    text(0.5, 0.5, 'No significant (channel, frequency) pairs', ...
        'Units','normalized', 'HorizontalAlignment','center', ...
        'Parent', ax, 'FontSize', 11);
    set(ax,'XTick',[],'YTick',[]); box(ax,'on');
    title(ax, sgtitle_str, 'FontWeight','bold');
    return
end

[sig_ch, sig_fr] = find(sig_mask);
nSig = numel(sig_ch);

% Pull out wrapped trajectories and build aligned/unwrapped copies.
phi_sig    = nan(nSig, nPos);    % wrapped, for circular mean
phi_indiv  = nan(nSig, nPos);    % unwrapped + aligned, for plotting
for k = 1:nSig
    p = squeeze(pref_phase(sig_ch(k), sig_fr(k), :)).';
    phi_sig(k, :)   = p;
    pu              = unwrap(p);
    phi_indiv(k, :) = pu - pu(1);
end

% Circular mean trajectory across the sig set, position by position.
zbar      = mean(exp(1i .* phi_sig), 1, 'omitnan');
mean_traj = unwrap(angle(zbar));
mean_traj = mean_traj - mean_traj(1);

% 95% bootstrap CI on the aligned circular mean.
nBoot = 500;
rng(2025);
boot = nan(nBoot, nPos);
for b = 1:nBoot
    idx = randi(nSig, nSig, 1);
    zb  = mean(exp(1i .* phi_sig(idx, :)), 1, 'omitnan');
    bm  = unwrap(angle(zb));
    boot(b, :) = bm - bm(1);
end
ci_lo = quantile(boot, 0.025, 1);
ci_hi = quantile(boot, 0.975, 1);

x = 1:nPos;

% CI band first so the lines sit on top.
fill(ax, [x, fliplr(x)], [ci_hi, fliplr(ci_lo)], ...
    [0.20 0.40 0.75], 'FaceAlpha', 0.20, 'EdgeColor','none', ...
    'HandleVisibility','off');

% Thin gray = individual (ch × freq) trajectories.
plot(ax, x, phi_indiv.', 'Color',[0.55 0.55 0.55 0.35], 'LineWidth', 0.6, ...
    'HandleVisibility','off');

% Thick blue = circular mean trajectory.
plot(ax, x, mean_traj, 'Color',[0.10 0.30 0.80], 'LineWidth', 2.4, ...
    'DisplayName','circular mean');

% Reference line at 0.
yline(ax, 0, '-', 'Color',[0.5 0.5 0.5], 'LineWidth', 0.5, ...
    'HandleVisibility','off');

xlim(ax, [0.5 nPos+0.5]);
ax.XTick      = x;
ax.XTickLabel = arrayfun(@(p) num2str(positions(p)), x, 'UniformOutput', false);
xlabel(ax, sprintf('Stimulus position (peripheral \\rightarrow foveal)  [n_{pos} = %d]', nPos));
ylabel(ax, 'Preferred phase, unwrapped & aligned (rad)');
yt = -2*pi:pi/2:2*pi;
ax.YTick      = yt;
ax.YTickLabel = arrayfun(@(v) tex_pi(v), yt, 'UniformOutput', false);
grid(ax,'on'); box(ax,'on');
legend(ax,'Location','best','FontSize',8);

title(ax, sprintf('%s   (n_{sig (ch,freq)} = %d)', sgtitle_str, nSig), ...
    'FontWeight','bold');
end


function s = tex_pi(v)
% Compact π-multiple tick label.
r = v / pi;
if r == 0
    s = '0';
elseif r == round(r)
    if r == 1,      s = '\pi';
    elseif r == -1, s = '-\pi';
    else            s = sprintf('%d\\pi', round(r));
    end
elseif r == 0.5,    s = '\pi/2';
elseif r == -0.5,   s = '-\pi/2';
else                s = sprintf('%.2f\\pi', r);
end
end


function fig = fig_systematicity_grid(R_obs, freq, thresh_chan_max, sgtitle_str)
[nCh, ~] = size(R_obs);
[rows, cols, fig] = make_grid_fig(nCh);

ymax = max([max(R_obs(:),[],'omitnan'); max(thresh_chan_max)]) * 1.05;
if ~isfinite(ymax) || ymax <= 0, ymax = 1; end

for ch = 1:nCh
    ax = subplot(rows, cols, ch); hold(ax,'on');
    R = R_obs(ch,:);
    t = thresh_chan_max(ch);
    if all(isnan(R))
        title(ax, sprintf('ch%d (no data)', ch), 'FontSize', 6);
        set(ax,'XTick',[],'YTick',[]); continue
    end
    if ~isnan(t)
        shade_sig(ax, freq, R(:)' >= t, ymax);
    end
    plot(ax, freq, R, 'Color',[0.20 0.40 0.75], 'LineWidth', 1.0);
    if ~isnan(t)
        yline(ax, t, '--', 'Color',[0.85 0.15 0.15], 'LineWidth', 0.9);
    end
    ylim(ax,[0 ymax]); xlim(ax,[min(freq) max(freq)]);
    title(ax, sprintf('ch%d  t=%.2f', ch, t), 'FontSize', 6);
    ax.FontSize = 5;
    if mod(ch-1, cols) ~= 0, ax.YTickLabel = {}; end
    if ch <= (rows-1)*cols,  ax.XTickLabel = {}; end
end
sgtitle(sgtitle_str, 'FontSize', 10, 'FontWeight', 'bold');
end


function fig = fig_pref_phase_grid(pref_phase, R_obs, sig_mask, freq, nPos, sgtitle_str)
% pref_phase : nCh × nFreq × nPos, radians in (-π, π]
[nCh, ~, ~] = size(pref_phase);
[rows, cols, fig] = make_grid_fig(nCh);

for ch = 1:nCh
    ax = subplot(rows, cols, ch);
    phi = squeeze(pref_phase(ch,:,:));
    if ~isempty(sig_mask)
        phi(~sig_mask(ch,:), :) = NaN;
    end
    h = imagesc(ax, 1:nPos, freq, phi);
    set(h, 'AlphaData', ~isnan(phi));
    set(ax, 'YDir','normal','Color','w');
    colormap(ax, hsv); caxis(ax,[-pi pi]);
    xlim(ax,[0.5 nPos+0.5]); ylim(ax,[min(freq) max(freq)]);

    if ~isempty(sig_mask)
        keep = sig_mask(ch,:);
        Rbar = mean(R_obs(ch, keep), 'omitnan');
    else
        Rbar = mean(R_obs(ch,:), 'omitnan');
    end
    if isnan(Rbar)
        title(ax, sprintf('ch%d (no sig)', ch), 'FontSize', 6);
    else
        title(ax, sprintf('ch%d  R=%.2f', ch, Rbar), 'FontSize', 6);
    end
    ax.XTick      = 1:nPos;
    ax.XTickLabel = arrayfun(@(p) sprintf('p%d',p), 1:nPos, 'UniformOutput', false);
    ax.FontSize   = 5;
    if mod(ch-1, cols) ~= 0, ax.YTickLabel = {}; end
    if ch <= (rows-1)*cols,  ax.XTickLabel = {}; end
end

cb = colorbar('Position',[0.93 0.10 0.012 0.82]);
cb.Ticks = [-pi -pi/2 0 pi/2 pi];
cb.TickLabels = {'-pi','-pi/2','0','pi/2','pi'};
cb.Label.String = 'Preferred phase (rad)';
cb.FontSize = 7;
sgtitle(sgtitle_str, 'FontSize', 10, 'FontWeight', 'bold');
end


function fig = fig_step_hist_grid(step_phase, sig_mask, nBins, sgtitle_str)
% step_phase : nCh × nFreq × (nPos-1), wrapped Δφ
[nCh, ~, ~] = size(step_phase);
[rows, cols, fig] = make_grid_fig(nCh);

for ch = 1:nCh
    tmp_ax = subplot(rows, cols, ch);
    pos    = tmp_ax.Position;
    delete(tmp_ax);
    pax    = polaraxes('Parent', fig, 'Position', pos);

    steps_ch = squeeze(step_phase(ch,:,:));
    if ~isempty(sig_mask)
        steps_ch = steps_ch(sig_mask(ch,:), :);
    end
    vals = steps_ch(:);
    vals = vals(~isnan(vals));

    if isempty(vals)
        title(pax, sprintf('ch%d (n=0)', ch), 'FontSize', 6);
        pax.ThetaTickLabel = {}; pax.RTickLabel = {};
        continue
    end

    polarhistogram(pax, vals, nBins, ...
        'FaceColor',[0.30 0.45 0.75], 'EdgeColor','none', 'FaceAlpha', 0.85);

    zbar     = mean(exp(1i .* vals));
    mean_dir = angle(zbar);
    Rval     = abs(zbar);
    rmax     = pax.RLim(2);
    hold(pax,'on');
    polarplot(pax, [mean_dir mean_dir], [0 rmax*Rval], ...
        'Color',[0.85 0.20 0.20], 'LineWidth', 1.3);
    hold(pax,'off');

    title(pax, sprintf('ch%d  R=%.2f  mean=%+d deg', ...
        ch, Rval, round(rad2deg(mean_dir))), 'FontSize', 6);
    pax.ThetaTickLabel = {}; pax.RTickLabel = {};
    pax.FontSize = 5;
end
sgtitle(sgtitle_str, 'FontSize', 10, 'FontWeight', 'bold');
end


function [rows, cols, fig] = make_grid_fig(nCh)
cols = ceil(sqrt(nCh));
rows = ceil(nCh / cols);
fig  = figure('Visible','off','Units','centimeters', ...
    'Position',[1 1 max(36, 2.4*cols) max(24, 2.4*rows)]);
set(fig,'PaperUnits','centimeters', ...
    'PaperSize',fig.Position(3:4), ...
    'PaperPosition',[0 0 fig.Position(3:4)]);
end


function shade_sig(ax, freq, sig, ymax)
if ~any(sig), return; end
starts = find(diff([0 sig]) ==  1);
ends   = find(diff([sig 0]) == -1);
for k = 1:numel(starts)
    xp = [freq(starts(k)) freq(ends(k)) freq(ends(k)) freq(starts(k))];
    yp = [0 0 ymax ymax];
    fill(ax, xp, yp, [0.85 0.15 0.15], ...
        'FaceAlpha', 0.12, 'EdgeColor','none', 'HandleVisibility','off');
end
end
