% =====================================================================
% Phase progression across stimulus positions, per channel
%
% Question: for each recording channel, does the preferred phase (angle
% of complex coherence) depend *systematically* on stimulus position?
%
% Per (channel × frequency):
%   - Complex coherence c_pos = mean(y · exp(i·phase)) across trials in
%     each stimulus position. y is the DV weight (MUA / LFP amplitude,
%     RT, or binary hit_miss). angle(c_pos) = preferred phase; abs(c_pos)
%     = coherence magnitude.
%   - Systematicity statistic: circular-linear correlation r_cl between
%     position index (1..nPos, peripheral→foveal) and the preferred-phase
%     vector (Mardia/Jupp). r_cl ∈ [0,1], unsigned.
%   - Null: shuffle the position label across trials, recompute r_cl,
%     nPerm times → per-channel p-value at each freq. R_obs / R_null hold
%     the r_cl values. The per-channel permutations are SYNCHRONISED
%     (shared perm_seed_base), so the channel-/monkey-average null below
%     is valid for spatially-correlated array LFP.
%
% This runs, per animal × DV, top to bottom:
%   1. submit one SLURM job per channel and aggregate the results,
%   2. save phase_progression.mat + channel_avg_results.mat,
%   3. render the multi-page scanning-figures PDF (pages 1–3 below),
% then after the loop a monkey-average summary and a cross-DV summary.
%
% Figure pages (one PDF per animal × DV):
%   1. systematicity per channel — R_obs(ch,f) vs freq + max-stat 95% line.
%   2. preferred-phase heatmaps per channel (all freqs).
%   3. preferred-phase heatmaps per channel (sig freqs only).
% Plus a cross-DV summary (# FWER-sig channels vs freq) and the
% monkey-average systematicity curve.
%
% Excluded channels: detected data-driven per (animal × DV) — a channel is
% dropped if its observed r_cl is entirely degenerate (all-zero, i.e. a
% dead/constant input, e.g. ch37 in hermes; or all-NaN, i.e. no valid
% trials) or its permutation null is degenerate (all-zero/all-NaN). Bad
% channels are forced non-sig and blanked in the figures AND excluded from
% the channel-/monkey-average curves, so a dead channel can neither be
% plotted nor dilute the pooled statistic. The excluded set is recorded as
% `bad_ch` in channel_avg_results.mat.
%
% NOTE on position ordering: positions are unique(trialinfo col 16) in
% sorted numerical order.
% =====================================================================
clearvars; close all; clc

%% Dependencies
addpath /opt/fieldtrip_github/; ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/Phase_coherence/functions
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/scanning/functions
addpath /mnt/hpc/projects/MWSampling/4Shivangi
clc

%% Settings
base     = '/mnt/hpc/projects/MWSampling/4Shivangi';
dv_types = {'mua','lfp','RT','hit_miss'};   % DV used as the coherence weight
animals  = {'hermes','klecks'};
nPerm    = 1000;
alpha    = 0.05;
rng(2025);

% Accumulator for the cross-DV summary figure built after the main loop.
summary = cell(numel(animals), numel(dv_types));

%% Per animal × DV: compute → save → figures
for a = 1:numel(animals)
    animalName = animals{a};
    fprintf('\n=== %s ===\n', animalName);

    data_load_folder = fullfile(base, ['results_' animalName], 'multi_lin_reg', 'cp10_till_100');
    output_folder    = fullfile(base, ['results_' animalName], 'scanning', 'phase_progression', 'cp10_till_100');
    if ~exist(output_folder,'dir'), mkdir(output_folder); end

    cd(data_load_folder); load('ph_all_sess.mat'); load('frequency.mat');
    freq = frequency;
    [~, nFreq, nCh] = size(ph_comb.phase_all);

    for d = 1:numel(dv_types)
        dv = dv_types{d};
        [~, trlInfo, keepIdx, ~] = get_dv_data(ph_comb, dv);
        fprintf('--- %s (%s) ---\n', animalName, upper(dv));

        % Positions taken from the trials this DV actually uses.
        positions = unique(trlInfo(keepIdx, 16));
        nPos      = numel(positions);
        if nPos < 2
            warning('%s/%s: <2 positions, skipping.', animalName, dv); continue
        end

        dv_dir = fullfile(output_folder, dv);
        if ~exist(dv_dir,'dir'), mkdir(dv_dir); end

        % ── Submit one SLURM job per channel ──────────────────────────
        % seed is per-channel, but perm_seed_base is SHARED so the null is
        % synchronised across channels (see phase_progression_chan.m) — this
        % is what makes the channel-/monkey-average null valid for
        % spatially-correlated array LFP.
        cfg = cell(1, nCh);
        for ch = 1:nCh
            cfg{ch}.ichan          = ch;
            cfg{ch}.nPerm          = nPerm;
            cfg{ch}.dv             = dv;
            cfg{ch}.infile         = data_load_folder;
            cfg{ch}.outfile        = dv_dir;
            cfg{ch}.seed           = 2025 + ch;
            cfg{ch}.perm_seed_base = 2025;      % SHARED across channels
        end
        slurmfun(@phase_progression_chan, cfg, ...
            'partition',     '8GB', ...
            'stopOnError',   false, ...
            'useUserPath',   true);

        % ── Aggregate per-channel results from disk ───────────────────
        pref_phase = nan(nCh, nFreq, nPos);     % φ(ch, f, pos)
        coh_mag    = nan(nCh, nFreq, nPos);     % |c|(ch, f, pos)
        R_obs      = nan(nCh, nFreq);           % observed systematicity (r_cl)
        R_null     = nan(nCh, nFreq, nPerm);    % permutation null
        p_val      = nan(nCh, nFreq);           % per-freq p-value
        n_pos      = zeros(nCh, nPos);          % trials per (ch, pos)

        for ch = 1:nCh
            chan_file = fullfile(dv_dir, num2str(ch), 'phase_progression_chan.mat');
            if ~isfile(chan_file)
                warning('Missing per-channel file for %s/%s ch %d.', animalName, dv, ch);
                continue
            end
            tmp = load(chan_file);
            R_obs(ch, :)         = tmp.R_obs;
            R_null(ch, :, :)     = tmp.R_null;
            pref_phase(ch, :, :) = reshape(tmp.pref_phase, [1 nFreq nPos]);
            coh_mag(ch, :, :)    = reshape(tmp.coh_mag,    [1 nFreq nPos]);
            p_val(ch, :)         = tmp.p_val;
            n_pos(ch, :)         = tmp.n_p;
        end

        save(fullfile(dv_dir, 'phase_progression.mat'), ...
            'pref_phase','coh_mag', ...
            'R_obs','R_null','p_val','positions','freq','n_pos','-v7.3');

        % ── Degenerate-channel detection (data-driven, per animal × DV) ─
        % Drop a channel if its observed r_cl is entirely degenerate
        % (all-zero → dead/constant input; all-NaN → no valid trials) or its
        % permutation null is degenerate (all-zero/all-NaN → no threshold).
        % Applied to BOTH the figures and the channel-/monkey-average below,
        % so a dead channel can neither be plotted nor dilute the pooled stat.
        bad_obs  = all(R_obs  == 0 | isnan(R_obs),  2);            % nCh × 1
        bad_null = all(all(R_null == 0 | isnan(R_null), 3), 2);    % nCh × 1
        bad_ch   = find(bad_obs | bad_null)';
        good_ch  = setdiff(1:nCh, bad_ch);
        if ~isempty(bad_ch)
            fprintf('  Excluding degenerate channels: %s\n', mat2str(bad_ch));
        end

        % ── Per-channel + channel-average significance ─────────────────
        % Per-channel max-stat threshold: max across freq per perm.
        tmax_per_chan   = squeeze(max(R_null, [], 2));         % nCh × nPerm
        thresh_per_chan = quantile(tmax_per_chan, 1-alpha, 2); % nCh × 1

        % Channel-average null: average the (synchronised) per-channel null
        % curves per permutation index over GOOD channels only, then
        % freq-wise max → single max-stat threshold that preserves
        % cross-channel dependence (and is not diluted by dead channels).
        R_null_chan_avg = squeeze(mean(R_null(good_ch,:,:), 1, 'omitnan'));  % nFreq × nPerm
        R_obs_chan_avg  = mean(R_obs(good_ch,:), 1, 'omitnan');              % 1 × nFreq
        tmax_chan_avg   = max(R_null_chan_avg, [], 1);          % 1 × nPerm
        thresh_chan_avg = quantile(tmax_chan_avg, 1-alpha);

        save(fullfile(dv_dir, 'channel_avg_results.mat'), ...
            'R_obs','R_obs_chan_avg','R_null_chan_avg','tmax_chan_avg', ...
            'thresh_chan_avg','thresh_per_chan','freq','positions','bad_ch','good_ch');

        fprintf('Saved %s (%s); channel-avg threshold = %.4f\n', ...
            upper(dv), animalName, thresh_chan_avg);

        %% ── Figures for this (animal, DV) ─────────────────────────────
        % sig_mask uses the per-channel FWER-corrected (max-stat) threshold;
        % degenerate channels (bad_ch, detected above) are forced non-sig.
        sig_mask = bsxfun(@ge, R_obs, thresh_per_chan);   % nCh × nFreq
        sig_mask(bad_ch, :) = false;

        thresh_chan_max = thresh_per_chan;      % FWER threshold for the figures
        thresh_chan_max(bad_ch) = NaN;

        % Blanked COPIES for plotting so the saved/averaged arrays stay intact.
        R_fig   = R_obs;      R_fig(bad_ch, :)      = NaN;
        phi_fig = pref_phase; phi_fig(bad_ch, :, :) = NaN;

        % Stash for the cross-DV summary plot.
        summary{a, d} = struct('animal', animalName, 'dv', dv, ...
            'freq', freq(:)', 'n_sig_f', sum(sig_mask, 1));

        plot_dir = fullfile(base, 'Plots','scanning','phase_progression','cp10_till_100', dv);
        if ~exist(plot_dir,'dir'), mkdir(plot_dir); end
        out_pdf = fullfile(plot_dir, sprintf('%s_scanning_figures.pdf', animalName));
        out_ps  = strrep(out_pdf, '.pdf', '.ps');
        if isfile(out_pdf), delete(out_pdf); end   % start fresh
        if isfile(out_ps),  delete(out_ps);  end

        tag = sprintf('%s — %s', animalName, upper(dv));
        % R2020 has exportgraphics but NOT its 'Append' option (R2022a+),
        % so pages are appended to a PS and converted to PDF at the end.

        % Page 1 — systematicity per channel
        hf = fig_systematicity_grid(R_fig, freq, thresh_chan_max, ...
            sprintf('%s :: R per channel (dashed = max-stat 95%%)', tag));
        print(hf, '-dpsc', '-append', '-painters', out_ps); close(hf);

        % Page 2 — preferred phase per channel, all freqs
        hf = fig_pref_phase_grid(phi_fig, R_fig, [], freq, nPos, ...
            sprintf('%s :: preferred phase per channel (all freqs)', tag));
        print(hf, '-dpsc', '-append', '-painters', out_ps); close(hf);

        % Page 3 — preferred phase per channel, sig freqs only
        hf = fig_pref_phase_grid(phi_fig, R_fig, sig_mask, freq, nPos, ...
            sprintf('%s :: preferred phase per channel (sig freqs, p<%.2f)', tag, alpha));
        print(hf, '-dpsc', '-append', '-painters', out_ps); close(hf);

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

%% Monkey-average systematicity (channel-avg curves across animals)
nAnimals = numel(animals);
for d = 1:numel(dv_types)
    dv = dv_types{d};
    R_animals = []; null_animals = []; thresh_animals = nan(1, nAnimals);
    for a = 1:nAnimals
        avg_file = fullfile(base, ['results_' animals{a}], 'scanning','phase_progression', ...
            'cp10_till_100', dv, 'channel_avg_results.mat');
        if ~isfile(avg_file)
            warning('Channel-avg not found for %s (%s).', dv, animals{a}); continue
        end
        tmp = load(avg_file);
        R_animals    = cat(1, R_animals,    tmp.R_obs_chan_avg);
        null_animals = cat(3, null_animals, tmp.R_null_chan_avg);  % nFreq × nPerm × nAnimals
        thresh_animals(a) = tmp.thresh_chan_avg;
    end
    if size(R_animals,1) < nAnimals
        warning('Skipping monkey-avg for %s (not all animals present).', dv); continue
    end

    R_monkey_avg  = mean(R_animals, 1);
    null_monkey   = mean(null_animals, 3);            % nFreq × nPerm
    tmax_monkey   = max(null_monkey, [], 1);
    thresh_monkey = quantile(tmax_monkey, 1-alpha);
    sig_freqs     = tmp.freq(R_monkey_avg >= thresh_monkey);

    save_dir = fullfile(base, 'results_combined','scanning','phase_progression','cp10_till_100', dv);
    if ~exist(save_dir,'dir'), mkdir(save_dir); end
    save(fullfile(save_dir, 'monkey_avg_results.mat'), ...
        'R_monkey_avg','null_monkey','tmax_monkey','thresh_monkey', ...
        'R_animals','animals','sig_freqs','-v7.3');

    fprintf('%s monkey-avg threshold = %.4f, n sig freqs = %d\n', ...
        upper(dv), thresh_monkey, numel(sig_freqs));

    % Plot: r_cl vs freq with monkey + per-animal thresholds.
    plot_dir = fullfile(base, 'Plots','scanning','phase_progression','cp10_till_100', dv);
    if ~exist(plot_dir,'dir'), mkdir(plot_dir); end
    hf = figure('Visible','off','Position',[100 100 700 380]); hold on;
    plot(tmp.freq, R_monkey_avg, 'k', 'LineWidth', 2, 'DisplayName','Monkey avg');
    cols = lines(nAnimals);
    for a = 1:nAnimals
        plot(tmp.freq, R_animals(a,:), 'Color', cols(a,:), 'LineWidth', 1, ...
            'DisplayName', animals{a});
    end
    yline(thresh_monkey,'-','Color',[0.85 0.15 0.15],'LineWidth',1.4, ...
        'DisplayName','Monkey 95% threshold');
    for a = 1:nAnimals
        if isnan(thresh_animals(a)), continue; end
        yline(thresh_animals(a),'--','Color',cols(a,:),'LineWidth',1.0, ...
            'DisplayName',sprintf('%s 95%% threshold',animals{a}));
    end
    xlabel('Frequency (Hz)'); ylabel('r_{cl} (position–phase circular-linear correlation)');
    title(sprintf('%s — phase progression systematicity (channel avg)', upper(dv)));
    legend('Location','best','FontSize',7); grid on;
    saveas(hf, fullfile(plot_dir,'systematicity_vs_freq.pdf'));
    close(hf);
end

%% Cross-DV summary: # sig channels vs frequency, one curve per (animal, DV)
summary_dir = fullfile(base, 'Plots','scanning','phase_progression','cp10_till_100');
if ~exist(summary_dir,'dir'), mkdir(summary_dir); end
sum_pdf = fullfile(summary_dir, 'summary_sig_count_vs_freq.pdf');
sum_ps  = strrep(sum_pdf, '.pdf', '.ps');
if isfile(sum_pdf), delete(sum_pdf); end
if isfile(sum_ps),  delete(sum_ps);  end

hf = fig_sig_count_vs_freq(summary, animals, dv_types);
print(hf, '-dpsc', '-painters', sum_ps); close(hf);
[conv_status, conv_msg] = system(sprintf('ps2pdf "%s" "%s"', sum_ps, sum_pdf));
if conv_status == 0 && isfile(sum_pdf)
    delete(sum_ps);
    fprintf('Saved %s\n', sum_pdf);
else
    warning('ps2pdf failed for summary (status %d): %s\nLeaving PS at %s', ...
        conv_status, strtrim(conv_msg), sum_ps);
end

fprintf('\nDone. Per-channel results in results_<animal>/scanning/phase_progression/.\n');

%% =====================================================================
%% Local functions
%% =====================================================================
function [Y, trlInfo, keepIdx, isPerCh] = get_dv_data(ph_comb, dv)
% Return the DV values, the trialinfo aligned with phase_all, a keep-mask,
% and whether Y is per-channel.
%   Y       — nTrials × nCh (mua, lfp, RT) or nTrials × 1 (hit_miss)
%   trlInfo — nTrials × nCols, stimulus position in col 16
%   keepIdx — logical(nTrials × 1), trials this DV is defined on
%   isPerCh — true if Y has a per-channel column
switch dv
    case 'mua'
        Y       = ph_comb.MUA_ERP_ampl_all;
        trlInfo = ph_comb.trialinfo;
        keepIdx = true(size(Y,1), 1);
        isPerCh = true;
    case 'lfp'
        Y       = ph_comb.LFP_ERP_ampl_all;
        trlInfo = ph_comb.trialinfo;
        keepIdx = true(size(Y,1), 1);
        isPerCh = true;
    case 'RT'
        Y       = ph_comb.RT;
        trlInfo = ph_comb.RT_trialinfo;
        keepIdx = trlInfo(:,20) == 1;       % hits only
        isPerCh = true;
    case 'hit_miss'
        Y           = ph_comb.trialinfo(:,20);
        Y(Y == 5)   = 0;                    % recode miss code 5 → 0
        trlInfo     = ph_comb.trialinfo;
        keepIdx     = true(size(Y,1), 1);
        isPerCh     = false;
    otherwise
        error('Unknown DV: %s', dv);
end
end


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
title(ax,  '# sig channels vs frequency — all (animal, DV)', 'FontWeight','bold');
grid(ax,'on'); box(ax,'on');
legend(ax,'Location','bestoutside','FontSize',8);
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
