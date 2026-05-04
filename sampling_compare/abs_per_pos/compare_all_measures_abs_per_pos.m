% Compare all phase-dependent analyses (H2: abs per stimulus position).
% Identical layout to compare_all_measures.m (complex/) but loads from
% *_abs_per_pos result folders. Values are already real magnitudes so no
% abs() wrapping is needed in the helper functions.

clear all;
close all;
clc

%% SETTINGS

animal = 'klecks';   % 'hermes' or 'klecks'

%% PATHS

addpath /opt/fieldtrip_github/
ft_defaults
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/Phase_coherence/functions
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/Correlation_analysis/functions

base_results = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animal]);

coh_root  = fullfile(base_results, 'phase_coherence',   'abs_per_pos', 'cp10_till_100');
corr_root = fullfile(base_results, 'phase_correlation', 'abs_per_pos', 'cp10_till_100');
reg_root  = fullfile(base_results, 'multi_lin_reg',     'abs_per_pos', 'cp10_till_100');
data_root = fullfile(base_results, 'multi_lin_reg', 'cp10_till_100');   % shared input data (frequency.mat, ph_all_sess.mat)

save_root = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/sampling_compare_abs_per_pos', animal);
if ~exist(save_root, 'dir'), mkdir(save_root); end

load(fullfile(data_root, 'frequency.mat'));
freq = frequency;
nCh  = 64;

%% HELPER: load and average across channels
% H2 values are already magnitudes — mean() instead of abs(mean()) here.

function [avg_val, avg_thr, n_valid, all_vals, all_perms] = load_coh_avg(base_folder, val_file, val_var, perm_file, perm_var, nCh, condition)
coh_all  = [];
perm_all = [];
n_valid  = 0;

for ch = 1:nCh
    ch_folder = fullfile(base_folder, condition, num2str(ch));
    if ~exist(ch_folder, 'dir'), continue; end

    val_path  = fullfile(ch_folder, val_file);
    perm_path = fullfile(ch_folder, perm_file);
    if ~exist(val_path, 'file') || ~exist(perm_path, 'file'), continue; end

    tmp = load(val_path,  val_var);  val = tmp.(val_var);
    tmp = load(perm_path, perm_var); prm = tmp.(perm_var);

    if any(isnan(val)) || any(isnan(prm(:))), continue; end

    coh_all  = [coh_all; val];
    perm_all = cat(3, perm_all, prm);
    n_valid  = n_valid + 1;
end

if isempty(coh_all)
    avg_val = NaN(1, 1); avg_thr = NaN;
    all_vals = coh_all; all_perms = perm_all;
    return
end

avg_val  = mean(coh_all, 1);      % already magnitudes
perm_avg = mean(perm_all, 3);     % average magnitudes across channels
tmax     = max(perm_avg, [], 2);
avg_thr  = quantile(tmax, 0.95);
all_vals = coh_all; all_perms = perm_all;
end

%% HELPER: build per-channel significance map
% H2: values are already magnitudes — no abs() needed.

function sig_map = build_sig_map(base_folder, val_file, val_var, perm_file, perm_var, nCh, condition, nFreq)
sig_map = false(nCh, nFreq);
for ch = 1:nCh
    ch_folder = fullfile(base_folder, condition, num2str(ch));
    if ~exist(ch_folder, 'dir'), continue; end

    val_path  = fullfile(ch_folder, val_file);
    perm_path = fullfile(ch_folder, perm_file);
    if ~exist(val_path, 'file') || ~exist(perm_path, 'file'), continue; end

    tmp = load(val_path,  val_var);  val = tmp.(val_var);
    tmp = load(perm_path, perm_var); prm = tmp.(perm_var);

    if any(isnan(val)) || any(isnan(prm(:))), continue; end

    tmax = max(prm, [], 2);       % already magnitudes
    thr  = quantile(tmax, 0.95);
    sig_map(ch,:) = val >= thr;
end
end

%% HELPER: plot curve with pastel significance shading

function plot_pretty(freq, values, limit, line_color, shade_color, thr_color, ylab)
    sig_freq  = (values >= limit);
    sigonset  = find(conv(sig_freq, [1 -1]) == 1);
    sigoffset = find(conv(sig_freq, [1 -1]) == -1) - 1;

    plot(freq, values, 'Color', line_color, 'LineWidth', 2.5); hold on;

    ylims = ylim;
    for i = 1:min(length(sigonset), length(sigoffset))
        v = [freq(sigonset(i))  ylims(1); ...
             freq(sigoffset(i)) ylims(1); ...
             freq(sigoffset(i)) ylims(2); ...
             freq(sigonset(i))  ylims(2)];
        patch('Faces', [1 2 3 4], 'Vertices', v, ...
            'FaceColor', shade_color, 'FaceAlpha', 0.30, 'EdgeColor', 'none');
    end

    yline(limit, '--', 'Color', thr_color, 'LineWidth', 1.5);
    ylim(ylims);
    xlabel('Frequency (Hz)'); ylabel(ylab);
    set(gca, 'FontSize', 8, 'Box', 'on');
end

%% LOAD ALL DATA

fprintf('Loading coherence data for %s...\n', animal);

[coh_mua, thr_coh_mua, n_coh_mua] = load_coh_avg(...
    fullfile(coh_root, 'mua'), 'coherence.mat', 'coh', ...
    'coh_perm_pos.mat', 'coh_perm_pos', nCh, 'all_loc_difflev');

[coh_lfp, thr_coh_lfp, n_coh_lfp] = load_coh_avg(...
    fullfile(coh_root, 'lfp'), 'coherence.mat', 'coh', ...
    'coh_perm_pos.mat', 'coh_perm_pos', nCh, 'all_loc_difflev');

[coh_rt, thr_coh_rt, n_coh_rt] = load_coh_avg(...
    fullfile(coh_root, 'RT'), 'coherence.mat', 'coh', ...
    'coh_perm_pos.mat', 'coh_perm_pos', nCh, 'all_loc_difflev');

[coh_hm, thr_coh_hm, n_coh_hm] = load_coh_avg(...
    fullfile(corr_root, 'hit_miss'), 'itc.mat', 'itc', ...
    'itc_perm_pos.mat', 'itc_perm_pos', nCh, 'all_loc_difflev');

fprintf('Loading correlation data for %s...\n', animal);

[corr_mua, thr_corr_mua, nc_mua] = load_coh_avg(...
    fullfile(corr_root, 'mua'), 'correlation.mat', 'correlation', ...
    'corr_perm_pos.mat', 'corr_perm_pos', nCh, 'all_loc_difflev');

[corr_lfp, thr_corr_lfp, nc_lfp] = load_coh_avg(...
    fullfile(corr_root, 'lfp'), 'correlation.mat', 'correlation', ...
    'corr_perm_pos.mat', 'corr_perm_pos', nCh, 'all_loc_difflev');

[corr_rt, thr_corr_rt, nc_rt] = load_coh_avg(...
    fullfile(corr_root, 'RT'), 'correlation.mat', 'correlation', ...
    'corr_perm_pos.mat', 'corr_perm_pos', nCh, 'all_loc_difflev');

[corr_hm, thr_corr_hm, nc_hm] = load_coh_avg(...
    fullfile(corr_root, 'hit_miss'), 'pos.mat', 'pos', ...
    'pos_perm_pos.mat', 'pos_perm_pos', nCh, 'all_loc_difflev');

fprintf('Loading regression data for %s...\n', animal);

reg_file = fullfile(reg_root, 'multi_regression_channelwise_R2_abs_per_pos.mat');
has_reg  = exist(reg_file, 'file');
if has_reg
    load(reg_file, 'reg_results');

    session_dirs = dir(fullfile(base_results, [animal '_*']));
    if ~isempty(session_dirs)
        freq_file = fullfile(base_results, session_dirs(1).name, ...
            'Phase_analysis', 'hit_miss', '100iter_cut@cp_m10', '1', 'freqpow.mat');
    else
        freq_file = '';
    end
    if ~isempty(freq_file) && exist(freq_file, 'file')
        tmp = load(freq_file); freqs_reg = tmp.freqpow.freq;
    else
        freqs_reg = freq;
    end

    reg_dvs = {'MUA_ERP_ampl_all', 'LFP_ERP_ampl_all', 'RT', 'hit_miss'};
else
    warning('Regression results file not found: %s', reg_file);
    freqs_reg = freq;
end

%% COLOR SCHEME

line_colors = [
    0.00 0.55 0.55;   % MUA:      dark teal
    0.55 0.35 0.70;   % LFP:      medium purple
    0.90 0.40 0.30;   % RT:       coral
    0.30 0.60 0.45;   % Hit/Miss: sage green
];

shade_colors = [
    0.60 0.90 0.88;   % MUA:      pastel teal
    0.80 0.70 0.95;   % LFP:      pastel lavender
    1.00 0.75 0.65;   % RT:       pastel peach
    0.70 0.92 0.78;   % Hit/Miss: pastel mint
];

thr_colors = repmat([0.40 0.40 0.40], 4, 1);

%% FIGURE 1: average curves with significance shading

row_labels = {'MUA', 'LFP', 'RT', 'Hit/Miss'};
col_labels = {'Coherence', 'Correlation', 'Regression R^2'};

coh_vals    = {coh_mua,   coh_lfp,   coh_rt,   coh_hm};
coh_thrs    = {thr_coh_mua, thr_coh_lfp, thr_coh_rt, thr_coh_hm};
coh_ns      = {n_coh_mua,  n_coh_lfp,  n_coh_rt,  n_coh_hm};
coh_ylabels = {'Coherence', 'Coherence', 'Coherence', 'ITC'};

corr_vals    = {corr_mua,   corr_lfp,   corr_rt,   corr_hm};
corr_thrs    = {thr_corr_mua, thr_corr_lfp, thr_corr_rt, thr_corr_hm};
corr_ns      = {nc_mua,     nc_lfp,     nc_rt,     nc_hm};
corr_ylabels = {'Correlation', 'Correlation', 'Correlation', 'POS'};

f1 = figure('Name', ['Phase Analysis (Abs-Per-Pos) - ' animal], ...
    'Units', 'centimeters', 'Position', [1 1 48 38]);
set(f1, 'PaperUnits', 'centimeters', 'PaperSize', [48 38], 'PaperPosition', [0 0 48 38]);

nRows = 4; nCols = 3;

for row = 1:nRows
    lc = line_colors(row,:); sc = shade_colors(row,:); tc = thr_colors(row,:);

    % --- Column 1: Coherence ---
    subplot(nRows, nCols, (row-1)*nCols + 1); hold on;
    plot_pretty(freq, coh_vals{row}, coh_thrs{row}, lc, sc, tc, coh_ylabels{row});
    title(sprintf('%s (n=%d)', row_labels{row}, coh_ns{row}), 'FontSize', 9);
    if row == 1
        text(0.5, 1.22, col_labels{1}, 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
    end

    % --- Column 2: Correlation ---
    subplot(nRows, nCols, (row-1)*nCols + 2); hold on;
    plot_pretty(freq, corr_vals{row}, corr_thrs{row}, lc, sc, tc, corr_ylabels{row});
    title(sprintf('%s (n=%d)', row_labels{row}, corr_ns{row}), 'FontSize', 9);
    if row == 1
        text(0.5, 1.22, col_labels{2}, 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
    end

    % --- Column 3: Regression R² ---
    subplot(nRows, nCols, (row-1)*nCols + 3); hold on;
    if has_reg
        dv = reg_dvs{row};
        if isfield(reg_results, dv) && isfield(reg_results.(dv), 'channel_avg_R') ...
                && isfield(reg_results.(dv).channel_avg_R, 'phase')
            avg_R2  = reg_results.(dv).channel_avg_R.phase;
            avg_thr = reg_results.(dv).channel_avg_thresh.phase;
            plot_pretty(freqs_reg, avg_R2, avg_thr, lc, sc, tc, 'R^2');
            title(row_labels{row}, 'FontSize', 9);
        else
            title([row_labels{row} ' (no data)'], 'FontSize', 9);
            xlabel('Frequency (Hz)'); ylabel('R^2');
        end
    else
        title([row_labels{row} ' (no data)'], 'FontSize', 9);
        xlabel('Frequency (Hz)'); ylabel('R^2');
    end
    set(gca, 'FontSize', 8, 'Box', 'on');
    if row == 1
        text(0.5, 1.22, col_labels{3}, 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
    end
end

print(f1, fullfile(save_root, 'comparison_all_measures.pdf'), '-dpdf');

%% FIGURE 2: SIGNIFICANCE HEATMAPS

fprintf('Building significance heatmaps...\n');
nFreq = numel(freq);

sig_coh_mua = build_sig_map(fullfile(coh_root,'mua'), 'coherence.mat','coh', 'coh_perm_pos.mat','coh_perm_pos', nCh, 'all_loc_difflev', nFreq);
sig_coh_lfp = build_sig_map(fullfile(coh_root,'lfp'), 'coherence.mat','coh', 'coh_perm_pos.mat','coh_perm_pos', nCh, 'all_loc_difflev', nFreq);
sig_coh_rt  = build_sig_map(fullfile(coh_root,'RT'),  'coherence.mat','coh', 'coh_perm_pos.mat','coh_perm_pos', nCh, 'all_loc_difflev', nFreq);
sig_coh_hm  = build_sig_map(fullfile(corr_root,'hit_miss'), 'itc.mat','itc', 'itc_perm_pos.mat','itc_perm_pos', nCh, 'all_loc_difflev', nFreq);

sig_corr_mua = build_sig_map(fullfile(corr_root,'mua'), 'correlation.mat','correlation', 'corr_perm_pos.mat','corr_perm_pos', nCh, 'all_loc_difflev', nFreq);
sig_corr_lfp = build_sig_map(fullfile(corr_root,'lfp'), 'correlation.mat','correlation', 'corr_perm_pos.mat','corr_perm_pos', nCh, 'all_loc_difflev', nFreq);
sig_corr_rt  = build_sig_map(fullfile(corr_root,'RT'),  'correlation.mat','correlation', 'corr_perm_pos.mat','corr_perm_pos', nCh, 'all_loc_difflev', nFreq);
sig_corr_hm  = build_sig_map(fullfile(corr_root,'hit_miss'), 'pos.mat','pos', 'pos_perm_pos.mat','pos_perm_pos', nCh, 'all_loc_difflev', nFreq);

sig_coh_all  = {sig_coh_mua, sig_coh_lfp, sig_coh_rt, sig_coh_hm};
sig_corr_all = {sig_corr_mua, sig_corr_lfp, sig_corr_rt, sig_corr_hm};

if has_reg
    nFreq_r    = length(freqs_reg);
    sig_reg_all = cell(1,4);
    for row = 1:4
        dv        = reg_dvs{row};
        sig_phase = false(nCh, nFreq_r);
        if isfield(reg_results, dv)
            numCh_r = size(reg_results.(dv).R2_phase, 1);
            for ch = 1:min(numCh_r, length(reg_results.(dv).thresholds))
                if isfield(reg_results.(dv).thresholds(ch), 'thresh_phase')
                    thr = reg_results.(dv).thresholds(ch).thresh_phase;
                    if ~isempty(thr)
                        sig_phase(ch,:) = reg_results.(dv).R2_phase(ch,:) > thr(1);
                    end
                end
            end
        end
        sig_reg_all{row} = sig_phase;
    end
end

f2 = figure('Name', ['Significance Heatmaps (Abs-Per-Pos) - ' animal], ...
    'Units', 'centimeters', 'Position', [1 1 48 38]);
set(f2, 'PaperUnits', 'centimeters', 'PaperSize', [48 38], 'PaperPosition', [0 0 48 38]);

heatmap_titles_coh  = {'Coh: MUA', 'Coh: LFP', 'Coh: RT', 'Coh: Hit/Miss (ITC)'};
heatmap_titles_corr = {'Corr: MUA', 'Corr: LFP', 'Corr: RT', 'Corr: Hit/Miss (POS)'};
heatmap_titles_reg  = {'Reg: MUA', 'Reg: LFP', 'Reg: RT', 'Reg: Hit/Miss'};

for row = 1:4
    lc = line_colors(row,:);

    subplot(4, 3, (row-1)*3 + 1);
    imagesc(freq, 1:nCh, sig_coh_all{row});
    set(gca, 'YDir', 'normal');
    colormap(gca, [0.95 0.95 0.95; lc]); caxis([0 1]);
    xlabel('Freq (Hz)'); ylabel('Channel');
    title(heatmap_titles_coh{row}, 'FontSize', 9);
    set(gca, 'FontSize', 8);
    if row == 1
        text(0.5, 1.22, 'Coherence', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
    end

    subplot(4, 3, (row-1)*3 + 2);
    imagesc(freq, 1:nCh, sig_corr_all{row});
    set(gca, 'YDir', 'normal');
    colormap(gca, [0.95 0.95 0.95; lc]); caxis([0 1]);
    xlabel('Freq (Hz)'); ylabel('Channel');
    title(heatmap_titles_corr{row}, 'FontSize', 9);
    set(gca, 'FontSize', 8);
    if row == 1
        text(0.5, 1.22, 'Correlation', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
    end

    subplot(4, 3, (row-1)*3 + 3);
    if has_reg && ~isempty(sig_reg_all{row})
        imagesc(freqs_reg, 1:size(sig_reg_all{row},1), sig_reg_all{row});
    end
    set(gca, 'YDir', 'normal');
    colormap(gca, [0.95 0.95 0.95; lc]); caxis([0 1]);
    xlabel('Freq (Hz)'); ylabel('Channel');
    title(heatmap_titles_reg{row}, 'FontSize', 9);
    set(gca, 'FontSize', 8);
    if row == 1
        text(0.5, 1.22, 'Regression R^2', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
    end
end

print(f2, fullfile(save_root, 'significance_heatmaps.pdf'), '-dpdf');

%% MONKEY-AVERAGE COMPARISON

fprintf('Loading monkey-average results...\n');

results_combined = '/mnt/hpc/projects/MWSampling/4Shivangi/results_combined';

coh_combined  = fullfile(results_combined, 'phase_coherence',   'abs_per_pos', 'cp10_till_100');
corr_combined = fullfile(results_combined, 'phase_correlation', 'abs_per_pos', 'cp10_till_100');
reg_combined  = fullfile(results_combined, 'multi_lin_reg',     'abs_per_pos', 'cp10_till_100');

monkey_save_root = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/sampling_compare_abs_per_pos/monkey_avg');
if ~exist(monkey_save_root, 'dir'), mkdir(monkey_save_root); end

% --- Coherence monkey-average ---
mk_coh_vals = cell(1,4); mk_coh_thrs = NaN(1,4);

coh_monkey_files = {
    fullfile(coh_combined,  'mua', 'all_loc_difflev', 'monkey_avg_results.mat');
    fullfile(coh_combined,  'lfp', 'all_loc_difflev', 'monkey_avg_results.mat');
    fullfile(coh_combined,  'RT',  'all_loc_difflev', 'monkey_avg_results.mat');
    fullfile(corr_combined, 'hit_miss_itc', 'all_loc_difflev', 'monkey_avg_results_itc.mat');
};

for row = 1:4
    if isfile(coh_monkey_files{row})
        tmp = load(coh_monkey_files{row});
        if row == 4
            mk_coh_vals{row} = tmp.itc_monkey_avg;
            mk_coh_thrs(row) = tmp.thresh_monkey_avg_itc;
        else
            mk_coh_vals{row} = tmp.coh_monkey_avg;
            mk_coh_thrs(row) = tmp.thresh_monkey_avg;
        end
        if ~exist('freq_monkey', 'var'), freq_monkey = tmp.freq; end
    else
        warning('Monkey-avg coherence file not found: %s', coh_monkey_files{row});
    end
end

% --- Correlation monkey-average ---
mk_corr_vals = cell(1,4); mk_corr_thrs = NaN(1,4);

corr_monkey_files = {
    fullfile(corr_combined, 'mua',      'all_loc_difflev', 'monkey_avg_results.mat');
    fullfile(corr_combined, 'lfp',      'all_loc_difflev', 'monkey_avg_results.mat');
    fullfile(corr_combined, 'RT',       'all_loc_difflev', 'monkey_avg_results.mat');
    fullfile(corr_combined, 'hit_miss', 'all_loc_difflev', 'monkey_avg_results_pos.mat');
};

for row = 1:4
    if isfile(corr_monkey_files{row})
        tmp = load(corr_monkey_files{row});
        if row == 4
            mk_corr_vals{row} = tmp.pos_monkey_avg;
            mk_corr_thrs(row) = tmp.thresh_monkey_avg_pos;
        else
            mk_corr_vals{row} = tmp.corr_monkey_avg;
            mk_corr_thrs(row) = tmp.thresh_monkey_avg;
        end
    else
        warning('Monkey-avg correlation file not found: %s', corr_monkey_files{row});
    end
end

% --- Regression monkey-average ---
mk_reg_vals = cell(1,4); mk_reg_thrs = NaN(1,4);
reg_monkey_dvs = {'MUA_ERP_ampl_all', 'LFP_ERP_ampl_all', 'RT', 'hit_miss'};

for row = 1:4
    reg_monkey_file = fullfile(reg_combined, reg_monkey_dvs{row}, 'monkey_avg_results.mat');
    if isfile(reg_monkey_file)
        tmp = load(reg_monkey_file);
        mk_reg_vals{row} = tmp.monkey_avg_obs.phase;
        mk_reg_thrs(row) = tmp.thresh_monkey.phase;
    else
        warning('Monkey-avg regression file not found: %s', reg_monkey_file);
    end
end

%% FIGURE 3: MONKEY-AVERAGE — curves with significance shading

mk_coh_ylabels  = {'Coherence', 'Coherence', 'Coherence', 'ITC'};
mk_corr_ylabels = {'Correlation', 'Correlation', 'Correlation', 'POS'};

f3 = figure('Name', 'Phase Analysis (Abs-Per-Pos) - Monkey Average', ...
    'Units', 'centimeters', 'Position', [1 1 48 38]);
set(f3, 'PaperUnits', 'centimeters', 'PaperSize', [48 38], 'PaperPosition', [0 0 48 38]);

for row = 1:nRows
    lc = line_colors(row,:); sc = shade_colors(row,:); tc = thr_colors(row,:);

    subplot(nRows, nCols, (row-1)*nCols + 1); hold on;
    if ~isempty(mk_coh_vals{row})
        plot_pretty(freq_monkey, mk_coh_vals{row}, mk_coh_thrs(row), lc, sc, tc, mk_coh_ylabels{row});
        title(row_labels{row}, 'FontSize', 9);
    else
        title([row_labels{row} ' (no data)'], 'FontSize', 9);
        xlabel('Frequency (Hz)'); ylabel(mk_coh_ylabels{row});
    end
    if row == 1
        text(0.5, 1.22, col_labels{1}, 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
    end

    subplot(nRows, nCols, (row-1)*nCols + 2); hold on;
    if ~isempty(mk_corr_vals{row})
        plot_pretty(freq_monkey, mk_corr_vals{row}, mk_corr_thrs(row), lc, sc, tc, mk_corr_ylabels{row});
        title(row_labels{row}, 'FontSize', 9);
    else
        title([row_labels{row} ' (no data)'], 'FontSize', 9);
        xlabel('Frequency (Hz)'); ylabel(mk_corr_ylabels{row});
    end
    if row == 1
        text(0.5, 1.22, col_labels{2}, 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
    end

    subplot(nRows, nCols, (row-1)*nCols + 3); hold on;
    if ~isempty(mk_reg_vals{row})
        plot_pretty(freqs_reg, mk_reg_vals{row}, mk_reg_thrs(row), lc, sc, tc, 'R^2');
        title(row_labels{row}, 'FontSize', 9);
    else
        title([row_labels{row} ' (no data)'], 'FontSize', 9);
        xlabel('Frequency (Hz)'); ylabel('R^2');
    end
    set(gca, 'FontSize', 8, 'Box', 'on');
    if row == 1
        text(0.5, 1.22, col_labels{3}, 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
    end
end

print(f3, fullfile(monkey_save_root, 'monkey_avg_comparison_all_measures.pdf'), '-dpdf');
