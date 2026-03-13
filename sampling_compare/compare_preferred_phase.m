%% Compare preferred phase across coherence and regression

clear all; close all; clc

%% SETTINGS 

animal = 'klecks';   % 'hermes' or 'klecks'

%% PATHS

addpath /opt/fieldtrip_github/
ft_defaults
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/Phase_coherence/functions
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/Correlation_analysis/functions
clc

base_results = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animal]);

coh_root  = fullfile(base_results, 'phase_coherence', 'cp10_till_100');
corr_root = fullfile(base_results, 'phase_correlation', 'cp10_till_100');
reg_root  = fullfile(base_results, 'multi_lin_reg', 'cp10_till_100');

save_root = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/sampling_compare', animal);
if ~exist(save_root, 'dir'), mkdir(save_root); end

load(fullfile(coh_root, 'frequency.mat'));
freq = frequency;
nFreq = numel(freq);
nCh = 64;

%% LOAD COHERENCE PREFERRED PHASE (phase_spec from coherence.mat)

coh_phase = struct();
coh_measures = {'mua', 'lfp', 'RT'};

for m = 1:length(coh_measures)
    measure = coh_measures{m};
    phase_map = NaN(nCh, nFreq);
    

    for ch = 1:nCh
        ch_folder = fullfile(coh_root, measure, 'all_loc_difflev', num2str(ch));
        coh_file = fullfile(ch_folder, 'coherence.mat');
        if ~exist(coh_file, 'file'), continue; end

        tmp = load(coh_file, 'phase_spec');
        if isfield(tmp, 'phase_spec') && ~any(isnan(tmp.phase_spec))
            phase_map(ch,:) = tmp.phase_spec;
        end
    end

    coh_phase.(measure) = phase_map;
end

%% LOAD HIT/MISS ITC PREFERRED PHASE (recompute from raw data)
% ITC files only store magnitude; recompute angle from raw trial phases

ph_data = load(fullfile(reg_root, 'ph_all_sess.mat'), 'ph_comb');
ph_comb = ph_data.ph_comb;

hit_idx  = find(ph_comb.trialinfo(:,20) == 1);
miss_idx = find(ph_comb.trialinfo(:,20) == 5);

itc_phase_map = NaN(nCh, nFreq);
for ch = 1:nCh
    phase_ch = ph_comb.phase_all(:,:,ch);  % trials x freq
    if all(isnan(phase_ch(:))), continue; end

    % ITC with inverted miss phases (same as Corr_RT_hitmiss script)
    phase_inv = phase_ch;
    phase_inv(miss_idx,:) = phase_inv(miss_idx,:) + pi;

    for f = 1:nFreq
        cavg = mean(exp(1i * phase_inv(:,f)));
        itc_phase_map(ch,f) = angle(cavg);
    end
end
coh_phase.hit_miss = itc_phase_map;

%% LOAD REGRESSION PREFERRED PHASE (phi_pref)

reg_file = fullfile(reg_root, 'multi_regression_channelwise_R2.mat');
has_reg = exist(reg_file, 'file');
reg_phase = struct();

if has_reg
    load(reg_file, 'reg_results');

    % Get regression frequency axis
    session_dirs = dir(fullfile(base_results, [animal '_*']));
    if ~isempty(session_dirs)
        freq_file = fullfile(base_results, session_dirs(1).name, ...
            'Phase_analysis', 'hit_miss', '100iter_cut@cp_m10', '1', 'freqpow.mat');
    else
        freq_file = '';
    end
    if ~isempty(freq_file) && exist(freq_file, 'file')
        tmp = load(freq_file);
        freqs_reg = tmp.freqpow.freq;
    else
        freqs_reg = freq;
    end

    reg_dvs    = {'MUA_ERP_ampl_all', 'LFP_ERP_ampl_all', 'RT', 'hit_miss'};
    reg_labels = {'mua', 'lfp', 'RT', 'hit_miss'};

    for m = 1:length(reg_dvs)
        dv = reg_dvs{m};
        if isfield(reg_results, dv) && isfield(reg_results.(dv), 'phi_pref')
            reg_phase.(reg_labels{m}) = reg_results.(dv).phi_pref;
        else
            reg_phase.(reg_labels{m}) = NaN(nCh, length(freqs_reg));
        end
    end
else
    warning('Regression file not found: %s', reg_file);
    freqs_reg = freq;
    for m = 1:4
        reg_phase.(reg_labels{m}) = NaN(nCh, nFreq);
    end
end

%% BUILD SIGNIFICANCE MASKS (per-channel permutation thresholds)

% Coherence sig masks for MUA, LFP, RT (per-channel coh_perm thresholds)
coh_sig = struct();
for m = 1:length(coh_measures)
    measure = coh_measures{m};
    sig_map = false(nCh, nFreq);
    for ch = 1:nCh
        ch_folder = fullfile(coh_root, measure, 'all_loc_difflev', num2str(ch));
        coh_file  = fullfile(ch_folder, 'coherence.mat');
        perm_file = fullfile(ch_folder, 'coh_perm.mat');
        if ~exist(coh_file, 'file') || ~exist(perm_file, 'file'), continue; end

        tmp = load(coh_file, 'coh');
        val = tmp.coh;
        tmp = load(perm_file, 'coh_perm');
        prm = tmp.coh_perm;
        if any(isnan(val)) || any(isnan(prm(:))), continue; end

        tmax = max(prm, [], 2);
        thr  = quantile(tmax, 0.95);
        sig_map(ch,:) = val >= thr;
    end
    coh_sig.(measure) = sig_map;
end

% Hit/Miss sig mask (ITC per-channel permutation thresholds)
sig_hm = false(nCh, nFreq);
for ch = 1:nCh
    ch_folder = fullfile(corr_root, 'hit_miss', 'all_loc_difflev', num2str(ch));
    itc_file  = fullfile(ch_folder, 'itc.mat');
    perm_file = fullfile(ch_folder, 'itc_perm.mat');
    if ~exist(itc_file, 'file') || ~exist(perm_file, 'file'), continue; end

    tmp = load(itc_file, 'itc');
    val = tmp.itc;
    tmp = load(perm_file, 'itc_perm');
    prm = tmp.itc_perm;
    if any(isnan(val)) || any(isnan(prm(:))), continue; end

    tmax = max(prm, [], 2);
    thr  = quantile(tmax, 0.95);
    sig_hm(ch,:) = val >= thr;
end
coh_sig.hit_miss = sig_hm;

% Regression sig masks (per-channel thresholds from reg_results)
reg_sig = struct();
if has_reg
    nFreq_r = length(freqs_reg);
    for m = 1:length(reg_dvs)
        dv = reg_dvs{m};
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
        reg_sig.(reg_labels{m}) = sig_phase;
    end
else
    for m = 1:4
        reg_sig.(reg_labels{m}) = false(nCh, nFreq);
    end
end

%% CIRCULAR PASTEL COLORMAP

% HSV-based pastel circular colormap (wraps around)
n_cmap = 256;
hue = linspace(0, 1, n_cmap+1)'; hue = hue(1:end-1);
sat = 0.45 * ones(n_cmap, 1);   % pastel saturation
val = 0.95 * ones(n_cmap, 1);   % high brightness
cmap_circ = hsv2rgb([hue, sat, val]);

%% FIGURE 1: preferred phase heatmaps

row_labels  = {'MUA', 'LFP', 'RT', 'Hit/Miss'};
row_keys    = {'mua', 'lfp', 'RT', 'hit_miss'};
col_labels  = {'Coherence (phase\_spec)', 'Regression (\phi_{pref})'};

% Per-row subtitle clarifying the measure used
coh_subtitles = {'MUA (amp-weighted phase)', 'LFP (amp-weighted phase)', ...
                 'RT (amp-weighted phase)', 'Hit/Miss (ITC angle)'};
reg_subtitles = {'MUA (atan2(\beta_{sin},\beta_{cos}))', ...
                 'LFP (atan2(\beta_{sin},\beta_{cos}))', ...
                 'RT (atan2(\beta_{sin},\beta_{cos}))', ...
                 'Hit/Miss (atan2(\beta_{sin},\beta_{cos}))'};

f1 = figure('Name', ['Preferred Phase - ' animal], ...
    'Units', 'centimeters', 'Position', [1 1 36 40]);
set(f1, 'PaperUnits', 'centimeters', 'PaperSize', [36 40], 'PaperPosition', [0 0 36 40]);

nRows = 4;
nCols = 2;

for row = 1:nRows
    key = row_keys{row};

    % --- Column 1: Coherence preferred phase ---
    subplot(nRows, nCols, (row-1)*nCols + 1);

    data_coh = coh_phase.(key);
    h_img = imagesc(freq, 1:nCh, data_coh);
    set(gca, 'YDir', 'normal', 'Color', [1 1 1]);  % white background
    colormap(gca, cmap_circ);
    caxis([-pi pi]);

    % Significant = full opacity, non-significant = 10% opacity
    alpha_coh = ones(nCh, nFreq) * 0.2;
    alpha_coh(coh_sig.(key)) = 1;
    set(h_img, 'AlphaData', alpha_coh);

    xlabel('Frequency (Hz)'); ylabel('Channel');
    title(coh_subtitles{row}, 'FontSize', 9);
    set(gca, 'FontSize', 8, 'Box', 'on');

    if row == 1
        text(0.5, 1.22, col_labels{1}, 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
    end

    % --- Column 2: Regression preferred phase ---
    subplot(nRows, nCols, (row-1)*nCols + 2);

    if has_reg
        data_reg = reg_phase.(key);
        nR = size(data_reg,1);
        nF = size(data_reg,2);
        h_img2 = imagesc(freqs_reg, 1:nR, data_reg);
        set(gca, 'YDir', 'normal', 'Color', [1 1 1]);
        colormap(gca, cmap_circ);
        caxis([-pi pi]);

        % Significant = full opacity, non-significant = 10% opacity
        alpha_reg = ones(nR, nF) * 0.2;
        alpha_reg(reg_sig.(key)(1:nR, 1:nF)) = 1;
        set(h_img2, 'AlphaData', alpha_reg);
    end
    xlabel('Frequency (Hz)'); ylabel('Channel');
    title(reg_subtitles{row}, 'FontSize', 9);
    set(gca, 'FontSize', 8, 'Box', 'on');

    if row == 1
        text(0.5, 1.22, col_labels{2}, 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
    end
end

% Single shared colorbar at the bottom
cb = colorbar('Location', 'southoutside');
cb.Ticks = [-pi -pi/2 0 pi/2 pi];
cb.TickLabels = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};
cb.Position = [0.25 0.02 0.5 0.015];  % centered below all subplots
cb.Label.String = 'Preferred Phase (rad)';

sgtitle(sprintf('Preferred Phase (%s)', animal), ...
    'FontSize', 14, 'FontWeight', 'bold');

print(f1, fullfile(save_root, 'preferred_phase_comparison.pdf'), '-dpdf');
fprintf('Figure 1 saved.\n');

%% FIGURE 2: POLAR HISTOGRAMS at key frequencies (sig channels only)
%  For each DV, pick frequencies with most significant channels.
%  Show coherence and regression side by side.

addpath /mnt/hpc/projects/MWSampling/4Shivangi/software_folder/CircStat2012a

% Pick top 3 frequencies per DV (most sig channels in coherence)
nTopFreq = 3;
all_keys = {'mua', 'lfp', 'RT', 'hit_miss'};

top_freq_idx = cell(4,1);
for row = 1:4
    key = row_keys{row};
    n_sig_per_freq = sum(coh_sig.(key), 1);  % 1 x nFreq
    [~, sorted_idx] = sort(n_sig_per_freq, 'descend');
    top_freq_idx{row} = sorted_idx(1:min(nTopFreq, length(sorted_idx)));
end

% 4 rows (DVs) x 6 columns (3 freq x 2 methods)
f2 = figure('Name', ['Polar Histograms - ' animal], ...
    'Units', 'centimeters', 'Position', [1 1 52 38]);
set(f2, 'PaperUnits', 'centimeters', 'PaperSize', [52 38], 'PaperPosition', [0 0 52 38]);

nBins = 18;  % 20-degree bins

for row = 1:4
    key = row_keys{row};

    for fi = 1:nTopFreq
        fidx = top_freq_idx{row}(fi);
        f_hz = freq(fidx);

        % Get significant channels at this frequency
        sig_ch_coh = find(coh_sig.(key)(:, fidx));

        % --- Coherence polar histogram ---
        sp_idx = (row-1)*nTopFreq*2 + (fi-1)*2 + 1;
        ax = subplot(4, nTopFreq*2, sp_idx, polaraxes);

        if ~isempty(sig_ch_coh)
            phases_coh = coh_phase.(key)(sig_ch_coh, fidx);
            phases_coh = phases_coh(~isnan(phases_coh));
            if ~isempty(phases_coh)
                polarhistogram(ax, phases_coh, nBins, ...
                    'FaceColor', [0.55 0.83 0.78], 'FaceAlpha', 0.7, ...
                    'EdgeColor', [0.3 0.6 0.55]);
                hold(ax, 'on');
                % Mean direction arrow
                mu = circ_mean(phases_coh);
                R  = circ_r(phases_coh);
                rlim_val = max(ax.RLim);
                polarplot(ax, [mu mu], [0 R*rlim_val], 'k-', 'LineWidth', 2.5);
            end
        end
        title(ax, sprintf('%s Coh %.0fHz (n=%d)', row_labels{row}, f_hz, length(sig_ch_coh)), ...
            'FontSize', 7);
        ax.FontSize = 6;

        % --- Regression polar histogram ---
        sp_idx2 = sp_idx + 1;
        ax2 = subplot(4, nTopFreq*2, sp_idx2, polaraxes);

        if has_reg
            % Find closest regression freq index
            [~, ridx] = min(abs(freqs_reg - f_hz));
            sig_ch_reg = find(reg_sig.(key)(:, ridx));

            if ~isempty(sig_ch_reg)
                phases_reg = reg_phase.(key)(sig_ch_reg, ridx);
                phases_reg = phases_reg(~isnan(phases_reg));
                if ~isempty(phases_reg)
                    polarhistogram(ax2, phases_reg, nBins, ...
                        'FaceColor', [0.78 0.68 0.90], 'FaceAlpha', 0.7, ...
                        'EdgeColor', [0.55 0.35 0.70]);
                    hold(ax2, 'on');
                    mu2 = circ_mean(phases_reg);
                    R2  = circ_r(phases_reg);
                    rlim_val2 = max(ax2.RLim);
                    polarplot(ax2, [mu2 mu2], [0 R2*rlim_val2], 'k-', 'LineWidth', 2.5);
                end
            end
            title(ax2, sprintf('%s Reg %.0fHz (n=%d)', row_labels{row}, f_hz, length(sig_ch_reg)), ...
                'FontSize', 7);
        else
            title(ax2, sprintf('%s Reg %.0fHz (no data)', row_labels{row}, f_hz), 'FontSize', 7);
        end
        ax2.FontSize = 6;
    end
end

sgtitle(sprintf('Preferred Phase Distribution at Key Frequencies (%s)\nRed arrow = circular mean, n = sig channels', animal), ...
    'FontSize', 12, 'FontWeight', 'bold');

print(f2, fullfile(save_root, 'preferred_phase_polar.pdf'), '-dpdf');

%% FIGURE 3: PHASE CONSISTENCY ACROSS MEASURES
%  Pairwise circular correlation of preferred phase across channels,
%  at each frequency, using only channels significant in BOTH measures.

measure_pairs = {
    'mua', 'lfp',      'MUA vs LFP';
    'mua', 'RT',       'MUA vs RT';
    'mua', 'hit_miss', 'MUA vs Hit/Miss';
    'lfp', 'RT',       'LFP vs RT';
    'lfp', 'hit_miss', 'LFP vs Hit/Miss';
    'RT',  'hit_miss', 'RT vs Hit/Miss';
};

nPairs = size(measure_pairs, 1);

% Compute pairwise circ_corrcc at each frequency (coherence phases)
rho_coh  = NaN(nPairs, nFreq);
pval_coh = NaN(nPairs, nFreq);

for p = 1:nPairs
    k1 = measure_pairs{p,1};
    k2 = measure_pairs{p,2};
    ph1 = coh_phase.(k1);   % nCh x nFreq
    ph2 = coh_phase.(k2);
    sig1 = coh_sig.(k1);
    sig2 = coh_sig.(k2);

    for fi = 1:nFreq
        % Channels significant in both measures
        valid = sig1(:,fi) & sig2(:,fi) & ~isnan(ph1(:,fi)) & ~isnan(ph2(:,fi));
        if sum(valid) < 4, continue; end  % need at least 4 for meaningful corr

        [rho_coh(p,fi), pval_coh(p,fi)] = circ_corrcc(ph1(valid,fi), ph2(valid,fi));
    end
end

% Same for regression phases (using regression sig masks)
rho_reg  = NaN(nPairs, nFreq);
pval_reg = NaN(nPairs, nFreq);

if has_reg
    nFreq_r = length(freqs_reg);
    rho_reg  = NaN(nPairs, nFreq_r);
    pval_reg = NaN(nPairs, nFreq_r);

    for p = 1:nPairs
        k1 = measure_pairs{p,1};
        k2 = measure_pairs{p,2};
        ph1 = reg_phase.(k1);
        ph2 = reg_phase.(k2);
        sig1 = reg_sig.(k1);
        sig2 = reg_sig.(k2);

        nR = min([size(ph1,1), size(ph2,1), size(sig1,1), size(sig2,1)]);

        for fi = 1:nFreq_r
            valid = sig1(1:nR,fi) & sig2(1:nR,fi) & ~isnan(ph1(1:nR,fi)) & ~isnan(ph2(1:nR,fi));
            if sum(valid) < 4, continue; end

            [rho_reg(p,fi), pval_reg(p,fi)] = circ_corrcc(ph1(valid,fi), ph2(valid,fi));
        end
    end
end

% Plot: 6 pairs x 2 columns (coherence, regression)
pair_colors = [
    0.00 0.55 0.55;   % teal
    0.55 0.35 0.70;   % purple
    0.90 0.40 0.30;   % coral
    0.30 0.60 0.45;   % sage
    0.85 0.60 0.20;   % amber
    0.40 0.40 0.70;   % steel blue
];

f3 = figure('Name', ['Phase Consistency - ' animal], ...
    'Units', 'centimeters', 'Position', [1 1 42 36]);
set(f3, 'PaperUnits', 'centimeters', 'PaperSize', [42 36], 'PaperPosition', [0 0 42 36]);

for p = 1:nPairs
    % Coherence
    subplot(nPairs, 2, (p-1)*2 + 1);
    hold on;
    plot(freq, rho_coh(p,:), 'Color', pair_colors(p,:), 'LineWidth', 2);
    % Mark significant correlations (p < 0.05)
    sig_f = find(pval_coh(p,:) < 0.05);
    if ~isempty(sig_f)
        plot(freq(sig_f), rho_coh(p,sig_f), '.', 'Color', pair_colors(p,:), ...
            'MarkerSize', 12);
    end
    yline(0, 'k--', 'LineWidth', 0.5);
    xlabel('Frequency (Hz)'); ylabel('\rho_{circ}');
    title(sprintf('Coh: %s', measure_pairs{p,3}), 'FontSize', 9);
    ylim([-1 1]);
    set(gca, 'FontSize', 8, 'Box', 'on');

    if p == 1
        text(0.5, 1.22, 'Coherence', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
    end

    % Regression
    subplot(nPairs, 2, (p-1)*2 + 2);
    hold on;
    if has_reg
        plot(freqs_reg, rho_reg(p,:), 'Color', pair_colors(p,:), 'LineWidth', 2);
        sig_f_r = find(pval_reg(p,:) < 0.05);
        if ~isempty(sig_f_r)
            plot(freqs_reg(sig_f_r), rho_reg(p,sig_f_r), '.', 'Color', pair_colors(p,:), ...
                'MarkerSize', 12);
        end
    end
    yline(0, 'k--', 'LineWidth', 0.5);
    xlabel('Frequency (Hz)'); ylabel('\rho_{circ}');
    title(sprintf('Reg: %s', measure_pairs{p,3}), 'FontSize', 9);
    ylim([-1 1]);
    set(gca, 'FontSize', 8, 'Box', 'on');

    if p == 1
        text(0.5, 1.22, 'Regression', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
    end
end

sgtitle(sprintf('Pairwise Phase Consistency Across Measures (%s)\nDots = p < 0.05, only sig channels in both measures', animal), ...
    'FontSize', 12, 'FontWeight', 'bold');

print(f3, fullfile(save_root, 'phase_consistency_across_measures.pdf'), '-dpdf');
