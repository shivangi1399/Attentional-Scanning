%% Compare Preferred Phase (H3: abs per stimulus position x difficulty bin)
%
% H3 does not save coh_complex per channel either — coherence preferred phase
% is computed here from ph_all_sess.mat using the same per-cell grouping
% (position x difficulty bin) as the H3 analysis scripts:
%
%   For each channel ch, frequency f, DV:
%     For each (position p, difficulty bin d):
%       coh_complex_pd = mean(DV_pd .* exp(i*phase_pd))   [MUA/LFP/RT]
%                      = mean(exp(i*phase_inv_pd))         [hit_miss]
%     phi(ch,f) = angle(circmean over (p,d) of coh_complex_pd)
%
% Significance uses saved per-channel coh_perm_pos_diff (magnitudes).
% Regression preferred phase loads phi_pref from
%   multi_regression_channelwise_R2_abs_per_pos_diff.mat (same shape as H1/H2).
%
% Figures:
%   Fig 1   — Preferred phase heatmaps (channels x freq), magnitude-masked
%   Fig 2   — Polar histograms at the most-significant frequencies
%   Fig 3   — Pairwise circular-circular consistency across DVs
%   Fig 4   — Monkey-average preferred phase heatmap

clearvars; close all; clc

%% 1. SETTINGS

animal       = 'hermes';   % 'hermes' or 'klecks'
nDiffBins    = 4;
nonsig_alpha = 0.2;

%% 2. PATHS & DEPENDENCIES

addpath /opt/fieldtrip_github/
ft_defaults
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/Phase_coherence/functions
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/Correlation_analysis/functions
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/multiple_linear_reg/functions
addpath /mnt/hpc/projects/MWSampling/4Shivangi/software_folder/CircStat2012a
clc

base_results = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animal]);

coh_root  = fullfile(base_results, 'phase_coherence',   'abs_per_pos_diff', 'cp10_till_100');
corr_root = fullfile(base_results, 'phase_correlation', 'abs_per_pos_diff', 'cp10_till_100');
reg_root  = fullfile(base_results, 'multi_lin_reg',     'abs_per_pos_diff', 'cp10_till_100');
data_root = fullfile(base_results, 'multi_lin_reg', 'cp10_till_100');   % shared input data (frequency.mat, ph_all_sess.mat)

save_root = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/sampling_compare_abs_per_pos_diff', animal);
if ~exist(save_root, 'dir'), mkdir(save_root); end

load(fullfile(data_root, 'frequency.mat'));
freq  = frequency;
nFreq = numel(freq);
nCh   = 64;

%% 3. COMPUTE COHERENCE PREFERRED PHASE FROM RAW DATA (H3 per (pos x diff))

fprintf('Computing coherence preferred phase for %s (H3)...\n', animal);

tmp_ph = load(fullfile(data_root, 'ph_all_sess.mat'), 'ph_comb');
a_ph   = tmp_ph.ph_comb;

coh_phase      = struct();   % preferred phase per (ch, freq) for each measure
coh_measures   = {'mua', 'lfp', 'RT', 'hit_miss'};

dv_specs = struct( ...
    'mua',      struct('field','MUA_ERP_ampl_all', 'ti','MUA_ERP_trialinfo'), ...
    'lfp',      struct('field','LFP_ERP_ampl_all', 'ti','LFP_ERP_trialinfo'), ...
    'RT',       struct('field','RT',               'ti','RT_trialinfo'), ...
    'hit_miss', struct('field','',                 'ti','trialinfo'));

% For each measure, build per-trial bin labels and (p,d) cell list, then
% compute the preferred phase as circmean(angle(coh_complex_cell)) across cells.
for m = 1:length(coh_measures)
    key = coh_measures{m};
    spec = dv_specs.(key);
    ti   = a_ph.(spec.ti);
    fprintf('  %s...\n', key);

    if strcmp(key,'RT')
        keep_idx  = find(ti(:,20) == 1);   % hit trials only
    elseif strcmp(key,'hit_miss')
        keep_idx  = find(ti(:,20) == 1 | ti(:,20) == 5);
    else
        keep_idx  = (1:size(ti,1))';
    end

    pos_labels = ti(keep_idx, 16);
    positions  = unique(pos_labels);
    nPos       = numel(positions);

    diff_bin = bin_difficulty_per_pos(ti(keep_idx,18), pos_labels, positions, nDiffBins);

    cell_pos = repmat(positions(:),1,nDiffBins); cell_pos = cell_pos(:);
    cell_dif = repmat(1:nDiffBins, nPos,1);      cell_dif = cell_dif(:);
    nCell    = numel(cell_pos);

    if strcmp(key,'hit_miss')
        hit_labels = (ti(keep_idx,20) == 1);
    end

    phase_map = NaN(nCh, nFreq);

    for ch = 1:nCh
        for f = 1:nFreq
            cplx_cell = NaN(nCell, 1);
            for c = 1:nCell
                mask  = (pos_labels == cell_pos(c)) & (diff_bin == cell_dif(c));
                if sum(mask) < 2, continue; end
                idx_g = keep_idx(mask);
                ph_c  = a_ph.phase_all(idx_g, f, ch);

                if strcmp(key,'hit_miss')
                    miss_c = ~hit_labels(mask);
                    ph_inv = ph_c;
                    ph_inv(miss_c) = mod(ph_c(miss_c)+pi, 2*pi) - pi;
                    cplx_cell(c) = mean(exp(1i * ph_inv));
                else
                    dv_c = a_ph.(spec.field)(idx_g, ch);
                    ok   = ~isnan(ph_c) & ~isnan(dv_c);
                    if sum(ok) < 2, continue; end
                    cplx_cell(c) = mean(dv_c(ok) .* exp(1i * ph_c(ok)));
                end
            end

            phi_cell  = angle(cplx_cell);
            valid     = ~isnan(phi_cell);
            if any(valid)
                phase_map(ch,f) = angle(mean(exp(1i * phi_cell(valid))));
            end
        end
    end
    coh_phase.(key) = phase_map;
end

%% 4. LOAD REGRESSION PREFERRED PHASE

reg_file  = fullfile(reg_root, 'multi_regression_channelwise_R2_abs_per_pos_diff.mat');
has_reg   = exist(reg_file, 'file');
reg_phase = struct();

reg_dvs    = {'MUA_ERP_ampl_all', 'LFP_ERP_ampl_all', 'RT', 'hit_miss'};
reg_labels = {'mua', 'lfp', 'RT', 'hit_miss'};

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
    for m = 1:length(reg_labels)
        reg_phase.(reg_labels{m}) = NaN(nCh, nFreq);
    end
end

%% 5. SIGNIFICANCE MASKS (per-channel coh_perm_pos_diff)

coh_sig = struct();
for m = 1:length(coh_measures)
    measure = coh_measures{m};
    sig_map = false(nCh, nFreq);

    if strcmp(measure, 'hit_miss')
        val_folder = fullfile(corr_root, 'hit_miss', 'all_loc_difflev');
        val_file   = 'itc.mat';   val_var  = 'itc';
        perm_file  = 'itc_perm_pos_diff.mat'; perm_var = 'itc_perm_pos_diff';
    else
        val_folder = fullfile(coh_root, measure, 'all_loc_difflev');
        val_file   = 'coherence.mat'; val_var  = 'coh';
        perm_file  = 'coh_perm_pos_diff.mat'; perm_var = 'coh_perm_pos_diff';
    end

    for ch = 1:nCh
        ch_folder = fullfile(val_folder, num2str(ch));
        vf = fullfile(ch_folder, val_file);
        pf = fullfile(ch_folder, perm_file);
        if ~exist(vf,'file') || ~exist(pf,'file'), continue; end

        tmp = load(vf, val_var);    val = tmp.(val_var);
        tmp = load(pf, perm_var);   prm = tmp.(perm_var);
        if any(isnan(val)) || any(isnan(prm(:))), continue; end

        thr = quantile(max(prm,[],2), 0.95);
        sig_map(ch,:) = val >= thr;
    end

    coh_sig.(measure) = sig_map;
end

reg_sig = struct();
if has_reg
    nFreq_r = length(freqs_reg);
    for m = 1:length(reg_dvs)
        dv        = reg_dvs{m};
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
    for m = 1:length(reg_labels)
        reg_sig.(reg_labels{m}) = false(nCh, nFreq);
    end
end

%% 6. CIRCULAR PASTEL COLORMAP

n_cmap   = 256;
hue      = linspace(0, 1, n_cmap+1)'; hue = hue(1:end-1);
sat      = 0.45 * ones(n_cmap, 1);
val_cmap = 0.95 * ones(n_cmap, 1);
cmap_circ = hsv2rgb([hue, sat, val_cmap]);

%% 7. SHARED LABELS

row_labels = {'MUA', 'LFP', 'RT', 'Hit/Miss'};
row_keys   = {'mua', 'lfp', 'RT', 'hit_miss'};
nDVs       = length(row_keys);

col_labels    = {'Coherence (H3 per pos x diff)', 'Regression (\phi_{pref})'};
coh_subtitles = {'MUA (per-cell weighted phase)', 'LFP (per-cell weighted phase)', ...
    'RT (per-cell weighted phase)', 'Hit/Miss (per-cell ITC)'};
reg_subtitles = {'MUA (atan2(\beta_{sin},\beta_{cos}))', ...
    'LFP (atan2(\beta_{sin},\beta_{cos}))', ...
    'RT (atan2(\beta_{sin},\beta_{cos}))', ...
    'Hit/Miss (atan2(\beta_{sin},\beta_{cos}))'};

%% FIGURE 1 — Preferred Phase Heatmaps (channels x frequency)

f1 = figure('Name', ['Preferred Phase (Abs-Per-Pos-Diff) - ' animal], ...
    'Units', 'centimeters', 'Position', [1 1 36 40]);
set(f1, 'PaperUnits', 'centimeters', 'PaperSize', [36 40], 'PaperPosition', [0 0 36 40]);

for row = 1:nDVs
    key = row_keys{row};

    subplot(nDVs, 2, (row-1)*2 + 1);
    data_coh = coh_phase.(key);
    h_img    = imagesc(freq, 1:nCh, data_coh);
    set(gca, 'YDir', 'normal', 'Color', [1 1 1]);
    colormap(gca, cmap_circ); caxis([-pi pi]);

    alpha_coh = ones(nCh, nFreq) * nonsig_alpha;
    alpha_coh(coh_sig.(key)) = 1;
    set(h_img, 'AlphaData', alpha_coh);

    xlabel('Frequency (Hz)'); ylabel('Channel');
    title(coh_subtitles{row}, 'FontSize', 9);
    set(gca, 'FontSize', 8, 'Box', 'on');
    if row == 1
        text(0.5, 1.22, col_labels{1}, 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
    end

    subplot(nDVs, 2, (row-1)*2 + 2);
    if has_reg
        data_reg = reg_phase.(key);
        nR = size(data_reg,1); nF = size(data_reg,2);
        h_img2 = imagesc(freqs_reg, 1:nR, data_reg);
        set(gca, 'YDir', 'normal', 'Color', [1 1 1]);
        colormap(gca, cmap_circ); caxis([-pi pi]);

        alpha_reg = ones(nR, nF) * nonsig_alpha;
        alpha_reg(reg_sig.(key)(1:nR,1:nF)) = 1;
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

cb = colorbar('Location', 'southoutside');
cb.Ticks = [-pi -pi/2 0 pi/2 pi];
cb.TickLabels = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};
cb.Position = [0.25 0.02 0.5 0.015];
cb.Label.String = 'Preferred Phase (rad)';
sgtitle(sprintf(['Preferred Phase — Abs-Per-Pos-Diff (%s)\n' ...
    'Significance = COHERENCE MAGNITUDE (per-cell arrow length)'], animal), ...
    'FontSize', 13, 'FontWeight', 'bold');
print(f1, fullfile(save_root, 'preferred_phase_comparison.pdf'), '-dpdf');
fprintf('Figure 1 saved.\n');

%% FIGURE 2 — Polar Histograms at Top-significant Frequencies

nTopFreq = 3;
top_freq_idx = cell(nDVs, 1);
for row = 1:nDVs
    key = row_keys{row};
    n_sig_per_freq = sum(coh_sig.(key), 1);
    [~, sorted_idx] = sort(n_sig_per_freq, 'descend');
    top_freq_idx{row} = sorted_idx(1:min(nTopFreq, length(sorted_idx)));
end

f2 = figure('Name', ['Polar Histograms (Abs-Per-Pos-Diff) - ' animal], ...
    'Units', 'centimeters', 'Position', [1 1 52 38]);
set(f2, 'PaperUnits', 'centimeters', 'PaperSize', [52 38], 'PaperPosition', [0 0 52 38]);

nBins = 18;

for row = 1:nDVs
    key = row_keys{row};

    for fi = 1:nTopFreq
        fidx       = top_freq_idx{row}(fi);
        f_hz       = freq(fidx);
        sig_ch_coh = find(coh_sig.(key)(:,fidx));

        sp_idx = (row-1)*nTopFreq*2 + (fi-1)*2 + 1;
        ax = subplot(nDVs, nTopFreq*2, sp_idx, polaraxes);
        if ~isempty(sig_ch_coh)
            phases_coh = coh_phase.(key)(sig_ch_coh, fidx);
            phases_coh = phases_coh(~isnan(phases_coh));
            if ~isempty(phases_coh)
                polarhistogram(ax, phases_coh, nBins, ...
                    'FaceColor', [0.55 0.83 0.78], 'FaceAlpha', 0.7, 'EdgeColor', [0.3 0.6 0.55]);
                hold(ax, 'on');
                mu = circ_mean(phases_coh); R = circ_r(phases_coh);
                polarplot(ax, [mu mu], [0 R*max(ax.RLim)], 'k-', 'LineWidth', 2.5);
            end
        end
        title(ax, sprintf('%s Coh %.0fHz (n=%d)', row_labels{row}, f_hz, length(sig_ch_coh)), 'FontSize', 7);
        ax.FontSize = 6;

        ax2 = subplot(nDVs, nTopFreq*2, sp_idx+1, polaraxes);
        if has_reg
            [~, ridx]  = min(abs(freqs_reg - f_hz));
            sig_ch_reg = find(reg_sig.(key)(:,ridx));
            if ~isempty(sig_ch_reg)
                phases_reg = reg_phase.(key)(sig_ch_reg, ridx);
                phases_reg = phases_reg(~isnan(phases_reg));
                if ~isempty(phases_reg)
                    polarhistogram(ax2, phases_reg, nBins, ...
                        'FaceColor', [0.78 0.68 0.90], 'FaceAlpha', 0.7, 'EdgeColor', [0.55 0.35 0.70]);
                    hold(ax2, 'on');
                    mu2 = circ_mean(phases_reg); R2 = circ_r(phases_reg);
                    polarplot(ax2, [mu2 mu2], [0 R2*max(ax2.RLim)], 'k-', 'LineWidth', 2.5);
                end
            end
            title(ax2, sprintf('%s Reg %.0fHz (n=%d)', row_labels{row}, f_hz, length(sig_ch_reg)), 'FontSize', 7);
        else
            title(ax2, sprintf('%s Reg %.0fHz (no data)', row_labels{row}, f_hz), 'FontSize', 7);
        end
        ax2.FontSize = 6;
    end
end

sgtitle(sprintf('Preferred Phase Distribution — Abs-Per-Pos-Diff (%s)\nBlack arrow = circular mean', animal), ...
    'FontSize', 12, 'FontWeight', 'bold');
print(f2, fullfile(save_root, 'preferred_phase_polar.pdf'), '-dpdf');
fprintf('Figure 2 saved.\n');

%% FIGURE 3 — Pairwise Phase Consistency Across DVs

measure_pairs = {
    'mua', 'lfp',      'MUA vs LFP';
    'mua', 'RT',       'MUA vs RT';
    'mua', 'hit_miss', 'MUA vs Hit/Miss';
    'lfp', 'RT',       'LFP vs RT';
    'lfp', 'hit_miss', 'LFP vs Hit/Miss';
    'RT',  'hit_miss', 'RT vs Hit/Miss';
    };
nPairs = size(measure_pairs, 1);

rho_coh  = NaN(nPairs, nFreq); pval_coh = NaN(nPairs, nFreq);
for p = 1:nPairs
    k1 = measure_pairs{p,1}; k2 = measure_pairs{p,2};
    ph1 = coh_phase.(k1);    ph2 = coh_phase.(k2);
    sg1 = coh_sig.(k1);      sg2 = coh_sig.(k2);
    for fi = 1:nFreq
        valid = sg1(:,fi) & sg2(:,fi) & ~isnan(ph1(:,fi)) & ~isnan(ph2(:,fi));
        if sum(valid) < 4, continue; end
        [rho_coh(p,fi), pval_coh(p,fi)] = circ_corrcc(ph1(valid,fi), ph2(valid,fi));
    end
end

rho_reg = NaN(nPairs, nFreq); pval_reg = NaN(nPairs, nFreq);
if has_reg
    nFreq_r = length(freqs_reg);
    rho_reg = NaN(nPairs,nFreq_r); pval_reg = NaN(nPairs,nFreq_r);
    for p = 1:nPairs
        k1 = measure_pairs{p,1}; k2 = measure_pairs{p,2};
        ph1 = reg_phase.(k1);    ph2 = reg_phase.(k2);
        sg1 = reg_sig.(k1);      sg2 = reg_sig.(k2);
        nR  = min([size(ph1,1),size(ph2,1),size(sg1,1),size(sg2,1)]);
        for fi = 1:nFreq_r
            valid = sg1(1:nR,fi) & sg2(1:nR,fi) & ~isnan(ph1(1:nR,fi)) & ~isnan(ph2(1:nR,fi));
            if sum(valid) < 4, continue; end
            [rho_reg(p,fi), pval_reg(p,fi)] = circ_corrcc(ph1(valid,fi), ph2(valid,fi));
        end
    end
end

pair_colors = [0.00 0.55 0.55; 0.55 0.35 0.70; 0.90 0.40 0.30; 0.30 0.60 0.45; 0.85 0.60 0.20; 0.40 0.40 0.70];

f3 = figure('Name', ['Phase Consistency (Abs-Per-Pos-Diff) - ' animal], ...
    'Units', 'centimeters', 'Position', [1 1 42 36]);
set(f3, 'PaperUnits', 'centimeters', 'PaperSize', [42 36], 'PaperPosition', [0 0 42 36]);

for p = 1:nPairs
    subplot(nPairs, 2, (p-1)*2+1); hold on;
    plot(freq, rho_coh(p,:), 'Color', pair_colors(p,:), 'LineWidth', 2);
    sig_f = find(pval_coh(p,:) < 0.05);
    if ~isempty(sig_f), plot(freq(sig_f), rho_coh(p,sig_f), '.', 'Color', pair_colors(p,:), 'MarkerSize', 12); end
    yline(0,'k--','LineWidth',0.5);
    xlabel('Frequency (Hz)'); ylabel('\rho_{circ}');
    title(sprintf('Coh: %s', measure_pairs{p,3}), 'FontSize', 9);
    ylim([-1 1]); set(gca,'FontSize',8,'Box','on');
    if p == 1
        text(0.5,1.22,'Coherence','Units','normalized','HorizontalAlignment','center','FontSize',13,'FontWeight','bold');
    end

    subplot(nPairs, 2, (p-1)*2+2); hold on;
    if has_reg
        plot(freqs_reg, rho_reg(p,:), 'Color', pair_colors(p,:), 'LineWidth', 2);
        sig_f_r = find(pval_reg(p,:) < 0.05);
        if ~isempty(sig_f_r), plot(freqs_reg(sig_f_r), rho_reg(p,sig_f_r), '.', 'Color', pair_colors(p,:), 'MarkerSize', 12); end
    end
    yline(0,'k--','LineWidth',0.5);
    xlabel('Frequency (Hz)'); ylabel('\rho_{circ}');
    title(sprintf('Reg: %s', measure_pairs{p,3}), 'FontSize', 9);
    ylim([-1 1]); set(gca,'FontSize',8,'Box','on');
    if p == 1
        text(0.5,1.22,'Regression','Units','normalized','HorizontalAlignment','center','FontSize',13,'FontWeight','bold');
    end
end

sgtitle(sprintf('Pairwise Phase Consistency — Abs-Per-Pos-Diff (%s)\nDots = p<0.05, only sig channels in both', animal), ...
    'FontSize', 12, 'FontWeight', 'bold');
print(f3, fullfile(save_root, 'phase_consistency_across_measures.pdf'), '-dpdf');
fprintf('Figure 3 saved.\n');

%% FIGURE 4 — Monkey-average preferred phase (channel circular mean per animal,
%             then circular mean across animals)

animals_all = {'hermes', 'klecks'};
nAnimals    = numel(animals_all);

monkey_save_root = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/sampling_compare_abs_per_pos_diff/monkey_avg');
if ~exist(monkey_save_root, 'dir'), mkdir(monkey_save_root); end

animal_coh_avg = cell(nAnimals,1);
animal_reg_avg = cell(nAnimals,1);
freqs_reg_a    = freqs_reg;

for a = 1:nAnimals
    animalName = animals_all{a};
    fprintf('Loading %s for monkey-average...\n', animalName);

    a_base      = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName]);
    a_data_root = fullfile(a_base, 'multi_lin_reg', 'cp10_till_100');   % shared input data
    a_reg_root  = fullfile(a_base, 'multi_lin_reg', 'abs_per_pos_diff', 'cp10_till_100');
    a_coh_root  = fullfile(a_base, 'phase_coherence',   'abs_per_pos_diff', 'cp10_till_100');

    tmp_d  = load(fullfile(a_data_root, 'ph_all_sess.mat'), 'ph_comb');
    a_phd  = tmp_d.ph_comb;
    tmp_f  = load(fullfile(a_data_root, 'frequency.mat'));
    a_freq = tmp_f.frequency; a_nFreq = numel(a_freq); a_nCh = 64;

    a_coh_phase = struct();
    for m = 1:length(coh_measures)
        key = coh_measures{m};
        spec = dv_specs.(key);
        ti   = a_phd.(spec.ti);
        if strcmp(key,'RT')
            keep_idx = find(ti(:,20) == 1);
        elseif strcmp(key,'hit_miss')
            keep_idx = find(ti(:,20) == 1 | ti(:,20) == 5);
        else
            keep_idx = (1:size(ti,1))';
        end
        pos_lab  = ti(keep_idx, 16);
        a_pos    = unique(pos_lab); a_nPos = numel(a_pos);
        d_bin    = bin_difficulty_per_pos(ti(keep_idx,18), pos_lab, a_pos, nDiffBins);
        c_pos    = repmat(a_pos(:),1,nDiffBins); c_pos = c_pos(:);
        c_dif    = repmat(1:nDiffBins, a_nPos,1); c_dif = c_dif(:);
        a_nCell  = numel(c_pos);

        if strcmp(key,'hit_miss'), hit_lab = (ti(keep_idx,20)==1); end

        ph_map = NaN(a_nCh, a_nFreq);
        for ch = 1:a_nCh
            for f = 1:a_nFreq
                cp = NaN(a_nCell,1);
                for c = 1:a_nCell
                    mk = (pos_lab == c_pos(c)) & (d_bin == c_dif(c));
                    if sum(mk) < 2, continue; end
                    idx_g = keep_idx(mk);
                    ph_c  = a_phd.phase_all(idx_g, f, ch);
                    if strcmp(key,'hit_miss')
                        miss_c = ~hit_lab(mk);
                        ph_inv = ph_c;
                        ph_inv(miss_c) = mod(ph_c(miss_c)+pi,2*pi)-pi;
                        cp(c) = mean(exp(1i*ph_inv));
                    else
                        dv_c = a_phd.(spec.field)(idx_g, ch);
                        ok   = ~isnan(ph_c) & ~isnan(dv_c);
                        if sum(ok) < 2, continue; end
                        cp(c) = mean(dv_c(ok) .* exp(1i*ph_c(ok)));
                    end
                end
                phi_c = angle(cp); v = ~isnan(phi_c);
                if any(v), ph_map(ch,f) = angle(mean(exp(1i*phi_c(v)))); end
            end
        end
        a_coh_phase.(key) = ph_map;
    end

    % Channel circular mean per animal
    a_coh_avg = struct();
    for m = 1:length(coh_measures)
        key = coh_measures{m};
        ph  = a_coh_phase.(key);
        avg_ph = NaN(1, a_nFreq);
        for f = 1:a_nFreq
            v = ~isnan(ph(:,f));
            if any(v), avg_ph(f) = angle(mean(exp(1i*ph(v,f)))); end
        end
        a_coh_avg.(key) = avg_ph;
    end
    animal_coh_avg{a} = a_coh_avg;

    % Regression preferred phase per animal -> channel circular mean
    a_reg_file = fullfile(a_reg_root, 'multi_regression_channelwise_R2_abs_per_pos_diff.mat');
    a_reg_avg  = struct();
    if exist(a_reg_file,'file')
        ar = load(a_reg_file,'reg_results'); arr = ar.reg_results;
        for m = 1:length(reg_dvs)
            dv = reg_dvs{m}; key = reg_labels{m};
            if isfield(arr,dv) && isfield(arr.(dv),'phi_pref')
                ph = arr.(dv).phi_pref; nF = size(ph,2);
                avg_ph = NaN(1,nF);
                for f = 1:nF
                    v = ~isnan(ph(:,f));
                    if any(v), avg_ph(f) = angle(mean(exp(1i*ph(v,f)))); end
                end
                a_reg_avg.(key) = avg_ph;
            else
                a_reg_avg.(key) = NaN(1, length(freqs_reg_a));
            end
        end
    else
        for m = 1:length(reg_labels)
            a_reg_avg.(reg_labels{m}) = NaN(1, length(freqs_reg_a));
        end
    end
    animal_reg_avg{a} = a_reg_avg;
end

% Stack to monkey-average (circular mean across animals)
mk_coh = struct(); mk_reg = struct();
for m = 1:length(coh_measures)
    key = coh_measures{m};
    M = NaN(nAnimals, nFreq);
    for a = 1:nAnimals
        v = animal_coh_avg{a}.(key);
        L = min(length(v), nFreq); M(a,1:L) = v(1:L);
    end
    avg = NaN(1,nFreq);
    for f = 1:nFreq
        vv = ~isnan(M(:,f));
        if any(vv), avg(f) = angle(mean(exp(1i*M(vv,f)))); end
    end
    mk_coh.(key) = avg;
end

nFreqR = length(freqs_reg_a);
for m = 1:length(reg_labels)
    key = reg_labels{m};
    M = NaN(nAnimals, nFreqR);
    for a = 1:nAnimals
        v = animal_reg_avg{a}.(key);
        L = min(length(v), nFreqR); M(a,1:L) = v(1:L);
    end
    avg = NaN(1,nFreqR);
    for f = 1:nFreqR
        vv = ~isnan(M(:,f));
        if any(vv), avg(f) = angle(mean(exp(1i*M(vv,f)))); end
    end
    mk_reg.(key) = avg;
end

f4 = figure('Name', 'Monkey-Avg Preferred Phase (Abs-Per-Pos-Diff)', ...
    'Units', 'centimeters', 'Position', [1 1 36 18]);
set(f4, 'PaperUnits', 'centimeters', 'PaperSize', [36 18], 'PaperPosition', [0 0 36 18]);

subplot(1,2,1); hold on;
ms = {'mua','lfp','RT','hit_miss'};
mc = [0.00 0.55 0.55; 0.55 0.35 0.70; 0.90 0.40 0.30; 0.30 0.60 0.45];
for m = 1:length(ms), plot(freq, mk_coh.(ms{m}), 'Color', mc(m,:), 'LineWidth',2); end
yline(0,'k:'); ylim([-pi pi]); xlabel('Frequency (Hz)'); ylabel('\phi (rad)');
legend(row_labels,'Location','best'); title('Coherence — Monkey Avg \phi');

subplot(1,2,2); hold on;
for m = 1:length(ms), plot(freqs_reg_a, mk_reg.(ms{m}), 'Color', mc(m,:), 'LineWidth',2); end
yline(0,'k:'); ylim([-pi pi]); xlabel('Frequency (Hz)'); ylabel('\phi (rad)');
legend(row_labels,'Location','best'); title('Regression — Monkey Avg \phi');

sgtitle('Monkey-average preferred phase (Abs-Per-Pos-Diff)','FontSize',12,'FontWeight','bold');
print(f4, fullfile(monkey_save_root, 'monkey_avg_preferred_phase.pdf'), '-dpdf');
fprintf('Figure 4 (monkey avg) saved.\n');
