%% Compare Preferred Phase (H2: abs per stimulus position)
%
% H2 does not save coh_complex per channel — coherence preferred phase is
% computed here from ph_all_sess.mat using the same per-position grouping
% as the H2 analysis scripts:
%
%   For each channel ch, frequency f, DV:
%     For each stimulus position p:
%       coh_complex_p = mean(DV_p .* exp(i*phase_p))   [MUA/LFP/RT]
%                     = mean(exp(i*phase_inv_p))         [hit_miss, inverted miss phases]
%     preferred_phase(ch,f) = angle(mean_over_positions(coh_complex_p))
%
% Significance uses saved per-channel coh_perm_pos (magnitudes).
% Regression preferred phase loads phi_pref from
%   multi_regression_channelwise_R2_abs_per_pos.mat (same structure as H1).
%
% Figures produced (same layout as compare_preferred_phase.m):
%   Fig 1  — Preferred phase heatmaps, channels x freq (per animal)
%   Fig 2  — Polar histograms at key frequencies (per animal)
%   Fig 3  — Pairwise phase consistency across DVs (per animal)
%   Fig 4  — Monkey-average preferred phase heatmap
%   Fig 5  — Monkey-average hit-only preferred phase
%   Fig 6  — Monkey-average miss-only preferred phase

clearvars; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

animal      = 'hermes';   % 'hermes' or 'klecks'
nonsig_alpha = 0.2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. PATHS & DEPENDENCIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath /opt/fieldtrip_github/
ft_defaults
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/Phase_coherence/functions
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/Correlation_analysis/functions
addpath /mnt/hpc/projects/MWSampling/4Shivangi/software_folder/CircStat2012a
clc

base_results = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animal]);

coh_root  = fullfile(base_results, 'phase_coherence',   'abs_per_pos', 'cp10_till_100');
corr_root = fullfile(base_results, 'phase_correlation', 'abs_per_pos', 'cp10_till_100');
reg_root  = fullfile(base_results, 'multi_lin_reg',     'abs_per_pos', 'cp10_till_100');
data_root = fullfile(base_results, 'multi_lin_reg', 'cp10_till_100');   % shared input data (ph_all_sess.mat, frequency.mat)

save_root = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/sampling_compare_abs_per_pos', animal);
if ~exist(save_root, 'dir'), mkdir(save_root); end

load(fullfile(data_root, 'frequency.mat'));
freq  = frequency;
nFreq = numel(freq);
nCh   = 64;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. COMPUTE COHERENCE PREFERRED PHASE FROM RAW DATA (H2 per-position)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% H2 analysis scripts do not save coh_complex per channel. We recompute
% the preferred phase here using the same per-position grouping.

fprintf('Computing coherence preferred phase from raw data for %s...\n', animal);

tmp_ph = load(fullfile(data_root, 'ph_all_sess.mat'), 'ph_comb');
a_ph   = tmp_ph.ph_comb;

coh_phase       = struct();
coh_measures    = {'mua', 'lfp', 'RT', 'hit_miss'};

% --- MUA ---
fprintf('  MUA...\n');
positions_mua = unique(a_ph.MUA_ERP_trialinfo(:,16));
nPos_mua      = numel(positions_mua);
phase_map     = NaN(nCh, nFreq);
for ch = 1:nCh
    for f = 1:nFreq
        cplx_pos = NaN(nPos_mua, 1);
        for p = 1:nPos_mua
            mask = a_ph.MUA_ERP_trialinfo(:,16) == positions_mua(p);
            if sum(mask) < 2, continue; end
            ph_p  = a_ph.phase_all(mask, f, ch);
            dv_p  = a_ph.MUA_ERP_ampl_all(mask, ch);
            ok    = ~isnan(ph_p) & ~isnan(dv_p);
            if sum(ok) < 2, continue; end
            cplx_pos(p) = mean(dv_p(ok) .* exp(1i * ph_p(ok)));
        end
        phi_pos   = angle(cplx_pos);
        valid_pos = ~isnan(phi_pos);
        if any(valid_pos)
            phase_map(ch,f) = angle(mean(exp(1i * phi_pos(valid_pos))));
        end
    end
end
coh_phase.mua = phase_map;

% --- LFP ---
fprintf('  LFP...\n');
positions_lfp = unique(a_ph.LFP_ERP_trialinfo(:,16));
nPos_lfp      = numel(positions_lfp);
phase_map     = NaN(nCh, nFreq);
for ch = 1:nCh
    for f = 1:nFreq
        cplx_pos = NaN(nPos_lfp, 1);
        for p = 1:nPos_lfp
            mask = a_ph.LFP_ERP_trialinfo(:,16) == positions_lfp(p);
            if sum(mask) < 2, continue; end
            ph_p  = a_ph.phase_all(mask, f, ch);
            dv_p  = a_ph.LFP_ERP_ampl_all(mask, ch);
            ok    = ~isnan(ph_p) & ~isnan(dv_p);
            if sum(ok) < 2, continue; end
            cplx_pos(p) = mean(dv_p(ok) .* exp(1i * ph_p(ok)));
        end
        phi_pos   = angle(cplx_pos);
        valid_pos = ~isnan(phi_pos);
        if any(valid_pos)
            phase_map(ch,f) = angle(mean(exp(1i * phi_pos(valid_pos))));
        end
    end
end
coh_phase.lfp = phase_map;

% --- RT (hit trials only) ---
fprintf('  RT...\n');
hit_idx_rt   = find(a_ph.RT_trialinfo(:,20) == 1);
positions_rt = unique(a_ph.RT_trialinfo(hit_idx_rt, 16));
nPos_rt      = numel(positions_rt);
phase_map    = NaN(nCh, nFreq);
for ch = 1:nCh
    for f = 1:nFreq
        cplx_pos = NaN(nPos_rt, 1);
        for p = 1:nPos_rt
            pos_mask   = a_ph.RT_trialinfo(hit_idx_rt,16) == positions_rt(p);
            pos_global = hit_idx_rt(pos_mask);
            if sum(pos_mask) < 2, continue; end
            ph_p  = a_ph.phase_all(pos_global, f, ch);
            dv_p  = a_ph.RT(pos_global, ch);
            ok    = ~isnan(ph_p) & ~isnan(dv_p);
            if sum(ok) < 2, continue; end
            cplx_pos(p) = mean(dv_p(ok) .* exp(1i * ph_p(ok)));
        end
        phi_pos   = angle(cplx_pos);
        valid_pos = ~isnan(phi_pos);
        if any(valid_pos)
            phase_map(ch,f) = angle(mean(exp(1i * phi_pos(valid_pos))));
        end
    end
end
coh_phase.RT = phase_map;

% --- Hit/Miss: ITC with inverted miss phases, per position ---
fprintf('  Hit/Miss...\n');
all_idx_hm   = find(a_ph.trialinfo(:,20) == 1 | a_ph.trialinfo(:,20) == 5);
hit_labels   = (a_ph.trialinfo(all_idx_hm, 20) == 1);
positions_hm = unique(a_ph.trialinfo(all_idx_hm, 16));
nPos_hm      = numel(positions_hm);
phase_map    = NaN(nCh, nFreq);
for ch = 1:nCh
    for f = 1:nFreq
        cplx_pos = NaN(nPos_hm, 1);
        for p = 1:nPos_hm
            mask = a_ph.trialinfo(all_idx_hm, 16) == positions_hm(p);
            if sum(mask) < 2, continue; end
            ph_p     = a_ph.phase_all(all_idx_hm(mask), f, ch);
            miss_p   = ~hit_labels(mask);
            ph_inv   = ph_p;
            ph_inv(miss_p) = mod(ph_p(miss_p) + pi, 2*pi) - pi;
            cplx_pos(p) = mean(exp(1i * ph_inv));
        end
        phi_pos   = angle(cplx_pos);
        valid_pos = ~isnan(phi_pos);
        if any(valid_pos)
            phase_map(ch,f) = angle(mean(exp(1i * phi_pos(valid_pos))));
        end
    end
end
coh_phase.hit_miss = phase_map;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. LOAD REGRESSION PREFERRED PHASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

reg_file  = fullfile(reg_root, 'multi_regression_channelwise_R2_abs_per_pos.mat');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 5. BUILD SIGNIFICANCE MASKS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coherence significance: per-channel coh_perm_pos (magnitudes).
% coh (magnitude) >= quantile(max(coh_perm_pos,[],2), 0.95).

coh_sig = struct();
for m = 1:length(coh_measures)
    measure = coh_measures{m};
    sig_map = false(nCh, nFreq);

    if strcmp(measure, 'hit_miss')
        val_folder  = fullfile(corr_root, 'hit_miss', 'all_loc_difflev');
        val_file    = 'itc.mat';   val_var  = 'itc';
        perm_file   = 'itc_perm_pos.mat'; perm_var = 'itc_perm_pos';
    elseif strcmp(measure, 'RT')
        val_folder  = fullfile(coh_root, 'RT', 'all_loc_difflev');
        val_file    = 'coherence.mat'; val_var  = 'coh';
        perm_file   = 'coh_perm_pos.mat'; perm_var = 'coh_perm_pos';
    else
        val_folder  = fullfile(coh_root, measure, 'all_loc_difflev');
        val_file    = 'coherence.mat'; val_var  = 'coh';
        perm_file   = 'coh_perm_pos.mat'; perm_var = 'coh_perm_pos';
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

% Regression significance (from per-channel thresholds in reg_results)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 6. CIRCULAR PASTEL COLORMAP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_cmap   = 256;
hue      = linspace(0, 1, n_cmap+1)'; hue = hue(1:end-1);
sat      = 0.45 * ones(n_cmap, 1);
val_cmap = 0.95 * ones(n_cmap, 1);
cmap_circ = hsv2rgb([hue, sat, val_cmap]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7. SHARED LABELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

row_labels = {'MUA', 'LFP', 'RT', 'Hit/Miss'};
row_keys   = {'mua', 'lfp', 'RT', 'hit_miss'};
nDVs       = length(row_keys);

col_labels    = {'Coherence (H2 per-pos)', 'Regression (\phi_{pref})'};
coh_subtitles = {'MUA (per-pos weighted phase)', 'LFP (per-pos weighted phase)', ...
    'RT (per-pos weighted phase)', 'Hit/Miss (per-pos ITC)'};
reg_subtitles = {'MUA (atan2(\beta_{sin},\beta_{cos}))', ...
    'LFP (atan2(\beta_{sin},\beta_{cos}))', ...
    'RT (atan2(\beta_{sin},\beta_{cos}))', ...
    'Hit/Miss (atan2(\beta_{sin},\beta_{cos}))'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 1 — Preferred Phase Heatmaps (channels x frequency)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = figure('Name', ['Preferred Phase (Abs-Per-Pos) - ' animal], ...
    'Units', 'centimeters', 'Position', [1 1 36 40]);
set(f1, 'PaperUnits', 'centimeters', 'PaperSize', [36 40], 'PaperPosition', [0 0 36 40]);

for row = 1:nDVs
    key = row_keys{row};

    % --- Column 1: Coherence ---
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

    % --- Column 2: Regression ---
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
sgtitle(sprintf(['Preferred Phase — Abs-Per-Pos (%s)\n' ...
    'Significance = COHERENCE MAGNITUDE (how big the per-position arrows are)'], animal), ...
    'FontSize', 13, 'FontWeight', 'bold');
print(f1, fullfile(save_root, 'preferred_phase_comparison.pdf'), '-dpdf');
fprintf('Figure 1 saved.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 2 — Polar Histograms at Key Frequencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nTopFreq = 3;
top_freq_idx = cell(nDVs, 1);
for row = 1:nDVs
    key = row_keys{row};
    n_sig_per_freq = sum(coh_sig.(key), 1);
    [~, sorted_idx] = sort(n_sig_per_freq, 'descend');
    top_freq_idx{row} = sorted_idx(1:min(nTopFreq, length(sorted_idx)));
end

f2 = figure('Name', ['Polar Histograms (Abs-Per-Pos) - ' animal], ...
    'Units', 'centimeters', 'Position', [1 1 52 38]);
set(f2, 'PaperUnits', 'centimeters', 'PaperSize', [52 38], 'PaperPosition', [0 0 52 38]);

nBins = 18;

for row = 1:nDVs
    key = row_keys{row};

    for fi = 1:nTopFreq
        fidx       = top_freq_idx{row}(fi);
        f_hz       = freq(fidx);
        sig_ch_coh = find(coh_sig.(key)(:,fidx));

        % Coherence polar histogram
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

        % Regression polar histogram
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

sgtitle(sprintf('Preferred Phase Distribution — Abs-Per-Pos (%s)\nBlack arrow = circular mean, n = sig channels', animal), ...
    'FontSize', 12, 'FontWeight', 'bold');
print(f2, fullfile(save_root, 'preferred_phase_polar.pdf'), '-dpdf');
fprintf('Figure 2 saved.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 3 — Pairwise Phase Consistency Across DVs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

f3 = figure('Name', ['Phase Consistency (Abs-Per-Pos) - ' animal], ...
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

sgtitle(sprintf('Pairwise Phase Consistency — Abs-Per-Pos (%s)\nDots = p<0.05, only sig channels in both', animal), ...
    'FontSize', 12, 'FontWeight', 'bold');
print(f3, fullfile(save_root, 'phase_consistency_across_measures.pdf'), '-dpdf');
fprintf('Figure 3 saved.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8. MONKEY-AVERAGE: COLLECT DATA FROM BOTH ANIMALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preferred phase is computed from raw data per animal.
% Level 2 (per animal): circular mean across all channels.
% Level 3 (monkey avg): circular mean across animals.

animals_all  = {'hermes', 'klecks'};
nAnimals     = numel(animals_all);

monkey_save_root = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/sampling_compare_abs_per_pos/monkey_avg');
if ~exist(monkey_save_root, 'dir'), mkdir(monkey_save_root); end

reg_dvs_all    = {'MUA_ERP_ampl_all', 'LFP_ERP_ampl_all', 'RT', 'hit_miss'};
reg_labels_all = {'mua', 'lfp', 'RT', 'hit_miss'};

animal_coh_avg     = cell(nAnimals, 1);
animal_coh_avg_sig = cell(nAnimals, 1);
animal_reg_avg     = cell(nAnimals, 1);
animal_reg_avg_sig = cell(nAnimals, 1);
animal_rho_coh_all = cell(nAnimals, 1);
animal_rho_reg_all = cell(nAnimals, 1);
animal_coh_phase_all = cell(nAnimals, 1);
animal_coh_sig_all   = cell(nAnimals, 1);
animal_reg_phase_all = cell(nAnimals, 1);
animal_reg_sig_all   = cell(nAnimals, 1);
animal_hit_coh_phase  = cell(nAnimals, 1);
animal_miss_coh_phase = cell(nAnimals, 1);
animal_hit_reg_phase  = cell(nAnimals, 1);
animal_miss_reg_phase = cell(nAnimals, 1);

hm_dvs       = {'MUA_ERP_ampl_all', 'LFP_ERP_ampl_all', 'RT', 'hit_miss'};
hm_ti_fields = {'MUA_ERP_trialinfo', 'LFP_ERP_trialinfo', 'RT_trialinfo', 'trialinfo'};
hm_labels    = {'mua', 'lfp', 'RT', 'hit_miss'};
nDV          = length(hm_dvs);

for a = 1:nAnimals
    animalName = animals_all{a};
    fprintf('\n=== Loading preferred phase for %s ===\n', animalName);

    a_base      = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName]);
    a_coh_root  = fullfile(a_base, 'phase_coherence',   'abs_per_pos', 'cp10_till_100');
    a_corr_root = fullfile(a_base, 'phase_correlation', 'abs_per_pos', 'cp10_till_100');
    a_reg_root  = fullfile(a_base, 'multi_lin_reg',     'abs_per_pos', 'cp10_till_100');
    a_data_root = fullfile(a_base, 'multi_lin_reg', 'cp10_till_100');   % shared input data

    tmp_f   = load(fullfile(a_data_root, 'frequency.mat'));
    a_freq  = tmp_f.frequency;
    a_nFreq = numel(a_freq);
    a_nCh   = 64;

    % --- Compute coherence preferred phase from raw data ---
    tmp_d  = load(fullfile(a_data_root, 'ph_all_sess.mat'), 'ph_comb');
    a_phd  = tmp_d.ph_comb;

    a_coh_phase = struct();
    a_coh_sig   = struct();

    % MUA
    pos_mua = unique(a_phd.MUA_ERP_trialinfo(:,16));
    nPm = numel(pos_mua);
    ph_map = NaN(a_nCh, a_nFreq);
    for ch = 1:a_nCh
        for f = 1:a_nFreq
            cp = NaN(nPm,1);
            for p = 1:nPm
                mk = a_phd.MUA_ERP_trialinfo(:,16)==pos_mua(p);
                if sum(mk)<2, continue; end
                ph_p=a_phd.phase_all(mk,f,ch); dv_p=a_phd.MUA_ERP_ampl_all(mk,ch);
                ok=~isnan(ph_p)&~isnan(dv_p); if sum(ok)<2, continue; end
                cp(p)=mean(dv_p(ok).*exp(1i*ph_p(ok)));
            end
            pp=angle(cp); vp=~isnan(pp);
            if any(vp), ph_map(ch,f)=angle(mean(exp(1i*pp(vp)))); end
        end
    end
    a_coh_phase.mua = ph_map;

    % LFP
    pos_lfp = unique(a_phd.LFP_ERP_trialinfo(:,16));
    nPl = numel(pos_lfp);
    ph_map = NaN(a_nCh, a_nFreq);
    for ch = 1:a_nCh
        for f = 1:a_nFreq
            cp = NaN(nPl,1);
            for p = 1:nPl
                mk = a_phd.LFP_ERP_trialinfo(:,16)==pos_lfp(p);
                if sum(mk)<2, continue; end
                ph_p=a_phd.phase_all(mk,f,ch); dv_p=a_phd.LFP_ERP_ampl_all(mk,ch);
                ok=~isnan(ph_p)&~isnan(dv_p); if sum(ok)<2, continue; end
                cp(p)=mean(dv_p(ok).*exp(1i*ph_p(ok)));
            end
            pp=angle(cp); vp=~isnan(pp);
            if any(vp), ph_map(ch,f)=angle(mean(exp(1i*pp(vp)))); end
        end
    end
    a_coh_phase.lfp = ph_map;

    % RT
    hidx_rt  = find(a_phd.RT_trialinfo(:,20)==1);
    pos_rt   = unique(a_phd.RT_trialinfo(hidx_rt,16));
    nPr      = numel(pos_rt);
    ph_map   = NaN(a_nCh, a_nFreq);
    for ch = 1:a_nCh
        for f = 1:a_nFreq
            cp = NaN(nPr,1);
            for p = 1:nPr
                pm=a_phd.RT_trialinfo(hidx_rt,16)==pos_rt(p); pg=hidx_rt(pm);
                if sum(pm)<2, continue; end
                ph_p=a_phd.phase_all(pg,f,ch); dv_p=a_phd.RT(pg,ch);
                ok=~isnan(ph_p)&~isnan(dv_p); if sum(ok)<2, continue; end
                cp(p)=mean(dv_p(ok).*exp(1i*ph_p(ok)));
            end
            pp=angle(cp); vp=~isnan(pp);
            if any(vp), ph_map(ch,f)=angle(mean(exp(1i*pp(vp)))); end
        end
    end
    a_coh_phase.RT = ph_map;

    % Hit/Miss
    aidx_hm  = find(a_phd.trialinfo(:,20)==1|a_phd.trialinfo(:,20)==5);
    hlbl_hm  = (a_phd.trialinfo(aidx_hm,20)==1);
    pos_hm   = unique(a_phd.trialinfo(aidx_hm,16));
    nPhm     = numel(pos_hm);
    ph_map   = NaN(a_nCh, a_nFreq);
    for ch = 1:a_nCh
        for f = 1:a_nFreq
            cp = NaN(nPhm,1);
            for p = 1:nPhm
                mk=a_phd.trialinfo(aidx_hm,16)==pos_hm(p);
                if sum(mk)<2, continue; end
                ph_p=a_phd.phase_all(aidx_hm(mk),f,ch);
                miss_p=~hlbl_hm(mk); ph_inv=ph_p;
                ph_inv(miss_p)=mod(ph_p(miss_p)+pi,2*pi)-pi;
                cp(p)=mean(exp(1i*ph_inv));
            end
            pp=angle(cp); vp=~isnan(pp);
            if any(vp), ph_map(ch,f)=angle(mean(exp(1i*pp(vp)))); end
        end
    end
    a_coh_phase.hit_miss = ph_map;

    % --- Significance masks ---
    for m = 1:length(coh_measures)
        measure = coh_measures{m};
        sig_map = false(a_nCh, a_nFreq);
        if strcmp(measure,'hit_miss')
            vfolder = fullfile(a_corr_root,'hit_miss','all_loc_difflev');
            vfile='itc.mat'; vvar='itc'; pfile='itc_perm_pos.mat'; pvar='itc_perm_pos';
        elseif strcmp(measure,'RT')
            vfolder = fullfile(a_coh_root,'RT','all_loc_difflev');
            vfile='coherence.mat'; vvar='coh'; pfile='coh_perm_pos.mat'; pvar='coh_perm_pos';
        else
            vfolder = fullfile(a_coh_root,measure,'all_loc_difflev');
            vfile='coherence.mat'; vvar='coh'; pfile='coh_perm_pos.mat'; pvar='coh_perm_pos';
        end
        for ch = 1:a_nCh
            vf=fullfile(vfolder,num2str(ch),vfile); pf=fullfile(vfolder,num2str(ch),pfile);
            if ~exist(vf,'file')||~exist(pf,'file'), continue; end
            tv=load(vf,vvar); val=tv.(vvar);
            tp=load(pf,pvar); prm=tp.(pvar);
            if any(isnan(val))||any(isnan(prm(:))), continue; end
            sig_map(ch,:) = val >= quantile(max(prm,[],2),0.95);
        end
        a_coh_sig.(measure) = sig_map;
    end

    % --- Level 2: circular mean across channels (coherence) ---
    a_coh_avg     = struct();
    a_coh_avg_sig = struct();
    for m = 1:length(coh_measures)
        key   = coh_measures{m};
        ph    = a_coh_phase.(key);
        avg_ph = NaN(1, a_nFreq);
        sig_f  = false(1, a_nFreq);
        for f = 1:a_nFreq
            valid = ~isnan(ph(:,f));
            if any(valid), avg_ph(f) = angle(mean(exp(1i*ph(valid,f)))); end
        end
        % Channel-average significance: mean magnitude >= channel_avg threshold
        if strcmp(key,'hit_miss')
            ca_file = fullfile(a_corr_root,'hit_miss','all_loc_difflev','channel_avg_results_itc.mat');
            if isfile(ca_file)
                ca = load(ca_file);
                sig_f = ca.itc_chan_avg >= ca.thresh_chan_avg_itc;
            end
        elseif strcmp(key,'RT')
            ca_file = fullfile(a_coh_root,'RT','all_loc_difflev','channel_avg_results.mat');
            if isfile(ca_file)
                ca = load(ca_file);
                sig_f = ca.coh_chan_avg >= ca.thresh_chan_avg;
            end
        else
            ca_file = fullfile(a_coh_root,key,'all_loc_difflev','channel_avg_results.mat');
            if isfile(ca_file)
                ca = load(ca_file);
                sig_f = ca.coh_chan_avg >= ca.thresh_chan_avg;
            end
        end
        a_coh_avg.(key)     = avg_ph;
        a_coh_avg_sig.(key) = sig_f;
    end

    animal_coh_avg{a}     = a_coh_avg;
    animal_coh_avg_sig{a} = a_coh_avg_sig;

    % --- Load regression preferred phase + level 2 avg ---
    a_reg_file  = fullfile(a_reg_root, 'multi_regression_channelwise_R2_abs_per_pos.mat');
    a_has_reg   = exist(a_reg_file,'file');
    a_reg_phase = struct(); a_reg_sig = struct();

    if a_has_reg
        a_reg = load(a_reg_file,'reg_results'); a_rr = a_reg.reg_results;
        sd = dir(fullfile(a_base,[animalName '_*']));
        if ~isempty(sd)
            ff = fullfile(a_base,sd(1).name,'Phase_analysis','hit_miss','100iter_cut@cp_m10','1','freqpow.mat');
        else, ff = ''; end
        if ~isempty(ff) && exist(ff,'file'), tmp=load(ff); a_freqs_reg=tmp.freqpow.freq;
        else, a_freqs_reg=a_freq; end
        a_nFreq_r = length(a_freqs_reg);

        for m = 1:length(reg_dvs_all)
            dv=reg_dvs_all{m}; lbl=reg_labels_all{m};
            if isfield(a_rr,dv) && isfield(a_rr.(dv),'phi_pref')
                a_reg_phase.(lbl) = a_rr.(dv).phi_pref;
            else
                a_reg_phase.(lbl) = NaN(a_nCh,a_nFreq_r);
            end
            sp = false(a_nCh,a_nFreq_r);
            if isfield(a_rr,dv)
                nr=size(a_rr.(dv).R2_phase,1);
                for ch=1:min(nr,length(a_rr.(dv).thresholds))
                    if isfield(a_rr.(dv).thresholds(ch),'thresh_phase')
                        thr=a_rr.(dv).thresholds(ch).thresh_phase;
                        if ~isempty(thr), sp(ch,:)=a_rr.(dv).R2_phase(ch,:)>thr(1); end
                    end
                end
            end
            a_reg_sig.(lbl) = sp;
        end
    else
        a_freqs_reg=a_freq; a_nFreq_r=a_nFreq;
        for m=1:length(reg_labels_all)
            a_reg_phase.(reg_labels_all{m})=NaN(a_nCh,a_nFreq);
            a_reg_sig.(reg_labels_all{m})=false(a_nCh,a_nFreq);
        end
    end

    a_reg_avg=struct(); a_reg_avg_sig=struct();
    for m=1:length(reg_labels_all)
        key=reg_labels_all{m}; dv=reg_dvs_all{m};
        ph=a_reg_phase.(key); nF=size(ph,2);
        avg_ph=NaN(1,nF); sig_f=false(1,nF);
        for f=1:nF
            v=~isnan(ph(:,f)); if any(v), avg_ph(f)=angle(mean(exp(1i*ph(v,f)))); end
        end
        if a_has_reg && isfield(a_rr,dv) && isfield(a_rr.(dv),'channel_avg_R') && isfield(a_rr.(dv),'channel_avg_thresh')
            obs=a_rr.(dv).channel_avg_R.phase; thr=a_rr.(dv).channel_avg_thresh.phase;
            nFs=min(length(obs),nF); sig_f(1:nFs)=obs(1:nFs)>thr;
        end
        a_reg_avg.(key)=avg_ph; a_reg_avg_sig.(key)=sig_f;
    end
    animal_reg_avg{a}=a_reg_avg; animal_reg_avg_sig{a}=a_reg_avg_sig;
    animal_coh_phase_all{a}=a_coh_phase; animal_coh_sig_all{a}=a_coh_sig;
    animal_reg_phase_all{a}=a_reg_phase; animal_reg_sig_all{a}=a_reg_sig;

    % --- Pairwise phase consistency (coherence) ---
    a_rho_coh = NaN(nPairs,a_nFreq);
    for p=1:nPairs
        k1=measure_pairs{p,1}; k2=measure_pairs{p,2};
        ph1=a_coh_phase.(k1); ph2=a_coh_phase.(k2);
        sg1=a_coh_sig.(k1);   sg2=a_coh_sig.(k2);
        for fi=1:a_nFreq
            vld=sg1(:,fi)&sg2(:,fi)&~isnan(ph1(:,fi))&~isnan(ph2(:,fi));
            if sum(vld)<4, continue; end
            a_rho_coh(p,fi)=circ_corrcc(ph1(vld,fi),ph2(vld,fi));
        end
    end
    animal_rho_coh_all{a}=a_rho_coh;

    % --- Hit-only and miss-only preferred phase ---
    fprintf('  Computing hit/miss preferred phase for %s...\n', animalName);

    hit_coh=struct(); miss_coh=struct();
    hit_reg=struct(); miss_reg=struct();

    for m=1:nDV
        dv_name=hm_dvs{m}; ti_field=hm_ti_fields{m}; lbl=hm_labels{m};
        hit_coh.(lbl)=NaN(a_nCh,a_nFreq); miss_coh.(lbl)=NaN(a_nCh,a_nFreq);
        hit_reg.(lbl)=NaN(a_nCh,a_nFreq); miss_reg.(lbl)=NaN(a_nCh,a_nFreq);

        if isfield(a_phd,ti_field), dv_ti=a_phd.(ti_field); else, dv_ti=a_phd.trialinfo; end
        nTp=size(a_phd.phase_all,1);

        if strcmp(lbl,'hit_miss')
            all_hm=find(a_phd.trialinfo(:,20)==1|a_phd.trialinfo(:,20)==5);
            pos_hm2=unique(a_phd.trialinfo(all_hm,16)); nPh2=numel(pos_hm2);
            hit_mask2=(a_phd.trialinfo(all_hm,20)==1);
            ph_hit_idx2=all_hm(hit_mask2); ph_miss_idx2=all_hm(~hit_mask2);
            for ch=1:a_nCh
                for f=1:a_nFreq
                    % hit-only ITC per position
                    cp_h=NaN(nPh2,1); cp_m=NaN(nPh2,1);
                    for p=1:nPh2
                        mk_h=a_phd.trialinfo(ph_hit_idx2,16)==pos_hm2(p);
                        if sum(mk_h)>=2, cp_h(p)=mean(exp(1i*a_phd.phase_all(ph_hit_idx2(mk_h),f,ch))); end
                        mk_m=a_phd.trialinfo(ph_miss_idx2,16)==pos_hm2(p);
                        if sum(mk_m)>=2, cp_m(p)=mean(exp(1i*a_phd.phase_all(ph_miss_idx2(mk_m),f,ch))); end
                    end
                    mch=mean(cp_h,'omitnan'); if ~isnan(mch), hit_coh.(lbl)(ch,f)=angle(mch); end
                    mcm=mean(cp_m,'omitnan'); if ~isnan(mcm), miss_coh.(lbl)(ch,f)=angle(mcm); end
                    hit_reg.(lbl)(ch,f)=hit_coh.(lbl)(ch,f);
                    miss_reg.(lbl)(ch,f)=miss_coh.(lbl)(ch,f);
                end
            end
        else
            if strcmp(lbl,'RT')
                hit_idx2=find(dv_ti(:,20)==1); pos2=unique(dv_ti(hit_idx2,16));
                hit_phase_idx=hit_idx2; miss_phase_idx=[];
                dv_mat=a_phd.RT;
            else
                hit_idx2=find(dv_ti(:,20)==1); miss_idx2=find(dv_ti(:,20)==5);
                pos2=unique(dv_ti(:,16));
                hit_phase_idx=hit_idx2; miss_phase_idx=miss_idx2;
                dv_mat=a_phd.(dv_name);
            end
            nP2=numel(pos2);

            for ch=1:a_nCh
                if size(dv_mat,2)<ch, continue; end
                dv_ch=dv_mat(:,ch);
                for f=1:a_nFreq
                    % hit
                    if ~isempty(hit_phase_idx)
                        cp_h=NaN(nP2,1);
                        for p=1:nP2
                            mk=dv_ti(hit_phase_idx,16)==pos2(p); pg=hit_phase_idx(mk);
                            if sum(mk)<2, continue; end
                            ph_p=a_phd.phase_all(pg,f,ch); dv_p=dv_ch(pg);
                            ok=~isnan(ph_p)&~isnan(dv_p); if sum(ok)<2, continue; end
                            cp_h(p)=mean(dv_p(ok).*exp(1i*ph_p(ok)));
                        end
                        mc=mean(cp_h,'omitnan'); if ~isnan(mc), hit_coh.(lbl)(ch,f)=angle(mc); end
                    end
                    % miss
                    if ~isempty(miss_phase_idx) && ~strcmp(lbl,'RT')
                        cp_m=NaN(nP2,1);
                        for p=1:nP2
                            mk=dv_ti(miss_phase_idx,16)==pos2(p); pg=miss_phase_idx(mk);
                            if sum(mk)<2, continue; end
                            ph_p=a_phd.phase_all(pg,f,ch); dv_p=dv_ch(pg);
                            ok=~isnan(ph_p)&~isnan(dv_p); if sum(ok)<2, continue; end
                            cp_m(p)=mean(dv_p(ok).*exp(1i*ph_p(ok)));
                        end
                        mc=mean(cp_m,'omitnan'); if ~isnan(mc), miss_coh.(lbl)(ch,f)=angle(mc); end
                    end
                    % regression phi_pref (from precomputed results)
                    if a_has_reg && isfield(a_reg_phase,lbl)
                        [~,ridx]=min(abs(a_freqs_reg-a_freq(f)));
                        hit_reg.(lbl)(ch,f)=a_reg_phase.(lbl)(ch,ridx);
                        miss_reg.(lbl)(ch,f)=a_reg_phase.(lbl)(ch,ridx);
                    end
                end
            end
        end
    end

    animal_hit_coh_phase{a}=hit_coh; animal_miss_coh_phase{a}=miss_coh;
    animal_hit_reg_phase{a}=hit_reg;  animal_miss_reg_phase{a}=miss_reg;

    %% Per-animal hit/miss figures
    a_save = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/sampling_compare_abs_per_pos', animalName);
    if ~exist(a_save,'dir'), mkdir(a_save); end

    coh_st_hm = {'MUA (per-pos weighted)', 'LFP (per-pos weighted)', 'RT (per-pos weighted)', 'Hit/Miss (per-pos ITC)'};
    reg_st_hm = {'MUA (\phi_{pref})', 'LFP (\phi_{pref})', 'RT (\phi_{pref})', 'Hit/Miss (\phi_{pref})'};

    for trial_type = {'hit','miss'}
        tt = trial_type{1};
        if strcmp(tt,'hit'), ph_coh=hit_coh; ph_reg=hit_reg;
        else, ph_coh=miss_coh; ph_reg=miss_reg; end

        fig_tt = figure('Name', sprintf('%s-Only Phase (Abs-Per-Pos) - %s', tt, animalName), ...
            'Units','centimeters','Position',[1 1 36 40]);
        set(fig_tt,'PaperUnits','centimeters','PaperSize',[36 40],'PaperPosition',[0 0 36 40]);

        for row=1:nDV
            lbl=hm_labels{row};

            subplot(nDV,2,(row-1)*2+1);
            data=ph_coh.(lbl);
            h_img=imagesc(a_freq,1:a_nCh,data);
            set(gca,'YDir','normal','Color',[1 1 1]);
            colormap(gca,cmap_circ); caxis([-pi pi]);
            set(h_img,'AlphaData',~isnan(data));
            xlabel('Frequency (Hz)'); ylabel('Channel');
            title(coh_st_hm{row},'FontSize',9); set(gca,'FontSize',8,'Box','on');
            if row==1
                text(0.5,1.22,'Coherence (per-pos)','Units','normalized','HorizontalAlignment','center','FontSize',13,'FontWeight','bold');
            end

            subplot(nDV,2,(row-1)*2+2);
            data=ph_reg.(lbl);
            h_img=imagesc(a_freq,1:a_nCh,data);
            set(gca,'YDir','normal','Color',[1 1 1]);
            colormap(gca,cmap_circ); caxis([-pi pi]);
            set(h_img,'AlphaData',~isnan(data));
            xlabel('Frequency (Hz)'); ylabel('Channel');
            title(reg_st_hm{row},'FontSize',9); set(gca,'FontSize',8,'Box','on');
            if row==1
                text(0.5,1.22,'Regression (\phi_{pref})','Units','normalized','HorizontalAlignment','center','FontSize',13,'FontWeight','bold');
            end
        end

        cb=colorbar('Location','southoutside');
        cb.Ticks=[-pi -pi/2 0 pi/2 pi];
        cb.TickLabels={'-\pi','-\pi/2','0','\pi/2','\pi'};
        cb.Position=[0.25 0.02 0.5 0.015]; cb.Label.String='Preferred Phase (rad)';
        sgtitle(sprintf('%s Trials Only — Abs-Per-Pos (%s)', [upper(tt(1)) tt(2:end)], animalName), 'FontSize',14,'FontWeight','bold');
        print(fig_tt, fullfile(a_save, sprintf('preferred_phase_%s_only.pdf',tt)), '-dpdf');
        fprintf('%s-only figure for %s saved.\n', tt, animalName);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 4 — Monkey-Average Preferred Phase Heatmap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coherence column: circular mean across channels per animal (level 2),
%   then circular mean across animals (level 3).
% Regression column: same approach using channel-average phi_pref.

f4 = figure('Name','Monkey-Avg Preferred Phase (Abs-Per-Pos)', ...
    'Units','centimeters','Position',[1 1 36 40]);
set(f4,'PaperUnits','centimeters','PaperSize',[36 40],'PaperPosition',[0 0 36 40]);

for row=1:nDVs
    key=row_keys{row};

    % --- Coherence column ---
    subplot(nDVs,2,(row-1)*2+1);
    phase_stack=NaN(3,nFreq); alpha_stack=ones(3,nFreq)*nonsig_alpha;
    for a=1:nAnimals
        ph_a=animal_coh_avg{a}.(key); nFa=min(length(ph_a),nFreq);
        phase_stack(a,1:nFa)=ph_a(1:nFa);
        sig_a=animal_coh_avg_sig{a}.(key);
        alpha_stack(a,find(sig_a(1:min(end,nFreq))))=1;
    end
    % monkey avg: circular mean across animals
    phase_stack(3,:)=angle(mean(exp(1i*phase_stack(1:nAnimals,:)),1,'omitnan'));
    alpha_stack(isnan(phase_stack))=0;

    h_img=imagesc(freq,1:3,phase_stack);
    set(gca,'YDir','normal','Color',[1 1 1]);
    colormap(gca,cmap_circ); caxis([-pi pi]);
    set(h_img,'AlphaData',alpha_stack);
    yticks(1:3); yticklabels([animals_all,{'Monkey avg'}]);
    xlabel('Frequency (Hz)');
    title(sprintf('%s — Coherence',row_labels{row}),'FontSize',9);
    set(gca,'FontSize',8,'Box','on');
    if row==1
        text(0.5,1.22,'Coherence (per-pos avg)','Units','normalized','HorizontalAlignment','center','FontSize',13,'FontWeight','bold');
    end

    % --- Regression column ---
    subplot(nDVs,2,(row-1)*2+2);
    phase_stack_r=NaN(3,length(freqs_reg)); alpha_stack_r=ones(3,length(freqs_reg))*nonsig_alpha;
    for a=1:nAnimals
        vals=animal_reg_avg{a}.(key); nF=min(length(vals),length(freqs_reg));
        phase_stack_r(a,1:nF)=vals(1:nF);
        sv=animal_reg_avg_sig{a}.(key);
        alpha_stack_r(a,find(sv(1:min(end,length(freqs_reg)))))=1;
    end
    phase_stack_r(3,:)=angle(mean(exp(1i*phase_stack_r(1:nAnimals,:)),1,'omitnan'));
    alpha_stack_r(isnan(phase_stack_r))=0;
    % monkey-avg regression significance from precomputed file
    dv_idx=find(strcmp(reg_labels_all,key),1);
    if ~isempty(dv_idx)
        mk_file=fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/results_combined', ...
            'multi_lin_reg','abs_per_pos','cp10_till_100',reg_dvs_all{dv_idx},'monkey_avg_results.mat');
        if isfile(mk_file)
            mk_r=load(mk_file,'monkey_avg_obs','thresh_monkey');
            nFr2=min(length(mk_r.monkey_avg_obs.phase),length(freqs_reg));
            mk_sig=mk_r.monkey_avg_obs.phase(1:nFr2)>mk_r.thresh_monkey.phase;
            alpha_stack_r(3,find(mk_sig))=1;
        end
    end

    h_img2=imagesc(freqs_reg,1:3,phase_stack_r);
    set(gca,'YDir','normal','Color',[1 1 1]);
    colormap(gca,cmap_circ); caxis([-pi pi]);
    set(h_img2,'AlphaData',alpha_stack_r);
    yticks(1:3); yticklabels([animals_all,{'Monkey avg'}]);
    xlabel('Frequency (Hz)');
    title(sprintf('%s — Regression',row_labels{row}),'FontSize',9);
    set(gca,'FontSize',8,'Box','on');
    if row==1
        text(0.5,1.22,'Regression (\phi_{pref})','Units','normalized','HorizontalAlignment','center','FontSize',13,'FontWeight','bold');
    end
end

cb=colorbar('Location','southoutside');
cb.Ticks=[-pi -pi/2 0 pi/2 pi]; cb.TickLabels={'-\pi','-\pi/2','0','\pi/2','\pi'};
cb.Position=[0.25 0.02 0.5 0.015]; cb.Label.String='Preferred Phase (rad)';
sgtitle({'Preferred Phase — Abs-Per-Pos', ...
    'Significance = COHERENCE MAGNITUDE (how big the per-position arrows are)'}, ...
    'FontSize',13,'FontWeight','bold');
print(f4,fullfile(monkey_save_root,'monkey_avg_preferred_phase.pdf'),'-dpdf');
fprintf('Figure 4 saved.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 5 — Monkey-Average Hit-Only Preferred Phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

coh_st_hm={'MUA (per-pos weighted)','LFP (per-pos weighted)','RT (per-pos weighted)','Hit/Miss (per-pos ITC)'};
reg_st_hm={'MUA (\phi_{pref})','LFP (\phi_{pref})','RT (\phi_{pref})','Hit/Miss (\phi_{pref})'};

f_hit_avg = figure('Name','Monkey-Avg Hit-Only Phase (Abs-Per-Pos)', ...
    'Units','centimeters','Position',[1 1 36 40]);
set(f_hit_avg,'PaperUnits','centimeters','PaperSize',[36 40],'PaperPosition',[0 0 36 40]);

for row=1:nDV
    lbl=hm_labels{row};

    subplot(nDV,2,(row-1)*2+1);
    phase_stack=NaN(3,nFreq); cplx_stack=complex(NaN(nAnimals,nFreq));
    for a=1:nAnimals
        ph_a=animal_hit_coh_phase{a}.(lbl);
        avg_cplx=complex(NaN(1,size(ph_a,2)));
        for f=1:size(ph_a,2)
            v=~isnan(ph_a(:,f));
            if any(v), avg_cplx(f)=mean(exp(1i*ph_a(v,f))); end
        end
        nFa=min(length(avg_cplx),nFreq);
        phase_stack(a,1:nFa)=angle(avg_cplx(1:nFa));
        cplx_stack(a,1:nFa)=avg_cplx(1:nFa);
    end
    phase_stack(3,:)=angle(mean(cplx_stack,1,'omitnan'));

    h_img=imagesc(freq,1:3,phase_stack);
    set(gca,'YDir','normal','Color',[1 1 1]);
    colormap(gca,cmap_circ); caxis([-pi pi]);
    set(h_img,'AlphaData',~isnan(phase_stack));
    yticks(1:3); yticklabels([animals_all,{'Monkey avg'}]);
    xlabel('Frequency (Hz)'); title(coh_st_hm{row},'FontSize',9);
    set(gca,'FontSize',8,'Box','on');
    if row==1
        text(0.5,1.22,'Coherence (per-pos)','Units','normalized','HorizontalAlignment','center','FontSize',13,'FontWeight','bold');
    end

    subplot(nDV,2,(row-1)*2+2);
    phase_stack=NaN(3,nFreq);
    for a=1:nAnimals
        ph=animal_hit_reg_phase{a}.(lbl);
        avg=NaN(1,size(ph,2));
        for f=1:size(ph,2), v=~isnan(ph(:,f)); if any(v), avg(f)=angle(mean(exp(1i*ph(v,f)))); end; end
        nFa=min(length(avg),nFreq); phase_stack(a,1:nFa)=avg(1:nFa);
    end
    phase_stack(3,:)=angle(mean(exp(1i*phase_stack(1:nAnimals,:)),1,'omitnan'));

    h_img=imagesc(freq,1:3,phase_stack);
    set(gca,'YDir','normal','Color',[1 1 1]);
    colormap(gca,cmap_circ); caxis([-pi pi]);
    set(h_img,'AlphaData',~isnan(phase_stack));
    yticks(1:3); yticklabels([animals_all,{'Monkey avg'}]);
    xlabel('Frequency (Hz)'); title(reg_st_hm{row},'FontSize',9);
    set(gca,'FontSize',8,'Box','on');
    if row==1
        text(0.5,1.22,'Regression (\phi_{pref})','Units','normalized','HorizontalAlignment','center','FontSize',13,'FontWeight','bold');
    end
end

cb=colorbar('Location','southoutside');
cb.Ticks=[-pi -pi/2 0 pi/2 pi]; cb.TickLabels={'-\pi','-\pi/2','0','\pi/2','\pi'};
cb.Position=[0.25 0.02 0.5 0.015]; cb.Label.String='Preferred Phase (rad)';
sgtitle('Hit Trials Only — Abs-Per-Pos','FontSize',14,'FontWeight','bold');
print(f_hit_avg,fullfile(monkey_save_root,'monkey_avg_preferred_phase_hit_only.pdf'),'-dpdf');
fprintf('Figure 5 (monkey-avg hit-only) saved.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 6 — Monkey-Average Miss-Only Preferred Phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f_miss_avg = figure('Name','Monkey-Avg Miss-Only Phase (Abs-Per-Pos)', ...
    'Units','centimeters','Position',[1 1 36 40]);
set(f_miss_avg,'PaperUnits','centimeters','PaperSize',[36 40],'PaperPosition',[0 0 36 40]);

for row=1:nDV
    lbl=hm_labels{row};

    subplot(nDV,2,(row-1)*2+1);
    phase_stack=NaN(3,nFreq); cplx_stack=complex(NaN(nAnimals,nFreq));
    for a=1:nAnimals
        ph_a=animal_miss_coh_phase{a}.(lbl);
        avg_cplx=complex(NaN(1,size(ph_a,2)));
        for f=1:size(ph_a,2)
            v=~isnan(ph_a(:,f));
            if any(v), avg_cplx(f)=mean(exp(1i*ph_a(v,f))); end
        end
        nFa=min(length(avg_cplx),nFreq);
        phase_stack(a,1:nFa)=angle(avg_cplx(1:nFa));
        cplx_stack(a,1:nFa)=avg_cplx(1:nFa);
    end
    phase_stack(3,:)=angle(mean(cplx_stack,1,'omitnan'));

    h_img=imagesc(freq,1:3,phase_stack);
    set(gca,'YDir','normal','Color',[1 1 1]);
    colormap(gca,cmap_circ); caxis([-pi pi]);
    set(h_img,'AlphaData',~isnan(phase_stack));
    yticks(1:3); yticklabels([animals_all,{'Monkey avg'}]);
    xlabel('Frequency (Hz)'); title(coh_st_hm{row},'FontSize',9);
    set(gca,'FontSize',8,'Box','on');
    if row==1
        text(0.5,1.22,'Coherence (per-pos)','Units','normalized','HorizontalAlignment','center','FontSize',13,'FontWeight','bold');
    end

    subplot(nDV,2,(row-1)*2+2);
    phase_stack=NaN(3,nFreq);
    for a=1:nAnimals
        ph=animal_miss_reg_phase{a}.(lbl);
        avg=NaN(1,size(ph,2));
        for f=1:size(ph,2), v=~isnan(ph(:,f)); if any(v), avg(f)=angle(mean(exp(1i*ph(v,f)))); end; end
        nFa=min(length(avg),nFreq); phase_stack(a,1:nFa)=avg(1:nFa);
    end
    phase_stack(3,:)=angle(mean(exp(1i*phase_stack(1:nAnimals,:)),1,'omitnan'));

    h_img=imagesc(freq,1:3,phase_stack);
    set(gca,'YDir','normal','Color',[1 1 1]);
    colormap(gca,cmap_circ); caxis([-pi pi]);
    set(h_img,'AlphaData',~isnan(phase_stack));
    yticks(1:3); yticklabels([animals_all,{'Monkey avg'}]);
    xlabel('Frequency (Hz)'); title(reg_st_hm{row},'FontSize',9);
    set(gca,'FontSize',8,'Box','on');
    if row==1
        text(0.5,1.22,'Regression (\phi_{pref})','Units','normalized','HorizontalAlignment','center','FontSize',13,'FontWeight','bold');
    end
end

cb=colorbar('Location','southoutside');
cb.Ticks=[-pi -pi/2 0 pi/2 pi]; cb.TickLabels={'-\pi','-\pi/2','0','\pi/2','\pi'};
cb.Position=[0.25 0.02 0.5 0.015]; cb.Label.String='Preferred Phase (rad)';
sgtitle('Miss Trials Only — Abs-Per-Pos','FontSize',14,'FontWeight','bold');
print(f_miss_avg,fullfile(monkey_save_root,'monkey_avg_preferred_phase_miss_only.pdf'),'-dpdf');
fprintf('Figure 6 (monkey-avg miss-only) saved.\n');
