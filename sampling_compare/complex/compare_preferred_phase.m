%% Compare Preferred Phase: Coherence vs Multiple Linear Regression
%
% This script compares the preferred oscillatory phase estimated by two
% independent methods — phase coherence and multiple linear regression —
% for four dependent variables (MUA amplitude, LFP amplitude, RT, and
% hit/miss outcome).
%
% METHOD 1 — Coherence: the preferred phase is extracted from precomputed
%   coherence results (phase_spec), which reflects the DV-amplitude-weighted
%   mean phase direction at each frequency.
%
% METHOD 2 — Regression: the preferred phase comes from the multiple linear
%   regression model DV ~ cos(phase) + sin(phase), yielding
%   phi_pref = atan2(beta_sin, beta_cos).
%
% -------------------------------------------------------------------------
% SIGNIFICANCE APPROACH — COHERENCE vs REGRESSION
% -------------------------------------------------------------------------
%
%   Per-channel level:
%     Coherence:   |coh_complex| >= 95th pctile of max(|coh_perm_complex|)
%                  across frequencies (max-stat corrected); from
%                  coh_perm_complex.mat per channel.
%     Regression:  R2_phase(ch,f) > thresh_phase(ch), where thresh_phase is
%                  the 95th pctile of max(null_R2_phase) across frequencies;
%                  stored in reg_results.(dv).thresholds(ch).thresh_phase.
%
%   Channel-average level (per-animal row in Fig 4):
%     Coherence:   coh_chan_avg >= thresh_chan_avg; both from
%                  channel_avg_results.mat, computed by Coh_lfp_mua.m as
%                  mean(|coh_complex|) vs permuted equivalent.
%     Regression:  channel_avg_R.phase(f) > channel_avg_thresh.phase;
%                  both from reg_results.(dv) in
%                  multi_regression_channelwise_R2.mat, computed by
%                  regress_stats_R2.m as mean(R2_phase) across channels vs
%                  mean(null_R2_phase) across channels (max-stat corrected).
%
%   Monkey-average level (monkey-avg row in Fig 4):
%     Coherence:   coh_monkey_avg >= thresh_monkey_avg; from
%                  monkey_avg_results.mat, computed by Coh_lfp_mua.m.
%     Regression:  monkey_avg_obs.phase(f) > thresh_monkey.phase; from
%                  monkey_avg_results.mat in results_combined, computed by
%                  regress_stats_R2.m as mean(channel_avg_R2) across animals
%                  vs mean(channel_avg_null) across animals (max-stat).
%
%   KEY DIFFERENCE: Coherence significance tests phase consistency (magnitude
%   of complex mean across channels/animals). Regression significance tests
%   whether phase explains DV variance (R² increment). These are related but
%   distinct: coherence requires consistent preferred direction; regression
%   only requires a significant phase–DV relationship per channel.
%
%   Hit-only / miss-only figures: no significance overlay (no hit/miss-
%   specific permutation tests have been run for regression).
%
% -------------------------------------------------------------------------
%
% The script produces the following outputs:
%
%   SINGLE-ANIMAL FIGURES (for the selected animal):
%     Fig 1 — Phase heatmaps: channels x frequency, coherence vs regression
%             (significant channels at full opacity, non-significant faded)
%     Fig 2 — Polar histograms at the 3 frequencies with most significant
%             channels, showing the distribution of preferred phases
%     Fig 3 — Pairwise phase consistency across DVs (circular correlation
%             of preferred phases between every pair of DVs, computed only
%             over channels significant in both)
%
%   MONKEY-AVERAGE FIGURES (both animals combined):
%     Fig 4 — Monkey-average preferred phase heatmap (per-animal circular
%             mean across all channels + cross-animal average); significance
%             as described above.
%
%   HIT-ONLY / MISS-ONLY FIGURES:
%     Per-animal and monkey-average heatmaps of preferred phase computed
%     separately from hit trials and miss trials. For each trial subset,
%     the coherence estimate uses DV-weighted mean phase and the regression
%     estimate fits DV ~ cos(phase) + sin(phase).
%     No significance overlay (hit/miss-specific permutation tests pending).

clearvars; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

animal = 'hermes';   % 'hermes' or 'klecks'
nonsig_alpha = 0.2;  % transparency for non-significant bins (0 = invisible, 1 = full color)

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

coh_root  = fullfile(base_results, 'phase_coherence',   'complex', 'cp10_till_100');
reg_root  = fullfile(base_results, 'multi_lin_reg',     'complex', 'cp10_till_100');
corr_root = fullfile(base_results, 'phase_correlation', 'complex', 'cp10_till_100');
data_root = fullfile(base_results, 'multi_lin_reg', 'cp10_till_100');   % shared input data (frequency.mat, ph_all_sess.mat)

save_root = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/sampling_compare', animal);
if ~exist(save_root, 'dir'), mkdir(save_root); end

% Frequency axis (shared across coherence channels)
load(fullfile(data_root, 'frequency.mat'));
freq  = frequency;
nFreq = numel(freq);
nCh   = 64;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. LOAD COHERENCE PREFERRED PHASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each DV and channel, load coh_complex from the precomputed coherence.mat
% file and extract the preferred phase as angle(coh_complex). coh_complex is the
% amplitude-weighted complex mean across trials (1 x nFreq); its angle gives the
% DV-weighted mean phase direction. Hit/miss uses angle(itc_complex) from itc.mat.

coh_phase    = struct();
coh_measures = {'mua', 'lfp', 'RT', 'hit_miss'};

for m = 1:length(coh_measures)
    measure   = coh_measures{m};
    phase_map = NaN(nCh, nFreq);

    for ch = 1:nCh
        ch_folder = fullfile(coh_root, measure, 'all_loc_difflev', num2str(ch));
        coh_file  = fullfile(ch_folder, 'coherence.mat');
        if ~exist(coh_file, 'file'), continue; end

        tmp = load(coh_file, 'coh_complex');
        if isfield(tmp, 'coh_complex') && ~any(isnan(tmp.coh_complex))
            phase_map(ch,:) = angle(tmp.coh_complex);
        end
    end

    coh_phase.(measure) = phase_map;
end

% --- Hit/miss: preferred phase = angle(itc_complex), loaded directly ---
hm_phase_map = NaN(nCh, nFreq);
for ch = 1:nCh
    itc_file = fullfile(corr_root, 'hit_miss', 'all_loc_difflev', num2str(ch), 'itc.mat');
    if ~exist(itc_file, 'file'), continue; end
    tmp = load(itc_file, 'itc_complex');
    if isfield(tmp, 'itc_complex') && ~any(isnan(tmp.itc_complex))
        hm_phase_map(ch,:) = angle(tmp.itc_complex);
    end
end
coh_phase.hit_miss = hm_phase_map;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. LOAD REGRESSION PREFERRED PHASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each DV, load phi_pref from the channelwise multiple regression
% results. phi_pref = atan2(beta_sin, beta_cos) from the full model..

reg_file  = fullfile(reg_root, 'multi_regression_channelwise_R2.mat');
has_reg   = exist(reg_file, 'file');
reg_phase = struct();

reg_dvs    = {'MUA_ERP_ampl_all', 'LFP_ERP_ampl_all', 'RT', 'hit_miss'};
reg_labels = {'mua', 'lfp', 'RT', 'hit_miss'};

if has_reg
    load(reg_file, 'reg_results');

    % The regression may use a slightly different frequency axis; retrieve it
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
% Significance is determined by permutation testing. A channel x frequency
% bin is significant if the observed value exceeds the 95th percentile of
% the distribution of per-iteration maxima (max-corrected across freq).

% --- Coherence significance (from coh_perm_complex.mat) ---
coh_sig = struct();
for m = 1:length(coh_measures)
    measure = coh_measures{m};
    sig_map = false(nCh, nFreq);

    for ch = 1:nCh
        ch_folder = fullfile(coh_root, measure, 'all_loc_difflev', num2str(ch));
        coh_file  = fullfile(ch_folder, 'coherence.mat');
        perm_file = fullfile(ch_folder, 'coh_perm_complex.mat');
        if ~exist(coh_file, 'file') || ~exist(perm_file, 'file'), continue; end

        tmp = load(coh_file, 'coh_complex');          val = abs(tmp.coh_complex);
        tmp = load(perm_file, 'coh_perm_complex');    prm = abs(tmp.coh_perm_complex);
        if any(isnan(val)) || any(isnan(prm(:))), continue; end

        tmax = max(prm, [], 2);           % max across frequencies per iteration
        thr  = quantile(tmax, 0.95);      % 95th percentile of max distribution
        sig_map(ch,:) = val >= thr;
    end

    coh_sig.(measure) = sig_map;
end

% --- Hit/miss significance from ITC permutation (itc.mat / itc_perm.mat) ---
% Consistent with mua/lfp/RT: load complex values, take abs() inline.
% observed |itc_complex| >= 95th percentile of max(|itc_perm_complex|) across freq.
hm_sig_map = false(nCh, nFreq);
for ch = 1:nCh
    ch_folder = fullfile(corr_root, 'hit_miss', 'all_loc_difflev', num2str(ch));
    itc_file  = fullfile(ch_folder, 'itc.mat');
    perm_file = fullfile(ch_folder, 'itc_perm.mat');
    if ~exist(itc_file, 'file') || ~exist(perm_file, 'file'), continue; end

    tmp = load(itc_file,  'itc_complex');       val = abs(tmp.itc_complex);
    tmp = load(perm_file, 'itc_perm_complex');  prm = abs(tmp.itc_perm_complex);
    if any(isnan(val)) || any(isnan(prm(:))), continue; end

    tmax = max(prm, [], 2);
    thr  = quantile(tmax, 0.95);
    hm_sig_map(ch,:) = val >= thr;
end
coh_sig.hit_miss = hm_sig_map;

% --- Regression significance (from per-channel thresholds in reg_results) ---
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
% Phase is circular (-pi to pi), so the colormap must wrap

n_cmap   = 256;
hue      = linspace(0, 1, n_cmap+1)'; hue = hue(1:end-1);
sat      = 0.45 * ones(n_cmap, 1);
val      = 0.95 * ones(n_cmap, 1);
cmap_circ = hsv2rgb([hue, sat, val]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 7. SHARED LABELS FOR ALL FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

row_labels = {'MUA', 'LFP', 'RT', 'Hit/Miss'};
row_keys   = {'mua', 'lfp', 'RT', 'hit_miss'};
nDVs       = length(row_keys);

col_labels = {'Coherence (phase\_spec)', 'Regression (\phi_{pref})'};

coh_subtitles = {'MUA (amp-weighted phase)', 'LFP (amp-weighted phase)', ...
    'RT (amp-weighted phase)', 'Hit/Miss (amp-weighted phase)'};
reg_subtitles = {'MUA (atan2(\beta_{sin},\beta_{cos}))', ...
    'LFP (atan2(\beta_{sin},\beta_{cos}))', ...
    'RT (atan2(\beta_{sin},\beta_{cos}))', ...
    'Hit/Miss (atan2(\beta_{sin},\beta_{cos}))'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 1 — Preferred Phase Heatmaps (channels x frequency)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Each row is a DV, each column is a method. Colour encodes the preferred
% phase angle. Significant channel-frequency bins are shown at full opacity;
% non-significant bins are faded.

f1 = figure('Name', ['Preferred Phase - ' animal], ...
    'Units', 'centimeters', 'Position', [1 1 36 40]);
set(f1, 'PaperUnits', 'centimeters', 'PaperSize', [36 40], 'PaperPosition', [0 0 36 40]);

for row = 1:nDVs
    key = row_keys{row};

    % --- Column 1: Coherence ---
    subplot(nDVs, 2, (row-1)*2 + 1);
    data_coh  = coh_phase.(key);
    h_img     = imagesc(freq, 1:nCh, data_coh);
    set(gca, 'YDir', 'normal', 'Color', [1 1 1]);
    colormap(gca, cmap_circ);  caxis([-pi pi]);

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
        nR = size(data_reg, 1);  nF = size(data_reg, 2);
        h_img2 = imagesc(freqs_reg, 1:nR, data_reg);
        set(gca, 'YDir', 'normal', 'Color', [1 1 1]);
        colormap(gca, cmap_circ);  caxis([-pi pi]);

        alpha_reg = ones(nR, nF) * nonsig_alpha;
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

cb = colorbar('Location', 'southoutside');
cb.Ticks      = [-pi -pi/2 0 pi/2 pi];
cb.TickLabels = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};
cb.Position   = [0.25 0.02 0.5 0.015];
cb.Label.String = 'Preferred Phase (rad)';

sgtitle(sprintf('Preferred Phase (%s)', animal), 'FontSize', 14, 'FontWeight', 'bold');
print(f1, fullfile(save_root, 'preferred_phase_comparison.pdf'), '-dpdf');
fprintf('Figure 1 saved.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 2 — Polar Histograms at Key Frequencies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each DV, select the 3 frequencies with the most significant channels
% (in coherence). At each frequency, plot a polar histogram of the
% preferred phase across significant channels. A black line marks the
% circular mean direction, scaled by the mean resultant length.

nTopFreq = 3;

% Identify the top frequencies per DV
top_freq_idx = cell(nDVs, 1);
for row = 1:nDVs
    key = row_keys{row};
    n_sig_per_freq = sum(coh_sig.(key), 1);
    [~, sorted_idx] = sort(n_sig_per_freq, 'descend');
    top_freq_idx{row} = sorted_idx(1:min(nTopFreq, length(sorted_idx)));
end

f2 = figure('Name', ['Polar Histograms - ' animal], ...
    'Units', 'centimeters', 'Position', [1 1 52 38]);
set(f2, 'PaperUnits', 'centimeters', 'PaperSize', [52 38], 'PaperPosition', [0 0 52 38]);

nBins = 18;  % 20-degree bins

for row = 1:nDVs
    key = row_keys{row};

    for fi = 1:nTopFreq
        fidx = top_freq_idx{row}(fi);
        f_hz = freq(fidx);
        sig_ch_coh = find(coh_sig.(key)(:, fidx));

        % --- Coherence polar histogram ---
        sp_idx = (row-1)*nTopFreq*2 + (fi-1)*2 + 1;
        ax = subplot(nDVs, nTopFreq*2, sp_idx, polaraxes);

        if ~isempty(sig_ch_coh)
            phases_coh = coh_phase.(key)(sig_ch_coh, fidx);
            phases_coh = phases_coh(~isnan(phases_coh));
            if ~isempty(phases_coh)
                polarhistogram(ax, phases_coh, nBins, ...
                    'FaceColor', [0.55 0.83 0.78], 'FaceAlpha', 0.7, ...
                    'EdgeColor', [0.3 0.6 0.55]);
                hold(ax, 'on');
                mu = circ_mean(phases_coh);
                R  = circ_r(phases_coh);
                polarplot(ax, [mu mu], [0 R * max(ax.RLim)], 'k-', 'LineWidth', 2.5);
            end
        end
        title(ax, sprintf('%s Coh %.0fHz (n=%d)', row_labels{row}, f_hz, length(sig_ch_coh)), ...
            'FontSize', 7);
        ax.FontSize = 6;

        % --- Regression polar histogram ---
        ax2 = subplot(nDVs, nTopFreq*2, sp_idx + 1, polaraxes);

        if has_reg
            [~, ridx]  = min(abs(freqs_reg - f_hz));  % closest regression freq
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
                    polarplot(ax2, [mu2 mu2], [0 R2 * max(ax2.RLim)], 'k-', 'LineWidth', 2.5);
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

sgtitle(sprintf('Preferred Phase Distribution at Key Frequencies (%s)\nBlack arrow = circular mean, n = sig channels', animal), ...
    'FontSize', 12, 'FontWeight', 'bold');
print(f2, fullfile(save_root, 'preferred_phase_polar.pdf'), '-dpdf');
fprintf('Figure 2 saved.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 3 — Pairwise Phase Consistency Across DVs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For every pair of DVs (6 pairs), compute the circular correlation
% (circ_corrcc) of preferred phases across channels at each frequency.
% Only channels significant in BOTH DVs are included (min 4 required).
% This tests whether the same channels tend to prefer similar phases
% regardless of which DV is used.

measure_pairs = {
    'mua', 'lfp',      'MUA vs LFP';
    'mua', 'RT',       'MUA vs RT';
    'mua', 'hit_miss', 'MUA vs Hit/Miss';
    'lfp', 'RT',       'LFP vs RT';
    'lfp', 'hit_miss', 'LFP vs Hit/Miss';
    'RT',  'hit_miss', 'RT vs Hit/Miss';
    };
nPairs = size(measure_pairs, 1);

% --- Coherence phases ---
rho_coh  = NaN(nPairs, nFreq);
pval_coh = NaN(nPairs, nFreq);

for p = 1:nPairs
    k1 = measure_pairs{p,1};  k2 = measure_pairs{p,2};
    ph1 = coh_phase.(k1);     ph2 = coh_phase.(k2);
    sg1 = coh_sig.(k1);       sg2 = coh_sig.(k2);

    for fi = 1:nFreq
        valid = sg1(:,fi) & sg2(:,fi) & ~isnan(ph1(:,fi)) & ~isnan(ph2(:,fi));
        if sum(valid) < 4, continue; end
        [rho_coh(p,fi), pval_coh(p,fi)] = circ_corrcc(ph1(valid,fi), ph2(valid,fi));
    end
end

% --- Regression phases ---
rho_reg  = NaN(nPairs, nFreq);
pval_reg = NaN(nPairs, nFreq);

if has_reg
    nFreq_r  = length(freqs_reg);
    rho_reg  = NaN(nPairs, nFreq_r);
    pval_reg = NaN(nPairs, nFreq_r);

    for p = 1:nPairs
        k1 = measure_pairs{p,1};  k2 = measure_pairs{p,2};
        ph1 = reg_phase.(k1);     ph2 = reg_phase.(k2);
        sg1 = reg_sig.(k1);       sg2 = reg_sig.(k2);
        nR = min([size(ph1,1), size(ph2,1), size(sg1,1), size(sg2,1)]);

        for fi = 1:nFreq_r
            valid = sg1(1:nR,fi) & sg2(1:nR,fi) & ~isnan(ph1(1:nR,fi)) & ~isnan(ph2(1:nR,fi));
            if sum(valid) < 4, continue; end
            [rho_reg(p,fi), pval_reg(p,fi)] = circ_corrcc(ph1(valid,fi), ph2(valid,fi));
        end
    end
end

% --- Plot ---
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
    % Coherence column
    subplot(nPairs, 2, (p-1)*2 + 1);
    hold on;
    plot(freq, rho_coh(p,:), 'Color', pair_colors(p,:), 'LineWidth', 2);
    sig_f = find(pval_coh(p,:) < 0.05);
    if ~isempty(sig_f)
        plot(freq(sig_f), rho_coh(p,sig_f), '.', 'Color', pair_colors(p,:), 'MarkerSize', 12);
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

    % Regression column
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
fprintf('Figure 3 saved.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 8. MONKEY-AVERAGE: COLLECT DATA FROM BOTH ANIMALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each animal, preferred phase is loaded from precomputed complex-valued
% summary files:
%
%   Level 1 (per channel): phase = angle(coh_complex_all_ch) loaded from
%     channel_avg_results.mat; significance from per-channel coh_perm_complex.mat.
%   Level 2 (per animal):  phase = angle(coh_complex_chan_avg) from
%     channel_avg_results.mat; significance: coh_chan_avg >= thresh_chan_avg.
%   Level 3 (monkey avg):  phase = angle(coh_complex_monkey_avg) from
%     monkey_avg_results.mat; significance: coh_monkey_avg >= thresh_monkey_avg.
%
% Hit/miss uses the ITC equivalent at each level (itc_complex / itc_perm).
% Per-animal level-2 summaries feed Figure 4; pairwise circular correlations
% across channels (level 1) feed Figure 5.

animals_all  = {'hermes', 'klecks'};
nAnimals     = numel(animals_all);
animal_colors = lines(nAnimals);

monkey_save_root = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/sampling_compare/monkey_avg');
if ~exist(monkey_save_root, 'dir'), mkdir(monkey_save_root); end

reg_dvs_all    = {'MUA_ERP_ampl_all', 'LFP_ERP_ampl_all', 'RT', 'hit_miss'};
reg_labels_all = {'mua', 'lfp', 'RT', 'hit_miss'};

% Preallocate per-animal storage
animal_coh_avg     = cell(nAnimals, 1);   % level-2 phase: angle(complex channel avg), coherence
animal_coh_avg_sig = cell(nAnimals, 1);   % level-2 significance mask for channel avg
animal_reg_avg     = cell(nAnimals, 1);   % level-2 phase: circular mean across all channels, regression
animal_reg_avg_sig = cell(nAnimals, 1);   % level-2 significance mask: channel_avg_R.phase > channel_avg_thresh.phase
animal_rho_coh_all = cell(nAnimals, 1);   % pairwise rho (coherence)
animal_rho_reg_all = cell(nAnimals, 1);   % pairwise rho (regression)
animal_coh_sig = cell(nAnimals, 1);   % per-channel significance masks (coherence)
animal_reg_sig = cell(nAnimals, 1);   % per-channel significance masks (regression)

for a = 1:nAnimals
    animalName = animals_all{a};
    fprintf('\n=== Loading preferred phase data for %s ===\n', animalName);

    a_base      = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName]);
    a_coh_root  = fullfile(a_base, 'phase_coherence', 'complex', 'cp10_till_100');
    a_reg_root  = fullfile(a_base, 'multi_lin_reg',   'complex', 'cp10_till_100');
    a_data_root = fullfile(a_base, 'multi_lin_reg', 'cp10_till_100');   % shared input data

    tmp_f   = load(fullfile(a_data_root, 'frequency.mat'));
    a_freq  = tmp_f.frequency;
    a_nFreq = numel(a_freq);
    a_nCh   = 64;

    % --- Load coherence preferred phase + significance masks ---
    % Level 1 phase:  angle(coh_complex_all_ch) from channel_avg_results.mat.
    % Level 1 sig:    |coh_complex_all_ch(ch,:)| >= quantile(max(|perm|),0.95)
    %                 from per-channel coh_perm_complex.mat (max-statistic test).
    % Level 2 phase:  angle(coh_complex_chan_avg) from channel_avg_results.mat.
    % Level 2 sig:    coh_chan_avg >= thresh_chan_avg (same file).
    a_coh_phase = struct();
    a_coh_sig   = struct();
    a_coh_avg   = struct();

    for m = 1:length(coh_measures)
        measure   = coh_measures{m};
        sig_map   = false(a_nCh, a_nFreq);

        % Level 1 phase & Level 2: load from precomputed channel_avg_results.mat.
        %   coh_complex_all_ch   [nCh x nFreq] — per-channel complex coherence
        %   coh_complex_chan_avg  [1 x nFreq]   — complex mean across channels (level 2)
        %   coh_chan_avg / thresh_chan_avg        — magnitude and perm threshold (level 2)
        avg_file = fullfile(a_coh_root, measure, 'all_loc_difflev', 'channel_avg_results.mat');
        if isfile(avg_file)
            tmp_avg   = load(avg_file, 'coh_complex_all_ch', 'coh_complex_chan_avg', ...
                             'coh_chan_avg', 'thresh_chan_avg');
            phase_map = angle(tmp_avg.coh_complex_all_ch);   % level 1 phase

            % Level 1 significance: per-channel max-statistic perm threshold
            for ch = 1:a_nCh
                perm_file = fullfile(a_coh_root, measure, 'all_loc_difflev', ...
                                     num2str(ch), 'coh_perm_complex.mat');
                if ~isfile(perm_file), continue; end
                tmp_p = load(perm_file, 'coh_perm_complex');
                prm   = tmp_p.coh_perm_complex;
                if any(isnan(prm(:))), continue; end
                thr_c = quantile(max(abs(prm), [], 2), 0.95);
                val   = abs(tmp_avg.coh_complex_all_ch(ch,:));
                if any(isnan(val)), continue; end
                sig_map(ch,:) = val >= thr_c;
            end

            % Level 2: angle of complex channel average + significance mask
            ph_avg  = angle(tmp_avg.coh_complex_chan_avg);
            sig_avg = tmp_avg.coh_chan_avg >= tmp_avg.thresh_chan_avg;
            a_coh_avg.(measure)     = ph_avg;    % raw phase, no NaN masking
            a_coh_avg_sig.(measure) = sig_avg;   % significance stored separately
        else
            phase_map = NaN(a_nCh, a_nFreq);
            a_coh_avg.(measure)     = NaN(1, a_nFreq);
            a_coh_avg_sig.(measure) = false(1, a_nFreq);
        end

        a_coh_phase.(measure) = phase_map;
        a_coh_sig.(measure)   = sig_map;
    end

    % --- Hit/miss coherence preferred phase + significance masks ---
    % Level 1 phase:  angle(itc_complex) per channel from itc.mat.
    % Level 1 sig:    itc >= quantile(max(itc_perm),0.95) from itc_perm.mat.
    % Level 2 phase:  angle(itc_complex_chan_avg) if saved, else phase_chan_avg_itc,
    %                 from channel_avg_results_itc.mat; masked by itc_chan_avg >= thresh.
    a_corr_root = fullfile(a_base, 'phase_correlation', 'complex', 'cp10_till_100');
    hm_phase = NaN(a_nCh, a_nFreq);
    for ch = 1:a_nCh
        itc_file = fullfile(a_corr_root, 'hit_miss', 'all_loc_difflev', num2str(ch), 'itc.mat');
        if ~exist(itc_file, 'file'), continue; end
        tmp = load(itc_file, 'itc_complex');
        if isfield(tmp, 'itc_complex') && ~any(isnan(tmp.itc_complex))
            hm_phase(ch,:) = angle(tmp.itc_complex);
        end
    end
    a_coh_phase.hit_miss = hm_phase;

    % Level 1 significance: per-channel max-statistic test on ITC
    hm_sig = false(a_nCh, a_nFreq);
    for ch = 1:a_nCh
        ch_folder = fullfile(a_corr_root, 'hit_miss', 'all_loc_difflev', num2str(ch));
        itc_file  = fullfile(ch_folder, 'itc.mat');
        perm_file = fullfile(ch_folder, 'itc_perm.mat');
        if ~exist(itc_file, 'file') || ~exist(perm_file, 'file'), continue; end

        tmp_i = load(itc_file,  'itc_complex');       val = abs(tmp_i.itc_complex);
        tmp_p = load(perm_file, 'itc_perm_complex');  prm = abs(tmp_p.itc_perm_complex);
        if any(isnan(val)) || any(isnan(prm(:))), continue; end

        tmax = max(prm, [], 2);
        thr  = quantile(tmax, 0.95);
        hm_sig(ch,:) = val >= thr;
    end
    a_coh_sig.hit_miss = hm_sig;

    % Level 2: load angle(itc_complex_chan_avg) if available, else phase_chan_avg_itc,
    % from channel_avg_results_itc.mat; NaN where itc_chan_avg < thresh_chan_avg_itc.
    hm_avg_file = fullfile(a_corr_root, 'hit_miss', 'all_loc_difflev', 'channel_avg_results_itc.mat');
    if isfile(hm_avg_file)
        tmp_hm = load(hm_avg_file);
        if isfield(tmp_hm, 'itc_complex_chan_avg')
            ph_hm = angle(tmp_hm.itc_complex_chan_avg);
        else
            ph_hm = tmp_hm.phase_chan_avg_itc;
        end
        sig_hm = tmp_hm.itc_chan_avg >= tmp_hm.thresh_chan_avg_itc;
        a_coh_avg.hit_miss     = ph_hm;    % raw phase, no NaN masking
        a_coh_avg_sig.hit_miss = sig_hm;   % significance stored separately
    else
        a_coh_avg.hit_miss     = NaN(1, a_nFreq);
        a_coh_avg_sig.hit_miss = false(1, a_nFreq);
    end

    animal_coh_avg{a}     = a_coh_avg;
    animal_coh_avg_sig{a} = a_coh_avg_sig;

    % --- Load regression preferred phase + significance masks ---
    a_reg_file  = fullfile(a_reg_root, 'multi_regression_channelwise_R2.mat');
    a_has_reg   = exist(a_reg_file, 'file');
    a_reg_phase = struct();
    a_reg_sig   = struct();

    if a_has_reg
        a_reg = load(a_reg_file, 'reg_results');
        a_rr  = a_reg.reg_results;

        session_dirs = dir(fullfile(a_base, [animalName '_*']));
        if ~isempty(session_dirs)
            freq_file_a = fullfile(a_base, session_dirs(1).name, ...
                'Phase_analysis', 'hit_miss', '100iter_cut@cp_m10', '1', 'freqpow.mat');
        else
            freq_file_a = '';
        end
        if ~isempty(freq_file_a) && exist(freq_file_a, 'file')
            tmp = load(freq_file_a);
            a_freqs_reg = tmp.freqpow.freq;
        else
            a_freqs_reg = a_freq;
        end
        a_nFreq_r = length(a_freqs_reg);

        for m = 1:length(reg_dvs_all)
            dv  = reg_dvs_all{m};
            lbl = reg_labels_all{m};

            if isfield(a_rr, dv) && isfield(a_rr.(dv), 'phi_pref')
                a_reg_phase.(lbl) = a_rr.(dv).phi_pref;
            else
                a_reg_phase.(lbl) = NaN(a_nCh, a_nFreq_r);
            end

            sig_ph = false(a_nCh, a_nFreq_r);
            if isfield(a_rr, dv)
                numCh_r = size(a_rr.(dv).R2_phase, 1);
                for ch = 1:min(numCh_r, length(a_rr.(dv).thresholds))
                    if isfield(a_rr.(dv).thresholds(ch), 'thresh_phase')
                        thr = a_rr.(dv).thresholds(ch).thresh_phase;
                        if ~isempty(thr)
                            sig_ph(ch,:) = a_rr.(dv).R2_phase(ch,:) > thr(1);
                        end
                    end
                end
            end
            a_reg_sig.(lbl) = sig_ph;
        end
    else
        a_freqs_reg = a_freq;
        a_nFreq_r   = a_nFreq;
        for m = 1:length(reg_labels_all)
            a_reg_phase.(reg_labels_all{m}) = NaN(a_nCh, a_nFreq);
            a_reg_sig.(reg_labels_all{m})   = false(a_nCh, a_nFreq);
        end
    end

    % --- Circular mean across channels (regression) ---
    % avg_ph: circular mean over ALL valid (non-NaN) channels.
    % sig_f:  channel_avg_R.phase(f) > channel_avg_thresh.phase, i.e. the
    %         observed mean R² across channels exceeds the 95th pctile of
    %         the max-stat null (precomputed in regress_stats_R2.m and stored
    %         in reg_results.(dv).channel_avg_R / .channel_avg_thresh).
    a_reg_avg     = struct();
    a_reg_avg_sig = struct();
    for m = 1:length(reg_labels_all)
        key = reg_labels_all{m};
        dv  = reg_dvs_all{m};
        ph  = a_reg_phase.(key);
        nF  = size(ph, 2);
        avg_ph = NaN(1, nF);
        sig_f  = false(1, nF);
        for f = 1:nF
            valid_all = ~isnan(ph(:,f));
            if any(valid_all)
                avg_ph(f) = angle(mean(exp(1i * ph(valid_all,f))));
            end
        end
        % Channel-average significance from precomputed R² null distribution
        if a_has_reg && isfield(a_rr, dv) && ...
                isfield(a_rr.(dv), 'channel_avg_R') && ...
                isfield(a_rr.(dv), 'channel_avg_thresh')
            obs_avg = a_rr.(dv).channel_avg_R.phase;
            thr_avg = a_rr.(dv).channel_avg_thresh.phase;
            nF_s = min(length(obs_avg), nF);
            sig_f(1:nF_s) = obs_avg(1:nF_s) > thr_avg;
        end
        a_reg_avg.(key)     = avg_ph;
        a_reg_avg_sig.(key) = sig_f;
    end
    animal_reg_avg{a}     = a_reg_avg;
    animal_reg_avg_sig{a} = a_reg_avg_sig;

    % --- Pairwise phase consistency (coherence) ---
    a_rho_coh = NaN(nPairs, a_nFreq);
    for p = 1:nPairs
        k1 = measure_pairs{p,1};  k2 = measure_pairs{p,2};
        ph1 = a_coh_phase.(k1);   ph2 = a_coh_phase.(k2);
        sg1 = a_coh_sig.(k1);     sg2 = a_coh_sig.(k2);
        for fi = 1:a_nFreq
            valid = sg1(:,fi) & sg2(:,fi) & ~isnan(ph1(:,fi)) & ~isnan(ph2(:,fi));
            if sum(valid) < 4, continue; end
            a_rho_coh(p,fi) = circ_corrcc(ph1(valid,fi), ph2(valid,fi));
        end
    end
    animal_rho_coh_all{a} = a_rho_coh;

    % --- Pairwise phase consistency (regression) ---
    a_rho_reg = NaN(nPairs, a_nFreq_r);
    if a_has_reg
        for p = 1:nPairs
            k1 = measure_pairs{p,1};  k2 = measure_pairs{p,2};
            ph1 = a_reg_phase.(k1);   ph2 = a_reg_phase.(k2);
            sg1 = a_reg_sig.(k1);     sg2 = a_reg_sig.(k2);
            nR = min([size(ph1,1), size(ph2,1), size(sg1,1), size(sg2,1)]);
            for fi = 1:a_nFreq_r
                valid = sg1(1:nR,fi) & sg2(1:nR,fi) & ~isnan(ph1(1:nR,fi)) & ~isnan(ph2(1:nR,fi));
                if sum(valid) < 4, continue; end
                a_rho_reg(p,fi) = circ_corrcc(ph1(valid,fi), ph2(valid,fi));
            end
        end
    end
    animal_rho_reg_all{a} = a_rho_reg;
    animal_coh_sig{a}     = a_coh_sig;
    animal_reg_sig{a}     = a_reg_sig;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 4 — Monkey-Average Preferred Phase Heatmap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each DV, 3 rows per column:
%   Row 1 = hermes  — level-2 phase: angle(coh_complex_chan_avg), NaN where not sig
%   Row 2 = klecks  — same
%   Row 3 = monkey avg — coherence: angle(coh_complex_monkey_avg) from saved file;
%                        regression: angle(mean(exp(i*[row1;row2]))) across animals
% Transparent where NaN (not significant or no data).

% Level 3 monkey-average preferred phase: loaded from precomputed files,
% phase = angle(complex monkey average), NaN where magnitude < perm threshold.
%   MUA/LFP/RT: angle(coh_complex_monkey_avg) from monkey_avg_results.mat,
%               masked by coh_monkey_avg >= thresh_monkey_avg.
%   Hit/miss:   angle(itc_complex_monkey_avg) if saved, else phase_monkey_avg_itc,
%               from monkey_avg_results_itc.mat; masked by itc_monkey_avg >= thresh.
results_combined_pre = '/mnt/hpc/projects/MWSampling/4Shivangi/results_combined';
coh_combined_pre     = fullfile(results_combined_pre, 'phase_coherence',  'complex', 'cp10_till_100');
corr_combined_pre    = fullfile(results_combined_pre, 'phase_correlation', 'complex', 'cp10_till_100');

monkey_avg_phase = struct();
monkey_avg_sig   = struct();
dv_keys_coh = {'mua','lfp','RT'};
for m = 1:length(dv_keys_coh)
    key     = dv_keys_coh{m};
    mk_file = fullfile(coh_combined_pre, key, 'all_loc_difflev', 'monkey_avg_results.mat');
    if isfile(mk_file)
        mk    = load(mk_file, 'coh_complex_monkey_avg', 'coh_monkey_avg', 'thresh_monkey_avg');
        monkey_avg_phase.(key) = angle(mk.coh_complex_monkey_avg);   % raw phase
        monkey_avg_sig.(key)   = mk.coh_monkey_avg >= mk.thresh_monkey_avg;
    else
        monkey_avg_phase.(key) = NaN(1, nFreq);
        monkey_avg_sig.(key)   = false(1, nFreq);
    end
end

hm_mk_file = fullfile(corr_combined_pre, 'hit_miss_itc', 'all_loc_difflev', 'monkey_avg_results_itc.mat');
if isfile(hm_mk_file)
    hm_mk = load(hm_mk_file);
    if isfield(hm_mk, 'itc_complex_monkey_avg')
        ph_hm_mk = angle(hm_mk.itc_complex_monkey_avg);
    else
        ph_hm_mk = hm_mk.phase_monkey_avg_itc;
    end
    monkey_avg_phase.hit_miss = ph_hm_mk;
    monkey_avg_sig.hit_miss   = hm_mk.itc_monkey_avg >= hm_mk.thresh_monkey_avg_itc;
else
    monkey_avg_phase.hit_miss = NaN(1, nFreq);
    monkey_avg_sig.hit_miss   = false(1, nFreq);
end

f4 = figure('Name', 'Monkey-Avg Preferred Phase', ...
    'Units', 'centimeters', 'Position', [1 1 36 40]);
set(f4, 'PaperUnits', 'centimeters', 'PaperSize', [36 40], 'PaperPosition', [0 0 36 40]);

for row = 1:nDVs
    key = row_keys{row};

    % --- Coherence column ---
    subplot(nDVs, 2, (row-1)*2 + 1);
    phase_stack = NaN(3, nFreq);
    for a = 1:nAnimals
        phase_stack(a,:) = animal_coh_avg{a}.(key);   % raw phase (no NaN masking)
    end
    phase_stack(3,:) = monkey_avg_phase.(key);         % level 3: raw phase

    % Build alpha: NaN (no data) = 0, non-sig = nonsig_alpha, sig = 1
    alpha_stack = ones(3, nFreq) * nonsig_alpha;
    alpha_stack(isnan(phase_stack)) = 0;
    for a = 1:nAnimals
        alpha_stack(a, animal_coh_avg_sig{a}.(key)) = 1;
    end
    alpha_stack(3, monkey_avg_sig.(key)) = 1;

    h_img = imagesc(freq, 1:3, phase_stack);
    set(gca, 'YDir', 'normal', 'Color', [1 1 1]);
    colormap(gca, cmap_circ);  caxis([-pi pi]);
    set(h_img, 'AlphaData', alpha_stack);
    yticks(1:3); yticklabels([animals_all, {'Monkey avg'}]);
    xlabel('Frequency (Hz)');
    title(sprintf('%s — Coherence', row_labels{row}), 'FontSize', 9);
    set(gca, 'FontSize', 8, 'Box', 'on');
    if row == 1
        text(0.5, 1.22, 'Coherence (coh\_complex)', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
    end

    % --- Regression column ---
    subplot(nDVs, 2, (row-1)*2 + 2);
    phase_stack_r = NaN(3, length(freqs_reg));
    for a = 1:nAnimals
        vals = animal_reg_avg{a}.(key);
        nF   = min(length(vals), length(freqs_reg));
        phase_stack_r(a, 1:nF) = vals(1:nF);
    end
    phase_stack_r(3,:) = angle(mean(exp(1i * phase_stack_r(1:nAnimals,:)), 1, 'omitnan'));

    % Build alpha: NaN (no data) = 0, non-sig = nonsig_alpha, sig = 1
    % Per-animal rows: channel_avg_R.phase > channel_avg_thresh.phase
    %   (mean R² across channels vs max-stat null; from regress_stats_R2.m)
    % Monkey-avg row:  monkey_avg_obs.phase > thresh_monkey.phase
    %   (mean channel_avg_R² across animals vs max-stat null; from
    %    regress_stats_R2.m, saved in results_combined/monkey_avg_results.mat)
    alpha_stack_r = ones(3, length(freqs_reg)) * nonsig_alpha;
    alpha_stack_r(isnan(phase_stack_r)) = 0;
    for a = 1:nAnimals
        sig_vals = animal_reg_avg_sig{a}.(key);
        nF = min(length(sig_vals), length(freqs_reg));
        alpha_stack_r(a, find(sig_vals(1:nF))) = 1;
    end
    % Monkey-average significance: load from precomputed monkey_avg_results.mat
    dv_idx = find(strcmp(reg_labels_all, key), 1);
    mk_sig_r = false(1, length(freqs_reg));
    if ~isempty(dv_idx)
        mk_avg_file = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/results_combined', ...
            'multi_lin_reg', 'complex', 'cp10_till_100', reg_dvs_all{dv_idx}, 'monkey_avg_results.mat');
        if isfile(mk_avg_file)
            mk_r = load(mk_avg_file, 'monkey_avg_obs', 'thresh_monkey');
            nF_r = min(length(mk_r.monkey_avg_obs.phase), length(freqs_reg));
            mk_sig_r(1:nF_r) = mk_r.monkey_avg_obs.phase(1:nF_r) > mk_r.thresh_monkey.phase;
        else
            warning('monkey_avg_results.mat not found for %s regression; monkey-avg row unfaded.', key);
        end
    end
    alpha_stack_r(3, mk_sig_r) = 1;

    h_img2 = imagesc(freqs_reg, 1:3, phase_stack_r); %check from here for the new significance
    set(gca, 'YDir', 'normal', 'Color', [1 1 1]);
    colormap(gca, cmap_circ);  caxis([-pi pi]);
    set(h_img2, 'AlphaData', alpha_stack_r);
    yticks(1:3); yticklabels([animals_all, {'Monkey avg'}]);
    xlabel('Frequency (Hz)');
    title(sprintf('%s — Regression', row_labels{row}), 'FontSize', 9);
    set(gca, 'FontSize', 8, 'Box', 'on');
    if row == 1
        text(0.5, 1.22, 'Regression (\phi_{pref})', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
    end
end

cb = colorbar('Location', 'southoutside');
cb.Ticks      = [-pi -pi/2 0 pi/2 pi];
cb.TickLabels = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};
cb.Position   = [0.25 0.02 0.5 0.015];
cb.Label.String = 'Preferred Phase (rad)';

sgtitle('Preferred Phase', ...
    'FontSize', 14, 'FontWeight', 'bold');
print(f4, fullfile(monkey_save_root, 'monkey_avg_preferred_phase.pdf'), '-dpdf');
fprintf('Figure 4 saved.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 9. HIT-ONLY AND MISS-ONLY PREFERRED PHASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recompute preferred phase separately for hit trials and miss trials.
% This uses the raw single-trial phase data (ph_all_sess.mat) rather than
% the precomputed coherence results. For each channel and frequency:
%
%   Coherence estimate:  angle( sum( DV .* exp(i*phase) ) )
%     = DV-weighted mean phase direction
%
%   Regression estimate: fit full model
%     phi_pref = atan2(beta_sin, beta_cos)
%
% Trial alignment: some DVs (e.g. MUA, LFP) have different trial counts
% than the phase data due to artefact rejection. Trials are matched via
% the unique trial identifier in column 14 of trialinfo.

hm_dvs       = {'MUA_ERP_ampl_all', 'LFP_ERP_ampl_all', 'RT', 'hit_miss'};
hm_ti_fields = {'MUA_ERP_trialinfo', 'LFP_ERP_trialinfo', 'RT_trialinfo', 'trialinfo'};
hm_labels    = {'mua', 'lfp', 'RT', 'hit_miss'};
hm_row_labels = {'MUA', 'LFP', 'RT', 'Hit/Miss'};
nDV = length(hm_dvs);

coh_subtitles_hm = {'MUA (DV-weighted phase)', 'LFP (DV-weighted phase)', ...
    'RT (DV-weighted phase)', 'Hit/Miss (circular mean phase)'};
reg_subtitles_hm = {'MUA (full model \phi_{pref})', ...
    'LFP (full model \phi_{pref})', ...
    'RT (full model \phi_{pref})', ...
    'Hit/Miss (circular mean phase)'};

animal_hit_coh_phase  = cell(nAnimals, 1);
animal_miss_coh_phase = cell(nAnimals, 1);
animal_hit_reg_phase  = cell(nAnimals, 1);
animal_miss_reg_phase = cell(nAnimals, 1);

for a = 1:nAnimals
    animalName = animals_all{a};
    fprintf('\n=== Computing hit/miss preferred phase for %s ===\n', animalName);

    a_base      = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName]);
    a_data_root = fullfile(a_base, 'multi_lin_reg', 'cp10_till_100');   % shared input data

    tmp_f   = load(fullfile(a_data_root, 'frequency.mat'));
    a_freq  = tmp_f.frequency;
    a_nFreq = numel(a_freq);
    a_nCh   = 64;

    % Load combined single-trial phase data
    a_ph = load(fullfile(a_data_root, 'ph_all_sess.mat'), 'ph_comb');
    a_ph = a_ph.ph_comb;
    nTrials_phase = size(a_ph.phase_all, 1);

    hit_coh  = struct();  miss_coh = struct();
    hit_reg  = struct();  miss_reg = struct();

    for m = 1:nDV
        dv_name  = hm_dvs{m};
        ti_field = hm_ti_fields{m};
        lbl      = hm_labels{m};

        hit_coh.(lbl)  = NaN(a_nCh, a_nFreq);
        miss_coh.(lbl) = NaN(a_nCh, a_nFreq);
        hit_reg.(lbl)  = NaN(a_nCh, a_nFreq);
        miss_reg.(lbl) = NaN(a_nCh, a_nFreq);

        if strcmp(lbl, 'hit_miss')
            dv_vals = repmat(a_ph.trialinfo(:,20),1,64);
        else
            dv_vals = a_ph.(dv_name);
        end

        % Use DV-specific trialinfo if available (trial counts may differ)
        if isfield(a_ph, ti_field)
            dv_ti = a_ph.(ti_field);
        else
            dv_ti = a_ph.trialinfo;
        end

        % Align DV trials to phase trials via unique trial ID (column 14)
        nTrials_dv = size(dv_vals, 1);
        if nTrials_dv == nTrials_phase
            phase_idx = (1:nTrials_phase)';
            dv_idx    = (1:nTrials_dv)';
            trial_ti  = a_ph.trialinfo;
        else
            [~, phase_idx, dv_idx] = intersect(a_ph.trialinfo(:,14), dv_ti(:,14));
            trial_ti = a_ph.trialinfo(phase_idx, :);
        end

        % Split into hit (code 1) and miss (code 5)
        hit_mask  = trial_ti(:,20) == 1;
        miss_mask = trial_ti(:,20) == 5;

        ph_hit_idx  = phase_idx(hit_mask);
        ph_miss_idx = phase_idx(miss_mask);
        dv_hit_idx  = dv_idx(hit_mask);
        dv_miss_idx = dv_idx(miss_mask);

        hit_coh_map  = NaN(a_nCh, a_nFreq);
        miss_coh_map = NaN(a_nCh, a_nFreq);
        hit_reg_map  = NaN(a_nCh, a_nFreq);
        miss_reg_map = NaN(a_nCh, a_nFreq);

        for ch = 1:a_nCh
            phase_ch  = a_ph.phase_all(:,:,ch);
            amp_ch    = a_ph.amp_all(:,:,ch);       % [nTrials x nFreq]
            pup_ch    = a_ph.pup_baseline(:,ch);    % [nTrials x 1]
            mua_bl_ch = a_ph.MUA_baseline(:,ch);    % [nTrials x 1]
            if size(dv_vals, 2) < ch, continue; end
            dv_ch = dv_vals(:, ch);

            for f = 1:a_nFreq
                % --- Hit trials ---
                ph_h = phase_ch(ph_hit_idx, f);
                dv_h = dv_ch(dv_hit_idx);
                ok_h = ~isnan(ph_h) & ~isnan(dv_h);

                if strcmp(lbl, 'hit_miss')
                    % For binary hit_miss, within hit-only trials the
                    % DV is constant (all 1s), so use the complex ITC
                    % (mean of unit phasors) — angle gives preferred phase.
                    if any(ok_h)
                        hit_coh_map(ch,f) = mean(exp(1i * ph_h(ok_h)));
                    end
                    hit_reg_map(ch,f) = angle(mean(exp(1i * ph_h(ok_h))));
                else
                    if any(ok_h)
                        hit_coh_map(ch,f) = sum(dv_h(ok_h) .* exp(1i * ph_h(ok_h)));
                    end

                    % Regression: full model matching regress_stats_R2.m
                    % DV ~ 1 + pup_baseline + MUA_baseline + amp + sin(phase) + cos(phase)
                    amp_h = amp_ch(ph_hit_idx, f);
                    pup_h = pup_ch(ph_hit_idx);
                    mua_h = mua_bl_ch(ph_hit_idx);
                    ok_reg = ok_h & ~isnan(amp_h) & ~isnan(pup_h) & ~isnan(mua_h);
                    if sum(ok_reg) >= 6
                        X_h = [ones(sum(ok_reg),1), pup_h(ok_reg), mua_h(ok_reg), ...
                            amp_h(ok_reg), sin(ph_h(ok_reg)), cos(ph_h(ok_reg))];
                        b_h = regress(dv_h(ok_reg), X_h);
                        hit_reg_map(ch,f) = atan2(b_h(end-1), b_h(end));
                    end
                end

                % --- Miss trials ---
                ph_m = phase_ch(ph_miss_idx, f);
                dv_m = dv_ch(dv_miss_idx);
                ok_m = ~isnan(ph_m) & ~isnan(dv_m);

                if strcmp(lbl, 'hit_miss')
                    % Same as above: DV is constant (all 0s for misses),
                    % so use the complex ITC — angle gives preferred phase.
                    if any(ok_m)
                        miss_coh_map(ch,f) = mean(exp(1i * ph_m(ok_m)));
                    end
                    miss_reg_map(ch,f) = angle(mean(exp(1i * ph_m(ok_m))));
                else
                    if any(ok_m)
                        miss_coh_map(ch,f) = sum(dv_m(ok_m) .* exp(1i * ph_m(ok_m)));
                    end

                    amp_m = amp_ch(ph_miss_idx, f);
                    pup_m = pup_ch(ph_miss_idx);
                    mua_m = mua_bl_ch(ph_miss_idx);
                    ok_reg = ok_m & ~isnan(amp_m) & ~isnan(pup_m) & ~isnan(mua_m);
                    if sum(ok_reg) >= 6
                        X_m = [ones(sum(ok_reg),1), pup_m(ok_reg), mua_m(ok_reg), ...
                            amp_m(ok_reg), sin(ph_m(ok_reg)), cos(ph_m(ok_reg))];
                        b_m = regress(dv_m(ok_reg), X_m);
                        miss_reg_map(ch,f) = atan2(b_m(end-1), b_m(end));
                    end
                end
            end
        end

        hit_coh.(lbl)  = hit_coh_map;
        miss_coh.(lbl) = miss_coh_map;
        hit_reg.(lbl)  = hit_reg_map;
        miss_reg.(lbl) = miss_reg_map;
    end

    animal_hit_coh_phase{a}  = hit_coh;
    animal_miss_coh_phase{a} = miss_coh;
    animal_hit_reg_phase{a}  = hit_reg;
    animal_miss_reg_phase{a} = miss_reg;

    %% Per-animal hit-only and miss-only figures
    a_save_root = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/sampling_compare', animalName);
    if ~exist(a_save_root, 'dir'), mkdir(a_save_root); end

    % --- Hit-only figure ---
    fh = figure('Name', sprintf('Hit-Only Phase - %s', animalName), ...
        'Units', 'centimeters', 'Position', [1 1 36 40]);
    set(fh, 'PaperUnits', 'centimeters', 'PaperSize', [36 40], 'PaperPosition', [0 0 36 40]);

    for row = 1:nDV
        lbl = hm_labels{row};

        subplot(nDV, 2, (row-1)*2 + 1);
        data = angle(hit_coh.(lbl));
        h_img = imagesc(freq, 1:a_nCh, data);
        set(gca, 'YDir', 'normal', 'Color', [1 1 1]);
        colormap(gca, cmap_circ);  caxis([-pi pi]);
        set(h_img, 'AlphaData', ~isnan(data));
        xlabel('Frequency (Hz)'); ylabel('Channel');
        title(coh_subtitles_hm{row}, 'FontSize', 9);
        set(gca, 'FontSize', 8, 'Box', 'on');
        if row == 1
            text(0.5, 1.22, 'Coherence (DV-weighted)', 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
        end

        subplot(nDV, 2, (row-1)*2 + 2);
        data = hit_reg.(lbl);
        h_img = imagesc(freq, 1:a_nCh, data);
        set(gca, 'YDir', 'normal', 'Color', [1 1 1]);
        colormap(gca, cmap_circ);  caxis([-pi pi]);
        set(h_img, 'AlphaData', ~isnan(data));
        xlabel('Frequency (Hz)'); ylabel('Channel');
        title(reg_subtitles_hm{row}, 'FontSize', 9);
        set(gca, 'FontSize', 8, 'Box', 'on');
        if row == 1
            text(0.5, 1.22, 'Regression (\phi_{pref})', 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
        end
    end

    cb = colorbar('Location', 'southoutside');
    cb.Ticks = [-pi -pi/2 0 pi/2 pi];
    cb.TickLabels = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};
    cb.Position = [0.25 0.02 0.5 0.015];
    cb.Label.String = 'Preferred Phase (rad)';
    sgtitle(sprintf('Hit Trials Only (%s)', animalName), ...
        'FontSize', 14, 'FontWeight', 'bold');
    print(fh, fullfile(a_save_root, 'preferred_phase_hit_only.pdf'), '-dpdf');
    fprintf('Hit-only figure for %s saved.\n', animalName);

    % --- Miss-only figure ---
    fm = figure('Name', sprintf('Miss-Only Phase - %s', animalName), ...
        'Units', 'centimeters', 'Position', [1 1 36 40]);
    set(fm, 'PaperUnits', 'centimeters', 'PaperSize', [36 40], 'PaperPosition', [0 0 36 40]);

    for row = 1:nDV
        lbl = hm_labels{row};

        subplot(nDV, 2, (row-1)*2 + 1);
        data = angle(miss_coh.(lbl));
        h_img = imagesc(freq, 1:a_nCh, data);
        set(gca, 'YDir', 'normal', 'Color', [1 1 1]);
        colormap(gca, cmap_circ);  caxis([-pi pi]);
        set(h_img, 'AlphaData', ~isnan(data));
        xlabel('Frequency (Hz)'); ylabel('Channel');
        title(coh_subtitles_hm{row}, 'FontSize', 9);
        set(gca, 'FontSize', 8, 'Box', 'on');
        if row == 1
            text(0.5, 1.22, 'Coherence (DV-weighted)', 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
        end

        subplot(nDV, 2, (row-1)*2 + 2);
        data = miss_reg.(lbl);
        h_img = imagesc(freq, 1:a_nCh, data);
        set(gca, 'YDir', 'normal', 'Color', [1 1 1]);
        colormap(gca, cmap_circ);  caxis([-pi pi]);
        set(h_img, 'AlphaData', ~isnan(data));
        xlabel('Frequency (Hz)'); ylabel('Channel');
        title(reg_subtitles_hm{row}, 'FontSize', 9);
        set(gca, 'FontSize', 8, 'Box', 'on');
        if row == 1
            text(0.5, 1.22, 'Regression (\phi_{pref})', 'Units', 'normalized', ...
                'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
        end
    end

    cb = colorbar('Location', 'southoutside');
    cb.Ticks = [-pi -pi/2 0 pi/2 pi];
    cb.TickLabels = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};
    cb.Position = [0.25 0.02 0.5 0.015];
    cb.Label.String = 'Preferred Phase (rad)';
    sgtitle(sprintf('Miss Trials Only (%s)', animalName), ...
        'FontSize', 14, 'FontWeight', 'bold');
    print(fm, fullfile(a_save_root, 'preferred_phase_miss_only.pdf'), '-dpdf');
    fprintf('Miss-only figure for %s saved.\n', animalName);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 5 — Monkey-Average Hit-Only Preferred Phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Same layout as Figure 4 but computed from hit trials only.
% Each animal's contribution is the circular mean across all channels.

f_hit_avg = figure('Name', 'Monkey-Avg Hit-Only Phase', ...
    'Units', 'centimeters', 'Position', [1 1 36 40]);
set(f_hit_avg, 'PaperUnits', 'centimeters', 'PaperSize', [36 40], ...
    'PaperPosition', [0 0 36 40]);

for row = 1:nDV
    lbl = hm_labels{row};

    % Coherence column — complex channel average, then complex animal average
    subplot(nDV, 2, (row-1)*2 + 1);
    phase_stack = NaN(3, nFreq);
    cplx_stack  = complex(NaN(nAnimals, nFreq));
    for a = 1:nAnimals
        ph_cplx = animal_hit_coh_phase{a}.(lbl);   % complex [nCh x nFreq_a]
        avg_cplx = complex(NaN(1, size(ph_cplx, 2)));
        for f = 1:size(ph_cplx, 2)
            valid = ~isnan(ph_cplx(:,f));
            if any(valid), avg_cplx(f) = mean(ph_cplx(valid, f)); end
        end
        nF = min(length(avg_cplx), nFreq);
        phase_stack(a, 1:nF) = angle(avg_cplx(1:nF));
        cplx_stack(a,  1:nF) = avg_cplx(1:nF);
    end
    phase_stack(3,:) = angle(mean(cplx_stack, 1, 'omitnan'));

    h_img = imagesc(freq, 1:3, phase_stack);
    set(gca, 'YDir', 'normal', 'Color', [1 1 1]);
    colormap(gca, cmap_circ);  caxis([-pi pi]);
    set(h_img, 'AlphaData', ~isnan(phase_stack));
    yticks(1:3); yticklabels([animals_all, {'Monkey avg'}]);
    xlabel('Frequency (Hz)');
    title(coh_subtitles_hm{row}, 'FontSize', 9);
    set(gca, 'FontSize', 8, 'Box', 'on');
    if row == 1
        text(0.5, 1.22, 'Coherence (DV-weighted)', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
    end

    % Regression column
    subplot(nDV, 2, (row-1)*2 + 2);
    phase_stack = NaN(3, nFreq);
    for a = 1:nAnimals
        ph  = animal_hit_reg_phase{a}.(lbl);
        avg = NaN(1, size(ph,2));
        for f = 1:size(ph,2)
            v = ~isnan(ph(:,f));
            if any(v), avg(f) = angle(mean(exp(1i * ph(v,f)))); end
        end
        nF = min(length(avg), nFreq);
        phase_stack(a, 1:nF) = avg(1:nF);
    end
    phase_stack(3,:) = angle(mean(exp(1i * phase_stack(1:nAnimals,:)), 1, 'omitnan'));

    h_img = imagesc(freq, 1:3, phase_stack);
    set(gca, 'YDir', 'normal', 'Color', [1 1 1]);
    colormap(gca, cmap_circ);  caxis([-pi pi]);
    set(h_img, 'AlphaData', ~isnan(phase_stack));
    yticks(1:3); yticklabels([animals_all, {'Monkey avg'}]);
    xlabel('Frequency (Hz)');
    title(reg_subtitles_hm{row}, 'FontSize', 9);
    set(gca, 'FontSize', 8, 'Box', 'on');
    if row == 1
        text(0.5, 1.22, 'Regression (\phi_{pref})', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
    end
end

cb = colorbar('Location', 'southoutside');
cb.Ticks      = [-pi -pi/2 0 pi/2 pi];
cb.TickLabels = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};
cb.Position   = [0.25 0.02 0.5 0.015];
cb.Label.String = 'Preferred Phase (rad)';
sgtitle('Hit Trials Only', ...
    'FontSize', 14, 'FontWeight', 'bold');
print(f_hit_avg, fullfile(monkey_save_root, 'monkey_avg_preferred_phase_hit_only.pdf'), '-dpdf');
fprintf('Monkey-average hit-only figure saved.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 6 — Monkey-Average Miss-Only Preferred Phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f_miss_avg = figure('Name', 'Monkey-Avg Miss-Only Phase', ...
    'Units', 'centimeters', 'Position', [1 1 36 40]);
set(f_miss_avg, 'PaperUnits', 'centimeters', 'PaperSize', [36 40], ...
    'PaperPosition', [0 0 36 40]);

for row = 1:nDV
    lbl = hm_labels{row};

    % Coherence column — complex channel average, then complex animal average
    subplot(nDV, 2, (row-1)*2 + 1);
    phase_stack = NaN(3, nFreq);
    cplx_stack  = complex(NaN(nAnimals, nFreq));
    for a = 1:nAnimals
        ph_cplx = animal_miss_coh_phase{a}.(lbl);   % complex [nCh x nFreq_a]
        avg_cplx = complex(NaN(1, size(ph_cplx, 2)));
        for f = 1:size(ph_cplx, 2)
            valid = ~isnan(ph_cplx(:,f));
            if any(valid), avg_cplx(f) = mean(ph_cplx(valid, f)); end
        end
        nF = min(length(avg_cplx), nFreq);
        phase_stack(a, 1:nF) = angle(avg_cplx(1:nF));
        cplx_stack(a,  1:nF) = avg_cplx(1:nF);
    end
    phase_stack(3,:) = angle(mean(cplx_stack, 1, 'omitnan'));

    h_img = imagesc(freq, 1:3, phase_stack);
    set(gca, 'YDir', 'normal', 'Color', [1 1 1]);
    colormap(gca, cmap_circ);  caxis([-pi pi]);
    set(h_img, 'AlphaData', ~isnan(phase_stack));
    yticks(1:3); yticklabels([animals_all, {'Monkey avg'}]);
    xlabel('Frequency (Hz)');
    title(coh_subtitles_hm{row}, 'FontSize', 9);
    set(gca, 'FontSize', 8, 'Box', 'on');
    if row == 1
        text(0.5, 1.22, 'Coherence (DV-weighted)', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
    end

    % Regression column
    subplot(nDV, 2, (row-1)*2 + 2);
    phase_stack = NaN(3, nFreq);
    for a = 1:nAnimals
        ph  = animal_miss_reg_phase{a}.(lbl);
        avg = NaN(1, size(ph,2));
        for f = 1:size(ph,2)
            v = ~isnan(ph(:,f));
            if any(v), avg(f) = angle(mean(exp(1i * ph(v,f)))); end
        end
        nF = min(length(avg), nFreq);
        phase_stack(a, 1:nF) = avg(1:nF);
    end
    phase_stack(3,:) = angle(mean(exp(1i * phase_stack(1:nAnimals,:)), 1, 'omitnan'));

    h_img = imagesc(freq, 1:3, phase_stack);
    set(gca, 'YDir', 'normal', 'Color', [1 1 1]);
    colormap(gca, cmap_circ);  caxis([-pi pi]);
    set(h_img, 'AlphaData', ~isnan(phase_stack));
    yticks(1:3); yticklabels([animals_all, {'Monkey avg'}]);
    xlabel('Frequency (Hz)');
    title(reg_subtitles_hm{row}, 'FontSize', 9);
    set(gca, 'FontSize', 8, 'Box', 'on');
    if row == 1
        text(0.5, 1.22, 'Regression (\phi_{pref})', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
    end
end

cb = colorbar('Location', 'southoutside');
cb.Ticks      = [-pi -pi/2 0 pi/2 pi];
cb.TickLabels = {'-\pi', '-\pi/2', '0', '\pi/2', '\pi'};
cb.Position   = [0.25 0.02 0.5 0.015];
cb.Label.String = 'Preferred Phase (rad)';
sgtitle('Miss Trials Only', ...
    'FontSize', 14, 'FontWeight', 'bold');
print(f_miss_avg, fullfile(monkey_save_root, 'monkey_avg_preferred_phase_miss_only.pdf'), '-dpdf');
fprintf('Monkey-average miss-only figure saved.\n');
