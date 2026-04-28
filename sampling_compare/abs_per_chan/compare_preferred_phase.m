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
%             mean across sig channels + cross-animal average)
%     Fig 5 — Monkey-average pairwise phase consistency
%
%   HIT-ONLY / MISS-ONLY FIGURES:
%     Per-animal and monkey-average heatmaps of preferred phase computed
%     separately from hit trials and miss trials. For each trial subset,
%     the coherence estimate uses DV-weighted mean phase and the regression
%     estimate fits DV ~ cos(phase) + sin(phase).

clear all; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. SETTINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

animal = 'klecks';   % 'hermes' or 'klecks'

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

coh_root  = fullfile(base_results, 'phase_coherence', 'abs_per_chan', 'cp10_till_100');
reg_root  = fullfile(base_results, 'multi_lin_reg', 'abs_per_chan', 'cp10_till_100');
corr_root = fullfile(base_results, 'phase_correlation', 'abs_per_chan', 'cp10_till_100');
data_root = fullfile(base_results, 'multi_lin_reg', 'cp10_till_100');  % shared input data (ph_all_sess.mat)

save_root = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/sampling_compare', animal);
if ~exist(save_root, 'dir'), mkdir(save_root); end

% Frequency axis (shared across coherence channels)
load(fullfile(coh_root, 'frequency.mat'));
freq  = frequency;
nFreq = numel(freq);
nCh   = 64;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. LOAD COHERENCE PREFERRED PHASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each DV and channel, load the phase_spec field from the precomputed
% coherence.mat file. phase_spec is the amplitude-weighted mean phase
% direction (1 x nFreq) stored during the coherence analysis.

coh_phase    = struct();
coh_measures = {'mua', 'lfp', 'RT', 'hit_miss'};

for m = 1:length(coh_measures)
    measure   = coh_measures{m};
    phase_map = NaN(nCh, nFreq);

    for ch = 1:nCh
        ch_folder = fullfile(coh_root, measure, 'all_loc_difflev', num2str(ch));
        coh_file  = fullfile(ch_folder, 'coherence.mat');
        if ~exist(coh_file, 'file'), continue; end

        tmp = load(coh_file, 'phase_spec');
        if isfield(tmp, 'phase_spec') && ~any(isnan(tmp.phase_spec))
            phase_map(ch,:) = tmp.phase_spec;
        end
    end

    coh_phase.(measure) = phase_map;
end

% --- Hit/miss: compute preferred phase from ITC with inverted miss phases ---
% Standard amplitude-weighted coherence is not meaningful for binary DVs:
% multiplying unit vectors by 0/1 would simply discard all miss trials.
% Instead we use the ITC-with-inverted-miss-phases approach (VanRullen 2016):
%   1. Flip miss trial phases by pi so that, if hits and misses sit at
%      opposite phases, all vectors now point in the same direction.
%   2. Take the complex mean across ALL trials (hits + flipped misses).
%   3. angle(complex mean) = preferred phase for hits.
% This is analogous to phase_spec for continuous DVs: the resultant angle
% tells you which phase is associated with hits (and the opposite with misses).

data_file = fullfile(data_root, 'ph_all_sess.mat');
if exist(data_file, 'file')
    tmp     = load(data_file, 'ph_comb');
    ph_comb = tmp.ph_comb;

    % Select all valid trials: hits (code 1) and misses (code 5)
    all_idx     = find(ph_comb.trialinfo(:,20) == 1 | ph_comb.trialinfo(:,20) == 5);
    hit_labels  = (ph_comb.trialinfo(all_idx, 20) == 1);
    miss_labels = ~hit_labels;

    hm_phase_map = NaN(nCh, nFreq);

    for ch = 1:nCh
        phase_ch = ph_comb.phase_all(all_idx, :, ch);   % [nTrials x nFreq]

        % Invert miss phases by pi
        phase_inv = phase_ch;
        phase_inv(miss_labels, :) = mod(phase_ch(miss_labels, :) + pi, 2*pi) - pi;

        % Complex mean — angle gives the preferred hit phase
        for f = 1:nFreq
            cavg = mean(exp(1i * phase_inv(:, f)));
            hm_phase_map(ch, f) = angle(cavg);
        end
    end

    coh_phase.hit_miss = hm_phase_map;
    clear ph_comb
else
    warning('ph_all_sess.mat not found — hit_miss coherence phase remains NaN');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. LOAD REGRESSION PREFERRED PHASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each DV, load phi_pref from the channelwise multiple regression
% results. phi_pref = atan2(beta_sin, beta_cos) from the model
% DV ~ cos(phase) + sin(phase) + intercept.

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
% the distribution of per-iteration maxima (cluster-corrected across freq).

% --- Coherence significance (from coh_perm.mat) ---
coh_sig = struct();
for m = 1:length(coh_measures)
    measure = coh_measures{m};
    sig_map = false(nCh, nFreq);

    for ch = 1:nCh
        ch_folder = fullfile(coh_root, measure, 'all_loc_difflev', num2str(ch));
        coh_file  = fullfile(ch_folder, 'coherence.mat');
        perm_file = fullfile(ch_folder, 'coh_perm.mat');
        if ~exist(coh_file, 'file') || ~exist(perm_file, 'file'), continue; end

        tmp = load(coh_file, 'coh');   val = tmp.coh;
        tmp = load(perm_file, 'coh_perm');  prm = tmp.coh_perm;
        if any(isnan(val)) || any(isnan(prm(:))), continue; end

        tmax = max(prm, [], 2);           % max across frequencies per iteration
        thr  = quantile(tmax, 0.95);      % 95th percentile of max distribution
        sig_map(ch,:) = val >= thr;
    end

    coh_sig.(measure) = sig_map;
end

% --- Hit/miss significance from ITC permutation (itc.mat / itc_perm.mat) ---
% The ITC magnitude and its permutation null live in the phase_correlation
% results (computed in Corr_RT_hitmiss.m). We apply the same cluster-
% corrected threshold: observed ITC >= 95th percentile of per-iteration
% frequency-maxima from the shuffled-label null distribution.
hm_sig_map = false(nCh, nFreq);
for ch = 1:nCh
    ch_folder = fullfile(corr_root, 'hit_miss', 'all_loc_difflev', num2str(ch));
    itc_file  = fullfile(ch_folder, 'itc.mat');
    perm_file = fullfile(ch_folder, 'itc_perm.mat');
    if ~exist(itc_file, 'file') || ~exist(perm_file, 'file'), continue; end

    tmp = load(itc_file, 'itc');        val = tmp.itc;
    tmp = load(perm_file, 'itc_perm');  prm = tmp.itc_perm;
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
% Phase is circular (-pi to pi), so the colormap must wrap. We use a
% pastel HSV ring (low saturation, high brightness) for readability.

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
% non-significant bins are faded to 20% so the overall pattern remains
% visible while significance is clear.

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

    % --- Column 2: Regression ---
    subplot(nDVs, 2, (row-1)*2 + 2);
    if has_reg
        data_reg = reg_phase.(key);
        nR = size(data_reg, 1);  nF = size(data_reg, 2);
        h_img2 = imagesc(freqs_reg, 1:nR, data_reg);
        set(gca, 'YDir', 'normal', 'Color', [1 1 1]);
        colormap(gca, cmap_circ);  caxis([-pi pi]);

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
% Repeat the data loading for both animals. For each animal, compute:
%   (a) circular mean of preferred phase across significant channels (1 x nFreq)
%   (b) pairwise circular correlation between DVs (nPairs x nFreq)
% These per-animal summaries are then averaged across animals.

animals_all  = {'hermes', 'klecks'};
nAnimals     = numel(animals_all);
animal_colors = lines(nAnimals);

monkey_save_root = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/sampling_compare/monkey_avg');
if ~exist(monkey_save_root, 'dir'), mkdir(monkey_save_root); end

reg_dvs_all    = {'MUA_ERP_ampl_all', 'LFP_ERP_ampl_all', 'RT', 'hit_miss'};
reg_labels_all = {'mua', 'lfp', 'RT', 'hit_miss'};

% Preallocate per-animal storage
animal_coh_avg     = cell(nAnimals, 1);   % circular mean phase (coherence)
animal_reg_avg     = cell(nAnimals, 1);   % circular mean phase (regression)
animal_rho_coh_all = cell(nAnimals, 1);   % pairwise rho (coherence)
animal_rho_reg_all = cell(nAnimals, 1);   % pairwise rho (regression)

for a = 1:nAnimals
    animalName = animals_all{a};
    fprintf('\n=== Loading preferred phase data for %s ===\n', animalName);

    a_base     = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName]);
    a_coh_root  = fullfile(a_base, 'phase_coherence', 'abs_per_chan', 'cp10_till_100');
    a_reg_root  = fullfile(a_base, 'multi_lin_reg', 'abs_per_chan', 'cp10_till_100');
    a_data_root = fullfile(a_base, 'multi_lin_reg', 'cp10_till_100');  % shared input data

    tmp_f   = load(fullfile(a_coh_root, 'frequency.mat'));
    a_freq  = tmp_f.frequency;
    a_nFreq = numel(a_freq);
    a_nCh   = 64;

    % --- Load coherence preferred phase + significance masks ---
    a_coh_phase = struct();
    a_coh_sig   = struct();

    for m = 1:length(coh_measures)
        measure   = coh_measures{m};
        phase_map = NaN(a_nCh, a_nFreq);
        sig_map   = false(a_nCh, a_nFreq);

        for ch = 1:a_nCh
            ch_folder = fullfile(a_coh_root, measure, 'all_loc_difflev', num2str(ch));
            coh_file  = fullfile(ch_folder, 'coherence.mat');
            perm_file = fullfile(ch_folder, 'coh_perm.mat');
            if ~exist(coh_file, 'file'), continue; end

            tmp = load(coh_file, 'phase_spec');
            if isfield(tmp, 'phase_spec') && ~any(isnan(tmp.phase_spec))
                phase_map(ch,:) = tmp.phase_spec;
            end

            if exist(perm_file, 'file')
                tmp_c = load(coh_file, 'coh');
                tmp_p = load(perm_file, 'coh_perm');
                if ~any(isnan(tmp_c.coh)) && ~any(isnan(tmp_p.coh_perm(:)))
                    tmax = max(tmp_p.coh_perm, [], 2);
                    thr  = quantile(tmax, 0.95);
                    sig_map(ch,:) = tmp_c.coh >= thr;
                end
            end
        end

        a_coh_phase.(measure) = phase_map;
        a_coh_sig.(measure)   = sig_map;
    end

    % --- Hit/miss: ITC with inverted miss phases (same as section 3) ---
    a_corr_root = fullfile(a_base, 'phase_correlation', 'abs_per_chan', 'cp10_till_100');
    a_data_file = fullfile(a_data_root, 'ph_all_sess.mat');
    if exist(a_data_file, 'file')
        tmp_ph  = load(a_data_file, 'ph_comb');
        a_ph    = tmp_ph.ph_comb;

        a_all_idx     = find(a_ph.trialinfo(:,20) == 1 | a_ph.trialinfo(:,20) == 5);
        a_hit_labels  = (a_ph.trialinfo(a_all_idx, 20) == 1);
        a_miss_labels = ~a_hit_labels;

        hm_phase = NaN(a_nCh, a_nFreq);
        for ch = 1:a_nCh
            phase_ch  = a_ph.phase_all(a_all_idx, :, ch);
            phase_inv = phase_ch;
            phase_inv(a_miss_labels, :) = mod(phase_ch(a_miss_labels, :) + pi, 2*pi) - pi;
            for f = 1:a_nFreq
                hm_phase(ch, f) = angle(mean(exp(1i * phase_inv(:, f))));
            end
        end
        a_coh_phase.hit_miss = hm_phase;
        clear a_ph
    end

    % Hit/miss significance from ITC permutation
    hm_sig = false(a_nCh, a_nFreq);
    for ch = 1:a_nCh
        ch_folder = fullfile(a_corr_root, 'hit_miss', 'all_loc_difflev', num2str(ch));
        itc_file  = fullfile(ch_folder, 'itc.mat');
        perm_file = fullfile(ch_folder, 'itc_perm.mat');
        if ~exist(itc_file, 'file') || ~exist(perm_file, 'file'), continue; end

        tmp_i = load(itc_file, 'itc');        val = tmp_i.itc;
        tmp_p = load(perm_file, 'itc_perm');  prm = tmp_p.itc_perm;
        if any(isnan(val)) || any(isnan(prm(:))), continue; end

        tmax = max(prm, [], 2);
        thr  = quantile(tmax, 0.95);
        hm_sig(ch,:) = val >= thr;
    end
    a_coh_sig.hit_miss = hm_sig;

    % --- Circular mean across significant channels (coherence) ---
    a_coh_avg = struct();
    for m = 1:nDVs
        key = row_keys{m};
        ph  = a_coh_phase.(key);
        sg  = a_coh_sig.(key);
        avg_ph = NaN(1, a_nFreq);
        for f = 1:a_nFreq
            valid = sg(:,f) & ~isnan(ph(:,f));
            if sum(valid) >= 1
                avg_ph(f) = angle(mean(exp(1i * ph(valid,f))));
            end
        end
        a_coh_avg.(key) = avg_ph;
    end
    animal_coh_avg{a} = a_coh_avg;

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

    % --- Circular mean across significant channels (regression) ---
    a_reg_avg = struct();
    for m = 1:length(reg_labels_all)
        key = reg_labels_all{m};
        ph  = a_reg_phase.(key);
        sg  = a_reg_sig.(key);
        nR  = min(size(ph,1), size(sg,1));
        nF  = min(size(ph,2), size(sg,2));
        avg_ph = NaN(1, nF);
        for f = 1:nF
            valid = sg(1:nR,f) & ~isnan(ph(1:nR,f));
            if sum(valid) >= 1
                avg_ph(f) = angle(mean(exp(1i * ph(valid,f))));
            end
        end
        a_reg_avg.(key) = avg_ph;
    end
    animal_reg_avg{a} = a_reg_avg;

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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 4 — Monkey-Average Preferred Phase Heatmap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each DV: row 1 = hermes, row 2 = klecks, row 3 = circular mean
% across the two animals. Transparent where no data.

f4 = figure('Name', 'Monkey-Avg Preferred Phase', ...
    'Units', 'centimeters', 'Position', [1 1 36 40]);
set(f4, 'PaperUnits', 'centimeters', 'PaperSize', [36 40], 'PaperPosition', [0 0 36 40]);

for row = 1:nDVs
    key = row_keys{row};

    % --- Coherence column ---
    subplot(nDVs, 2, (row-1)*2 + 1);
    phase_stack = NaN(3, nFreq);
    for a = 1:nAnimals
        phase_stack(a,:) = animal_coh_avg{a}.(key);
    end
    phase_stack(3,:) = angle(mean(exp(1i * phase_stack(1:nAnimals,:)), 1, 'omitnan'));

    h_img = imagesc(freq, 1:3, phase_stack);
    set(gca, 'YDir', 'normal', 'Color', [1 1 1]);
    colormap(gca, cmap_circ);  caxis([-pi pi]);
    set(h_img, 'AlphaData', ~isnan(phase_stack));
    yticks(1:3); yticklabels([animals_all, {'Monkey avg'}]);
    xlabel('Frequency (Hz)');
    title(sprintf('%s — Coherence', row_labels{row}), 'FontSize', 9);
    set(gca, 'FontSize', 8, 'Box', 'on');
    if row == 1
        text(0.5, 1.22, 'Coherence (phase\_spec)', 'Units', 'normalized', ...
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

    h_img2 = imagesc(freqs_reg, 1:3, phase_stack_r);
    set(gca, 'YDir', 'normal', 'Color', [1 1 1]);
    colormap(gca, cmap_circ);  caxis([-pi pi]);
    set(h_img2, 'AlphaData', ~isnan(phase_stack_r));
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

sgtitle('Monkey-Average Preferred Phase (circular mean across sig channels)', ...
    'FontSize', 14, 'FontWeight', 'bold');
print(f4, fullfile(monkey_save_root, 'monkey_avg_preferred_phase.pdf'), '-dpdf');
fprintf('Figure 4 saved.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 5 — Monkey-Average Phase Consistency Across Measures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thin lines = individual animals, thick line = average across animals.

rho_coh_monkey = mean(cat(3, animal_rho_coh_all{:}), 3, 'omitnan');
rho_reg_monkey = mean(cat(3, animal_rho_reg_all{:}), 3, 'omitnan');

f5 = figure('Name', 'Monkey-Avg Phase Consistency', ...
    'Units', 'centimeters', 'Position', [1 1 42 36]);
set(f5, 'PaperUnits', 'centimeters', 'PaperSize', [42 36], 'PaperPosition', [0 0 42 36]);

for p = 1:nPairs
    % Coherence column
    subplot(nPairs, 2, (p-1)*2 + 1);
    hold on;
    for a = 1:nAnimals
        plot(freq, animal_rho_coh_all{a}(p,:), 'Color', [animal_colors(a,:) 0.4], 'LineWidth', 1);
    end
    plot(freq, rho_coh_monkey(p,:), 'Color', pair_colors(p,:), 'LineWidth', 2.5);
    yline(0, 'k--', 'LineWidth', 0.5);
    xlabel('Frequency (Hz)'); ylabel('\rho_{circ}');
    title(sprintf('Coh: %s', measure_pairs{p,3}), 'FontSize', 9);
    ylim([-1 1]);
    set(gca, 'FontSize', 8, 'Box', 'on');
    if p == 1
        text(0.5, 1.22, 'Coherence', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'FontSize', 13, 'FontWeight', 'bold');
        legend([animals_all, {'Monkey avg'}], 'Location', 'best', 'FontSize', 6);
    end

    % Regression column
    subplot(nPairs, 2, (p-1)*2 + 2);
    hold on;
    for a = 1:nAnimals
        plot(freqs_reg, animal_rho_reg_all{a}(p,:), 'Color', [animal_colors(a,:) 0.4], 'LineWidth', 1);
    end
    plot(freqs_reg, rho_reg_monkey(p,:), 'Color', pair_colors(p,:), 'LineWidth', 2.5);
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

sgtitle('Monkey-Average Pairwise Phase Consistency (thin = per animal, thick = average)', ...
    'FontSize', 12, 'FontWeight', 'bold');
print(f5, fullfile(monkey_save_root, 'monkey_avg_phase_consistency.pdf'), '-dpdf');
fprintf('Figure 5 saved.\n');

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
%   Regression estimate: fit DV ~ cos(phase) + sin(phase) + 1
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

    a_base     = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName]);
    a_coh_root = fullfile(a_base, 'multi_lin_reg', 'cp10_till_100');
    a_reg_root = fullfile(a_base, 'multi_lin_reg', 'cp10_till_100');

    tmp_f   = load(fullfile(a_coh_root, 'frequency.mat'));
    a_freq  = tmp_f.frequency;
    a_nFreq = numel(a_freq);
    a_nCh   = 64;

    % Load combined single-trial phase data
    a_ph = load(fullfile(a_reg_root, 'ph_all_sess.mat'), 'ph_comb');
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
                    % DV is constant (all 1s), so DV-weighted coherence
                    % and regression are degenerate. Use the circular
                    % mean phase (ITC direction) instead — this is the
                    % preferred phase at which hits occur.
                    hit_coh_map(ch,f) = angle(mean(exp(1i * ph_h(ok_h))));
                    hit_reg_map(ch,f) = angle(mean(exp(1i * ph_h(ok_h))));
                else
                    hit_coh_map(ch,f) = angle(sum(dv_h(ok_h) .* exp(1i * ph_h(ok_h))));

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
                    % so use circular mean phase — the preferred phase
                    % at which misses occur.
                    miss_coh_map(ch,f) = angle(mean(exp(1i * ph_m(ok_m))));
                    miss_reg_map(ch,f) = angle(mean(exp(1i * ph_m(ok_m))));
                else
                    miss_coh_map(ch,f) = angle(sum(dv_m(ok_m) .* exp(1i * ph_m(ok_m))));

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
        data = hit_coh.(lbl);
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
    sgtitle(sprintf('Preferred Phase — Hit Trials Only (%s)', animalName), ...
        'FontSize', 14, 'FontWeight', 'bold');
    annotation(fh, 'textbox', [0.01 0.93 0.98 0.04], 'String', ...
        ['MUA/LFP/RT: coherence = DV-weighted mean phase, regression = full model '...
        '(DV ~ pup + MUA_{bl} + amp + sin + cos) \phi_{pref} = atan2(\beta_{sin},\beta_{cos}). ' ...
        'Hit/Miss: both columns show circular mean phase angle(mean(e^{i\theta}))'], ...
        'EdgeColor', 'none', 'FontSize', 6, 'FitBoxToText', 'off', ...
        'HorizontalAlignment', 'center', 'Interpreter', 'tex');
    print(fh, fullfile(a_save_root, 'preferred_phase_hit_only.pdf'), '-dpdf');
    fprintf('Hit-only figure for %s saved.\n', animalName);

    % --- Miss-only figure ---
    fm = figure('Name', sprintf('Miss-Only Phase - %s', animalName), ...
        'Units', 'centimeters', 'Position', [1 1 36 40]);
    set(fm, 'PaperUnits', 'centimeters', 'PaperSize', [36 40], 'PaperPosition', [0 0 36 40]);

    for row = 1:nDV
        lbl = hm_labels{row};

        subplot(nDV, 2, (row-1)*2 + 1);
        data = miss_coh.(lbl);
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
    sgtitle(sprintf('Preferred Phase — Miss Trials Only (%s)', animalName), ...
        'FontSize', 14, 'FontWeight', 'bold');
    annotation(fm, 'textbox', [0.01 0.93 0.98 0.04], 'String', ...
        ['MUA/LFP/RT: coherence = DV-weighted mean phase, regression = full model '...
        '(DV ~ pup + MUA_{bl} + amp + sin + cos) \phi_{pref} = atan2(\beta_{sin},\beta_{cos}). ' ...
        'Hit/Miss: both columns show circular mean phase angle(mean(e^{i\theta}))'], ...
        'EdgeColor', 'none', 'FontSize', 6, 'FitBoxToText', 'off', ...
        'HorizontalAlignment', 'center', 'Interpreter', 'tex');
    print(fm, fullfile(a_save_root, 'preferred_phase_miss_only.pdf'), '-dpdf');
    fprintf('Miss-only figure for %s saved.\n', animalName);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 6 — Monkey-Average Hit-Only Preferred Phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Same layout as Figure 4 but computed from hit trials only.
% Each animal's contribution is the circular mean across all channels.

f_hit_avg = figure('Name', 'Monkey-Avg Hit-Only Phase', ...
    'Units', 'centimeters', 'Position', [1 1 36 40]);
set(f_hit_avg, 'PaperUnits', 'centimeters', 'PaperSize', [36 40], ...
    'PaperPosition', [0 0 36 40]);

for row = 1:nDV
    lbl = hm_labels{row};

    % Coherence column
    subplot(nDV, 2, (row-1)*2 + 1);
    phase_stack = NaN(3, nFreq);
    for a = 1:nAnimals
        ph  = animal_hit_coh_phase{a}.(lbl);
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
sgtitle('Monkey-Average Preferred Phase — Hit Trials Only', ...
    'FontSize', 14, 'FontWeight', 'bold');
annotation(f_hit_avg, 'textbox', [0.01 0.93 0.98 0.04], 'String', ...
    ['MUA/LFP/RT: coherence = DV-weighted mean phase, regression = full model '...
    '(DV ~ pup + MUA_{bl} + amp + sin + cos) \phi_{pref} = atan2(\beta_{sin},\beta_{cos}). ' ...
    'Hit/Miss: both columns show circular mean phase angle(mean(e^{i\theta}))'],...
    'EdgeColor', 'none', 'FontSize', 6, 'FitBoxToText', 'off', ...
    'HorizontalAlignment', 'center', 'Interpreter', 'tex');
print(f_hit_avg, fullfile(monkey_save_root, 'monkey_avg_preferred_phase_hit_only.pdf'), '-dpdf');
fprintf('Monkey-average hit-only figure saved.\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 7 — Monkey-Average Miss-Only Preferred Phase
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f_miss_avg = figure('Name', 'Monkey-Avg Miss-Only Phase', ...
    'Units', 'centimeters', 'Position', [1 1 36 40]);
set(f_miss_avg, 'PaperUnits', 'centimeters', 'PaperSize', [36 40], ...
    'PaperPosition', [0 0 36 40]);

for row = 1:nDV
    lbl = hm_labels{row};

    % Coherence column
    subplot(nDV, 2, (row-1)*2 + 1);
    phase_stack = NaN(3, nFreq);
    for a = 1:nAnimals
        ph  = animal_miss_coh_phase{a}.(lbl);
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
sgtitle('Monkey-Average Preferred Phase — Miss Trials Only', ...
    'FontSize', 14, 'FontWeight', 'bold');
annotation(f_miss_avg, 'textbox', [0.01 0.93 0.98 0.04], 'String', ...
    ['MUA/LFP/RT: coherence = DV-weighted mean phase, regression = full model '...
    '(DV ~ pup + MUA_{bl} + amp + sin + cos) \phi_{pref} = atan2(\beta_{sin},\beta_{cos}). ' ...
    'Hit/Miss: both columns show circular mean phase angle(mean(e^{i\theta}))'], ...
    'EdgeColor', 'none', 'FontSize', 6, 'FitBoxToText', 'off', ...
    'HorizontalAlignment', 'center', 'Interpreter', 'tex');
print(f_miss_avg, fullfile(monkey_save_root, 'monkey_avg_preferred_phase_miss_only.pdf'), '-dpdf');
fprintf('Monkey-average miss-only figure saved.\n');
