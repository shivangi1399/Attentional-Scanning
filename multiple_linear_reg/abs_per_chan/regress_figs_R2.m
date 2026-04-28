clearvars;
close all;
clc

% Description:
% -------------
% This script creates a set of visualizations for R² (variance explained)
% from the channel-wise regression analysis in regress_stats_R.m.
% Loads: multi_regression_channelwise_R2.mat

%% Load data

info_folder = '/mnt/hpc/projects/MWSampling/4Shivangi/results_hermes/multi_lin_reg/abs_per_chan/cp10_till_100';
cd(info_folder)
load('multi_regression_channelwise_R2.mat', 'reg_results');

% Load frequency vector
cd('/mnt/hpc/projects/MWSampling/4Shivangi/results_hermes/hermes_20170824_attentional-sampling_1/Phase_analysis/hit_miss/100iter_cut@cp_m10/1')
load('freqpow.mat')
freqs = freqpow.freq;

% Dependent variables
Y_vars   = {'RT','MUA_ERP_ampl_all','LFP_ERP_ampl_all','hit_miss'};
Y_labels = {'RT', 'MUA ERP', 'LFP ERP', 'Hit/Miss'};

p_thresh = 0.05;

%% Color schemes

% Define colors for predictors: Phase, MUA, Amp, Amp+Phase, Full
colors_pred = [
    0 0.4470 0.7410;   % Phase:     blue
    0.4660 0.6740 0.1882; % MUA:    green
    0.0 0.6 0.6;       % Amp:       teal
    0.6 0.3 0.6;       % Amp+Phase: purple
    0.85 0.33 0.1      % Full model: orange
];

% Custom colormap for R² heatmaps (white -> dark teal)
r2_cmap = [linspace(1,0,64)', linspace(1,0.5,64)', linspace(1,0.5,64)'];

% Pink colormap for phase strength (R)
pink_cmap = [linspace(1,0.8,64)', linspace(1,0.3,64)', linspace(1,0.6,64)'];

%% Helper: extract p-value from stats struct
% stats is an array of structs: reg_results.(dv).stats(ch).p_phase etc.

%% FIGURE 1: OVERVIEW HEATMAPS (R² per DV, all predictors)

nDV = length(Y_vars);
nPred = 5;  % Phase, MUA, Amp, Amp+Phase, Full

figure('Name','R² Heatmaps Overview','Position',[50 50 2400 1000]);

pred_fields  = {'R2_phase','R2_MUA','R2_Amp','R2_AmpPhase','R2_any'};
pred_labels  = {'Phase','MUA','Amp','Amp+Phase','Full Model'};
thresh_fields = {'thresh_phase','thresh_MUA','thresh_Amp','thresh_AmpPhase','thresh_any'};
p_fields      = {'p_phase','p_MUA','p_Amp','p_AmpPhase','p_any'};

for d = 1:nDV

    depVarName = Y_vars{d};
    numCh   = size(reg_results.(depVarName).R2_phase, 1);
    numFreq = length(freqs);

    % Build significance masks from per-channel thresholds
    sig_masks = false(nPred, numCh, numFreq);

    for ch = 1:min(numCh, length(reg_results.(depVarName).thresholds))
        for pm = 1:nPred
            tf = thresh_fields{pm};
            if isfield(reg_results.(depVarName).thresholds(ch), tf)
                thr = reg_results.(depVarName).thresholds(ch).(tf);
                if ~isempty(thr)
                    R2 = reg_results.(depVarName).(pred_fields{pm})(ch,:);
                    sig_masks(pm,ch,:) = R2 > thr(1);
                end
            end
        end
    end

    % Column 1: Phase strength (R)
    subplot(nDV, nPred+1, (d-1)*(nPred+1) + 1);
    imagesc(freqs, 1:numCh, reg_results.(depVarName).R_phase);
    colormap(gca, pink_cmap); colorbar;
    xlabel('Frequency (Hz)'); ylabel('Channel');
    title([Y_labels{d} ': Phase Strength (R)']);
    set(gca,'YDir','normal');
    ylim([1 numCh]); yticks(1:8:numCh);

    % Columns 2-6: R² heatmaps masked by significance
    for pm = 1:nPred
        subplot(nDV, nPred+1, (d-1)*(nPred+1) + 1 + pm);
        R2_mat = reg_results.(depVarName).(pred_fields{pm});
        R2_masked = R2_mat .* squeeze(sig_masks(pm,:,:));
        imagesc(freqs, 1:numCh, R2_masked);
        colormap(gca, r2_cmap); colorbar;
        clim([0 max(R2_masked(:)+eps)]);
        xlabel('Frequency (Hz)'); ylabel('Channel');
        title([Y_labels{d} ': ' pred_labels{pm} ' R²']);
        set(gca,'YDir','normal');
        ylim([1 numCh]); yticks(1:8:numCh);
    end

end

sgtitle('R² Heatmaps per Predictor (FWER-masked, non-sig set to 0)',...
        'FontSize',16,'FontWeight','bold');

%% FIGURE 2: SIGNIFICANCE HEATMAPS (binary, all DVs)

figure('Name','Significant Frequencies per Channel (R²)','Position',[50 50 2400 1000]);

for d = 1:nDV

    depVarName = Y_vars{d};
    numCh   = size(reg_results.(depVarName).R2_phase, 1);
    numFreq = length(freqs);

    sig_masks = false(nPred, numCh, numFreq);
    for ch = 1:min(numCh, length(reg_results.(depVarName).thresholds))
        for pm = 1:nPred
            tf = thresh_fields{pm};
            if isfield(reg_results.(depVarName).thresholds(ch), tf)
                thr = reg_results.(depVarName).thresholds(ch).(tf);
                if ~isempty(thr)
                    R2 = reg_results.(depVarName).(pred_fields{pm})(ch,:);
                    sig_masks(pm,ch,:) = R2 > thr(1);
                end
            end
        end
    end

    for pm = 1:nPred
        subplot(nDV, nPred, (d-1)*nPred + pm);
        imagesc(freqs, 1:numCh, squeeze(sig_masks(pm,:,:)));
        colormap(gca, [0.9 0.9 0.9; colors_pred(pm,:)]);
        clim([0 1]);
        xlabel('Frequency (Hz)'); ylabel('Channel');
        title([Y_labels{d} ': ' pred_labels{pm} ' Sig']);
        set(gca,'YDir','normal');
        ylim([1 numCh]); yticks(1:8:numCh);
    end

end

sgtitle('Significant Frequencies per Channel – R² (FWER corrected)',...
        'FontSize',16,'FontWeight','bold');

%% FIGURE 3: CHANNEL-AVERAGE R² WITH THRESHOLD

figure('Name','Channel-Average R²','Position',[50 50 1800 1000]);

avg_fields   = {'phase','MUA','Amp','AmpPhase','any'};  % fields in channel_avg_R / channel_avg_thresh

for d = 1:nDV

    depVarName = Y_vars{d};

    for pm = 1:nPred
        subplot(nDV, nPred, (d-1)*nPred + pm);
        hold on;

        af = avg_fields{pm};

        if ~isfield(reg_results.(depVarName), 'channel_avg_R') || ...
           ~isfield(reg_results.(depVarName).channel_avg_R, af)
            title([Y_labels{d} ': ' pred_labels{pm} ' (no data)']);
            continue
        end

        avg_R2  = reg_results.(depVarName).channel_avg_R.(af);
        avg_thr = reg_results.(depVarName).channel_avg_thresh.(af);

        % Find significant frequencies
        sig_idx = avg_R2 > avg_thr;

        % Plot R² curve
        plot(freqs, avg_R2, 'Color', colors_pred(pm,:), 'LineWidth', 2.5);

        % Threshold line
        yline(avg_thr, 'k--', 'LineWidth', 2);

        % Shade significant regions
        if any(sig_idx)
            y_fill_bot = zeros(1, length(freqs));
            y_fill_top = avg_R2;
            y_fill_top(~sig_idx) = 0;
            fill([freqs fliplr(freqs)], [y_fill_top fliplr(y_fill_bot)], ...
                colors_pred(pm,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        end

        xlabel('Frequency (Hz)'); ylabel('R² (channel avg)');
        title([Y_labels{d} ': ' pred_labels{pm}]);
        grid on;
    end

end

sgtitle('Channel-Average R² with FWER Threshold (shaded = significant)',...
        'FontSize',16,'FontWeight','bold');

%% FIGURE 4: PER-CHANNEL R² FOR EACH SIGNIFICANT CHANNEL (one figure per DV x predictor)

for d = 1:nDV

    depVarName = Y_vars{d};
    numCh = size(reg_results.(depVarName).R2_phase, 1);

    for pm = 1:nPred

        pvals  = ones(numCh, 1);
        thresh = nan(numCh, 1);

        for ch = 1:min(numCh, length(reg_results.(depVarName).stats))
            pf = p_fields{pm};
            tf = thresh_fields{pm};
            if isfield(reg_results.(depVarName).stats(ch), pf)
                tmp = reg_results.(depVarName).stats(ch).(pf);
                if ~isempty(tmp), pvals(ch) = tmp(1); end
            end
            if isfield(reg_results.(depVarName).thresholds(ch), tf)
                tmp = reg_results.(depVarName).thresholds(ch).(tf);
                if ~isempty(tmp), thresh(ch) = tmp(1); end
            end
        end

        sig_ch = find(pvals < p_thresh);

        if isempty(sig_ch)
            fprintf('No significant channels for %s - %s\n', Y_labels{d}, pred_labels{pm});
            continue
        end

        nSig  = length(sig_ch);
        nCols = ceil(sqrt(nSig));
        nRows = ceil(nSig / nCols);

        figure('Name', [Y_labels{d} ' - ' pred_labels{pm} ' R²'], ...
               'Position', [100 100 1800 900]);

        R2_mat = reg_results.(depVarName).(pred_fields{pm});

        % Common y-limit
        maxR2 = max(R2_mat(sig_ch,:), [], 'all');
        maxR2 = max(maxR2, max(thresh(sig_ch), [], 'omitnan')) * 1.1;

        for i = 1:nSig
            ch = sig_ch(i);
            subplot(nRows, nCols, i); hold on;

            plot(freqs, R2_mat(ch,:), 'Color', colors_pred(pm,:), 'LineWidth', 2);

            if ~isnan(thresh(ch))
                yline(thresh(ch), 'k--', 'LineWidth', 1.5);
            end

            xlabel('Frequency (Hz)'); ylabel('R²');
            title(sprintf('Ch %d', ch));
            ylim([0 maxR2]);
            xlim([freqs(1) freqs(end)]);
            grid on;
        end

        sgtitle(sprintf('%s    %s  (n = %d)', Y_labels{d}, pred_labels{pm}, nSig), ...
                'FontSize', 16, 'FontWeight', 'bold');
    end
end

%% FIGURE 5: R² HEATMAP (FWER-masked) – hot colormap version

figure('Name','R² Heatmaps (hot, FWER masked)','Position',[50 50 2200 1000]);

for d = 1:nDV

    depVarName = Y_vars{d};
    numCh   = size(reg_results.(depVarName).R2_phase, 1);
    numFreq = length(freqs);

    % Rebuild sig masks
    sig_masks_local = false(nPred, numCh, numFreq);
    for ch = 1:min(numCh, length(reg_results.(depVarName).thresholds))
        for pm = 1:nPred
            tf = thresh_fields{pm};
            if isfield(reg_results.(depVarName).thresholds(ch), tf)
                thr = reg_results.(depVarName).thresholds(ch).(tf);
                if ~isempty(thr)
                    R2 = reg_results.(depVarName).(pred_fields{pm})(ch,:);
                    sig_masks_local(pm,ch,:) = R2 > thr(1);
                end
            end
        end
    end

    for pm = 1:nPred
        subplot(nDV, nPred, (d-1)*nPred + pm);

        R2_mat = reg_results.(depVarName).(pred_fields{pm});
        imagesc(freqs, 1:numCh, R2_mat);
        hold on;
        colormap(gca, hot); colorbar;

        % Fade non-significant regions
        mask = squeeze(sig_masks_local(pm,:,:));
        for ch = 1:numCh
            nonsig_idx = find(~mask(ch,:));
            for fi = nonsig_idx
                if fi < numFreq
                    patch([freqs(fi) freqs(fi+1) freqs(fi+1) freqs(fi)], ...
                          [ch-0.5 ch-0.5 ch+0.5 ch+0.5], ...
                          [1 1 1], 'EdgeColor','none','FaceAlpha',0.55);
                end
            end
        end

        xlabel('Frequency (Hz)'); ylabel('Channel');
        title([Y_labels{d} ': ' pred_labels{pm} ' R²']);
        set(gca,'YDir','normal');
    end

end

sgtitle('R² Heatmaps (FWER-corrected, non-sig faded)',...
        'FontSize', 14, 'FontWeight', 'bold');

%% FIGURE 6: PREFERRED PHASE (reuses phi_pref and R_phase from reg_results)

figure('Name','Preferred Phase - All DVs','Position',[200 100 900 800]);

nRows = 2; nCols = 2;

for d = 1:nDV

    depVarName   = Y_vars{d};
    phi_pref     = reg_results.(depVarName).phi_pref;   % channels x freqs
    numCh_result = size(phi_pref, 1);

    % Extract per-channel phase p-values
    p_phase_ch = ones(numCh_result, 1);
    for ch = 1:min(length(reg_results.(depVarName).stats), numCh_result)
        if isfield(reg_results.(depVarName).stats(ch), 'p_phase') && ...
                ~isempty(reg_results.(depVarName).stats(ch).p_phase)
            p_phase_ch(ch) = reg_results.(depVarName).stats(ch).p_phase(1);
        end
    end

    sig_ch_phase = find(p_phase_ch < p_thresh);

    % Sig frequency mask for phase (per channel)
    numFreq = length(freqs);
    sig_freq_phase = false(numCh_result, numFreq);
    for ch = 1:min(numCh_result, length(reg_results.(depVarName).thresholds))
        if isfield(reg_results.(depVarName).thresholds(ch), 'thresh_phase')
            thr = reg_results.(depVarName).thresholds(ch).thresh_phase;
            if ~isempty(thr)
                sig_freq_phase(ch,:) = reg_results.(depVarName).R2_phase(ch,:) > thr(1);
            end
        end
    end

    % Create polar subplot
    sp = subplot(nRows, nCols, d);
    ax = polaraxes('Position', sp.Position);
    delete(sp);
    hold(ax, 'on');

    if isempty(sig_ch_phase)
        title(ax, [Y_labels{d} ' (No significant channels)']);
        continue
    end

    % Circular mean across significant frequencies per channel
    phi_ch_mean = NaN(length(sig_ch_phase), 1);

    for idx = 1:length(sig_ch_phase)
        ch = sig_ch_phase(idx);
        sig_fq = sig_freq_phase(ch,:);
        if any(sig_fq)
            phi_ch_mean(idx) = angle(mean(exp(1i * phi_pref(ch, sig_fq))));
        end
    end

    valid = ~isnan(phi_ch_mean);
    phi_ch_mean  = phi_ch_mean(valid);
    sig_ch_valid = sig_ch_phase(valid);

    mean_vector = mean(exp(1i * phi_ch_mean));
    phi_global  = angle(mean_vector);
    R_circ      = abs(mean_vector);

    % Individual channels
    polarscatter(ax, phi_ch_mean, ones(size(phi_ch_mean)), ...
        50, colors_pred(1,:), 'filled', 'MarkerFaceAlpha', 0.4);

    % Mean vector
    polarplot(ax, [phi_global phi_global], [0 R_circ], ...
        'Color', colors_pred(1,:), 'LineWidth', 4);

    ax.ThetaZeroLocation = 'top';
    ax.ThetaDir = 'clockwise';
    rlim(ax, [0 1]);

    title(ax, [Y_labels{d} ...
        ' (n=' num2str(length(sig_ch_valid)) ...
        ', R=' num2str(R_circ, 2) ')'], 'FontWeight', 'bold');

end

sgtitle('Preferred Phase Across Dependent Variables',...
        'FontSize', 16, 'FontWeight', 'bold');

%% FIGURE 7: PREFERRED PHASE vs FREQUENCY HEATMAP (significant bins only)

figure('Name','Preferred Phase vs Frequency','Position',[200 100 1200 900]);

% HSV colormap with NaN shown as white
hsv_cmap = hsv(256);

for d = 1:nDV

    depVarName = Y_vars{d};
    phi_pref_mat = reg_results.(depVarName).phi_pref;   % channels x freqs
    numCh_result = size(phi_pref_mat, 1);
    numFreq      = length(freqs);

    % Build phase significance mask (same logic as Figures 1-2)
    sig_phase_mask = false(numCh_result, numFreq);
    for ch = 1:min(numCh_result, length(reg_results.(depVarName).thresholds))
        if isfield(reg_results.(depVarName).thresholds(ch), 'thresh_phase')
            thr = reg_results.(depVarName).thresholds(ch).thresh_phase;
            if ~isempty(thr)
                sig_phase_mask(ch,:) = reg_results.(depVarName).R2_phase(ch,:) > thr(1);
            end
        end
    end

    % Mask non-significant bins with NaN
    phi_masked = phi_pref_mat;
    phi_masked(~sig_phase_mask) = NaN;

    subplot(2, 2, d);
    imagesc(freqs, 1:numCh_result, phi_masked, 'AlphaData', ~isnan(phi_masked));
    set(gca, 'Color', [1 1 1]);  % white background for NaN
    colormap(gca, hsv_cmap);
    clim([-pi pi]);
    cb = colorbar;
    cb.Ticks = [-pi -pi/2 0 pi/2 pi];
    cb.TickLabels = {'-\pi','-\pi/2','0','\pi/2','\pi'};
    xlabel('Frequency (Hz)'); ylabel('Channel');
    title([Y_labels{d} ': Preferred Phase']);
    set(gca, 'YDir', 'normal');
    ylim([1 numCh_result]); yticks(1:8:numCh_result);

end

sgtitle('Preferred Phase vs Frequency (significant bins only)',...
        'FontSize', 16, 'FontWeight', 'bold');

%% FIGURE 8: PROPORTION OF SIGNIFICANT CHANNELS PER FREQUENCY

figure('Name','Proportion Significant Channels','Position',[200 100 1400 900]);

for d = 1:nDV

    depVarName = Y_vars{d};
    numCh   = size(reg_results.(depVarName).R2_phase, 1);
    numFreq = length(freqs);

    subplot(2, 2, d);
    hold on;

    for pm = 1:nPred
        % Build significance mask for this predictor
        sig_mask_pm = false(numCh, numFreq);
        for ch = 1:min(numCh, length(reg_results.(depVarName).thresholds))
            tf = thresh_fields{pm};
            if isfield(reg_results.(depVarName).thresholds(ch), tf)
                thr = reg_results.(depVarName).thresholds(ch).(tf);
                if ~isempty(thr)
                    R2 = reg_results.(depVarName).(pred_fields{pm})(ch,:);
                    sig_mask_pm(ch,:) = R2 > thr(1);
                end
            end
        end

        prop_sig = sum(sig_mask_pm, 1) / numCh;
        plot(freqs, prop_sig, 'Color', colors_pred(pm,:), 'LineWidth', 2);
    end

    xlabel('Frequency (Hz)'); ylabel('Proportion significant');
    title(Y_labels{d});
    ylim([0 1]);
    legend(pred_labels, 'Location', 'best', 'FontSize', 8);
    grid on;

end

sgtitle('Proportion of Significant Channels per Frequency (FWER corrected)',...
        'FontSize', 16, 'FontWeight', 'bold');

%% FIGURE 9: MONKEY-AVERAGE R² WITH THRESHOLD

animals = {'hermes', 'klecks'};
nAnimals = numel(animals);
animal_colors = lines(nAnimals);

figure('Name','Monkey-Average R²','Position',[50 50 1800 1000]);

for d = 1:nDV

    depVarName = Y_vars{d};
    monkey_file = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/results_combined/multi_lin_reg/abs_per_chan/cp10_till_100', ...
        depVarName, 'monkey_avg_results.mat');

    if ~isfile(monkey_file)
        for pm = 1:nPred
            subplot(nDV, nPred, (d-1)*nPred + pm);
            title([Y_labels{d} ': ' pred_labels{pm} ' (no data)']);
        end
        continue
    end

    mk = load(monkey_file);

    for pm = 1:nPred
        subplot(nDV, nPred, (d-1)*nPred + pm);
        hold on;

        af = avg_fields{pm};
        avg_R2  = mk.monkey_avg_obs.(af);
        avg_thr = mk.thresh_monkey.(af);

        % Find significant frequencies
        sig_idx = avg_R2 > avg_thr;

        % Plot R² curve
        plot(freqs, avg_R2, 'Color', colors_pred(pm,:), 'LineWidth', 2.5);

        % Threshold line
        yline(avg_thr, 'k--', 'LineWidth', 2);

        % Shade significant regions
        if any(sig_idx)
            y_fill_bot = zeros(1, length(freqs));
            y_fill_top = avg_R2;
            y_fill_top(~sig_idx) = 0;
            fill([freqs fliplr(freqs)], [y_fill_top fliplr(y_fill_bot)], ...
                colors_pred(pm,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        end

        xlabel('Frequency (Hz)'); ylabel('R² (monkey avg)');
        title([Y_labels{d} ': ' pred_labels{pm}]);
        grid on;
    end

end

sgtitle('Monkey-Average R² with FWER Threshold (shaded = significant)',...
        'FontSize', 16, 'FontWeight', 'bold');

%% FIGURE 10: PER-ANIMAL + MONKEY-AVERAGE OVERLAY

figure('Name','Per-Animal + Monkey-Average R²','Position',[50 50 1800 1000]);

for d = 1:nDV

    depVarName = Y_vars{d};
    monkey_file = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/results_combined/multi_lin_reg/abs_per_chan/cp10_till_100', ...
        depVarName, 'monkey_avg_results.mat');

    if ~isfile(monkey_file)
        for pm = 1:nPred
            subplot(nDV, nPred, (d-1)*nPred + pm);
            title([Y_labels{d} ': ' pred_labels{pm} ' (no data)']);
        end
        continue
    end

    mk = load(monkey_file);

    for pm = 1:nPred
        subplot(nDV, nPred, (d-1)*nPred + pm);
        hold on;

        af = avg_fields{pm};

        % Per-animal curves
        for a = 1:nAnimals
            plot(freqs, mk.obs_monkey.(af)(a,:), 'Color', animal_colors(a,:), 'LineWidth', 1.5);
        end

        % Monkey average
        plot(freqs, mk.monkey_avg_obs.(af), 'k', 'LineWidth', 2.5);

        % Threshold
        yline(mk.thresh_monkey.(af), '--r', 'LineWidth', 1.5);

        xlabel('Frequency (Hz)'); ylabel('R²');
        title([Y_labels{d} ': ' pred_labels{pm}]);
        grid on;

        if d == 1 && pm == 1
            legend([animals, {'Monkey avg', 'Threshold'}], 'Location', 'best', 'FontSize', 7);
        end
    end

end

sgtitle('Per-Animal + Monkey-Average R² with Threshold',...
        'FontSize', 16, 'FontWeight', 'bold');
