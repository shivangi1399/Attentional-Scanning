% =====================================================================
% Figures: regression F-statistic (H1 — complex/)
%
% Loads permutation-corrected F-statistic results from
% regress_stats_Fstat.m and produces visualisations of the
% significance of the sin(φ) + cos(φ) phase regressors across
% frequencies and channels.
%
% See sampling_compare/README.md for the Way-1 / Way-2 framing.
% =====================================================================
clearvars;
close all;
clc

%% Load data

data_folder    = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks/multi_lin_reg/cp10_till_100';
results_folder = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks/multi_lin_reg/complex/cp10_till_100';
cd(results_folder)
load('multi_regression_perm_maxstat.mat', 'reg_results');
cd(data_folder)
load('ph_all_sess.mat') % loads ph_comb

ph_all  = ph_comb.phase_all; % trials x freqs x channels
[numTrials, numFreq, numCh] = size(ph_all);

% Load frequency vector 
cd('/mnt/hpc/projects/MWSampling/4Shivangi/results_hermes/hermes_20170824_attentional-sampling_1/Phase_analysis/hit_miss/100iter_cut@cp_m10/1')
load('freqpow.mat')
freqs = freqpow.freq;

% Dependent variables
Y_vars = {'RT','MUA_ERP_ampl_all','LFP_ERP_ampl_all','hit_miss'};
Y_labels = {'RT', 'MUA ERP', 'LFP ERP', 'Hit/Miss'};

p_thresh = 0.05;

%% Color schemes

% DV colors for dv
colors_dv = [
    0.0 0.6 0.6;  % Teal
];

% Define colors for predictors
colors_pred = [
    0 0.4470 0.7410;   % Phase: blue
    0.0 0.6 0.6;       % Amp: Teal
    0.6 0.3 0.6;       % Amp+Phase: purple
    0.85 0.33 0.1      % Full model: orange
];

% Custom colormaps
pink_cmap = [linspace(1,0.8,64)', linspace(1,0.3,64)', linspace(1,0.6,64)'];

%% FIGURE 1: OVERVIEW HEATMAPS FOR ALL DVs

figure('Name','Significant Frequencies per Channel','Position',[50 50 2200 1000]);

for d = 1:length(Y_vars)

    depVarName = Y_vars{d};

    % Extract observed stats
    obs_F_phase    = reg_results.(depVarName).obs_F_phase;
    obs_F_Amp      = reg_results.(depVarName).obs_F_Amp;
    obs_F_AmpPhase = reg_results.(depVarName).obs_F_AmpPhase;
    obs_F_any      = reg_results.(depVarName).obs_F_any;
    R_phase        = reg_results.(depVarName).R_phase;

    numCh = size(obs_F_phase,1);
    numFreq = length(freqs);

    % Initialize significance matrices
    sig_freq_phase    = false(numCh, numFreq);
    sig_freq_Amp      = false(numCh, numFreq);
    sig_freq_AmpPhase = false(numCh, numFreq);
    sig_freq_any      = false(numCh, numFreq);

    % ------------------------------------
    % Apply permutation thresholds
    % ------------------------------------
    for ch = 1:numCh

        if ch > length(reg_results.(depVarName).FWER_thresholds)
            continue
        end

        % ----- PHASE -----
        if isfield(reg_results.(depVarName).FWER_thresholds(ch),'thresh_phase')
            tmp = reg_results.(depVarName).FWER_thresholds(ch).thresh_phase;
            if ~isempty(tmp)
                sig_freq_phase(ch,:) = ...
                    obs_F_phase(ch,:) > tmp(1);
            end
        end

        % ----- AMP -----
        if isfield(reg_results.(depVarName).FWER_thresholds(ch),'thresh_Amp')
            tmp = reg_results.(depVarName).FWER_thresholds(ch).thresh_Amp;
            if ~isempty(tmp)
                sig_freq_Amp(ch,:) = ...
                    obs_F_Amp(ch,:) > tmp(1);
            end
        end

        % ----- AMP+PHASE -----
        if isfield(reg_results.(depVarName).FWER_thresholds(ch),'thresh_AmpPhase')
            tmp = reg_results.(depVarName).FWER_thresholds(ch).thresh_AmpPhase;
            if ~isempty(tmp)
                sig_freq_AmpPhase(ch,:) = ...
                    obs_F_AmpPhase(ch,:) > tmp(1);
            end
        end

        % ----- FULL MODEL (ANY vs NULL) -----
        if isfield(reg_results.(depVarName).FWER_thresholds(ch),'thresh_any')
            tmp = reg_results.(depVarName).FWER_thresholds(ch).thresh_any;
            if ~isempty(tmp)
                sig_freq_any(ch,:) = ...
                    obs_F_any(ch,:) > tmp(1);
            end
        end

    end

    % ------------------------------------
    % Plotting
    % ------------------------------------

    % Phase strength (R)
    subplot(length(Y_vars),5,(d-1)*5 + 1);
    imagesc(freqs, 1:numCh, R_phase);
    colormap(gca, pink_cmap); colorbar;
    xlabel('Frequency (Hz)');
    ylabel('Channel');
    title([Y_labels{d} ': Phase Strength (R)']);
    set(gca,'YDir','normal');
    ylim([1 numCh]); yticks(1:8:numCh);

    % Phase significance
    subplot(length(Y_vars),5,(d-1)*5 + 2);
    imagesc(freqs, 1:numCh, sig_freq_phase);
    colormap(gca,[0.9 0.9 0.9; colors_pred(1,:)]);
    clim([0 1]);
    xlabel('Frequency (Hz)'); ylabel('Channel');
    title([Y_labels{d} ': Phase Sig']);
    set(gca,'YDir','normal');
    ylim([1 numCh]); yticks(1:8:numCh);

    % Amp significance
    subplot(length(Y_vars),5,(d-1)*5 + 3);
    imagesc(freqs, 1:numCh, sig_freq_Amp);
    colormap(gca,[0.9 0.9 0.9; colors_pred(2,:)]);
    clim([0 1]);
    xlabel('Frequency (Hz)'); ylabel('Channel');
    title([Y_labels{d} ': Amp Sig']);
    set(gca,'YDir','normal');
    ylim([1 numCh]); yticks(1:8:numCh);

    % Amp+Phase significance
    subplot(length(Y_vars),5,(d-1)*5 + 4);
    imagesc(freqs, 1:numCh, sig_freq_AmpPhase);
    colormap(gca,[0.9 0.9 0.9; colors_pred(3,:)]);
    clim([0 1]);
    xlabel('Frequency (Hz)'); ylabel('Channel');
    title([Y_labels{d} ': Amp+Phase Sig']);
    set(gca,'YDir','normal');
    ylim([1 numCh]); yticks(1:8:numCh);

    % Full model significance
    subplot(length(Y_vars),5,(d-1)*5 + 5);
    imagesc(freqs, 1:numCh, sig_freq_any);
    colormap(gca,[0.9 0.9 0.9; colors_pred(4,:)]);
    clim([0 1]);
    xlabel('Frequency (Hz)'); ylabel('Channel');
    title([Y_labels{d} ': Full Model Sig']);
    set(gca,'YDir','normal');
    ylim([1 numCh]); yticks(1:8:numCh);

end

sgtitle('Significant Frequencies per Channel (FWER Corrected across freq)','FontSize',16,'FontWeight','bold');

%% FIGURE 2: F-STATISTICS ACROSS FREQUENCIES

figure('Name','F-Statistics & Full Model','Position',[50 50 2400 1200]);

for d = 1:length(Y_vars)

    depVarName = Y_vars{d};

    % Extract observed statistics
    obs_F_phase    = reg_results.(depVarName).obs_F_phase;
    obs_F_Amp      = reg_results.(depVarName).obs_F_Amp;
    obs_F_AmpPhase = reg_results.(depVarName).obs_F_AmpPhase;
    obs_F_any      = reg_results.(depVarName).obs_F_any;

    numCh = size(obs_F_phase,1);

    % Preallocate
    p_phase = ones(numCh,1);
    p_Amp = ones(numCh,1);
    p_AmpPhase = ones(numCh,1);
    p_any = ones(numCh,1);

    thresh_phase = zeros(numCh,1);
    thresh_Amp = zeros(numCh,1);
    thresh_AmpPhase = zeros(numCh,1);
    thresh_any = zeros(numCh,1);

    % -----------------------------
    % Extract thresholds + p-values
    % -----------------------------
    for ch = 1:numCh

        if ch > length(reg_results.(depVarName).pvals)
            continue
        end

        % ----- Thresholds -----
        if isfield(reg_results.(depVarName).FWER_thresholds(ch),'thresh_phase')
            tmp = reg_results.(depVarName).FWER_thresholds(ch).thresh_phase;
            if ~isempty(tmp), thresh_phase(ch) = tmp(1); end
        end

        if isfield(reg_results.(depVarName).FWER_thresholds(ch),'thresh_Amp')
            tmp = reg_results.(depVarName).FWER_thresholds(ch).thresh_Amp;
            if ~isempty(tmp), thresh_Amp(ch) = tmp(1); end
        end

        if isfield(reg_results.(depVarName).FWER_thresholds(ch),'thresh_AmpPhase')
            tmp = reg_results.(depVarName).FWER_thresholds(ch).thresh_AmpPhase;
            if ~isempty(tmp), thresh_AmpPhase(ch) = tmp(1); end
        end

        if isfield(reg_results.(depVarName).FWER_thresholds(ch),'thresh_any')
            tmp = reg_results.(depVarName).FWER_thresholds(ch).thresh_any;
            if ~isempty(tmp), thresh_any(ch) = tmp(1); end
        end

        % ----- P-values -----
        if isfield(reg_results.(depVarName).pvals(ch),'p_phase')
            tmp = reg_results.(depVarName).pvals(ch).p_phase;
            if ~isempty(tmp), p_phase(ch) = tmp(1); end
        end

        if isfield(reg_results.(depVarName).pvals(ch),'p_Amp')
            tmp = reg_results.(depVarName).pvals(ch).p_Amp;
            if ~isempty(tmp), p_Amp(ch) = tmp(1); end
        end

        if isfield(reg_results.(depVarName).pvals(ch),'p_AmpPhase')
            tmp = reg_results.(depVarName).pvals(ch).p_AmpPhase;
            if ~isempty(tmp), p_AmpPhase(ch) = tmp(1); end
        end

        if isfield(reg_results.(depVarName).pvals(ch),'p_any')
            tmp = reg_results.(depVarName).pvals(ch).p_any;
            if ~isempty(tmp), p_any(ch) = tmp(1); end
        end
    end

    % Significant channels
    sig_ch_phase    = find(p_phase < p_thresh);
    sig_ch_Amp      = find(p_Amp < p_thresh);
    sig_ch_AmpPhase = find(p_AmpPhase < p_thresh);
    sig_ch_any      = find(p_any < p_thresh);

    % ---- PHASE ----
    subplot(length(Y_vars),4,(d-1)*4 + 1); hold on;
    for ch = sig_ch_phase'
        plot(freqs, obs_F_phase(ch,:), 'Color',[colors_pred(1,:) 0.3],'LineWidth',1);
    end
    if ~isempty(sig_ch_phase)
        plot(freqs, mean(obs_F_phase(sig_ch_phase,:),1), ...
            'Color', colors_pred(1,:), 'LineWidth',3);
        plot([freqs(1) freqs(end)], ...
            [mean(thresh_phase(sig_ch_phase)) mean(thresh_phase(sig_ch_phase))], ...
            'k--','LineWidth',2);
    end
    xlabel('Frequency (Hz)'); ylabel('F-stat'); grid on;
    title([Y_labels{d} ': Phase (n=' num2str(length(sig_ch_phase)) ')']);

    % ---- AMP ----
    subplot(length(Y_vars),4,(d-1)*4 + 2); hold on;
    for ch = sig_ch_Amp'
        plot(freqs, obs_F_Amp(ch,:), 'Color',[colors_pred(2,:) 0.3],'LineWidth',1);
    end
    if ~isempty(sig_ch_Amp)
        plot(freqs, mean(obs_F_Amp(sig_ch_Amp,:),1), ...
            'Color', colors_pred(2,:), 'LineWidth',3);
        plot([freqs(1) freqs(end)], ...
            [mean(thresh_Amp(sig_ch_Amp)) mean(thresh_Amp(sig_ch_Amp))], ...
            'k--','LineWidth',2);
    end
    xlabel('Frequency (Hz)'); ylabel('F-stat'); grid on;
    title([Y_labels{d} ': Amp (n=' num2str(length(sig_ch_Amp)) ')']);

    % ---- AMP+PHASE ----
    subplot(length(Y_vars),4,(d-1)*4 + 3); hold on;
    for ch = sig_ch_AmpPhase'
        plot(freqs, obs_F_AmpPhase(ch,:), 'Color',[colors_pred(3,:) 0.3],'LineWidth',1);
    end
    if ~isempty(sig_ch_AmpPhase)
        plot(freqs, mean(obs_F_AmpPhase(sig_ch_AmpPhase,:),1), ...
            'Color', colors_pred(3,:), 'LineWidth',3);
        plot([freqs(1) freqs(end)], ...
            [mean(thresh_AmpPhase(sig_ch_AmpPhase)) mean(thresh_AmpPhase(sig_ch_AmpPhase))], ...
            'k--','LineWidth',2);
    end
    xlabel('Frequency (Hz)'); ylabel('F-stat'); grid on;
    title([Y_labels{d} ': Amp+Phase (n=' num2str(length(sig_ch_AmpPhase)) ')']);

    % FULL MODEL 
    subplot(length(Y_vars),4,(d-1)*4 + 4); hold on;
    for ch = sig_ch_any'
        plot(freqs, obs_F_any(ch,:), 'Color',[colors_pred(4,:) 0.3],'LineWidth',1);
    end
    if ~isempty(sig_ch_any)
        plot(freqs, mean(obs_F_any(sig_ch_any,:),1), ...
            'Color', colors_pred(4,:), 'LineWidth',3);
        plot([freqs(1) freqs(end)], ...
            [mean(thresh_any(sig_ch_any)) mean(thresh_any(sig_ch_any))], ...
            'k--','LineWidth',2);
    end
    xlabel('Frequency (Hz)'); ylabel('F-stat'); grid on;
    title([Y_labels{d} ': Full Model (n=' num2str(length(sig_ch_any)) ')']);

end

sgtitle('F-Statistics Across Frequency Spectrum','FontSize',16,'FontWeight','bold');

%% FIGURE 3: F-STATISTICS ACROSS FREQUENCIES PER CHANNEL

models = { ...
    'Phase',     'obs_F_phase',    'p_phase',    'thresh_phase',    colors_pred(1,:); ...
    'Amp',       'obs_F_Amp',      'p_Amp',      'thresh_Amp',      colors_pred(2,:); ...
    'Amp+Phase', 'obs_F_AmpPhase', 'p_AmpPhase', 'thresh_AmpPhase', colors_pred(3,:); ...
    'FullModel', 'obs_F_any',      'p_any',      'thresh_any',      colors_pred(4,:) ...
    };

for d = 1:length(Y_vars)

    depVarName = Y_vars{d};
    nCh = size(reg_results.(depVarName).obs_F_phase,1);

    for m = 1:size(models,1)

        model_label  = models{m,1};
        F_field      = models{m,2};
        p_field      = models{m,3};
        thresh_field = models{m,4};
        col          = models{m,5};

        % Extract observed F values
        obs_F = reg_results.(depVarName).(F_field);

        % Preallocate
        pvals  = ones(nCh,1);
        thresh = nan(nCh,1);

        % --------------------------------------
        % Extract p-values and thresholds safely
        % --------------------------------------
        for ch = 1:nCh

            if ch <= length(reg_results.(depVarName).pvals)

                % p-values
                if isfield(reg_results.(depVarName).pvals(ch), p_field)
                    tmp = reg_results.(depVarName).pvals(ch).(p_field);
                    if ~isempty(tmp)
                        pvals(ch) = tmp(1);
                    end
                end

                % thresholds
                if isfield(reg_results.(depVarName).FWER_thresholds(ch), thresh_field)
                    tmp = reg_results.(depVarName).FWER_thresholds(ch).(thresh_field);
                    if ~isempty(tmp)
                        thresh(ch) = tmp(1);
                    end
                end
            end
        end

        % Significant channels
        sig_ch = find(pvals < p_thresh);

        if isempty(sig_ch)
            fprintf('No significant channels for %s - %s\n', ...
                Y_labels{d}, model_label);
            continue
        end

        % --------------------------------------
        % Create figure
        % --------------------------------------
        figure('Name',[Y_labels{d} ' - ' model_label], ...
               'Position',[100 100 1800 900]);

        nSig = length(sig_ch);

        % Automatic grid layout
        nCols = ceil(sqrt(nSig));
        nRows = ceil(nSig / nCols);

        % Determine common y-limit for this figure
        maxF = 0;
        for ch = sig_ch'
            maxF = max(maxF, max(obs_F(ch,:)));
            if ~isnan(thresh(ch))
                maxF = max(maxF, thresh(ch));
            end
        end
        yMax = maxF * 1.1;

        % --------------------------------------
        % Plot each significant channel
        % --------------------------------------
        for i = 1:nSig

            ch = sig_ch(i);

            subplot(nRows,nCols,i); hold on;

            % F-stat curve
            plot(freqs, obs_F(ch,:), ...
                'Color', col, ...
                'LineWidth', 2);

            % Channel-specific threshold
            if ~isnan(thresh(ch))
                yline(thresh(ch), 'k--', 'LineWidth', 1.5);
            end

            xlabel('Frequency (Hz)');
            ylabel('F-stat');
            title(sprintf('Ch %d ', ch));

            ylim([0 yMax]);
            xlim([freqs(1) freqs(end)]);
            grid on;
        end

        % Super title
        sgtitle(sprintf('%s    %s  (n = %d)', ...
            Y_labels{d}, model_label, nSig), ...
            'FontSize',16,'FontWeight','bold');

    end
end

%% FIGURE 4: Heatmap of Model Component F-stats per Channel (FWER masked)

figure('Name','Model Component F-stats Heatmap',...
       'Position',[50 50 2200 1000]);

for d = 1:length(Y_vars)

    depVarName = Y_vars{d};

    obs_F_phase    = reg_results.(depVarName).obs_F_phase;
    obs_F_Amp      = reg_results.(depVarName).obs_F_Amp;
    obs_F_AmpPhase = reg_results.(depVarName).obs_F_AmpPhase;
    obs_F_any      = reg_results.(depVarName).obs_F_any;

    numCh_result = size(obs_F_phase,1);
    numFreq = length(freqs);

    % Initialize significance masks
    sig_mask_phase    = false(numCh_result, numFreq);
    sig_mask_Amp      = false(numCh_result, numFreq);
    sig_mask_AmpPhase = false(numCh_result, numFreq);
    sig_mask_any      = false(numCh_result, numFreq);

    % ---------------------------------------------------
    % Apply stored FWER thresholds 
    % ---------------------------------------------------
    for ch = 1:numCh_result

        if ch > length(reg_results.(depVarName).FWER_thresholds)
            continue
        end

        % ----- PHASE -----
        if isfield(reg_results.(depVarName).FWER_thresholds(ch),'thresh_phase')
            tmp = reg_results.(depVarName).FWER_thresholds(ch).thresh_phase;
            if ~isempty(tmp)
                sig_mask_phase(ch,:) = obs_F_phase(ch,:) > tmp(1);
            end
        end

        % ----- AMP -----
        if isfield(reg_results.(depVarName).FWER_thresholds(ch),'thresh_Amp')
            tmp = reg_results.(depVarName).FWER_thresholds(ch).thresh_Amp;
            if ~isempty(tmp)
                sig_mask_Amp(ch,:) = obs_F_Amp(ch,:) > tmp(1);
            end
        end

        % ----- AMP+PHASE -----
        if isfield(reg_results.(depVarName).FWER_thresholds(ch),'thresh_AmpPhase')
            tmp = reg_results.(depVarName).FWER_thresholds(ch).thresh_AmpPhase;
            if ~isempty(tmp)
                sig_mask_AmpPhase(ch,:) = obs_F_AmpPhase(ch,:) > tmp(1);
            end
        end

        % ----- FULL MODEL -----
        if isfield(reg_results.(depVarName).FWER_thresholds(ch),'thresh_any')
            tmp = reg_results.(depVarName).FWER_thresholds(ch).thresh_any;
            if ~isempty(tmp)
                sig_mask_any(ch,:) = obs_F_any(ch,:) > tmp(1);
            end
        end

    end

    % ---------------------------------------------------
    % PLOTTING
    % ---------------------------------------------------

    components = {obs_F_phase, obs_F_Amp, obs_F_AmpPhase, obs_F_any};
    masks      = {sig_mask_phase, sig_mask_Amp, sig_mask_AmpPhase, sig_mask_any};
    titles     = {'Phase','Amp','Amp+Phase','Full Model'};

    for c = 1:4

        subplot(length(Y_vars),4,(d-1)*4 + c);

        imagesc(freqs, 1:numCh_result, components{c});
        hold on;
        colormap(gca, hot);
        colorbar;

        % Fade non-significant regions
        mask = masks{c};

        for ch = 1:numCh_result
            nonsig_idx = find(~mask(ch,:));
            for f = nonsig_idx
                if f < numFreq
                    patch([freqs(f) freqs(f+1) freqs(f+1) freqs(f)], ...
                          [ch-0.5 ch-0.5 ch+0.5 ch+0.5], ...
                          [1 1 1], 'EdgeColor','none','FaceAlpha',0.5);
                end
            end
        end

        xlabel('Frequency (Hz)');
        ylabel('Channel');
        title([Y_labels{d} ': ' titles{c} ' F-stat']);
        set(gca,'YDir','normal');

    end

end

sgtitle('F-Statistics Heatmaps (FWER-corrected across freq)','FontSize',14,'FontWeight','bold');

%% FIGURE 5: Preferred Phase 

figure('Name','Preferred Phase - All DVs',...
       'Position',[200 100 900 800]);

nDV = length(Y_vars);
nRows = 2;  
nCols = 2;

for d = 1:nDV

    depVarName = Y_vars{d};
    phi_pref   = reg_results.(depVarName).phi_pref;   % channels x freqs
    numCh_result = size(phi_pref,1);

    % Extract phase p-values
    p_phase = ones(numCh_result,1);

    for ch = 1:min(length(reg_results.(depVarName).pvals), numCh_result)
        if isfield(reg_results.(depVarName).pvals(ch),'p_phase') && ...
                ~isempty(reg_results.(depVarName).pvals(ch).p_phase)
            p_phase(ch) = reg_results.(depVarName).pvals(ch).p_phase(1);
        end
    end

    sig_ch_phase = find(p_phase < p_thresh);

    % Create polar subplot 
    sp = subplot(nRows,nCols,d);       % create subplot
    ax = polaraxes('Position', sp.Position);  % polar axes in same position
    delete(sp);                        % remove the original Cartesian axes
    hold(ax,'on');

    if isempty(sig_ch_phase)
        title(ax, [Y_labels{d} ' (No significant channels)']);
        continue
    end

    % Circular mean across significant frequencies per channel
    phi_ch_mean = zeros(length(sig_ch_phase),1);  % preallocate

    for idx = 1:length(sig_ch_phase)
        ch = sig_ch_phase(idx);

        % Identify significant frequencies for this channel
        sig_freq_ch = sig_freq_phase(ch,:);  % sig_freq_phase defined in Figure 1 section

        if any(sig_freq_ch)
            phi_sig_ch = phi_pref(ch, sig_freq_ch);  % take only significant freqs
            phi_ch_mean(idx) = angle(mean(exp(1i*phi_sig_ch)));  % circular mean
        else
            phi_ch_mean(idx) = NaN;  % if no significant freq, set NaN
        end
    end

    % Remove channels with NaN (no sig freq)
    phi_ch_mean = phi_ch_mean(~isnan(phi_ch_mean));
    sig_ch_phase = sig_ch_phase(~isnan(phi_ch_mean));


    % Global circular mean across channels
    mean_vector = mean(exp(1i*phi_ch_mean));
    phi_global  = angle(mean_vector);
    R           = abs(mean_vector);

    % Plot individual channels
    polarscatter(ax, phi_ch_mean, ones(size(phi_ch_mean)), ...
        50, colors_dv(1,:), ...
        'filled', ...
        'MarkerFaceAlpha', 0.4);

    % Plot mean vector
    polarplot(ax, [phi_global phi_global], [0 R], ...
        'Color', colors_dv(1,:), ...
        'LineWidth',4);

    % Formatting
    ax.ThetaZeroLocation = 'top';
    ax.ThetaDir = 'clockwise';
    rlim(ax,[0 1]);

    title(ax, [Y_labels{d} ...
          ' (n=' num2str(length(sig_ch_phase)) ...
          ', R=' num2str(R,2) ')'], ...
          'FontWeight','bold');

end

sgtitle('Preferred Phase Across Dependent Variables',...
        'FontSize',16,'FontWeight','bold');

%% FIGURE 6: Preferred Phase - DVs with Frequency Bands

figure('Name','Preferred Phase by Frequency Band - All DVs',...
       'Position',[200 100 1200 800]);

% Define frequency bands (Hz) and colors
freq_bands = [1 4; 4 8; 8 13; 13 30; 30 100]; % Delta, Theta, Alpha, Beta, Gamma
band_labels = {'Delta','Theta','Alpha','Beta','Gamma'};
band_colors = [
    148, 103, 189;  % Purple
    0, 169, 164;    % Teal
    106, 90, 205;   % Slate Blue
    255, 179, 186;  % Pastel Pink
    173, 255, 173   % Pastel Green 
] / 255; 

nDV = length(Y_vars);
nRows = 2;  
nCols = 2;

for d = 1:nDV
    depVarName = Y_vars{d};
    phi_pref   = reg_results.(depVarName).phi_pref;   % channels x freqs
    numCh_result = size(phi_pref,1);

    % Extract phase p-values
    p_phase = ones(numCh_result,1);
    for ch = 1:min(length(reg_results.(depVarName).pvals), numCh_result)
        if isfield(reg_results.(depVarName).pvals(ch),'p_phase') && ...
                ~isempty(reg_results.(depVarName).pvals(ch).p_phase)
            p_phase(ch) = reg_results.(depVarName).pvals(ch).p_phase(1);
        end
    end
    sig_ch_phase = find(p_phase < p_thresh);

    % Create polar subplot
    sp = subplot(nRows,nCols,d);
    ax = polaraxes('Position', sp.Position);
    delete(sp);
    hold(ax,'on');

    if isempty(sig_ch_phase)
        title(ax, [Y_labels{d} ' (No significant channels)']);
        continue
    end

    legend_handles = [];  % for legend entries
    legend_entries = {};

    % For each frequency band, compute circular mean across sig freqs
    for b = 1:length(band_labels)
        band_freqs = freqs >= freq_bands(b,1) & freqs < freq_bands(b,2);

        phi_ch_mean_band = zeros(length(sig_ch_phase),1); % per channel
        for idx = 1:length(sig_ch_phase)
            ch = sig_ch_phase(idx);
            sig_freq_ch = sig_freq_phase(ch,:) & band_freqs;  % only sig freqs in band
            if any(sig_freq_ch)
                phi_sig_ch = phi_pref(ch, sig_freq_ch);
                phi_ch_mean_band(idx) = angle(mean(exp(1i*phi_sig_ch)));  % circular mean
            else
                phi_ch_mean_band(idx) = NaN;  % no sig freqs in this band
            end
        end

        % Remove NaNs
        phi_ch_mean_band = phi_ch_mean_band(~isnan(phi_ch_mean_band));
        if isempty(phi_ch_mean_band)
            continue  % skip this band if no sig frequencies
        end

        % Global mean vector across channels for this band
        mean_vector = mean(exp(1i*phi_ch_mean_band));
        phi_global  = angle(mean_vector);
        R           = abs(mean_vector);

        % Plot individual channels
        polarscatter(ax, phi_ch_mean_band, ones(size(phi_ch_mean_band)), ...
            50, band_colors(b,:), 'filled', 'MarkerFaceAlpha', 0.4);

        % Plot mean vector for the band
        h = polarplot(ax, [phi_global phi_global], [0 R], ...
            'Color', band_colors(b,:), 'LineWidth', 4);

        % Collect handles for legend
        legend_handles = [legend_handles h];
        legend_entries{end+1} = band_labels{b};
    end

    % Formatting
    ax.ThetaZeroLocation = 'top';
    ax.ThetaDir = 'clockwise';
    rlim(ax,[0 1]);
    title(ax, [Y_labels{d} ' (n=' num2str(length(sig_ch_phase)) ')'], ...
          'FontWeight','bold');

    % Add legend for frequency bands
    if ~isempty(legend_handles)
        legend(ax, legend_handles, legend_entries, 'Location', 'eastoutside');
    end
end

sgtitle('Preferred Phase by Frequency Band Across Dependent Variables',...
        'FontSize',16,'FontWeight','bold');

%% FIGURE 7 : Grid of Channel-Specific model comparison - we can't directly compare that

for dv_idx = 1:length(Y_vars)
    
    depVarName = Y_vars{dv_idx};
    
    obs_F_phase    = reg_results.(depVarName).obs_F_phase;
    obs_F_Amp      = reg_results.(depVarName).obs_F_Amp;
    obs_F_AmpPhase = reg_results.(depVarName).obs_F_AmpPhase;
    obs_F_any      = reg_results.(depVarName).obs_F_any;
    
    % Extract p-values for all channels
    numCh_result = size(obs_F_phase,1);
    p_phase    = ones(numCh_result,1);
    p_Amp      = ones(numCh_result,1);
    p_AmpPhase = ones(numCh_result,1);
    
    for ch = 1:min(length(reg_results.(depVarName).pvals), numCh_result)
        if isfield(reg_results.(depVarName).pvals(ch),'p_phase') && ~isempty(reg_results.(depVarName).pvals(ch).p_phase)
            p_phase(ch) = reg_results.(depVarName).pvals(ch).p_phase;
        end
        if isfield(reg_results.(depVarName).pvals(ch),'p_Amp') && ~isempty(reg_results.(depVarName).pvals(ch).p_Amp)
            p_Amp(ch) = reg_results.(depVarName).pvals(ch).p_Amp;
        end
        if isfield(reg_results.(depVarName).pvals(ch),'p_AmpPhase') && ~isempty(reg_results.(depVarName).pvals(ch).p_AmpPhase)
            p_AmpPhase(ch) = reg_results.(depVarName).pvals(ch).p_AmpPhase;
        end
    end
    
    % Find significant channels (at least one predictor significant)
    sig_channels = find(p_phase < p_thresh | p_Amp < p_thresh | p_AmpPhase < p_thresh);
    
    if ~isempty(sig_channels)
        
        % Determine grid layout
        n_sig = length(sig_channels);
        ncols = ceil(sqrt(n_sig));
        nrows = ceil(n_sig / ncols);
        
        figure('Name',sprintf('Channel-Specific Model Comparison - %s',Y_labels{dv_idx}),...
               'Position',[50 50 400*ncols 300*nrows]);
        
        for idx = 1:n_sig
            ch = sig_channels(idx);
            
            subplot(nrows, ncols, idx);
            hold on;
            
            % Extract thresholds for this channel
            thresh_phase = 0; thresh_Amp = 0; thresh_AmpPhase = 0;
            if ch <= length(reg_results.(depVarName).FWER_thresholds)
                if ~isempty(reg_results.(depVarName).FWER_thresholds(ch).thresh_phase)
                    thresh_phase = reg_results.(depVarName).FWER_thresholds(ch).thresh_phase;
                end
                if ~isempty(reg_results.(depVarName).FWER_thresholds(ch).thresh_Amp)
                    thresh_Amp = reg_results.(depVarName).FWER_thresholds(ch).thresh_Amp;
                end
                if ~isempty(reg_results.(depVarName).FWER_thresholds(ch).thresh_AmpPhase)
                    thresh_AmpPhase = reg_results.(depVarName).FWER_thresholds(ch).thresh_AmpPhase;
                end
            end
            
            % Plot F-statistics
            h1 = plot(freqs, obs_F_phase(ch,:), 'Color', colors_pred(1,:), 'LineWidth', 1.5, 'DisplayName', 'Phase');
            h2 = plot(freqs, obs_F_Amp(ch,:), 'Color', colors_pred(2,:), 'LineWidth', 1.5, 'DisplayName', 'Amp');
            h3 = plot(freqs, obs_F_AmpPhase(ch,:), 'Color', colors_pred(3,:), 'LineWidth', 1.5, 'DisplayName', 'Amp+Phase');
            h4 = plot(freqs, obs_F_any(ch,:), 'Color', colors_pred(4,:), 'LineWidth', 1.5, 'DisplayName', 'Full Model');
            
            % Plot thresholds
            h5 = []; h6 = []; h7 = [];
            if thresh_phase > 0
                h5 = plot([freqs(1) freqs(end)], [thresh_phase thresh_phase], '--', 'Color', colors_pred(1,:), 'LineWidth', 1, 'DisplayName', 'Phase Thresh');
            end
            if thresh_Amp > 0
                h6 = plot([freqs(1) freqs(end)], [thresh_Amp thresh_Amp], '--', 'Color', colors_pred(2,:), 'LineWidth', 1, 'DisplayName', 'Amp Thresh');
            end
            if thresh_AmpPhase > 0
                h7 = plot([freqs(1) freqs(end)], [thresh_AmpPhase thresh_AmpPhase], '--', 'Color', colors_pred(3,:), 'LineWidth', 1, 'DisplayName', 'Amp+Phase Thresh');
            end
            
            xlabel('Frequency (Hz)');
            ylabel('F-statistic');
            title(sprintf('Ch %d', ch));
            grid on;
            
            if idx == 1
                % Build legend handles dynamically based on what thresholds exist
                legend_handles = [h1, h2, h3, h4];
                if ~isempty(h5), legend_handles = [legend_handles, h5]; end
                if ~isempty(h6), legend_handles = [legend_handles, h6]; end
                if ~isempty(h7), legend_handles = [legend_handles, h7]; end
                legend(legend_handles, 'Location','best','FontSize',8);
            end
        end
        
        sgtitle(sprintf('Model Comparison per Channel - %s (n=%d significant)', Y_labels{dv_idx}, n_sig),...
                'FontSize',14,'FontWeight','bold');
    else
        fprintf('No significant channels for %s\n', Y_labels{dv_idx});
    end
end

%% FIGURE 8: Model Comparison - Variance Explained

figure('Name','Model Comparison','Position',[50 50 1600 1000]);

for d = 1:length(Y_vars)
    depVarName = Y_vars{d};
    
    % You'll need to compute R² for each model during regression
    % For now, we can use F-stats as a proxy for relative improvement but
    % it's not entirely correct
    obs_F_phase    = reg_results.(depVarName).obs_F_phase;
    obs_F_Amp      = reg_results.(depVarName).obs_F_Amp;
    obs_F_AmpPhase = reg_results.(depVarName).obs_F_AmpPhase;
    obs_F_any      = reg_results.(depVarName).obs_F_any;  % full vs null
    
    % Get significant channels
    p_phase = ones(size(obs_F_phase,1),1);
    for ch = 1:length(reg_results.(depVarName).pvals)
        if ~isempty(reg_results.(depVarName).pvals(ch).p_phase)
            p_phase(ch) = reg_results.(depVarName).pvals(ch).p_phase;
        end
    end
    sig_ch = find(p_phase < p_thresh);
    
    if isempty(sig_ch), continue; end
    
    % Plot for significant channels
    subplot(2,2,d); hold on;
    
    % Phase-only contribution (relative to full model)
    plot(freqs, mean(obs_F_phase(sig_ch,:),1), ...
        'Color', colors_pred(1,:), 'LineWidth', 2.5, 'DisplayName', 'Phase');
    
    % Amplitude-only contribution
    plot(freqs, mean(obs_F_Amp(sig_ch,:),1), ...
        'Color', colors_pred(2,:), 'LineWidth', 2.5, 'DisplayName', 'Amplitude');
    
    % Combined Amp+Phase contribution
    plot(freqs, mean(obs_F_AmpPhase(sig_ch,:),1), ...
        'Color', colors_pred(3,:), 'LineWidth', 2.5, 'DisplayName', 'Amp+Phase');
    
    % Full model (vs null)
    plot(freqs, mean(obs_F_any(sig_ch,:),1), ...
        'k--', 'LineWidth', 2.5, 'DisplayName', 'Full Model');
    
    xlabel('Frequency (Hz)'); 
    ylabel('F-statistic (contribution)');
    title([Y_labels{d} ' - Model Components']);
    legend('Location','best');
    grid on;
end

sgtitle('Relative Contribution of Model Components','FontSize',16,'FontWeight','bold');

