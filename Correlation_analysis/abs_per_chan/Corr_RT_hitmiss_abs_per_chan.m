% =====================================================================
% Circular-linear correlation: pre-stimulus phase vs. reaction time
% AND phase-opposition / inverted-miss-ITC for hit/miss
% Hypothesis H1+H4 (abs_per_chan/)
%
% Claim: trials within a channel share a preferred phase (H1 at trial
% level), but channels are NOT required to share a preferred phase.
%
% Recipe: per-channel computation on all trials together
%   RT:  circ_corrcl(phase, RT)                    (non-negative magnitude)
%   POS: ITC_hits + ITC_misses                     (non-negative magnitude)
%   ITC: abs(mean(exp(i·phase_inverted)))          (non-negative magnitude)
% Way 2 across channels (arithmetic mean of channel magnitudes); Way 2
% across animals.
%
% See sampling_compare/README.md for the Way-1 / Way-2 framing.
% =====================================================================
clear all
close all
clc

%% Settings
animals = {'hermes', 'klecks'};
permut_n = 1000;
nCh = 64;

%% Specify paths

addpath /opt/fieldtrip_github/
ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/Correlation_analysis
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/Correlation_analysis/functions
addpath /mnt/hpc/projects/MWSampling/4Shivangi
addpath /mnt/hpc/projects/MWSampling/4Shivangi/software_folder/CircStat2012a
clc

for a = 1:numel(animals)
    animalName = animals{a};
    fprintf('\n=== Processing %s ===\n', animalName);

    data_folder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName], 'multi_lin_reg', 'cp10_till_100');
    corr_folder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName], 'phase_correlation', 'abs_per_chan', 'cp10_till_100');
    if ~exist(corr_folder,'dir'), mkdir(corr_folder); end

    %% =====================================================================
    %  SECTION A: Phase vs RT correlation (all locations and difficulty levels)
    %  =====================================================================

    cd(data_folder)
    load('ph_all_sess.mat')

    % Select hit trials only (misses have NaN RT)
    hit_idx = find(ph_comb.RT_trialinfo(:,20) == 1);

    output_RT = fullfile(corr_folder, 'RT');

    %% A1. Real data — circ_corrcl(phase, RT)

    for ichan = 1:nCh
        ichan
        phase = ph_comb.phase_all(hit_idx, :, ichan);  % [nHits x nFreq]
        rt    = ph_comb.RT(hit_idx, ichan);             % [nHits x 1]

        % Remove NaN RT (safety)
        valid = ~isnan(rt);
        phase_clean = phase(valid, :);
        rt_clean    = rt(valid);

        % Skip channel if no valid RT values
        if isempty(rt_clean)
            fprintf('Channel %d skipped (no valid RT)\n', ichan);
            continue
        end

        nFreq = size(phase_clean, 2);
        correlation = nan(1, nFreq);
        pvalue      = nan(1, nFreq);

        for foi = 1:nFreq
            [correlation(1, foi), pvalue(1, foi)] = circ_corrcl(phase_clean(:, foi), rt_clean);
        end

        chan_folder = fullfile(output_RT, 'all_loc_difflev', num2str(ichan));
        if ~exist(chan_folder, 'dir'), mkdir(chan_folder); end
        cd(chan_folder)

        save correlation correlation
        save pvalue pvalue
    end

    %% A2. Permutation (SLURM)

    cd(data_folder)
    load('ph_all_sess.mat')

    hit_idx  = find(ph_comb.RT_trialinfo(:,20) == 1);
    nTrials  = length(hit_idx);

    rng(2025)
    perm_indices = arrayfun(@(x) randperm(nTrials), 1:permut_n, 'UniformOutput', false);

    cfg = cell(1, nCh);
    for ichan = 1:nCh
        cfg{ichan}.ichan        = ichan;
        cfg{ichan}.permut_n     = permut_n;
        cfg{ichan}.infile       = data_folder;
        cfg{ichan}.outfile      = fullfile(output_RT, 'all_loc_difflev');
        cfg{ichan}.perm_indices = perm_indices;
        cfg{ichan}.trial_idx    = hit_idx;
    end

    slurmfun(@circlin_corr_RT_perm, cfg, ...
        'partition',   '8GB', ...
        'stopOnError', false, ...
        'useUserPath', true);

    %% =====================================================================
    %  SECTION B: Phase vs Hit/Miss correlation (all locations and difficulty levels)
    %  =====================================================================

    cd(data_folder)
    load('ph_all_sess.mat')

    % All trials (hits + misses)
    all_idx    = find(ph_comb.trialinfo(:,20) == 1 | ph_comb.trialinfo(:,20) == 5);
    hit_labels = (ph_comb.trialinfo(all_idx, 20) == 1);  % logical

    output_HM = fullfile(corr_folder, 'hit_miss');

    %% B1. Real data — POS + ITC with inverted miss phases

    for ichan = 1:nCh
        ichan

        phase = ph_comb.phase_all(all_idx, :, ichan);  % [nTrials x nFreq]
        nFreq = size(phase, 2);

        pos = nan(1, nFreq);  % Phase Opposition Sum
        itc = nan(1, nFreq);  % ITC with inverted miss phases

        hit_idx  = hit_labels;
        miss_idx = ~hit_labels;

        % POS: ITC_hits + ITC_misses (VanRullen 2016) - more sensitive than PBI
        % Tests whether hits and misses each have distinct preferred phases
        for foi = 1:nFreq
            itc_h = abs(mean(exp(1i * phase(hit_idx, foi))));
            itc_m = abs(mean(exp(1i * phase(miss_idx, foi))));
            pos(1, foi) = itc_h + itc_m;
        end

        % ITC with inverted miss phases
        % Tests specifically whether hits and misses are at opposite phases
        phase_inv = phase;
        % + pi only. mod(x+pi,2*pi)-pi is the WRAP-into-[-pi,pi) idiom,
        % and the phases are already wrapped, so it added pi and took it
        % straight back off - the inversion was a no-op and this ITC was
        % plain ITC over hits and misses pooled, blind to the labels.
        % exp(1i*.) does not care that the result leaves [-pi,pi).
        phase_inv(miss_idx, :) = phase(miss_idx, :) + pi;

        for foi = 1:nFreq
            itc(1, foi) = abs(mean(exp(1i * phase_inv(:, foi))));
        end

        chan_folder = fullfile(output_HM, 'all_loc_difflev', num2str(ichan));
        if ~exist(chan_folder, 'dir'), mkdir(chan_folder); end
        cd(chan_folder)

        save pos pos
        save itc itc
    end

    %% B2. Permutation (SLURM)

    cd(data_folder)
    load('ph_all_sess.mat')

    all_idx    = find(ph_comb.trialinfo(:,20) == 1 | ph_comb.trialinfo(:,20) == 5);
    hit_labels = (ph_comb.trialinfo(all_idx, 20) == 1);
    nTrials    = length(all_idx);

    rng(2025)
    perm_indices = arrayfun(@(x) randperm(nTrials), 1:permut_n, 'UniformOutput', false);

    cfg = cell(1, nCh);
    for ichan = 1:nCh
        cfg{ichan}.ichan        = ichan;
        cfg{ichan}.permut_n     = permut_n;
        cfg{ichan}.infile       = data_folder;
        cfg{ichan}.outfile      = fullfile(output_HM, 'all_loc_difflev');
        cfg{ichan}.perm_indices = perm_indices;
        cfg{ichan}.trial_idx    = all_idx;
        cfg{ichan}.hit_labels   = hit_labels;
    end

    slurmfun(@circlin_corr_hitmiss_perm, cfg, ...
        'partition',   '8GB', ...
        'stopOnError', false, ...
        'useUserPath', true);

    %% =====================================================================
    %  PLOTTING — Phase vs RT correlation (per animal)
    %  =====================================================================

    save_root_RT = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/phase_correlation/abs_per_chan', animalName, 'cp10_till_100', 'RT', 'all_loc_difflev');
    if ~exist(save_root_RT, 'dir'), mkdir(save_root_RT); end

    cd(data_folder)
    load('frequency.mat')
    freq = frequency;
    nFreq = numel(freq);

    limit_max = nan(1, nCh);
    Corr_chan  = false(nCh, nFreq);

    f1 = figure(1);
    set(f1, 'Units', 'centimeters', 'Position', [1 1 50 40]);
    set(f1, 'PaperUnits', 'centimeters', 'PaperSize', [50 40], 'PaperPosition', [0 0 50 40]);
    for ch = 1:nCh
        ch_folder = fullfile(output_RT, 'all_loc_difflev', num2str(ch));
        if ~exist(ch_folder, 'dir')
            warning('Skipping channel %d (folder missing)', ch);
            continue
        end
        cd(ch_folder);

        if ~exist('correlation.mat','file') || ~exist('corr_perm.mat','file')
            warning('Skipping channel %d (missing files)', ch);
            continue
        end

        load correlation
        load corr_perm

        if any(isnan(correlation)) || any(isnan(corr_perm(:)))
            warning('Skipping channel %d (NaN values)', ch);
            continue
        end

        tmax = max(corr_perm, [], 2);
        limit_max(ch) = quantile(tmax, 0.95);

        if isnan(limit_max(ch)), continue; end

        subplot(8, 8, ch);
        plot_sigfreq(freq, correlation, limit_max(ch));
        title(['Ch ' num2str(ch)])

        Corr_chan(ch,:) = correlation >= limit_max(ch);
    end

    print(f1, fullfile(save_root_RT, 'all_channels_RT_corr.pdf'), '-dpdf');

    % Summary heatmap
    figure;
    imagesc(freq, 1:nCh, Corr_chan);
    set(gca, 'YDir', 'normal');
    xlabel('Frequency (Hz)'); ylabel('Channels');
    title(sprintf('%s — Significant Phase-RT Correlation per Channel', animalName));
    caxis([0 1]); colorbar;
    saveas(gcf, fullfile(save_root_RT, 'summary_RT_corr.pdf'));

    % Combined across channels
    valid_idx = ~isnan(limit_max);
    corr_all = [];
    corr_perm_all = [];

    for ch = find(valid_idx)
        cd(fullfile(output_RT, 'all_loc_difflev', num2str(ch)));
        load correlation
        load corr_perm
        corr_all = [corr_all; correlation];
        corr_perm_all = cat(3, corr_perm_all, corr_perm);
    end

    if ~isempty(corr_all)
        corr_avg = nanmean(corr_all, 1);
        corr_perm_avg = nanmean(corr_perm_all, 3);
        tmax_all = nanmax(corr_perm_avg, [], 2);
        limit_avg = quantile(tmax_all, 0.95);

        figure;
        plot_sigfreq(freq, corr_avg, limit_avg);
        ylim([0 0.1]); title(sprintf('%s — All Channels Combined - Phase vs RT Correlation', animalName));
        xlabel('Frequency (Hz)'); ylabel('Correlation');
        saveas(gcf, fullfile(save_root_RT, 'combined_RT_corr.pdf'));
    end

    %% Channel-average null distribution — RT correlation

    fprintf('Computing RT correlation channel-average null for %s...\n', animalName);

    corr_all_ch = NaN(nCh, nFreq);
    corr_perm_all_ch = [];

    for ch = 1:nCh
        ch_folder = fullfile(output_RT, 'all_loc_difflev', num2str(ch));
        if ~exist(ch_folder, 'dir'), continue; end

        corr_file = fullfile(ch_folder, 'correlation.mat');
        perm_file = fullfile(ch_folder, 'corr_perm.mat');
        if ~isfile(corr_file) || ~isfile(perm_file), continue; end

        load(corr_file, 'correlation');
        load(perm_file, 'corr_perm');

        if any(isnan(correlation)) || any(isnan(corr_perm(:))), continue; end

        corr_all_ch(ch,:) = correlation;
        corr_perm_all_ch = cat(3, corr_perm_all_ch, corr_perm);
    end

    corr_chan_avg = mean(corr_all_ch, 1, 'omitnan');

    corr_perm_chan_avg = mean(corr_perm_all_ch, 3);
    tmax_chan_avg = max(corr_perm_chan_avg, [], 2);
    thresh_chan_avg = quantile(tmax_chan_avg, 0.95);

    save_file = fullfile(output_RT, 'all_loc_difflev', 'channel_avg_results.mat');
    save(save_file, 'corr_chan_avg', 'corr_perm_chan_avg', 'tmax_chan_avg', 'thresh_chan_avg', 'corr_all_ch', 'freq');
    fprintf('Saved RT correlation channel-average results for %s\n', animalName);

    %% =====================================================================
    %  PLOTTING — Phase vs Hit/Miss (POS) (per animal)
    %  =====================================================================

    save_root_HM = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/phase_correlation/abs_per_chan', animalName, 'cp10_till_100', 'hit_miss', 'all_loc_difflev');
    if ~exist(save_root_HM, 'dir'), mkdir(save_root_HM); end

    limit_max_pos = nan(1, nCh);
    POS_chan       = false(nCh, nFreq);

    f2 = figure(2);
    set(f2, 'Units', 'centimeters', 'Position', [1 1 50 40]);
    set(f2, 'PaperUnits', 'centimeters', 'PaperSize', [50 40], 'PaperPosition', [0 0 50 40]);
    for ch = 1:nCh
        ch_folder = fullfile(output_HM, 'all_loc_difflev', num2str(ch));
        if ~exist(ch_folder, 'dir'), continue; end
        cd(ch_folder);

        if ~exist('pos.mat','file') || ~exist('pos_perm.mat','file')
            continue
        end

        load pos
        load pos_perm

        if any(isnan(pos)) || any(isnan(pos_perm(:)))
            continue
        end

        tmax = max(pos_perm, [], 2);
        limit_max_pos(ch) = quantile(tmax, 0.95);

        if isnan(limit_max_pos(ch)), continue; end

        subplot(8, 8, ch);
        plot_sigfreq(freq, pos, limit_max_pos(ch));
        title(['Ch ' num2str(ch)])

        POS_chan(ch,:) = pos >= limit_max_pos(ch);
    end

    print(f2, fullfile(save_root_HM, 'all_channels_hitmiss_pos.pdf'), '-dpdf');

    % Summary heatmap
    figure;
    imagesc(freq, 1:nCh, POS_chan);
    set(gca, 'YDir', 'normal');
    xlabel('Frequency (Hz)'); ylabel('Channels');
    title(sprintf('%s — Significant POS (Phase Opposition Sum) per Channel', animalName));
    caxis([0 1]); colorbar;
    saveas(gcf, fullfile(save_root_HM, 'summary_hitmiss_pos.pdf'));

    % Combined across channels
    valid_idx = ~isnan(limit_max_pos);
    pos_all = [];
    pos_perm_all = [];

    for ch = find(valid_idx)
        cd(fullfile(output_HM, 'all_loc_difflev', num2str(ch)));
        load pos
        load pos_perm
        pos_all = [pos_all; pos];
        pos_perm_all = cat(3, pos_perm_all, pos_perm);
    end

    if ~isempty(pos_all)
        pos_avg = nanmean(pos_all, 1);
        pos_perm_avg = nanmean(pos_perm_all, 3);
        tmax_all = nanmax(pos_perm_avg, [], 2);
        limit_avg = quantile(tmax_all, 0.95);

        figure;
        plot_sigfreq(freq, pos_avg, limit_avg);
        title(sprintf('%s — All Channels Combined - POS (Phase Opposition Sum)', animalName));
        xlabel('Frequency (Hz)'); ylabel('POS');
        saveas(gcf, fullfile(save_root_HM, 'combined_hitmiss_pos.pdf'));
    end

    %% Channel-average null distribution — POS

    fprintf('Computing POS channel-average null for %s...\n', animalName);

    pos_all_ch = NaN(nCh, nFreq);
    pos_perm_all_ch = [];

    for ch = 1:nCh
        ch_folder = fullfile(output_HM, 'all_loc_difflev', num2str(ch));
        if ~exist(ch_folder, 'dir'), continue; end

        pos_file = fullfile(ch_folder, 'pos.mat');
        perm_file = fullfile(ch_folder, 'pos_perm.mat');
        if ~isfile(pos_file) || ~isfile(perm_file), continue; end

        load(pos_file, 'pos');
        load(perm_file, 'pos_perm');

        if any(isnan(pos)) || any(isnan(pos_perm(:))), continue; end

        pos_all_ch(ch,:) = pos;
        pos_perm_all_ch = cat(3, pos_perm_all_ch, pos_perm);
    end

    pos_chan_avg = mean(pos_all_ch, 1, 'omitnan');

    pos_perm_chan_avg = mean(pos_perm_all_ch, 3);
    tmax_chan_avg_pos = max(pos_perm_chan_avg, [], 2);
    thresh_chan_avg_pos = quantile(tmax_chan_avg_pos, 0.95);

    save_file = fullfile(output_HM, 'all_loc_difflev', 'channel_avg_results_pos.mat');
    save(save_file, 'pos_chan_avg', 'pos_perm_chan_avg', 'tmax_chan_avg_pos', 'thresh_chan_avg_pos', 'pos_all_ch', 'freq');
    fprintf('Saved POS channel-average results for %s\n', animalName);

    %% =====================================================================
    %  PLOTTING — Phase vs Hit/Miss (ITC) (per animal)
    %  =====================================================================

    save_root_ITC = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/phase_correlation/abs_per_chan', animalName, 'cp10_till_100', 'hit_miss_itc', 'all_loc_difflev');
    if ~exist(save_root_ITC, 'dir'), mkdir(save_root_ITC); end

    limit_max_itc = nan(1, nCh);
    ITC_chan       = false(nCh, nFreq);

    f3 = figure(3);
    set(f3, 'Units', 'centimeters', 'Position', [1 1 50 40]);
    set(f3, 'PaperUnits', 'centimeters', 'PaperSize', [50 40], 'PaperPosition', [0 0 50 40]);
    for ch = 1:nCh
        ch_folder = fullfile(output_HM, 'all_loc_difflev', num2str(ch));
        if ~exist(ch_folder, 'dir'), continue; end
        cd(ch_folder);

        if ~exist('itc.mat','file') || ~exist('itc_perm.mat','file')
            continue
        end

        load itc
        load itc_perm

        if any(isnan(itc)) || any(isnan(itc_perm(:)))
            continue
        end

        tmax = max(itc_perm, [], 2);
        limit_max_itc(ch) = quantile(tmax, 0.95);

        if isnan(limit_max_itc(ch)), continue; end

        subplot(8, 8, ch);
        plot_sigfreq(freq, itc, limit_max_itc(ch));
        title(['Ch ' num2str(ch)])

        ITC_chan(ch,:) = itc >= limit_max_itc(ch);
    end

    print(f3, fullfile(save_root_ITC, 'all_channels_hitmiss_itc.pdf'), '-dpdf');

    % Summary heatmap
    figure;
    imagesc(freq, 1:nCh, ITC_chan);
    set(gca, 'YDir', 'normal');
    xlabel('Frequency (Hz)'); ylabel('Channels');
    title(sprintf('%s — Significant ITC (Inverted Miss Phases) per Channel', animalName));
    caxis([0 1]); colorbar;
    saveas(gcf, fullfile(save_root_ITC, 'summary_hitmiss_itc.pdf'));

    % Combined across channels
    valid_idx = ~isnan(limit_max_itc);
    itc_all = [];
    itc_perm_all = [];

    for ch = find(valid_idx)
        cd(fullfile(output_HM, 'all_loc_difflev', num2str(ch)));
        load itc
        load itc_perm
        itc_all = [itc_all; itc];
        itc_perm_all = cat(3, itc_perm_all, itc_perm);
    end

    if ~isempty(itc_all)
        itc_avg = nanmean(itc_all, 1);
        itc_perm_avg = nanmean(itc_perm_all, 3);
        tmax_all = nanmax(itc_perm_avg, [], 2);
        limit_avg = quantile(tmax_all, 0.95);

        figure;
        plot_sigfreq(freq, itc_avg, limit_avg);
        ylim([0 0.15]); title(sprintf('%s — All Channels Combined - ITC (Inverted Miss Phases)', animalName));
        xlabel('Frequency (Hz)'); ylabel('ITC');
        saveas(gcf, fullfile(save_root_ITC, 'combined_hitmiss_itc.pdf'));
    end

    %% Channel-average null distribution — ITC

    fprintf('Computing ITC channel-average null for %s...\n', animalName);

    itc_all_ch = NaN(nCh, nFreq);
    itc_perm_all_ch = [];

    for ch = 1:nCh
        ch_folder = fullfile(output_HM, 'all_loc_difflev', num2str(ch));
        if ~exist(ch_folder, 'dir'), continue; end

        itc_file = fullfile(ch_folder, 'itc.mat');
        perm_file = fullfile(ch_folder, 'itc_perm.mat');
        if ~isfile(itc_file) || ~isfile(perm_file), continue; end

        load(itc_file, 'itc');
        load(perm_file, 'itc_perm');

        if any(isnan(itc)) || any(isnan(itc_perm(:))), continue; end

        itc_all_ch(ch,:) = itc;
        itc_perm_all_ch = cat(3, itc_perm_all_ch, itc_perm);
    end

    itc_chan_avg = mean(itc_all_ch, 1, 'omitnan');

    itc_perm_chan_avg = mean(itc_perm_all_ch, 3);
    tmax_chan_avg_itc = max(itc_perm_chan_avg, [], 2);
    thresh_chan_avg_itc = quantile(tmax_chan_avg_itc, 0.95);

    save_file = fullfile(output_HM, 'all_loc_difflev', 'channel_avg_results_itc.mat');
    save(save_file, 'itc_chan_avg', 'itc_perm_chan_avg', 'tmax_chan_avg_itc', 'thresh_chan_avg_itc', 'itc_all_ch', 'freq');
    fprintf('Saved ITC channel-average results for %s\n', animalName);

    close all

end  % animal loop

%% =====================================================================
%  MONKEY-AVERAGE: combine across animals
%  =====================================================================

nAnimals = numel(animals);

%% Monkey-average — RT correlation

corr_monkey = [];
perm_monkey = [];

for a = 1:nAnimals
    animal_corr = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ...
        ['results_' animals{a}], 'phase_correlation', 'abs_per_chan', 'cp10_till_100');
    avg_file = fullfile(animal_corr, 'RT', 'all_loc_difflev', 'channel_avg_results.mat');

    if ~isfile(avg_file)
        warning('RT correlation channel-average results not found for %s.', animals{a});
        continue
    end

    tmp = load(avg_file);
    corr_monkey = cat(1, corr_monkey, tmp.corr_chan_avg);
    perm_monkey = cat(3, perm_monkey, tmp.corr_perm_chan_avg);
end

if size(corr_monkey, 1) == nAnimals
    freq = tmp.freq;

    corr_monkey_avg = mean(corr_monkey, 1);

    perm_monkey_avg = mean(perm_monkey, 3);
    tmax_monkey_avg = max(perm_monkey_avg, [], 2);
    thresh_monkey_avg = quantile(tmax_monkey_avg, 0.95);

    monkey_save_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/results_combined/phase_correlation/abs_per_chan/cp10_till_100', 'RT', 'all_loc_difflev');
    if ~exist(monkey_save_dir, 'dir'), mkdir(monkey_save_dir); end
    save(fullfile(monkey_save_dir, 'monkey_avg_results.mat'), ...
        'corr_monkey_avg', 'perm_monkey_avg', 'tmax_monkey_avg', ...
        'thresh_monkey_avg', 'corr_monkey', 'animals', 'freq');

    fprintf('RT Correlation Monkey-average threshold: %.4f\n', thresh_monkey_avg);

    monkey_plot_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/phase_correlation/abs_per_chan/monkey_avg/cp10_till_100/RT/all_loc_difflev');
    if ~exist(monkey_plot_dir, 'dir'), mkdir(monkey_plot_dir); end

    figure;
    plot_sigfreq(freq, corr_monkey_avg, thresh_monkey_avg);
    ylim([0 0.1]);
    title('Monkey-Average (Channel-Avg) Phase-RT Correlation');
    xlabel('Frequency (Hz)'); ylabel('Correlation');
    saveas(gcf, fullfile(monkey_plot_dir, 'monkey_avg_RT_corr.pdf'));

    figure; hold on;
    colors = lines(nAnimals);
    for a = 1:nAnimals
        plot(freq, corr_monkey(a,:), 'Color', colors(a,:), 'LineWidth', 1.5);
    end
    plot(freq, corr_monkey_avg, 'k', 'LineWidth', 2);
    yline(thresh_monkey_avg, '--r', 'LineWidth', 1.5);
    legend([animals, {'Monkey avg', 'Threshold'}], 'Location', 'best');
    xlabel('Frequency (Hz)'); ylabel('Correlation');
    title('Phase-RT Correlation: Per Animal + Monkey Average');
    saveas(gcf, fullfile(monkey_plot_dir, 'monkey_avg_RT_corr_overlay.pdf'));

    figure;
    histogram(tmax_monkey_avg, 50, 'FaceColor', [0.7 0.7 0.7]);
    xline(thresh_monkey_avg, '--r', 'LineWidth', 2);
    xlabel('Max correlation (permutation)'); ylabel('Count');
    title('Monkey-Average Null Distribution (tmax) — Phase-RT Correlation');
    legend({'Null distribution', '95th percentile'}, 'Location', 'best');
    saveas(gcf, fullfile(monkey_plot_dir, 'monkey_avg_RT_null_distribution.pdf'));
else
    warning('Not all animals have RT correlation channel-average results.');
end

%% Monkey-average — POS (Phase Opposition Sum)

pos_monkey = [];
perm_monkey_pos = [];

for a = 1:nAnimals
    animal_corr = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ...
        ['results_' animals{a}], 'phase_correlation', 'abs_per_chan', 'cp10_till_100');
    avg_file = fullfile(animal_corr, 'hit_miss', 'all_loc_difflev', 'channel_avg_results_pos.mat');

    if ~isfile(avg_file)
        warning('POS channel-average results not found for %s.', animals{a});
        continue
    end

    tmp = load(avg_file);
    pos_monkey = cat(1, pos_monkey, tmp.pos_chan_avg);
    perm_monkey_pos = cat(3, perm_monkey_pos, tmp.pos_perm_chan_avg);
end

if size(pos_monkey, 1) == nAnimals
    freq = tmp.freq;

    pos_monkey_avg = mean(pos_monkey, 1);

    perm_monkey_avg_pos = mean(perm_monkey_pos, 3);
    tmax_monkey_avg_pos = max(perm_monkey_avg_pos, [], 2);
    thresh_monkey_avg_pos = quantile(tmax_monkey_avg_pos, 0.95);

    monkey_save_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/results_combined/phase_correlation/abs_per_chan/cp10_till_100', 'hit_miss', 'all_loc_difflev');
    if ~exist(monkey_save_dir, 'dir'), mkdir(monkey_save_dir); end
    save(fullfile(monkey_save_dir, 'monkey_avg_results_pos.mat'), ...
        'pos_monkey_avg', 'perm_monkey_avg_pos', 'tmax_monkey_avg_pos', ...
        'thresh_monkey_avg_pos', 'pos_monkey', 'animals', 'freq');

    fprintf('POS Monkey-average threshold: %.4f\n', thresh_monkey_avg_pos);

    monkey_plot_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/phase_correlation/abs_per_chan/monkey_avg/cp10_till_100/hit_miss/all_loc_difflev');
    if ~exist(monkey_plot_dir, 'dir'), mkdir(monkey_plot_dir); end

    figure;
    plot_sigfreq(freq, pos_monkey_avg, thresh_monkey_avg_pos);
    title('Monkey-Average (Channel-Avg) POS (Phase Opposition Sum)');
    xlabel('Frequency (Hz)'); ylabel('POS');
    saveas(gcf, fullfile(monkey_plot_dir, 'monkey_avg_hitmiss_pos.pdf'));

    figure; hold on;
    colors = lines(nAnimals);
    for a = 1:nAnimals
        plot(freq, pos_monkey(a,:), 'Color', colors(a,:), 'LineWidth', 1.5);
    end
    plot(freq, pos_monkey_avg, 'k', 'LineWidth', 2);
    yline(thresh_monkey_avg_pos, '--r', 'LineWidth', 1.5);
    legend([animals, {'Monkey avg', 'Threshold'}], 'Location', 'best');
    xlabel('Frequency (Hz)'); ylabel('POS');
    title('POS: Per Animal + Monkey Average');
    saveas(gcf, fullfile(monkey_plot_dir, 'monkey_avg_hitmiss_pos_overlay.pdf'));

    figure;
    histogram(tmax_monkey_avg_pos, 50, 'FaceColor', [0.7 0.7 0.7]);
    xline(thresh_monkey_avg_pos, '--r', 'LineWidth', 2);
    xlabel('Max POS (permutation)'); ylabel('Count');
    title('Monkey-Average Null Distribution (tmax) — POS');
    legend({'Null distribution', '95th percentile'}, 'Location', 'best');
    saveas(gcf, fullfile(monkey_plot_dir, 'monkey_avg_hitmiss_pos_null_distribution.pdf'));
else
    warning('Not all animals have POS channel-average results.');
end

%% Monkey-average — ITC (Inverted Miss Phases)

itc_monkey = [];
perm_monkey_itc = [];

for a = 1:nAnimals
    animal_corr = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ...
        ['results_' animals{a}], 'phase_correlation', 'abs_per_chan', 'cp10_till_100');
    avg_file = fullfile(animal_corr, 'hit_miss', 'all_loc_difflev', 'channel_avg_results_itc.mat');

    if ~isfile(avg_file)
        warning('ITC channel-average results not found for %s.', animals{a});
        continue
    end

    tmp = load(avg_file);
    itc_monkey = cat(1, itc_monkey, tmp.itc_chan_avg);
    perm_monkey_itc = cat(3, perm_monkey_itc, tmp.itc_perm_chan_avg);
end

if size(itc_monkey, 1) == nAnimals
    freq = tmp.freq;

    itc_monkey_avg = mean(itc_monkey, 1);

    perm_monkey_avg_itc = mean(perm_monkey_itc, 3);
    tmax_monkey_avg_itc = max(perm_monkey_avg_itc, [], 2);
    thresh_monkey_avg_itc = quantile(tmax_monkey_avg_itc, 0.95);

    monkey_save_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/results_combined/phase_correlation/abs_per_chan/cp10_till_100', 'hit_miss_itc', 'all_loc_difflev');
    if ~exist(monkey_save_dir, 'dir'), mkdir(monkey_save_dir); end
    save(fullfile(monkey_save_dir, 'monkey_avg_results_itc.mat'), ...
        'itc_monkey_avg', 'perm_monkey_avg_itc', 'tmax_monkey_avg_itc', ...
        'thresh_monkey_avg_itc', 'itc_monkey', 'animals', 'freq');

    fprintf('ITC Monkey-average threshold: %.4f\n', thresh_monkey_avg_itc);

    monkey_plot_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/phase_correlation/abs_per_chan/monkey_avg/cp10_till_100/hit_miss_itc/all_loc_difflev');
    if ~exist(monkey_plot_dir, 'dir'), mkdir(monkey_plot_dir); end

    figure;
    plot_sigfreq(freq, itc_monkey_avg, thresh_monkey_avg_itc);
    ylim([0 0.15]);
    title('Monkey-Average (Channel-Avg) ITC (Inverted Miss Phases)');
    xlabel('Frequency (Hz)'); ylabel('ITC');
    saveas(gcf, fullfile(monkey_plot_dir, 'monkey_avg_hitmiss_itc.pdf'));

    figure; hold on;
    colors = lines(nAnimals);
    for a = 1:nAnimals
        plot(freq, itc_monkey(a,:), 'Color', colors(a,:), 'LineWidth', 1.5);
    end
    plot(freq, itc_monkey_avg, 'k', 'LineWidth', 2);
    yline(thresh_monkey_avg_itc, '--r', 'LineWidth', 1.5);
    legend([animals, {'Monkey avg', 'Threshold'}], 'Location', 'best');
    xlabel('Frequency (Hz)'); ylabel('ITC');
    title('ITC (Inverted Miss Phases): Per Animal + Monkey Average');
    saveas(gcf, fullfile(monkey_plot_dir, 'monkey_avg_hitmiss_itc_overlay.pdf'));

    figure;
    histogram(tmax_monkey_avg_itc, 50, 'FaceColor', [0.7 0.7 0.7]);
    xline(thresh_monkey_avg_itc, '--r', 'LineWidth', 2);
    xlabel('Max ITC (permutation)'); ylabel('Count');
    title('Monkey-Average Null Distribution (tmax) — ITC');
    legend({'Null distribution', '95th percentile'}, 'Location', 'best');
    saveas(gcf, fullfile(monkey_plot_dir, 'monkey_avg_hitmiss_itc_null_distribution.pdf'));
else
    warning('Not all animals have ITC channel-average results.');
end
