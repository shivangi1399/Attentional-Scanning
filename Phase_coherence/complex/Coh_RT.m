% =====================================================================
% Phase coherence: pre-stimulus phase vs. reaction time (hit trials only)
% Hypothesis H1 (complex/)
%
% Claim: a single optimal phase is shared across all trials, positions,
% difficulty levels, channels, and animals.
%
% Recipe: pool all trials in complex space within each channel; average
% complex resultants across channels and animals; take abs() only at
% the very end. Way 1 at every level.
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
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/Phase_coherence
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/Phase_coherence/functions
addpath /mnt/hpc/projects/MWSampling/4Shivangi
clc

for a = 1:numel(animals)
    animalName = animals{a};
    fprintf('\n=== Processing %s ===\n', animalName);

    data_folder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName], 'multi_lin_reg', 'cp10_till_100');
    coh_folder  = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName], 'phase_coherence', 'complex', 'cp10_till_100');

    %% =====================================================================
    %  Phase coherence with RT (all locations and difficulty levels)
    %  =====================================================================

    cd(data_folder)
    load('ph_all_sess.mat')

    % Select hit trials only (misses have NaN RT)
    hit_idx = find(ph_comb.RT_trialinfo(:,20) == 1);

    output_coh_RT = fullfile(coh_folder, 'RT');

    %% Real data — phase_coherence(phase, RT)

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

        [coh, phase_spec, coh_complex] = phase_coherence(phase_clean, rt_clean);

        chan_folder = fullfile(output_coh_RT, 'all_loc_difflev', num2str(ichan));
        if ~exist(chan_folder, 'dir'), mkdir(chan_folder); end
        save(fullfile(chan_folder, 'coherence.mat'), 'coh', 'phase_spec', 'coh_complex');
    end

    %% Permutation (SLURM)

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
        cfg{ichan}.outfile      = fullfile(output_coh_RT, 'all_loc_difflev');
        cfg{ichan}.perm_indices = perm_indices;
        cfg{ichan}.trial_idx    = hit_idx;
    end

    slurmfun(@phase_coherence_RT_perm, cfg, ...
        'partition',   '8GB', ...
        'stopOnError', false, ...
        'useUserPath', true);

    %% =====================================================================
    %  PLOTTING — Phase coherence with RT (per animal)
    %  =====================================================================

    save_root_RT = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/phase_coherence/complex', animalName, 'cp10_till_100', 'RT', 'all_loc_difflev');
    if ~exist(save_root_RT, 'dir'), mkdir(save_root_RT); end

    cd(coh_folder)
    load('frequency.mat')
    freq = frequency;
    nFreq = numel(freq);

    limit_maxc = nan(1, nCh);
    Coh_chan   = false(nCh, nFreq);

    f1 = figure(1);
    set(f1, 'Units', 'centimeters', 'Position', [1 1 50 40]);
    set(f1, 'PaperUnits', 'centimeters', 'PaperSize', [50 40], 'PaperPosition', [0 0 50 40]);
    f2 = figure(2);
    set(f2, 'Units', 'centimeters', 'Position', [1 1 50 40]);
    set(f2, 'PaperUnits', 'centimeters', 'PaperSize', [50 40], 'PaperPosition', [0 0 50 40]);

    for ch = 1:nCh
        ch_folder = fullfile(output_coh_RT, 'all_loc_difflev', num2str(ch));
        if ~exist(ch_folder, 'dir'), continue; end
        cd(ch_folder);

        if ~exist('coherence.mat','file') || ~exist('coh_perm_complex.mat','file')
            continue
        end

        load coherence
        load coh_perm_complex

        coh     = abs(coh_complex);
        coh_perm = abs(coh_perm_complex);

        if any(isnan(coh)) || any(isnan(coh_perm(:))) || any(isnan(phase_spec))
            continue
        end

        tmaxc = max(coh_perm, [], 2);
        limit_maxc(ch) = quantile(tmaxc, 0.95);

        if isnan(limit_maxc(ch)), continue; end

        figure(f1);
        subplot(8, 8, ch);
        plot_sig(freq, coh, limit_maxc(ch), 'Frequency', 'Coherence');
        title(['Ch ' num2str(ch)])

        figure(f2);
        subplot(8, 8, ch);
        plot_sig(freq, phase_spec, [], 'Frequency', 'Phase spec');
        title(['Ch ' num2str(ch)])

        Coh_chan(ch,:) = coh >= limit_maxc(ch);
    end

    print(f1, fullfile(save_root_RT, 'all_channels_RT_coherence.pdf'), '-dpdf');
    print(f2, fullfile(save_root_RT, 'all_channels_RT_phase.pdf'), '-dpdf');

    % Summary heatmaps
    f3 = figure;
    imagesc(freq, 1:nCh, Coh_chan);
    set(gca, 'YDir', 'normal');
    xlabel('Frequency (Hz)'); ylabel('Channels');
    title(sprintf('%s — Significant Phase-RT Coherence per Channel', animalName));
    caxis([0 1]); colorbar;
    saveas(f3, fullfile(save_root_RT, 'summary_RT_coherence.pdf'));

    % Combined across channels — average in complex space
    valid_idx = ~isnan(limit_maxc);
    coh_complex_all = [];
    coh_perm_complex_all = [];

    for ch = find(valid_idx)
        tmp_c = load(fullfile(output_coh_RT, 'all_loc_difflev', num2str(ch), 'coherence.mat'), 'coh_complex');
        tmp_p = load(fullfile(output_coh_RT, 'all_loc_difflev', num2str(ch), 'coh_perm_complex.mat'), 'coh_perm_complex');
        coh_complex_all      = [coh_complex_all; tmp_c.coh_complex];
        coh_perm_complex_all = cat(3, coh_perm_complex_all, tmp_p.coh_perm_complex);
    end

    if ~isempty(coh_complex_all)
        coh_avg      = abs(mean(coh_complex_all, 1));
        perm_avg_cplx = mean(coh_perm_complex_all, 3);
        coh_perm_avg = abs(perm_avg_cplx);
        tmax_all     = max(coh_perm_avg, [], 2);
        limit_avg    = quantile(tmax_all, 0.95);

        f5 = figure;
        plot_sig(freq, coh_avg, limit_avg, 'Frequency', 'Coherence');
        title(sprintf('%s — All Channels Combined - Phase vs RT Coherence', animalName));
        saveas(f5, fullfile(save_root_RT, 'combined_RT_coherence.pdf'));
    end

    %% Channel-average null distribution (per animal)

    fprintf('Computing RT channel-average permutation null for %s...\n', animalName);

    coh_complex_all_ch      = NaN(nCh, nFreq);   % complex coherence per channel
    coh_perm_complex_all_ch = [];                 % complex perm null [permut_n x nFreq x nCh]

    for ch = 1:nCh
        ch_folder = fullfile(output_coh_RT, 'all_loc_difflev', num2str(ch));
        if ~exist(ch_folder, 'dir'), continue; end

        coh_file  = fullfile(ch_folder, 'coherence.mat');
        perm_file = fullfile(ch_folder, 'coh_perm_complex.mat');
        if ~isfile(coh_file) || ~isfile(perm_file), continue; end

        tmp_c = load(coh_file,  'coh_complex');
        tmp_p = load(perm_file, 'coh_perm_complex');

        if any(isnan(tmp_c.coh_complex)) || any(isnan(tmp_p.coh_perm_complex(:))), continue; end

        coh_complex_all_ch(ch,:)    = tmp_c.coh_complex;
        coh_perm_complex_all_ch     = cat(3, coh_perm_complex_all_ch, tmp_p.coh_perm_complex);
    end

    % Average in complex space, decompose only at the end
    coh_complex_chan_avg      = mean(coh_complex_all_ch, 1, 'omitnan');
    coh_chan_avg              = abs(coh_complex_chan_avg);
    phase_chan_avg            = angle(coh_complex_chan_avg);

    coh_perm_chan_avg_complex = mean(coh_perm_complex_all_ch, 3);   % [permut_n x nFreq]
    coh_perm_chan_avg         = abs(coh_perm_chan_avg_complex);
    tmax_chan_avg             = max(coh_perm_chan_avg, [], 2);
    thresh_chan_avg           = quantile(tmax_chan_avg, 0.95);

    save_file = fullfile(output_coh_RT, 'all_loc_difflev', 'channel_avg_results.mat');
    save(save_file, 'coh_chan_avg', 'coh_complex_chan_avg', 'phase_chan_avg', ...
        'coh_perm_chan_avg_complex', 'coh_perm_chan_avg', ...
        'tmax_chan_avg', 'thresh_chan_avg', 'coh_complex_all_ch', 'freq');
    fprintf('Saved RT channel-average results for %s\n', animalName);

    close all

end  % animal loop

%% Monkey-average: combine across animals

nAnimals = numel(animals);

coh_complex_monkey  = [];   % [nAnimals x nFreq] complex
perm_complex_monkey = [];   % [permut_n x nFreq x nAnimals] complex

for a = 1:nAnimals
    animal_coh = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ...
        ['results_' animals{a}], 'phase_coherence', 'complex', 'cp10_till_100');
    avg_file = fullfile(animal_coh, 'RT', 'all_loc_difflev', 'channel_avg_results.mat');

    if ~isfile(avg_file)
        warning('RT channel-average results not found for %s. Run per-animal section first.', animals{a});
        continue
    end

    tmp = load(avg_file);
    coh_complex_monkey  = cat(1, coh_complex_monkey,  tmp.coh_complex_chan_avg);
    perm_complex_monkey = cat(3, perm_complex_monkey, tmp.coh_perm_chan_avg_complex);
end

if size(coh_complex_monkey, 1) == nAnimals
    freq = tmp.freq;

    % Average complex across animals, decompose only at the end
    coh_complex_monkey_avg = mean(coh_complex_monkey, 1);
    coh_monkey_avg         = abs(coh_complex_monkey_avg);
    phase_monkey_avg       = angle(coh_complex_monkey_avg);

    % Per-animal magnitudes (for overlay plot)
    coh_monkey = abs(coh_complex_monkey);

    perm_monkey_avg_complex = mean(perm_complex_monkey, 3);   % [permut_n x nFreq]
    perm_monkey_avg   = abs(perm_monkey_avg_complex);
    tmax_monkey_avg   = max(perm_monkey_avg, [], 2);
    thresh_monkey_avg = quantile(tmax_monkey_avg, 0.95);

    monkey_save_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/results_combined/phase_coherence/complex/cp10_till_100', 'RT', 'all_loc_difflev');
    if ~exist(monkey_save_dir, 'dir'), mkdir(monkey_save_dir); end
    save(fullfile(monkey_save_dir, 'monkey_avg_results.mat'), ...
        'coh_monkey_avg', 'coh_complex_monkey_avg', 'phase_monkey_avg', ...
        'perm_monkey_avg', 'tmax_monkey_avg', ...
        'thresh_monkey_avg', 'coh_monkey', 'coh_complex_monkey', 'animals', 'freq');

    fprintf('RT Monkey-average threshold: %.4f\n', thresh_monkey_avg);

    %% Monkey-level figures

    monkey_plot_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/phase_coherence/complex/monkey_avg/cp10_till_100/RT/all_loc_difflev');
    if ~exist(monkey_plot_dir, 'dir'), mkdir(monkey_plot_dir); end

    % Coherence: monkey average with threshold
    f7 = figure;
    plot_sig(freq, coh_monkey_avg, thresh_monkey_avg, 'Frequency (Hz)', 'Coherence');
    title('Monkey-Average (Channel-Avg) Phase-RT Coherence');
    saveas(f7, fullfile(monkey_plot_dir, 'monkey_avg_RT_coherence.pdf'));

    % Per-animal overlay
    f8 = figure; hold on;
    colors = lines(nAnimals);
    for a = 1:nAnimals
        plot(freq, coh_monkey(a,:), 'Color', colors(a,:), 'LineWidth', 1.5);
    end
    plot(freq, coh_monkey_avg, 'k', 'LineWidth', 2);
    yline(thresh_monkey_avg, '--r', 'LineWidth', 1.5);
    legend([animals, {'Monkey avg', 'Threshold'}], 'Location', 'best');
    xlabel('Frequency (Hz)'); ylabel('Coherence');
    title('Phase-RT Coherence: Per Animal + Monkey Average');
    saveas(f8, fullfile(monkey_plot_dir, 'monkey_avg_RT_coherence_overlay.pdf'));

    % Null distribution histogram
    f9 = figure;
    histogram(tmax_monkey_avg, 50, 'FaceColor', [0.7 0.7 0.7]);
    xline(thresh_monkey_avg, '--r', 'LineWidth', 2);
    xlabel('Max coherence (permutation)'); ylabel('Count');
    title('Monkey-Average Null Distribution (tmax) — Phase-RT Coherence');
    legend({'Null distribution', '95th percentile'}, 'Location', 'best');
    saveas(f9, fullfile(monkey_plot_dir, 'monkey_avg_RT_null_distribution.pdf'));

else
    warning('Not all animals have RT channel-average results. Run per-animal section for each animal first.');
end
