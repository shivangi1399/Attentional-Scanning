% Coherence between pre stimulus phase and reaction time - only absolute saved
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
    coh_folder  = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName], 'phase_coherence', 'abs_per_chan', 'cp10_till_100');
    if ~exist(coh_folder,'dir'), mkdir(coh_folder); end

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

        [coh, phase_spec] = phase_coherence(phase_clean, rt_clean);

        chan_folder = fullfile(output_coh_RT, 'all_loc_difflev', num2str(ichan));
        if ~exist(chan_folder, 'dir'), mkdir(chan_folder); end
        save(fullfile(chan_folder, 'coherence.mat'), 'coh', 'phase_spec');
    end

    %% Permutation (SLURM)

    cd(data_folder)
    load('ph_all_sess.mat')

    hit_idx  = find(ph_comb.RT_trialinfo(:,20) == 1);
    nTrials  = length(hit_idx);

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

    save_root_RT = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/phase_coherence/abs_per_chan', animalName, 'cp10_till_100', 'RT', 'all_loc_difflev');
    if ~exist(save_root_RT, 'dir'), mkdir(save_root_RT); end

    cd(coh_folder)
    load('frequency.mat')
    freq = frequency;
    nFreq = numel(freq);

    limit_maxc = nan(1, nCh);
    limit_maxp = nan(1, nCh);
    Coh_chan    = false(nCh, nFreq);
    Phase_chan  = false(nCh, nFreq);

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

        if ~exist('coherence.mat','file') || ~exist('coh_perm.mat','file') || ...
           ~exist('phase_spec_perm.mat','file')
            continue
        end

        load coherence
        load coh_perm
        load phase_spec_perm

        if any(isnan(coh)) || any(isnan(coh_perm(:))) || ...
           any(isnan(phase_spec)) || any(isnan(phase_spec_perm(:)))
            continue
        end

        tmaxc = max(coh_perm, [], 2);
        limit_maxc(ch) = quantile(tmaxc, 0.95);

        tmaxp = max(phase_spec_perm, [], 2);
        limit_maxp(ch) = quantile(tmaxp, 0.95);

        if isnan(limit_maxc(ch)) || isnan(limit_maxp(ch)), continue; end

        figure(f1);
        subplot(8, 8, ch);
        plot_sig(freq, coh, limit_maxc(ch), 'Frequency', 'Coherence');
        title(['Ch ' num2str(ch)])

        figure(f2);
        subplot(8, 8, ch);
        plot_sig(freq, phase_spec, limit_maxp(ch), 'Frequency', 'Phase spec');
        title(['Ch ' num2str(ch)])

        Coh_chan(ch,:)   = coh >= limit_maxc(ch);
        Phase_chan(ch,:) = phase_spec >= limit_maxp(ch);
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

    f4 = figure;
    imagesc(freq, 1:nCh, Phase_chan);
    set(gca, 'YDir', 'normal');
    xlabel('Frequency (Hz)'); ylabel('Channels');
    title(sprintf('%s — Significant Phase-RT Phase Spectrum per Channel', animalName));
    caxis([0 1]); colorbar;
    saveas(f4, fullfile(save_root_RT, 'summary_RT_phase.pdf'));

    % Combined across channels
    valid_idx = ~isnan(limit_maxc) & ~isnan(limit_maxp);
    coh_all = []; phase_all = [];
    coh_perm_all = []; phase_perm_all = [];

    for ch = find(valid_idx)
        cd(fullfile(output_coh_RT, 'all_loc_difflev', num2str(ch)));
        load coherence
        load coh_perm
        load phase_spec_perm
        coh_all   = [coh_all; coh];
        phase_all = [phase_all; phase_spec];
        coh_perm_all   = cat(3, coh_perm_all, coh_perm);
        phase_perm_all = cat(3, phase_perm_all, phase_spec_perm);
    end

    if ~isempty(coh_all)
        coh_avg = nanmean(coh_all, 1);
        coh_perm_avg = nanmean(coh_perm_all, 3);
        tmax_all = nanmax(coh_perm_avg, [], 2);
        limit_avg = quantile(tmax_all, 0.95);

        f5 = figure;
        plot_sig(freq, coh_avg, limit_avg, 'Frequency', 'Coherence');
        title(sprintf('%s — All Channels Combined - Phase vs RT Coherence', animalName));
        saveas(f5, fullfile(save_root_RT, 'combined_RT_coherence.pdf'));

        phase_avg = nanmean(phase_all, 1);
        phase_perm_avg = nanmean(phase_perm_all, 3);
        tmaxp_all = nanmax(phase_perm_avg, [], 2);
        limit_avgp = quantile(tmaxp_all, 0.95);

        f6 = figure;
        plot_sig(freq, phase_avg, limit_avgp, 'Frequency', 'Phase spec');
        title(sprintf('%s — All Channels Combined - Phase vs RT Phase Spec', animalName));
        saveas(f6, fullfile(save_root_RT, 'combined_RT_phase.pdf'));
    end

    %% Channel-average null distribution (per animal)

    fprintf('Computing RT channel-average permutation null for %s...\n', animalName);

    coh_all_ch = NaN(nCh, nFreq);
    coh_perm_all_ch = [];

    for ch = 1:nCh
        ch_folder = fullfile(output_coh_RT, 'all_loc_difflev', num2str(ch));
        if ~exist(ch_folder, 'dir'), continue; end

        coh_file = fullfile(ch_folder, 'coherence.mat');
        perm_file = fullfile(ch_folder, 'coh_perm.mat');
        if ~isfile(coh_file) || ~isfile(perm_file), continue; end

        load(coh_file, 'coh');
        load(perm_file, 'coh_perm');

        if any(isnan(coh)) || any(isnan(coh_perm(:))), continue; end

        coh_all_ch(ch,:) = coh;
        coh_perm_all_ch = cat(3, coh_perm_all_ch, coh_perm);
    end

    coh_chan_avg = mean(coh_all_ch, 1, 'omitnan');

    coh_perm_chan_avg = mean(coh_perm_all_ch, 3);
    tmax_chan_avg = max(coh_perm_chan_avg, [], 2);
    thresh_chan_avg = quantile(tmax_chan_avg, 0.95);

    save_file = fullfile(output_coh_RT, 'all_loc_difflev', 'channel_avg_results.mat');
    save(save_file, 'coh_chan_avg', 'coh_perm_chan_avg', 'tmax_chan_avg', 'thresh_chan_avg', 'coh_all_ch', 'freq');
    fprintf('Saved RT channel-average results for %s\n', animalName);

    close all

end  % animal loop

%% Monkey-average: combine across animals

nAnimals = numel(animals);

coh_monkey = [];
perm_monkey = [];

for a = 1:nAnimals
    animal_coh = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ...
        ['results_' animals{a}], 'phase_coherence', 'abs_per_chan', 'cp10_till_100');
    avg_file = fullfile(animal_coh, 'RT', 'all_loc_difflev', 'channel_avg_results.mat');

    if ~isfile(avg_file)
        warning('RT channel-average results not found for %s. Run per-animal section first.', animals{a});
        continue
    end

    tmp = load(avg_file);
    coh_monkey = cat(1, coh_monkey, tmp.coh_chan_avg);
    perm_monkey = cat(3, perm_monkey, tmp.coh_perm_chan_avg);
end

if size(coh_monkey, 1) == nAnimals
    freq = tmp.freq;

    coh_monkey_avg = mean(coh_monkey, 1);

    perm_monkey_avg = mean(perm_monkey, 3);
    tmax_monkey_avg = max(perm_monkey_avg, [], 2);
    thresh_monkey_avg = quantile(tmax_monkey_avg, 0.95);

    monkey_save_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/results_combined/phase_coherence/abs_per_chan/cp10_till_100', 'RT', 'all_loc_difflev');
    if ~exist(monkey_save_dir, 'dir'), mkdir(monkey_save_dir); end
    save(fullfile(monkey_save_dir, 'monkey_avg_results.mat'), ...
        'coh_monkey_avg', 'perm_monkey_avg', 'tmax_monkey_avg', ...
        'thresh_monkey_avg', 'coh_monkey', 'animals', 'freq');

    fprintf('RT Monkey-average threshold: %.4f\n', thresh_monkey_avg);

    %% Monkey-level figures

    monkey_plot_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/phase_coherence/abs_per_chan/monkey_avg/cp10_till_100/RT/all_loc_difflev');
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
