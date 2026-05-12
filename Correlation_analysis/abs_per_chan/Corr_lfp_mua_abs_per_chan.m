% =====================================================================
% Circular-linear correlation: pre-stimulus phase vs. post-stimulus
%                              LFP / MUA amplitude
% Hypothesis H1+H4 (abs_per_chan/)
%
% Claim: trials within a channel share a preferred phase (H1 at trial
% level), but channels are NOT required to share a preferred phase.
%
% Recipe: circ_corrcl on all trials per channel (the function returns a
% non-negative magnitude — abs is implicit at channel level, Way 2
% across channels); arithmetic mean of magnitudes across channels and
% animals.
%
% See sampling_compare/README.md for the Way-1 / Way-2 framing.
% =====================================================================
clear all
close all
clc

%% Settings
signal_types = {'lfp', 'mua'};  % choose: {'lfp'}, {'mua'}, or {'lfp', 'mua'}
animals = {'hermes', 'klecks'};
permut_n = 1000;

%% Specify paths

addpath /opt/fieldtrip_github/
ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/Correlation_analysis
addpath /mnt/hpc/projects/MWSampling/4Shivangi
addpath /mnt/hpc/projects/MWSampling/4Shivangi/software_folder/CircStat2012a
clc

for a = 1:numel(animals)
    animalName = animals{a};
    fprintf('\n=== Processing %s ===\n', animalName);

    %% Create data paths
    datafolder   = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName]);

    cd(datafolder),
    temp = dir;
    session_names = [];
    ii = 0;
    for i = 1:length(temp)
        if strfind(temp(i).name,animalName)
            ii = ii+1;
            session_names{ii,1} = temp(i).name;
        end
    end

    session_paths_files = [];
    session_paths_files = cellfun(@(x) fullfile(datafolder,x, 'clean_lfp.mat'), session_names, 'uniform',0);

    phase_paths = cellfun(@(x) fullfile(datafolder, x,'Phase_analysis/hit_miss'),session_names, 'uniform',0);
    output_folder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName], 'phase_correlation', 'abs_per_chan', 'cp10_till_100');
    if ~exist(output_folder,'dir'), mkdir(output_folder); end
    data_load_folder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName], 'multi_lin_reg', 'cp10_till_100');

    for s = 1:numel(signal_types)
        sig = signal_types{s};

        switch sig
            case 'lfp'
                amp_field = 'LFP_ERP_ampl_all';
                perm_func = @circlin_correlation_lfp;
            case 'mua'
                amp_field = 'MUA_ERP_ampl_all';
                perm_func = @circlin_correlation_mua;
        end

        fprintf('\n--- %s: %s ---\n', animalName, upper(sig));

        %% Correlation all locations and difficulty levels

        % circular-linear correlation (real)
        cd(data_load_folder)
        load('ph_all_sess.mat')

        for ichan = 1:64
            ichan

            phase = ph_comb.phase_all(:,:,ichan);
            erp_amp = ph_comb.(amp_field)(:,ichan);

            for foi = 1:length(ph_comb.phase_all(1,:,1))
                [correlation(1,foi),pvalue(1,foi)] = circ_corrcl(phase(:,foi), erp_amp);
            end

            chan_folder = fullfile(output_folder, sig, 'all_loc_difflev', num2str(ichan));
            if ~exist(chan_folder, 'dir'), mkdir(chan_folder); end
            cd(chan_folder)

            save correlation correlation
            save pvalue pvalue
        end

        % circular-linear correlation (permutation)

        cd(data_load_folder)
        load('ph_all_sess.mat')
        nTrials = size(ph_comb.phase_all, 1);

        rng(2025)
        perm_indices = arrayfun(@(x) randperm(nTrials), 1:permut_n, 'UniformOutput', false);

        cfg = cell(1,64);
        for i = 1:64
            cfg{i}.ichan = i;
            cfg{i}.permut_n = permut_n;
            cfg{i}.infile = fullfile(data_load_folder);
            cfg{i}.outfile = fullfile(output_folder, sig, 'all_loc_difflev');
            cfg{i}.perm_indices = perm_indices;
        end

        slurmfun(perm_func, cfg, ...
            'partition',     '8GB', ...
            'stopOnError',   false, ...
            'useUserPath',   true);

        %% Channel-average null distribution (per animal)

        fprintf('Computing %s channel-average permutation null for %s...\n', upper(sig), animalName);

        cd(data_load_folder)
        load('frequency.mat')
        freq = frequency;
        nCh = 64;
        nFreq = numel(freq);

        corr_all_ch = NaN(nCh, nFreq);
        corr_perm_all_ch = [];

        for ch = 1:nCh
            ch_folder = fullfile(output_folder, sig, 'all_loc_difflev', num2str(ch));
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

        save_file = fullfile(output_folder, sig, 'all_loc_difflev', 'channel_avg_results.mat');
        save(save_file, 'corr_chan_avg', 'corr_perm_chan_avg', 'tmax_chan_avg', 'thresh_chan_avg', 'corr_all_ch', 'freq');
        fprintf('Saved %s channel-average results for %s\n', upper(sig), animalName);

    end  % signal type loop

end  % animal loop

%% Monkey-average: combine across animals

nAnimals = numel(animals);

for s = 1:numel(signal_types)
    sig = signal_types{s};

    corr_monkey = [];
    perm_monkey = [];

    for a = 1:nAnimals
        animal_output = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ...
            ['results_' animals{a}], 'phase_correlation', 'abs_per_chan', 'cp10_till_100');
        avg_file = fullfile(animal_output, sig, 'all_loc_difflev', 'channel_avg_results.mat');

        if ~isfile(avg_file)
            warning('%s channel-average results not found for %s. Run per-animal section first.', upper(sig), animals{a});
            continue
        end

        tmp = load(avg_file);
        corr_monkey = cat(1, corr_monkey, tmp.corr_chan_avg);
        perm_monkey = cat(3, perm_monkey, tmp.corr_perm_chan_avg);
    end

    if size(corr_monkey, 1) == nAnimals
        corr_monkey_avg = mean(corr_monkey, 1);

        perm_monkey_avg = mean(perm_monkey, 3);
        tmax_monkey_avg = max(perm_monkey_avg, [], 2);
        thresh_monkey_avg = quantile(tmax_monkey_avg, 0.95);

        monkey_save_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/results_combined/phase_correlation/abs_per_chan/cp10_till_100', sig, 'all_loc_difflev');
        if ~exist(monkey_save_dir, 'dir'), mkdir(monkey_save_dir); end
        save(fullfile(monkey_save_dir, 'monkey_avg_results.mat'), ...
            'corr_monkey_avg', 'perm_monkey_avg', 'tmax_monkey_avg', ...
            'thresh_monkey_avg', 'corr_monkey', 'animals', 'freq');

        fprintf('%s Monkey-average threshold: %.4f\n', upper(sig), thresh_monkey_avg);

        monkey_plot_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/phase_correlation/abs_per_chan/monkey_avg/cp10_till_100', sig, 'all_loc_difflev');
        if ~exist(monkey_plot_dir, 'dir'), mkdir(monkey_plot_dir); end

        figure;
        plot_sigfreq(freq, corr_monkey_avg, thresh_monkey_avg);
        ylim([0 0.1]);
        title(sprintf('Monkey-Average (Channel-Avg) %s Correlation', upper(sig)));
        xlabel('Frequency (Hz)'); ylabel('Correlation');
        saveas(gcf, fullfile(monkey_plot_dir, 'monkey_avg_corr.pdf'));

        figure; hold on;
        colors = lines(nAnimals);
        for a = 1:nAnimals
            plot(freq, corr_monkey(a,:), 'Color', colors(a,:), 'LineWidth', 1.5);
        end
        plot(freq, corr_monkey_avg, 'k', 'LineWidth', 2);
        yline(thresh_monkey_avg, '--r', 'LineWidth', 1.5);
        legend([animals, {'Monkey avg', 'Threshold'}], 'Location', 'best');
        xlabel('Frequency (Hz)'); ylabel('Correlation');
        title(sprintf('%s Correlation: Per Animal + Monkey Average', upper(sig)));
        saveas(gcf, fullfile(monkey_plot_dir, 'monkey_avg_corr_overlay.pdf'));

        figure;
        histogram(tmax_monkey_avg, 50, 'FaceColor', [0.7 0.7 0.7]);
        xline(thresh_monkey_avg, '--r', 'LineWidth', 2);
        xlabel('Max correlation (permutation)'); ylabel('Count');
        title(sprintf('Monkey-Average Null Distribution (tmax) — %s Correlation', upper(sig)));
        legend({'Null distribution', '95th percentile'}, 'Location', 'best');
        saveas(gcf, fullfile(monkey_plot_dir, 'monkey_avg_null_distribution.pdf'));
    else
        warning('Not all animals have %s channel-average results. Run per-animal section for each animal first.', upper(sig));
    end

end  % signal type loop
