% Correlation between pre-stimulus phase and post-stimulus LFP/MUA amplitude
% H2 (abs_per_pos): group trials by stimulus position (trialinfo col 16),
% compute circ_corrcl within each position, average correlations across positions,
% then across channels and animals.
clear all; close all; clc

%% Settings
signal_types = {'lfp', 'mua'};
animals      = {'hermes', 'klecks'};
permut_n     = 1000;

%% Paths
addpath /opt/fieldtrip_github/; ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/Correlation_analysis/functions
addpath /mnt/hpc/projects/MWSampling/4Shivangi
addpath /mnt/hpc/projects/MWSampling/4Shivangi/software_folder/CircStat2012a
clc

for a = 1:numel(animals)
    animalName = animals{a};
    fprintf('\n=== Processing %s ===\n', animalName);

    output_folder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName], ...
        'phase_correlation', 'abs_per_pos', 'cp10_till_100');

    for s = 1:numel(signal_types)
        sig = signal_types{s};
        switch sig
            case 'lfp'
                amp_field = 'LFP_ERP_ampl_all';
                perm_func = @circlin_correlation_lfp_pos;
            case 'mua'
                amp_field = 'MUA_ERP_ampl_all';
                perm_func = @circlin_correlation_mua_pos;
        end
        data_load_folder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ...
            ['results_' animalName], 'multi_lin_reg', 'cp10_till_100');

        fprintf('\n--- %s: %s ---\n', animalName, upper(sig));
        cd(data_load_folder); load('ph_all_sess.mat')

        positions = unique(ph_comb.trialinfo(:,16));
        nPos      = numel(positions);
        nFreq_data = size(ph_comb.phase_all, 2);

        %% Real data — circ_corrcl per position, average across positions

        for ichan = 1:64
            ichan
            corr_pos = NaN(nPos, nFreq_data);

            for p = 1:nPos
                mask = ph_comb.trialinfo(:,16) == positions(p);
                if sum(mask) < 2, continue; end
                phase_p = ph_comb.phase_all(mask, :, ichan);
                amp_p   = ph_comb.(amp_field)(mask, ichan);
                for foi = 1:nFreq_data
                    corr_pos(p,foi) = circ_corrcl(phase_p(:,foi), amp_p);
                end
            end

            correlation = mean(corr_pos, 1, 'omitnan');
            pvalue      = nan(1, nFreq_data);   % p-values at position level not computed here

            chan_folder = fullfile(output_folder, sig, 'all_loc_difflev', num2str(ichan));
            if ~exist(chan_folder,'dir'), mkdir(chan_folder); end
            cd(chan_folder)
            save correlation correlation
            save pvalue pvalue
        end

        %% Permutation (SLURM)

        nTrials      = size(ph_comb.phase_all, 1);
        perm_indices = arrayfun(@(x) randperm(nTrials), 1:permut_n, 'UniformOutput', false);

        cfg = cell(1, 64);
        for i = 1:64
            cfg{i}.ichan        = i;
            cfg{i}.permut_n     = permut_n;
            cfg{i}.infile       = data_load_folder;
            cfg{i}.outfile      = fullfile(output_folder, sig, 'all_loc_difflev');
            cfg{i}.perm_indices = perm_indices;
            cfg{i}.positions    = positions;
        end

        slurmfun(perm_func, cfg, 'partition','8GB','stopOnError',false,'useUserPath',true);

        %% Channel-average null distribution (per animal)

        fprintf('Computing %s channel-average null for %s...\n', upper(sig), animalName);
        cd(output_folder); load('frequency.mat'); freq = frequency; nFreq = numel(freq);

        nCh = 64;
        corr_all_ch      = NaN(nCh, nFreq);
        corr_perm_all_ch = [];

        for ch = 1:nCh
            ch_folder = fullfile(output_folder, sig, 'all_loc_difflev', num2str(ch));
            if ~exist(ch_folder,'dir'), continue; end
            corr_file = fullfile(ch_folder, 'correlation.mat');
            perm_file = fullfile(ch_folder, 'corr_perm_pos.mat');
            if ~isfile(corr_file) || ~isfile(perm_file), continue; end
            load(corr_file, 'correlation');
            load(perm_file, 'corr_perm_pos');
            if any(isnan(correlation)) || any(isnan(corr_perm_pos(:))), continue; end
            corr_all_ch(ch,:) = correlation;
            corr_perm_all_ch  = cat(3, corr_perm_all_ch, corr_perm_pos);
        end

        corr_chan_avg      = mean(corr_all_ch, 1, 'omitnan');
        corr_perm_chan_avg = mean(corr_perm_all_ch, 3);
        tmax_chan_avg      = max(corr_perm_chan_avg, [], 2);
        thresh_chan_avg    = quantile(tmax_chan_avg, 0.95);

        save(fullfile(output_folder, sig, 'all_loc_difflev', 'channel_avg_results.mat'), ...
            'corr_chan_avg','corr_all_ch','corr_perm_chan_avg','tmax_chan_avg','thresh_chan_avg','freq');
        fprintf('Saved %s channel-average results for %s\n', upper(sig), animalName);

    end  % signal type loop
end  % animal loop

%% Monkey-average

nAnimals = numel(animals);
for s = 1:numel(signal_types)
    sig = signal_types{s};
    corr_monkey = []; perm_monkey = [];

    for a = 1:nAnimals
        avg_file = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animals{a}], ...
            'phase_correlation', 'abs_per_pos','cp10_till_100',sig,'all_loc_difflev','channel_avg_results.mat');
        if ~isfile(avg_file)
            warning('%s channel-average not found for %s.', upper(sig), animals{a}); continue
        end
        tmp = load(avg_file);
        corr_monkey = cat(1, corr_monkey, tmp.corr_chan_avg);
        perm_monkey = cat(3, perm_monkey, tmp.corr_perm_chan_avg);
    end

    if size(corr_monkey,1) == nAnimals
        freq = tmp.freq;
        corr_monkey_avg   = mean(corr_monkey, 1);
        perm_monkey_avg   = mean(perm_monkey, 3);
        tmax_monkey_avg   = max(perm_monkey_avg, [], 2);
        thresh_monkey_avg = quantile(tmax_monkey_avg, 0.95);

        monkey_save_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/results_combined', ...
            'phase_correlation', 'abs_per_pos','cp10_till_100',sig,'all_loc_difflev');
        if ~exist(monkey_save_dir,'dir'), mkdir(monkey_save_dir); end
        save(fullfile(monkey_save_dir,'monkey_avg_results.mat'), ...
            'corr_monkey_avg','perm_monkey_avg','tmax_monkey_avg', ...
            'thresh_monkey_avg','corr_monkey','animals','freq');
        fprintf('%s Monkey-average threshold: %.4f\n', upper(sig), thresh_monkey_avg);

        monkey_plot_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots', ...
            'phase_correlation', 'abs_per_pos','monkey_avg','cp10_till_100',sig,'all_loc_difflev');
        if ~exist(monkey_plot_dir,'dir'), mkdir(monkey_plot_dir); end

        figure; plot_sigfreq(freq, corr_monkey_avg, thresh_monkey_avg); ylim([0 0.1]);
        title(sprintf('Monkey-Avg (Abs-Per-Pos) %s Correlation', upper(sig)));
        xlabel('Frequency (Hz)'); ylabel('Correlation');
        saveas(gcf, fullfile(monkey_plot_dir,'monkey_avg_corr.pdf'));

        figure; hold on; colors = lines(nAnimals);
        for a = 1:nAnimals, plot(freq,corr_monkey(a,:),'Color',colors(a,:),'LineWidth',1.5); end
        plot(freq,corr_monkey_avg,'k','LineWidth',2);
        yline(thresh_monkey_avg,'--r','LineWidth',1.5);
        legend([animals,{'Monkey avg','Threshold'}],'Location','best');
        xlabel('Frequency (Hz)'); ylabel('Correlation');
        title(sprintf('%s Correlation: Per Animal + Monkey Avg (Abs-Per-Pos)', upper(sig)));
        saveas(gcf, fullfile(monkey_plot_dir,'monkey_avg_corr_overlay.pdf'));

        figure; histogram(tmax_monkey_avg,50,'FaceColor',[0.7 0.7 0.7]);
        xline(thresh_monkey_avg,'--r','LineWidth',2);
        xlabel('Max correlation (permutation)'); ylabel('Count');
        title(sprintf('Monkey-Avg Null — %s Correlation (Abs-Per-Pos)', upper(sig)));
        legend({'Null distribution','95th percentile'},'Location','best');
        saveas(gcf, fullfile(monkey_plot_dir,'monkey_avg_null_distribution.pdf'));
    else
        warning('Not all animals have %s correlation results.', upper(sig));
    end
end
