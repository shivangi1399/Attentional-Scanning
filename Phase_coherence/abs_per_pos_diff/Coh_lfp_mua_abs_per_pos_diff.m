% Coherence between pre-stimulus phase and post-stimulus LFP/MUA amplitude
% H3 (abs_per_pos_diff): group trials by (stimulus position x difficulty bin)
% cells. Difficulty (trialinfo col 18) is binned into nDiffBins quantile bins
% within each position. Compute complex coherence within each cell, abs()
% per cell, average magnitudes across cells, then across channels and animals.
clear all; close all; clc

%% Settings
signal_types = {'lfp', 'mua'};
animals      = {'hermes', 'klecks'};
permut_n     = 1000;
nDiffBins    = 4;

%% Paths
addpath /opt/fieldtrip_github/; ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/coherence_analysis
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/Phase_coherence/functions
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/multiple_linear_reg/functions
addpath /mnt/hpc/projects/MWSampling/4Shivangi
clc

for a = 1:numel(animals)
    animalName = animals{a};
    fprintf('\n=== Processing %s ===\n', animalName);

    output_folder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName], ...
        'phase_coherence', 'abs_per_pos_diff', 'cp10_till_100');

    for s = 1:numel(signal_types)
        sig = signal_types{s};
        switch sig
            case 'lfp'
                amp_field = 'LFP_ERP_ampl_all';
                perm_func = @phase_coherence_perm_lfp_per_pos_diff;
            case 'mua'
                amp_field = 'MUA_ERP_ampl_all';
                perm_func = @phase_coherence_perm_mua_per_pos_diff;
        end
        data_load_folder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ...
            ['results_' animalName], 'multi_lin_reg', 'cp10_till_100');

        fprintf('\n--- %s: %s ---\n', animalName, upper(sig));
        cd(data_load_folder); load('ph_all_sess.mat')

        positions  = unique(ph_comb.trialinfo(:,16));
        nPos       = numel(positions);
        nCh        = 64;
        trial_idx  = 1:size(ph_comb.trialinfo, 1);

        % Per-trial difficulty bin (within-position quantile bins of col 18)
        diff_bin = bin_difficulty_per_pos(ph_comb.trialinfo(:,18), ...
                                          ph_comb.trialinfo(:,16), positions, nDiffBins);

        % Enumerate (position x difficulty) cells
        cell_pos = repmat(positions(:), 1, nDiffBins); cell_pos = cell_pos(:);
        cell_dif = repmat(1:nDiffBins, nPos, 1);       cell_dif = cell_dif(:);
        nCell    = numel(cell_pos);

        %% Real data — abs per (position x difficulty) cell, average across cells

        for ichan = 1:nCh
            ichan
            coh_cell = NaN(nCell, size(ph_comb.phase_all, 2));

            for c = 1:nCell
                mask = (ph_comb.trialinfo(:,16) == cell_pos(c)) & (diff_bin == cell_dif(c));
                if sum(mask) < 2, continue; end
                phase_c = ph_comb.phase_all(mask, :, ichan);
                amp_c   = ph_comb.(amp_field)(mask, ichan);
                [~, ~, coh_complex_c] = phase_coherence(phase_c, amp_c);
                coh_cell(c,:) = abs(coh_complex_c);
            end

            coh = mean(coh_cell, 1, 'omitnan');

            chan_folder = fullfile(output_folder, sig, 'all_loc_difflev', num2str(ichan));
            if ~exist(chan_folder,'dir'), mkdir(chan_folder); end
            save(fullfile(chan_folder,'coherence.mat'), 'coh', 'coh_cell', 'cell_pos', 'cell_dif');
        end

        %% Permutation (SLURM) — same per-cell abs logic

        nTrials      = length(trial_idx);
        perm_indices = arrayfun(@(x) randperm(nTrials), 1:permut_n, 'UniformOutput', false);

        cfg = cell(1, nCh);
        for ichan = 1:nCh
            cfg{ichan}.ichan        = ichan;
            cfg{ichan}.permut_n     = permut_n;
            cfg{ichan}.infile       = data_load_folder;
            cfg{ichan}.outfile      = fullfile(output_folder, sig, 'all_loc_difflev');
            cfg{ichan}.perm_indices = perm_indices;
            cfg{ichan}.trial_idx    = trial_idx;
            cfg{ichan}.cell_pos     = cell_pos;
            cfg{ichan}.cell_dif     = cell_dif;
            cfg{ichan}.diff_bin     = diff_bin;
        end

        slurmfun(perm_func, cfg, 'partition','8GB','stopOnError',false,'useUserPath',true);

        %% Channel-average null distribution (per animal)

        fprintf('Computing %s channel-average null for %s...\n', upper(sig), animalName);
        cd(output_folder); load('frequency.mat'); freq = frequency; nFreq = numel(freq);

        coh_all_ch      = NaN(nCh, nFreq);
        coh_perm_all_ch = [];

        for ch = 1:nCh
            ch_folder = fullfile(output_folder, sig, 'all_loc_difflev', num2str(ch));
            if ~exist(ch_folder,'dir'), continue; end
            coh_file  = fullfile(ch_folder, 'coherence.mat');
            perm_file = fullfile(ch_folder, 'coh_perm_pos_diff.mat');
            if ~isfile(coh_file) || ~isfile(perm_file), continue; end
            load(coh_file,  'coh');
            load(perm_file, 'coh_perm_pos_diff');
            if any(isnan(coh)) || any(isnan(coh_perm_pos_diff(:))), continue; end
            coh_all_ch(ch,:) = coh;
            coh_perm_all_ch  = cat(3, coh_perm_all_ch, coh_perm_pos_diff);
        end

        coh_chan_avg      = mean(coh_all_ch, 1, 'omitnan');
        coh_perm_chan_avg = mean(coh_perm_all_ch, 3);
        tmax_chan_avg     = max(coh_perm_chan_avg, [], 2);
        thresh_chan_avg   = quantile(tmax_chan_avg, 0.95);

        save(fullfile(output_folder, sig, 'all_loc_difflev', 'channel_avg_results.mat'), ...
            'coh_chan_avg','coh_all_ch','coh_perm_chan_avg','tmax_chan_avg','thresh_chan_avg','freq');
        fprintf('Saved %s channel-average results for %s\n', upper(sig), animalName);

    end  % signal type loop
end  % animal loop

%% Monkey-average

nAnimals = numel(animals);
for s = 1:numel(signal_types)
    sig = signal_types{s};
    coh_monkey = []; perm_monkey = [];

    for a = 1:nAnimals
        avg_file = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animals{a}], ...
            'phase_coherence', 'abs_per_pos_diff','cp10_till_100',sig,'all_loc_difflev','channel_avg_results.mat');
        if ~isfile(avg_file)
            warning('%s channel-average not found for %s.', upper(sig), animals{a}); continue
        end
        tmp = load(avg_file);
        coh_monkey  = cat(1, coh_monkey,  tmp.coh_chan_avg);
        perm_monkey = cat(3, perm_monkey, tmp.coh_perm_chan_avg);
    end

    if size(coh_monkey,1) == nAnimals
        freq = tmp.freq;
        coh_monkey_avg    = mean(coh_monkey, 1);
        perm_monkey_avg   = mean(perm_monkey, 3);
        tmax_monkey_avg   = max(perm_monkey_avg, [], 2);
        thresh_monkey_avg = quantile(tmax_monkey_avg, 0.95);

        monkey_save_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/results_combined', ...
            'phase_coherence', 'abs_per_pos_diff','cp10_till_100',sig,'all_loc_difflev');
        if ~exist(monkey_save_dir,'dir'), mkdir(monkey_save_dir); end
        save(fullfile(monkey_save_dir,'monkey_avg_results.mat'), ...
            'coh_monkey_avg','perm_monkey_avg','tmax_monkey_avg', ...
            'thresh_monkey_avg','coh_monkey','animals','freq');
        fprintf('%s Monkey-average threshold: %.4f\n', upper(sig), thresh_monkey_avg);

        monkey_plot_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots', ...
            'phase_coherence', 'abs_per_pos_diff','monkey_avg','cp10_till_100',sig,'all_loc_difflev');
        if ~exist(monkey_plot_dir,'dir'), mkdir(monkey_plot_dir); end

        figure; plot_sigfreq(freq, coh_monkey_avg, thresh_monkey_avg); ylim([0 0.1]);
        title(sprintf('Monkey-Avg (Abs-Per-Pos-Diff) %s Coherence', upper(sig)));
        xlabel('Frequency (Hz)'); ylabel('Coherence');
        saveas(gcf, fullfile(monkey_plot_dir,'monkey_avg_coherence.pdf'));

        figure; hold on; colors = lines(nAnimals);
        for a = 1:nAnimals, plot(freq,coh_monkey(a,:),'Color',colors(a,:),'LineWidth',1.5); end
        plot(freq,coh_monkey_avg,'k','LineWidth',2);
        yline(thresh_monkey_avg,'--r','LineWidth',1.5);
        legend([animals,{'Monkey avg','Threshold'}],'Location','best');
        xlabel('Frequency (Hz)'); ylabel('Coherence');
        title(sprintf('%s Coherence: Per Animal + Monkey Avg (Abs-Per-Pos-Diff)', upper(sig)));
        saveas(gcf, fullfile(monkey_plot_dir,'monkey_avg_coherence_overlay.pdf'));

        figure; histogram(tmax_monkey_avg,50,'FaceColor',[0.7 0.7 0.7]);
        xline(thresh_monkey_avg,'--r','LineWidth',2);
        xlabel('Max coherence (permutation)'); ylabel('Count');
        title(sprintf('Monkey-Avg Null — %s Coherence (Abs-Per-Pos-Diff)', upper(sig)));
        legend({'Null distribution','95th percentile'},'Location','best');
        saveas(gcf, fullfile(monkey_plot_dir,'monkey_avg_null_distribution.pdf'));
    else
        warning('Not all animals have %s channel-average results.', upper(sig));
    end
end
