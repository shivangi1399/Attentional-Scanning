% =====================================================================
% Phase coherence: pre-stimulus phase vs. reaction time (hit trials only)
% Hypothesis H2 (abs_per_pos/)
%
% Claim: each stimulus position has its own phase-DV relationship;
% positions are NOT required to share a preferred phase.
%
% Recipe: per-position complex resultant; abs() per position
% (Way 2 across positions); arithmetic mean of magnitudes across
% positions, channels, and animals.
%
% Stimulus position read from RT_trialinfo column 16.
%
% See sampling_compare/README.md for the Way-1 / Way-2 framing.
% =====================================================================
clear all; close all; clc

%% Settings
animals  = {'hermes', 'klecks'};
permut_n = 1000;
nCh      = 64;

%% Paths
addpath /opt/fieldtrip_github/; ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/Phase_coherence/functions
addpath /mnt/hpc/projects/MWSampling/4Shivangi
clc

for a = 1:numel(animals)
    animalName = animals{a};
    fprintf('\n=== Processing %s ===\n', animalName);

    data_folder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName], 'multi_lin_reg', 'cp10_till_100');
    output_folder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName], ...
        'phase_coherence', 'abs_per_pos', 'cp10_till_100');
    output_coh_RT = fullfile(output_folder, 'RT');

    cd(data_folder); load('ph_all_sess.mat')

    hit_idx   = find(ph_comb.RT_trialinfo(:,20) == 1);
    positions = unique(ph_comb.RT_trialinfo(hit_idx, 16));
    nPos      = numel(positions);

    %% Real data — abs per stimulus position, average across positions

    for ichan = 1:nCh
        ichan
        coh_pos = NaN(nPos, size(ph_comb.phase_all, 2));

        for p = 1:nPos
            pos_mask   = ph_comb.RT_trialinfo(hit_idx, 16) == positions(p);
            pos_global = hit_idx(pos_mask);
            if sum(pos_mask) < 2, continue; end
            phase_p = ph_comb.phase_all(pos_global, :, ichan);
            rt_p    = ph_comb.RT(pos_global, ichan);
            valid   = ~isnan(rt_p);
            if sum(valid) < 2, continue; end
            [~, ~, coh_complex_p] = phase_coherence(phase_p(valid,:), rt_p(valid));
            coh_pos(p,:) = abs(coh_complex_p);
        end

        coh = mean(coh_pos, 1, 'omitnan');

        chan_folder = fullfile(output_coh_RT, 'all_loc_difflev', num2str(ichan));
        if ~exist(chan_folder,'dir'), mkdir(chan_folder); end
        save(fullfile(chan_folder,'coherence.mat'), 'coh', 'coh_pos', 'positions');
    end

    %% Permutation (SLURM)

    nTrials      = length(hit_idx);
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
        cfg{ichan}.positions    = positions;
    end

    slurmfun(@phase_coherence_RT_perm_pos, cfg, 'partition','8GB','stopOnError',false,'useUserPath',true);

    %% Channel-average null distribution (per animal)

    fprintf('Computing RT channel-average null for %s...\n', animalName);

    cd(output_folder); load('frequency.mat'); freq = frequency; nFreq = numel(freq);

    coh_all_ch      = NaN(nCh, nFreq);
    coh_perm_all_ch = [];

    for ch = 1:nCh
        ch_folder = fullfile(output_coh_RT, 'all_loc_difflev', num2str(ch));
        if ~exist(ch_folder,'dir'), continue; end
        coh_file  = fullfile(ch_folder, 'coherence.mat');
        perm_file = fullfile(ch_folder, 'coh_perm_pos.mat');
        if ~isfile(coh_file) || ~isfile(perm_file), continue; end
        load(coh_file,  'coh');
        load(perm_file, 'coh_perm_pos');
        if any(isnan(coh)) || any(isnan(coh_perm_pos(:))), continue; end
        coh_all_ch(ch,:) = coh;
        coh_perm_all_ch  = cat(3, coh_perm_all_ch, coh_perm_pos);
    end

    coh_chan_avg      = mean(coh_all_ch, 1, 'omitnan');
    coh_perm_chan_avg = mean(coh_perm_all_ch, 3);
    tmax_chan_avg     = max(coh_perm_chan_avg, [], 2);
    thresh_chan_avg   = quantile(tmax_chan_avg, 0.95);

    save(fullfile(output_coh_RT,'all_loc_difflev','channel_avg_results.mat'), ...
        'coh_chan_avg','coh_all_ch','coh_perm_chan_avg','tmax_chan_avg','thresh_chan_avg','freq');
    fprintf('Saved RT channel-average results for %s\n', animalName);

    close all
end  % animal loop

%% Monkey-average

nAnimals = numel(animals);
coh_monkey = []; perm_monkey = [];

for a = 1:nAnimals
    avg_file = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animals{a}], ...
        'phase_coherence', 'abs_per_pos','cp10_till_100','RT','all_loc_difflev','channel_avg_results.mat');
    if ~isfile(avg_file)
        warning('RT channel-average not found for %s.', animals{a}); continue
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
        'phase_coherence', 'abs_per_pos','cp10_till_100','RT','all_loc_difflev');
    if ~exist(monkey_save_dir,'dir'), mkdir(monkey_save_dir); end
    save(fullfile(monkey_save_dir,'monkey_avg_results.mat'), ...
        'coh_monkey_avg','perm_monkey_avg','tmax_monkey_avg', ...
        'thresh_monkey_avg','coh_monkey','animals','freq');
    fprintf('RT Monkey-average threshold: %.4f\n', thresh_monkey_avg);

    monkey_plot_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots', ...
        'phase_coherence', 'abs_per_pos','monkey_avg','cp10_till_100','RT','all_loc_difflev');
    if ~exist(monkey_plot_dir,'dir'), mkdir(monkey_plot_dir); end

    figure; plot_sigfreq(freq, coh_monkey_avg, thresh_monkey_avg); ylim([0 0.1]);
    title('Monkey-Avg (Abs-Per-Pos) Phase-RT Coherence');
    xlabel('Frequency (Hz)'); ylabel('Coherence');
    saveas(gcf, fullfile(monkey_plot_dir,'monkey_avg_RT_coherence.pdf'));

    figure; hold on; colors = lines(nAnimals);
    for a = 1:nAnimals, plot(freq,coh_monkey(a,:),'Color',colors(a,:),'LineWidth',1.5); end
    plot(freq,coh_monkey_avg,'k','LineWidth',2);
    yline(thresh_monkey_avg,'--r','LineWidth',1.5);
    legend([animals,{'Monkey avg','Threshold'}],'Location','best');
    xlabel('Frequency (Hz)'); ylabel('Coherence');
    title('Phase-RT Coherence: Per Animal + Monkey Avg (Abs-Per-Pos)');
    saveas(gcf, fullfile(monkey_plot_dir,'monkey_avg_RT_coherence_overlay.pdf'));

    figure; histogram(tmax_monkey_avg,50,'FaceColor',[0.7 0.7 0.7]);
    xline(thresh_monkey_avg,'--r','LineWidth',2);
    xlabel('Max coherence (permutation)'); ylabel('Count');
    title('Monkey-Avg Null — Phase-RT Coherence (Abs-Per-Pos)');
    legend({'Null distribution','95th percentile'},'Location','best');
    saveas(gcf, fullfile(monkey_plot_dir,'monkey_avg_RT_null_distribution.pdf'));
else
    warning('Not all animals have RT channel-average results.');
end
