% =====================================================================
% Circular-linear correlation: pre-stimulus phase vs. reaction time
% AND phase-opposition / inverted-miss-ITC for hit/miss
% Hypothesis H3 (abs_per_pos_diff/)
%
% Claim: each (position × difficulty bin) cell has its own phase-DV
% relationship; cells are NOT required to share a preferred phase.
%
% Recipe per cell (Way 2 across cells afterwards):
%   RT:  circ_corrcl(phase, RT)              per cell → mean
%   POS: ITC_hits + ITC_misses               per cell → mean
%   ITC: abs(mean(exp(i·phase_inverted)))    per cell → mean
%        (miss phases flipped by pi)
%
% Then arithmetic mean across channels and animals.
%
% Stimulus position read from trialinfo column 16. Difficulty
% (col 18) is binned into nDiffBins quantile bins WITHIN each position.
%
% See sampling_compare/README.md for the Way-1 / Way-2 framing.
% =====================================================================
clear all; close all; clc

%% Settings
animals   = {'hermes', 'klecks'};
permut_n  = 1000;
nCh       = 64;
nDiffBins = 4;

%% Paths
addpath /opt/fieldtrip_github/; ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/Correlation_analysis/functions
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/multiple_linear_reg/functions
addpath /mnt/hpc/projects/MWSampling/4Shivangi
addpath /mnt/hpc/projects/MWSampling/4Shivangi/software_folder/CircStat2012a
clc

for a = 1:numel(animals)
    animalName = animals{a};
    fprintf('\n=== Processing %s ===\n', animalName);

    data_folder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName], 'multi_lin_reg', 'cp10_till_100');
    corr_folder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName], 'phase_correlation', 'abs_per_pos_diff', 'cp10_till_100');

    %% =========================================================
    %  SECTION A: Phase vs RT — circ_corrcl per (pos x diff) cell
    %  =========================================================

    cd(data_folder); load('ph_all_sess.mat')

    hit_idx        = find(ph_comb.RT_trialinfo(:,20) == 1);
    positions_RT   = unique(ph_comb.RT_trialinfo(hit_idx, 16));
    nPos_RT        = numel(positions_RT);
    nFreq_data     = size(ph_comb.phase_all, 2);
    output_RT      = fullfile(corr_folder, 'RT');

    diff_bin_RT = bin_difficulty_per_pos( ...
        ph_comb.RT_trialinfo(hit_idx, 18), ...
        ph_comb.RT_trialinfo(hit_idx, 16), ...
        positions_RT, nDiffBins);

    cell_pos_RT = repmat(positions_RT(:), 1, nDiffBins); cell_pos_RT = cell_pos_RT(:);
    cell_dif_RT = repmat(1:nDiffBins, nPos_RT, 1);       cell_dif_RT = cell_dif_RT(:);
    nCell_RT    = numel(cell_pos_RT);

    %% A1. Real data

    for ichan = 1:nCh
        ichan
        corr_cell = NaN(nCell_RT, nFreq_data);

        for c = 1:nCell_RT
            cell_mask   = (ph_comb.RT_trialinfo(hit_idx, 16) == cell_pos_RT(c)) & ...
                          (diff_bin_RT == cell_dif_RT(c));
            cell_global = hit_idx(cell_mask);
            if sum(cell_mask) < 2, continue; end
            phase_c = ph_comb.phase_all(cell_global, :, ichan);
            rt_c    = ph_comb.RT(cell_global, ichan);
            valid   = ~isnan(rt_c);
            if sum(valid) < 2, continue; end
            phase_c = phase_c(valid,:); rt_c = rt_c(valid);
            for foi = 1:nFreq_data
                corr_cell(c,foi) = circ_corrcl(phase_c(:,foi), rt_c);
            end
        end

        correlation = mean(corr_cell, 1, 'omitnan');
        pvalue      = nan(1, nFreq_data);

        chan_folder = fullfile(output_RT, 'all_loc_difflev', num2str(ichan));
        if ~exist(chan_folder,'dir'), mkdir(chan_folder); end
        cd(chan_folder); save correlation correlation; save pvalue pvalue
        save corr_cell corr_cell cell_pos_RT cell_dif_RT
    end

    %% A2. Permutation (SLURM)

    cd(data_folder); load('ph_all_sess.mat')
    nTrials      = length(hit_idx);
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
        cfg{ichan}.cell_pos     = cell_pos_RT;
        cfg{ichan}.cell_dif     = cell_dif_RT;
        cfg{ichan}.diff_bin     = diff_bin_RT;
    end

    slurmfun(@circlin_corr_RT_perm_per_pos_diff, cfg, 'partition','8GB','stopOnError',false,'useUserPath',true);

    %% =========================================================
    %  SECTION B: Phase vs Hit/Miss — POS and ITC per (pos x diff) cell
    %  =========================================================

    cd(data_folder); load('ph_all_sess.mat')

    all_idx      = find(ph_comb.trialinfo(:,20) == 1 | ph_comb.trialinfo(:,20) == 5);
    hit_labels   = (ph_comb.trialinfo(all_idx, 20) == 1);
    positions_HM = unique(ph_comb.trialinfo(all_idx, 16));
    nPos_HM      = numel(positions_HM);
    output_HM    = fullfile(corr_folder, 'hit_miss');

    diff_bin_HM = bin_difficulty_per_pos( ...
        ph_comb.trialinfo(all_idx, 18), ...
        ph_comb.trialinfo(all_idx, 16), ...
        positions_HM, nDiffBins);

    cell_pos_HM = repmat(positions_HM(:), 1, nDiffBins); cell_pos_HM = cell_pos_HM(:);
    cell_dif_HM = repmat(1:nDiffBins, nPos_HM, 1);       cell_dif_HM = cell_dif_HM(:);
    nCell_HM    = numel(cell_pos_HM);

    %% B1. Real data — POS and ITC per (pos x diff) cell

    for ichan = 1:nCh
        ichan
        pos_cell = NaN(nCell_HM, nFreq_data);
        itc_cell = NaN(nCell_HM, nFreq_data);

        for c = 1:nCell_HM
            mask     = (ph_comb.trialinfo(all_idx, 16) == cell_pos_HM(c)) & ...
                       (diff_bin_HM == cell_dif_HM(c));
            if sum(mask) < 2, continue; end
            ph_c    = ph_comb.phase_all(all_idx(mask), :, ichan);
            hits_c  = hit_labels(mask);
            misses_c = ~hits_c;
            if sum(hits_c) < 1 || sum(misses_c) < 1, continue; end

            for foi = 1:nFreq_data
                itc_h = abs(mean(exp(1i * ph_c(hits_c,   foi))));
                itc_m = abs(mean(exp(1i * ph_c(misses_c, foi))));
                pos_cell(c,foi) = itc_h + itc_m;
            end

            phase_inv_c = ph_c;
            phase_inv_c(misses_c, :) = mod(ph_c(misses_c,:) + pi, 2*pi) - pi;
            for foi = 1:nFreq_data
                itc_cell(c,foi) = abs(mean(exp(1i * phase_inv_c(:,foi))));
            end
        end

        pos = mean(pos_cell, 1, 'omitnan');
        itc = mean(itc_cell, 1, 'omitnan');

        chan_folder = fullfile(output_HM, 'all_loc_difflev', num2str(ichan));
        if ~exist(chan_folder,'dir'), mkdir(chan_folder); end
        cd(chan_folder)
        save pos pos pos_cell cell_pos_HM cell_dif_HM
        save itc itc itc_cell
    end

    %% B2. Permutation (SLURM)

    cd(data_folder); load('ph_all_sess.mat')
    nTrials      = length(all_idx);
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
        cfg{ichan}.cell_pos     = cell_pos_HM;
        cfg{ichan}.cell_dif     = cell_dif_HM;
        cfg{ichan}.diff_bin     = diff_bin_HM;
    end

    slurmfun(@circlin_corr_hitmiss_perm_per_pos_diff, cfg, 'partition','8GB','stopOnError',false,'useUserPath',true);

    %% =========================================================
    %  Channel-average null — RT
    %  =========================================================

    fprintf('Computing RT channel-average null for %s...\n', animalName);
    cd(corr_folder); load('frequency.mat'); freq = frequency; nFreq = numel(freq);

    corr_all_ch = NaN(nCh, nFreq); corr_perm_all_ch = [];

    for ch = 1:nCh
        ch_folder = fullfile(output_RT, 'all_loc_difflev', num2str(ch));
        if ~exist(ch_folder,'dir'), continue; end
        corr_file = fullfile(ch_folder,'correlation.mat');
        perm_file = fullfile(ch_folder,'corr_perm_pos_diff.mat');
        if ~isfile(corr_file) || ~isfile(perm_file), continue; end
        load(corr_file,'correlation'); load(perm_file,'corr_perm_pos_diff');
        if any(isnan(correlation)) || any(isnan(corr_perm_pos_diff(:))), continue; end
        corr_all_ch(ch,:) = correlation;
        corr_perm_all_ch  = cat(3, corr_perm_all_ch, corr_perm_pos_diff);
    end

    corr_chan_avg      = mean(corr_all_ch, 1, 'omitnan');
    corr_perm_chan_avg = mean(corr_perm_all_ch, 3);
    tmax_chan_avg      = max(corr_perm_chan_avg, [], 2);
    thresh_chan_avg    = quantile(tmax_chan_avg, 0.95);
    save(fullfile(output_RT,'all_loc_difflev','channel_avg_results.mat'), ...
        'corr_chan_avg','corr_all_ch','corr_perm_chan_avg','tmax_chan_avg','thresh_chan_avg','freq');
    fprintf('Saved RT channel-average results for %s\n', animalName);

    %% =========================================================
    %  Channel-average null — POS
    %  =========================================================

    fprintf('Computing POS channel-average null for %s...\n', animalName);
    pos_all_ch = NaN(nCh, nFreq); pos_perm_all_ch = [];

    for ch = 1:nCh
        ch_folder = fullfile(output_HM, 'all_loc_difflev', num2str(ch));
        if ~exist(ch_folder,'dir'), continue; end
        pos_file  = fullfile(ch_folder,'pos.mat');
        perm_file = fullfile(ch_folder,'pos_perm_pos_diff.mat');
        if ~isfile(pos_file) || ~isfile(perm_file), continue; end
        load(pos_file,'pos'); load(perm_file,'pos_perm_pos_diff');
        if any(isnan(pos)) || any(isnan(pos_perm_pos_diff(:))), continue; end
        pos_all_ch(ch,:) = pos;
        pos_perm_all_ch  = cat(3, pos_perm_all_ch, pos_perm_pos_diff);
    end

    pos_chan_avg      = mean(pos_all_ch, 1, 'omitnan');
    pos_perm_chan_avg = mean(pos_perm_all_ch, 3);
    tmax_chan_avg_pos = max(pos_perm_chan_avg, [], 2);
    thresh_chan_avg_pos = quantile(tmax_chan_avg_pos, 0.95);
    save(fullfile(output_HM,'all_loc_difflev','channel_avg_results_pos.mat'), ...
        'pos_chan_avg','pos_all_ch','pos_perm_chan_avg','tmax_chan_avg_pos','thresh_chan_avg_pos','freq');
    fprintf('Saved POS channel-average results for %s\n', animalName);

    %% =========================================================
    %  Channel-average null — ITC
    %  =========================================================

    fprintf('Computing ITC channel-average null for %s...\n', animalName);
    itc_all_ch = NaN(nCh, nFreq); itc_perm_all_ch = [];

    for ch = 1:nCh
        ch_folder = fullfile(output_HM, 'all_loc_difflev', num2str(ch));
        if ~exist(ch_folder,'dir'), continue; end
        itc_file  = fullfile(ch_folder,'itc.mat');
        perm_file = fullfile(ch_folder,'itc_perm_pos_diff.mat');
        if ~isfile(itc_file) || ~isfile(perm_file), continue; end
        load(itc_file,'itc'); load(perm_file,'itc_perm_pos_diff');
        if any(isnan(itc)) || any(isnan(itc_perm_pos_diff(:))), continue; end
        itc_all_ch(ch,:) = itc;
        itc_perm_all_ch  = cat(3, itc_perm_all_ch, itc_perm_pos_diff);
    end

    itc_chan_avg      = mean(itc_all_ch, 1, 'omitnan');
    itc_perm_chan_avg = mean(itc_perm_all_ch, 3);
    tmax_chan_avg_itc = max(itc_perm_chan_avg, [], 2);
    thresh_chan_avg_itc = quantile(tmax_chan_avg_itc, 0.95);
    save(fullfile(output_HM,'all_loc_difflev','channel_avg_results_itc.mat'), ...
        'itc_chan_avg','itc_all_ch','itc_perm_chan_avg','tmax_chan_avg_itc','thresh_chan_avg_itc','freq');
    fprintf('Saved ITC channel-average results for %s\n', animalName);

    close all
end  % animal loop

%% =========================================================
%  Monkey-average — RT
%  =========================================================

nAnimals = numel(animals);
corr_monkey = []; perm_monkey = [];

for a = 1:nAnimals
    avg_file = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animals{a}], ...
        'phase_correlation', 'abs_per_pos_diff','cp10_till_100','RT','all_loc_difflev','channel_avg_results.mat');
    if ~isfile(avg_file), warning('RT results not found for %s.', animals{a}); continue; end
    tmp = load(avg_file);
    corr_monkey = cat(1,corr_monkey,tmp.corr_chan_avg);
    perm_monkey = cat(3,perm_monkey,tmp.corr_perm_chan_avg);
end

if size(corr_monkey,1) == nAnimals
    freq = tmp.freq;
    corr_monkey_avg   = mean(corr_monkey,1);
    perm_monkey_avg   = mean(perm_monkey,3);
    tmax_monkey_avg   = max(perm_monkey_avg,[],2);
    thresh_monkey_avg = quantile(tmax_monkey_avg,0.95);

    monkey_save_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/results_combined', ...
        'phase_correlation', 'abs_per_pos_diff','cp10_till_100','RT','all_loc_difflev');
    if ~exist(monkey_save_dir,'dir'), mkdir(monkey_save_dir); end
    save(fullfile(monkey_save_dir,'monkey_avg_results.mat'), ...
        'corr_monkey_avg','perm_monkey_avg','tmax_monkey_avg','thresh_monkey_avg','corr_monkey','animals','freq');
    fprintf('RT Monkey-average threshold: %.4f\n', thresh_monkey_avg);

    monkey_plot_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots', ...
        'phase_correlation', 'abs_per_pos_diff','monkey_avg','cp10_till_100','RT','all_loc_difflev');
    if ~exist(monkey_plot_dir,'dir'), mkdir(monkey_plot_dir); end

    figure; plot_sigfreq(freq,corr_monkey_avg,thresh_monkey_avg); ylim([0 0.1]);
    title('Monkey-Avg (Abs-Per-Pos-Diff) Phase-RT Correlation');
    xlabel('Frequency (Hz)'); ylabel('Correlation');
    saveas(gcf,fullfile(monkey_plot_dir,'monkey_avg_RT_corr.pdf'));

    figure; hold on; colors = lines(nAnimals);
    for a=1:nAnimals, plot(freq,corr_monkey(a,:),'Color',colors(a,:),'LineWidth',1.5); end
    plot(freq,corr_monkey_avg,'k','LineWidth',2); yline(thresh_monkey_avg,'--r','LineWidth',1.5);
    legend([animals,{'Monkey avg','Threshold'}],'Location','best');
    xlabel('Frequency (Hz)'); ylabel('Correlation');
    title('Phase-RT Correlation: Per Animal + Monkey Avg (Abs-Per-Pos-Diff)');
    saveas(gcf,fullfile(monkey_plot_dir,'monkey_avg_RT_corr_overlay.pdf'));
end

%% Monkey-average — POS

pos_monkey = []; perm_monkey_pos = [];

for a = 1:nAnimals
    avg_file = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animals{a}], ...
        'phase_correlation', 'abs_per_pos_diff','cp10_till_100','hit_miss','all_loc_difflev','channel_avg_results_pos.mat');
    if ~isfile(avg_file), warning('POS results not found for %s.', animals{a}); continue; end
    tmp = load(avg_file);
    pos_monkey      = cat(1,pos_monkey,tmp.pos_chan_avg);
    perm_monkey_pos = cat(3,perm_monkey_pos,tmp.pos_perm_chan_avg);
end

if size(pos_monkey,1) == nAnimals
    freq = tmp.freq;
    pos_monkey_avg      = mean(pos_monkey,1);
    perm_monkey_avg_pos = mean(perm_monkey_pos,3);
    tmax_monkey_avg_pos = max(perm_monkey_avg_pos,[],2);
    thresh_monkey_avg_pos = quantile(tmax_monkey_avg_pos,0.95);

    monkey_save_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/results_combined', ...
        'phase_correlation', 'abs_per_pos_diff','cp10_till_100','hit_miss','all_loc_difflev');
    if ~exist(monkey_save_dir,'dir'), mkdir(monkey_save_dir); end
    save(fullfile(monkey_save_dir,'monkey_avg_results_pos.mat'), ...
        'pos_monkey_avg','perm_monkey_avg_pos','tmax_monkey_avg_pos','thresh_monkey_avg_pos','pos_monkey','animals','freq');
    fprintf('POS Monkey-average threshold: %.4f\n', thresh_monkey_avg_pos);

    monkey_plot_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots', ...
        'phase_correlation', 'abs_per_pos_diff','monkey_avg','cp10_till_100','hit_miss','all_loc_difflev');
    if ~exist(monkey_plot_dir,'dir'), mkdir(monkey_plot_dir); end

    figure; plot_sigfreq(freq,pos_monkey_avg,thresh_monkey_avg_pos);
    title('Monkey-Avg (Abs-Per-Pos-Diff) POS'); xlabel('Frequency (Hz)'); ylabel('POS');
    saveas(gcf,fullfile(monkey_plot_dir,'monkey_avg_hitmiss_pos.pdf'));
end

%% Monkey-average — ITC

itc_monkey = []; perm_monkey_itc = [];

for a = 1:nAnimals
    avg_file = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animals{a}], ...
        'phase_correlation', 'abs_per_pos_diff','cp10_till_100','hit_miss','all_loc_difflev','channel_avg_results_itc.mat');
    if ~isfile(avg_file), warning('ITC results not found for %s.', animals{a}); continue; end
    tmp = load(avg_file);
    itc_monkey      = cat(1,itc_monkey,tmp.itc_chan_avg);
    perm_monkey_itc = cat(3,perm_monkey_itc,tmp.itc_perm_chan_avg);
end

if size(itc_monkey,1) == nAnimals
    freq = tmp.freq;
    itc_monkey_avg      = mean(itc_monkey,1);
    perm_monkey_avg_itc = mean(perm_monkey_itc,3);
    tmax_monkey_avg_itc = max(perm_monkey_avg_itc,[],2);
    thresh_monkey_avg_itc = quantile(tmax_monkey_avg_itc,0.95);

    monkey_save_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/results_combined', ...
        'phase_correlation', 'abs_per_pos_diff','cp10_till_100','hit_miss_itc','all_loc_difflev');
    if ~exist(monkey_save_dir,'dir'), mkdir(monkey_save_dir); end
    save(fullfile(monkey_save_dir,'monkey_avg_results_itc.mat'), ...
        'itc_monkey_avg','perm_monkey_avg_itc','tmax_monkey_avg_itc','thresh_monkey_avg_itc','itc_monkey','animals','freq');
    fprintf('ITC Monkey-average threshold: %.4f\n', thresh_monkey_avg_itc);

    monkey_plot_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots', ...
        'phase_correlation', 'abs_per_pos_diff','monkey_avg','cp10_till_100','hit_miss_itc','all_loc_difflev');
    if ~exist(monkey_plot_dir,'dir'), mkdir(monkey_plot_dir); end

    figure; plot_sigfreq(freq,itc_monkey_avg,thresh_monkey_avg_itc); ylim([0 0.15]);
    title('Monkey-Avg (Abs-Per-Pos-Diff) ITC (Inverted Miss Phases)');
    xlabel('Frequency (Hz)'); ylabel('ITC');
    saveas(gcf,fullfile(monkey_plot_dir,'monkey_avg_hitmiss_itc.pdf'));

    figure; hold on; colors = lines(nAnimals);
    for a=1:nAnimals, plot(freq,itc_monkey(a,:),'Color',colors(a,:),'LineWidth',1.5); end
    plot(freq,itc_monkey_avg,'k','LineWidth',2); yline(thresh_monkey_avg_itc,'--r','LineWidth',1.5);
    legend([animals,{'Monkey avg','Threshold'}],'Location','best');
    xlabel('Frequency (Hz)'); ylabel('ITC');
    title('ITC (Abs-Per-Pos-Diff): Per Animal + Monkey Average');
    saveas(gcf,fullfile(monkey_plot_dir,'monkey_avg_hitmiss_itc_overlay.pdf'));

    figure; histogram(tmax_monkey_avg_itc,50,'FaceColor',[0.7 0.7 0.7]);
    xline(thresh_monkey_avg_itc,'--r','LineWidth',2);
    xlabel('Max ITC (permutation)'); ylabel('Count');
    title('Monkey-Avg Null — ITC (Abs-Per-Pos-Diff)');
    legend({'Null distribution','95th percentile'},'Location','best');
    saveas(gcf,fullfile(monkey_plot_dir,'monkey_avg_hitmiss_itc_null_distribution.pdf'));
end
