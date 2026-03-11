clear all
close all
clc

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

%% Settings

data_folder = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks/multi_lin_reg/cp10_till_100';
corr_folder = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks/phase_correlation/cp10_till_100';
permut_n = 1000;
nCh = 64;

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
    phase_inv(miss_idx, :) = mod(phase(miss_idx, :) + pi, 2*pi) - pi;

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
%  PLOTTING — Phase vs RT correlation
%  =====================================================================

save_root_RT = '/mnt/hpc/projects/MWSampling/4Shivangi/Plots/correlation/klecks/cp10_till_100/RT/all_loc_difflev';
if ~exist(save_root_RT, 'dir'), mkdir(save_root_RT); end

cd(corr_folder)
load('frequency.mat')
freq = frequency;

limit_max = nan(1, nCh);
Corr_chan  = false(nCh, numel(freq));

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
title('Significant Phase-RT Correlation per Channel');
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
    ylim([0 0.1]); title('All Channels Combined - Phase vs RT Correlation');
    xlabel('Frequency (Hz)'); ylabel('Correlation');
    saveas(gcf, fullfile(save_root_RT, 'combined_RT_corr.pdf'));
end

%% =====================================================================
%  PLOTTING — Phase vs Hit/Miss (POS: Phase Opposition Sum)
%  =====================================================================

save_root_HM = '/mnt/hpc/projects/MWSampling/4Shivangi/Plots/correlation/klecks/cp10_till_100/hit_miss/all_loc_difflev';
if ~exist(save_root_HM, 'dir'), mkdir(save_root_HM); end

limit_max_pos = nan(1, nCh);
POS_chan       = false(nCh, numel(freq));

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
title('Significant POS (Phase Opposition Sum) per Channel');
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
    title('All Channels Combined - POS (Phase Opposition Sum)');
    xlabel('Frequency (Hz)'); ylabel('POS');
    saveas(gcf, fullfile(save_root_HM, 'combined_hitmiss_pos.pdf'));
end

%% =====================================================================
%  PLOTTING — Phase vs Hit/Miss (ITC with inverted miss phases)
%  =====================================================================

save_root_ITC = '/mnt/hpc/projects/MWSampling/4Shivangi/Plots/correlation/klecks/cp10_till_100/hit_miss_itc/all_loc_difflev';
if ~exist(save_root_ITC, 'dir'), mkdir(save_root_ITC); end

limit_max_itc = nan(1, nCh);
ITC_chan       = false(nCh, numel(freq));

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
title('Significant ITC (Inverted Miss Phases) per Channel');
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
    ylim([0 0.15]); title('All Channels Combined - ITC (Inverted Miss Phases)');
    xlabel('Frequency (Hz)'); ylabel('ITC');
    saveas(gcf, fullfile(save_root_ITC, 'combined_hitmiss_itc.pdf'));
end
