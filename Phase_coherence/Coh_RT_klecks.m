clear all
close all
clc

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

%% Settings

data_folder = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks/multi_lin_reg/cp10_till_100';
coh_folder  = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks/phase_coherence/cp10_till_100';
permut_n = 1000;
nCh = 64;

%% =====================================================================
%  SECTION A: Phase coherence with RT (all locations and difficulty levels)
%  =====================================================================

cd(data_folder)
load('ph_all_sess.mat')

% Select hit trials only (misses have NaN RT)
hit_idx = find(ph_comb.RT_trialinfo(:,20) == 1);

output_coh_RT = fullfile(coh_folder, 'RT');

%% A1. Real data — phase_coherence(phase, RT)

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
    cfg{ichan}.outfile      = fullfile(output_coh_RT, 'all_loc_difflev');
    cfg{ichan}.perm_indices = perm_indices;
    cfg{ichan}.trial_idx    = hit_idx;
end

slurmfun(@phase_coherence_RT_perm, cfg, ...
    'partition',   '8GB', ...
    'stopOnError', false, ...
    'useUserPath', true);

%% =====================================================================
%  PLOTTING — Phase coherence with RT
%  =====================================================================

save_root_RT = '/mnt/hpc/projects/MWSampling/4Shivangi/Plots/phase_coherence/klecks/cp10_till_100/RT/all_loc_difflev';
if ~exist(save_root_RT, 'dir'), mkdir(save_root_RT); end

cd(coh_folder)
load('frequency.mat')
freq = frequency;

limit_maxc = nan(1, nCh);
limit_maxp = nan(1, nCh);
Coh_chan    = false(nCh, numel(freq));
Phase_chan  = false(nCh, numel(freq));

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
title('Significant Phase-RT Coherence per Channel');
caxis([0 1]); colorbar;
saveas(f3, fullfile(save_root_RT, 'summary_RT_coherence.pdf'));

f4 = figure;
imagesc(freq, 1:nCh, Phase_chan);
set(gca, 'YDir', 'normal');
xlabel('Frequency (Hz)'); ylabel('Channels');
title('Significant Phase-RT Phase Spectrum per Channel');
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
    title('All Channels Combined - Phase vs RT Coherence');
    saveas(f5, fullfile(save_root_RT, 'combined_RT_coherence.pdf'));

    phase_avg = nanmean(phase_all, 1);
    phase_perm_avg = nanmean(phase_perm_all, 3);
    tmaxp_all = nanmax(phase_perm_avg, [], 2);
    limit_avgp = quantile(tmaxp_all, 0.95);

    f6 = figure;
    plot_sig(freq, phase_avg, limit_avgp, 'Frequency', 'Phase spec');
    title('All Channels Combined - Phase vs RT Phase Spec');
    saveas(f6, fullfile(save_root_RT, 'combined_RT_phase.pdf'));
end

