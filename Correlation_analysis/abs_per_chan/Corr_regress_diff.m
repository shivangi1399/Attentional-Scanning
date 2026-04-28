% Corr_regression.m
% Circular-linear correlation between phase and ERP amplitude residuals
% (after regressing out difficulty level) across all channels.
% Includes permutation testing (SLURM) and plotting of per-channel and
% combined results.

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
addpath /mnt/hpc/projects/MWSampling/4Shivangi
addpath /mnt/hpc/projects/MWSampling/4Shivangi/software_folder/CircStat2012a
clc

%% Create data paths
animalName = 'hermes';  % Change this to switch animal (e.g. 'klecks')

datafolder   = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName]);

cd(datafolder)
temp = dir;
session_names = [];
ii = 0;
for i = 1:length(temp)
    if contains(temp(i).name, animalName)
        ii = ii+1;
        session_names{ii,1} = temp(i).name;
    end
end

session_paths_files = cellfun(@(x) fullfile(datafolder,x,'clean_data.mat'), session_names, 'uniform',0);
phase_paths = cellfun(@(x) fullfile(datafolder,x,'Phase_analysis/hit_miss'), session_names, 'uniform',0);
output_folder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName], 'phase_correlation', 'abs_per_chan', 'cp10_till_100');
if ~exist(output_folder,'dir'), mkdir(output_folder); end
data_load_folder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName], 'multi_lin_reg', 'cp10_till_100');

%% Correlation for all locations and difficulty levels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(data_load_folder)
load('ph_all_sess.mat')

nChan = 64;
nFreq = length(ph_comb.phase_all(1,:,1));

hits2use = []; hits2use = find(ph_comb.ERP_trialinfo(:,20)==1);
miss2use = []; miss2use = find(ph_comb.ERP_trialinfo(:,20)==5);

ERP_hits = []; ERP_hits = ph_comb.ERP_ampl_all(hits2use,:);
DE_hits = []; DE_hits = ph_comb.ERP_trialinfo(hits2use,18);
ERP_miss = []; ERP_miss = ph_comb.ERP_ampl_all(miss2use,:);
DE_miss = []; DE_miss = ph_comb.ERP_trialinfo(miss2use,18);

for ichan = 1:nChan
    ichan
    
    % hits reg
    tbl = [];
    tbl = table(DE_hits,ERP_hits(:,ichan));
    model_hits  = fitlm(tbl);
    
    % misses reg
    tbl = [];
    tbl = table(DE_miss,ERP_miss(:,ichan));
    model_miss  = fitlm(tbl);
    
    %% correlation (real)
    
    temp = [];
    temp = table2array(model_hits.Residuals(:,1));
    hits2use = find(~isnan(temp));
    ampl_resdl = temp(hits2use);
    
    temp = [];
    temp = table2array(model_miss.Residuals(:,1));
    miss2use = find(~isnan(temp));
    ampl_resdl = [ampl_resdl; temp(miss2use)];
    
    for ifreq = 1:nFreq
        phase2use = [ph_comb.phase_all(hits2use,ifreq,ichan); ph_comb.phase_all(miss2use,ifreq,ichan)];
        
        if~isempty(ampl_resdl)
            [correlation(1,ifreq),pvalue(1,ifreq)] = circ_corrcl(phase2use, ampl_resdl);
        end
    end
    
    chan_folder = fullfile(output_folder, 'regression', 'all_loc_difflev', num2str(ichan));
    if ~exist(chan_folder, 'dir'), mkdir(chan_folder); end
    cd(chan_folder)
    
    save correlation correlation
    save pvalue pvalue
    
end

%% Permutation

cd(data_load_folder)
load('ph_all_sess.mat')

% selecting trials for regression
hits2use = []; hits2use = find(ph_comb.ERP_trialinfo(:,20)==1);
miss2use = []; miss2use = find(ph_comb.ERP_trialinfo(:,20)==5);

% permutation for correlation after regression
permut_n = 1000;
perm_indices_hits = arrayfun(@(x) randperm(length(hits2use)), 1:permut_n, 'UniformOutput', false);
perm_indices_miss = arrayfun(@(x) randperm(length(miss2use)), 1:permut_n, 'UniformOutput', false);

cfg = cell(1,64);
for ichan = 1:64
    cfg{ichan}.ichan         = ichan;
    cfg{ichan}.permut_n      = permut_n;
    cfg{ichan}.infile        = fullfile(data_load_folder);
    cfg{ichan}.outfile       = fullfile(output_folder,'regression', 'all_loc_difflev');
    cfg{ichan}.perm_ind_hits = perm_indices_hits;
    cfg{ichan}.perm_ind_miss = perm_indices_miss;
    cfg{ichan}.hits2use      = hits2use;
    cfg{ichan}.miss2use      = miss2use;
end

slurmfun(@circlin_regression_perm, cfg, ...
    'partition','8GB', ...
    'stopOnError',false, ...
    'useUserPath',true);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting all locations and difficulty levels

cd(output_folder)
load('frequency.mat')
freq = frequency;
nCh = 64;

limit_max = nan(1, nCh);   % initialize with NaN
Corr_chan = false(nCh, numel(freq));

figure;
for ch = 1:nCh
    ch_folder = fullfile(output_folder, 'regression','all_loc_difflev', num2str(ch));
    if ~exist(ch_folder, 'dir')
        warning(['Skipping channel ' num2str(ch) ' (folder missing)']);
        continue
    end
    
    cd(ch_folder);
    
    if ~exist('correlation.mat','file') || ~exist('corr_perm.mat','file')
        warning(['Skipping channel ' num2str(ch) ' (missing correlation or corr_perm)']);
        continue
    end
    
    load correlation
    load corr_perm
    
    % Skip if correlation or permutation data has NaN
    if any(isnan(correlation)) || any(isnan(corr_perm(:)))
        warning(['Skipping channel ' num2str(ch) ' (NaN values present)']);
        continue
    end
    
    tmax = max(corr_perm, [], 2);
    limit_max(ch) = quantile(tmax, 0.95);
    
    % Skip if threshold is NaN
    if isnan(limit_max(ch))
        warning(['Skipping channel ' num2str(ch) ' (threshold is NaN)']);
        continue
    end
    
    subplot(8, 8, ch);
    plot_sigfreq(freq, correlation, limit_max(ch));
    title(['Ch ' num2str(ch)])
    
    % Store significant frequencies
    Corr_chan(ch,:) = correlation >= limit_max(ch);
end

% Summary image of significant frequencies
figure;
imagesc(freq, 1:nCh, Corr_chan);
set(gca, 'YDir', 'normal');
xlabel('Frequency (Hz)'); ylabel('Channels');
title('Significant Frequencies per Channel');
caxis([0 1]); colorbar;

% Combine Across Channels

valid_idx = ~isnan(limit_max);
corr_all = [];
corr_perm_all = [];

for ch = find(valid_idx)
    cd(fullfile(output_folder, 'regression','all_loc_difflev', num2str(ch)));
    load correlation
    load corr_perm
    
    corr_all = [corr_all; correlation];
    corr_perm_all = cat(3, corr_perm_all, corr_perm);
end

if isempty(corr_all)
    warning('No valid channels for combined plot');
else
    corr_avg = mean(corr_all, 1);
    corr_perm_avg = mean(corr_perm_all, 3);
    tmax_all = max(corr_perm_avg, [], 2);
    limit_avg = quantile(tmax_all, 0.95);
    
    figure;
    plot_sigfreq(freq, corr_avg, limit_avg);
    ylim([0 0.1]); title('All Channels Combined');
    xlabel('Frequency (Hz)'); ylabel('Correlation');
end

