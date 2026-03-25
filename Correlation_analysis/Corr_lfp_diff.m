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
output_folder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName], 'phase_correlation', 'cp10_till_100');
permut_n = 1000;

%% Correlation all locations and difficulty levels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% circular-linear correlation (real)
cd(output_folder)
load('ph_all_sess.mat')

for ichan = 1:64
    ichan
    
    phase = ph_comb.phase_all(:,:,ichan);
    erp_amp = ph_comb.LFP_ERP_ampl_all(:,ichan);
    
    for foi = 1:length(ph_comb.phase_all(1,:,1))
        [correlation(1,foi),pvalue(1,foi)] = circ_corrcl(phase(:,foi), erp_amp);
    end
    
    chan_folder = fullfile(output_folder,'lfp', 'all_loc_difflev', num2str(ichan));
    if ~exist(chan_folder, 'dir'), mkdir(chan_folder); end
    cd(chan_folder)
    
    save correlation correlation
    save pvalue pvalue
end

%% circular-linear correlation (permutation)

permut_n = 1000;

cd(output_folder)
load('ph_all_sess.mat')
nTrials = size(ph_comb.phase_all, 1);

perm_indices = arrayfun(@(x) randperm(nTrials), 1:permut_n, 'UniformOutput', false);

cfg = cell(1,64);
for i = 1:64
    cfg{i}.ichan = i;
    cfg{i}.permut_n = permut_n;
    cfg{i}.infile = fullfile(output_folder);
    cfg{i}.outfile = fullfile(output_folder, 'lfp','all_loc_difflev');
    cfg{i}.perm_indices = perm_indices;
end

slurmfun(@circlin_correlation_lfp, cfg, ...
    'partition',     '8GB', ...
    'stopOnError',   false, ...
    'useUserPath',   true);

%% Correlation particular locations and difficulty levels - remove 1 percent %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(output_folder)
load('ph_all_sess.mat')

% Find unique locations and difficulty levels
targ_loc = unique(ph_comb.trialinfo(:,16));
diff_levels = unique(ph_comb.trialinfo(:,18));

min_difflev = [];
max_difflev = [];

for iloc = 1:length(targ_loc)
    hits = find(ph_comb.trialinfo(:,20) == 1 & ph_comb.trialinfo(:,16) == targ_loc(iloc));
    miss = find(ph_comb.trialinfo(:,20) == 5 & ph_comb.trialinfo(:,16) == targ_loc(iloc));
    de = ph_comb.trialinfo([hits; miss], 18);
    
    min_thres_q = quantile(de, 0.01); % Compute quantile thresholds
    max_thres_q = quantile(de, 0.99);
    
    [~, min_idx] = min(abs(diff_levels - min_thres_q));
    [~, max_idx] = min(abs(diff_levels - max_thres_q));
    
    min_difflev = [min_difflev; diff_levels(min_idx)];
    max_difflev = [max_difflev; diff_levels(max_idx)];
end

%% Correlation (real data)

for iloc = 1:length(targ_loc)
    iloc
    loc = targ_loc(iloc);
    min_th = min_difflev(iloc);
    max_th = max_difflev(iloc);
    
    % Select trials for this location and difficulty
    
    trial_idx = find(ph_comb.trialinfo(:,16) == loc & ...
                     ph_comb.trialinfo(:,18) >= min_th & ...
                     ph_comb.trialinfo(:,18) <= max_th);
    
    phase_data = ph_comb.phase_all(trial_idx,:,:);  % trials x freq x channels
    erp_amp    = ph_comb.LFP_ERP_ampl_all(trial_idx,:); % trials x channels
    
    for ichan = 1:64
        phase = squeeze(phase_data(:,:,ichan));  % trials x freq
        erp   = erp_amp(:,ichan);
        
        correlation = nan(1, size(phase,2));
        pvalue      = nan(1, size(phase,2));
        for foi = 1:size(phase,2)
            [correlation(foi), pvalue(foi)] = circ_corrcl(phase(:,foi), erp);
        end
        
        chan_folder = fullfile(output_folder, 'lfp',...
            'loc_difflev_1perc', ...
            sprintf('loc%d', loc), ...
            sprintf('%d_%d', min_th, max_th), ...
            num2str(ichan));
        if ~exist(chan_folder, 'dir'), mkdir(chan_folder); end
        save(fullfile(chan_folder, 'correlation.mat'), 'correlation')
        save(fullfile(chan_folder, 'pvalue.mat'), 'pvalue')
    end
    
    %% Permutation
    
    permut_n = 1000;
    nTrials  = length(trial_idx);
    perm_indices = arrayfun(@(x) randperm(nTrials), 1:permut_n, 'UniformOutput', false);
    
    cfg = cell(1,64);
    for ichan = 1:64
        cfg{ichan}.ichan        = ichan;
        cfg{ichan}.permut_n     = permut_n;
        cfg{ichan}.infile       = fullfile(output_folder);
        cfg{ichan}.outfile      = fullfile(output_folder, 'lfp', ...
            'loc_difflev_1perc', ...
            sprintf('loc%d', loc), ...
            sprintf('%d_%d', min_th, max_th));
        cfg{ichan}.perm_indices = perm_indices;
        cfg{ichan}.trial_idx    = trial_idx;
    end
    
    % Launch jobs
    slurmfun(@circlin_corr_select_trl_lfp, cfg, ...
        'partition',   '8GB', ...
        'stopOnError', false, ...
        'useUserPath', true);
end

%% Correlation particular locations and difficulty levels - keep all Dlev %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(output_folder)
load('ph_all_sess.mat')

% Find unique locations and difficulty levels
targ_loc = unique(ph_comb.trialinfo(:,16));
diff_levels = unique(ph_comb.trialinfo(:,18));

min_difflev = [];
max_difflev = [];

for iloc = 1:length(targ_loc)
    hits = find(ph_comb.trialinfo(:,20) == 1 & ph_comb.trialinfo(:,16) == targ_loc(iloc));
    miss = find(ph_comb.trialinfo(:,20) == 5 & ph_comb.trialinfo(:,16) == targ_loc(iloc));
    de = ph_comb.trialinfo([hits; miss], 18);
    
    min_difflev = [min_difflev; min(de)];
    max_difflev = [max_difflev; max(de)];
end

%% Correlation (real data)

for iloc = 1:length(targ_loc)
    iloc
    loc = targ_loc(iloc);
    min_th = min_difflev(iloc);
    max_th = max_difflev(iloc);
    
    % Select trials for this location and difficulty
    trial_idx = find(ph_comb.trialinfo(:,16) == loc & ...
                     ph_comb.trialinfo(:,18) >= min_th & ...
                     ph_comb.trialinfo(:,18) <= max_th);
    
    phase_data = ph_comb.phase_all(trial_idx,:,:);  % trials x freq x channels
    erp_amp    = ph_comb.LFP_ERP_ampl_all(trial_idx,:); % trials x channels
    
    for ichan = 1:64
        phase = squeeze(phase_data(:,:,ichan));  % trials x freq
        erp   = erp_amp(:,ichan);
        
        correlation = nan(1, size(phase,2));
        pvalue      = nan(1, size(phase,2));
        for foi = 1:size(phase,2)
            [correlation(foi), pvalue(foi)] = circ_corrcl(phase(:,foi), erp);
        end
        
        chan_folder = fullfile(output_folder, 'lfp',...
            'loc_difflev_all', ...
            sprintf('loc%d', loc), ...
            sprintf('%d_%d', min_th, max_th), ...
            num2str(ichan));
        if ~exist(chan_folder, 'dir'), mkdir(chan_folder); end
        save(fullfile(chan_folder, 'correlation.mat'), 'correlation')
        save(fullfile(chan_folder, 'pvalue.mat'), 'pvalue')
    end
    
    %% Permutation
    
    permut_n = 1000;
    nTrials  = length(trial_idx);
    perm_indices = arrayfun(@(x) randperm(nTrials), 1:permut_n, 'UniformOutput', false);
    
    cfg = cell(1,64);
    for ichan = 1:64
        cfg{ichan}.ichan        = ichan;
        cfg{ichan}.permut_n     = permut_n;
        cfg{ichan}.infile       = fullfile(output_folder);
        cfg{ichan}.outfile      = fullfile(output_folder, 'lfp',...
            'loc_difflev_all', ...
            sprintf('loc%d', loc), ...
            sprintf('%d_%d', min_th, max_th));
        cfg{ichan}.perm_indices = perm_indices;
        cfg{ichan}.trial_idx    = trial_idx;
    end
    
    % Launch jobs
    slurmfun(@circlin_corr_select_trl_lfp, cfg, ...
        'partition',   '8GB', ...
        'stopOnError', false, ...
        'useUserPath', true);
end

%% Correlation locations and difficulty levels - keep all Dlev - separate hit/miss  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(output_folder)
load('ph_all_sess.mat')
output_folder_s = fullfile(output_folder,'lfp','only_miss');

% Find unique locations and difficulty levels
targ_loc = unique(ph_comb.trialinfo(:,16));
diff_levels = unique(ph_comb.trialinfo(:,18));

min_difflev = [];
max_difflev = [];

for iloc = 1:length(targ_loc)
    %hits = find(ph_comb.trialinfo(:,20) == 1 & ph_comb.trialinfo(:,16) == targ_loc(iloc));
    miss = find(ph_comb.trialinfo(:,20) == 5 & ph_comb.trialinfo(:,16) == targ_loc(iloc));
    %de = ph_comb.trialinfo([hits], 18);
    de = ph_comb.trialinfo([miss], 18);
    
    min_difflev = [min_difflev; min(de)];
    max_difflev = [max_difflev; max(de)];
end

%% Correlation (real data)

for iloc = 1:length(targ_loc)
    iloc
    loc = targ_loc(iloc);
    min_th = min_difflev(iloc);
    max_th = max_difflev(iloc);
    
    % Select trials for this location and difficulty
    trial_idx = find(ph_comb.trialinfo(:,16) == loc & ...
                     ph_comb.trialinfo(:,20) >= 1 & ... %hit or miss
                     ph_comb.trialinfo(:,18) >= min_th & ...
                     ph_comb.trialinfo(:,18) <= max_th);
    
    phase_data = ph_comb.phase_all(trial_idx,:,:);  % trials x freq x channels
    erp_amp    = ph_comb.LFP_ERP_ampl_all(trial_idx,:); % trials x channels
    
    for ichan = 1:64
        phase = squeeze(phase_data(:,:,ichan));  % trials x freq
        erp   = erp_amp(:,ichan);
        
        correlation = nan(1, size(phase,2));
        pvalue      = nan(1, size(phase,2));
        for foi = 1:size(phase,2)
            [correlation(foi), pvalue(foi)] = circ_corrcl(phase(:,foi), erp);
        end
        
        chan_folder = fullfile(output_folder_s, ...
            'loc_difflev_all', ...
            sprintf('loc%d', loc), ...
            sprintf('%d_%d', min_th, max_th), ...
            num2str(ichan));
        if ~exist(chan_folder, 'dir'), mkdir(chan_folder); end
        save(fullfile(chan_folder, 'correlation.mat'), 'correlation')
        save(fullfile(chan_folder, 'pvalue.mat'), 'pvalue')
    end
    
    %% Permutation
    
    permut_n = 1000;
    nTrials  = length(trial_idx);
    perm_indices = arrayfun(@(x) randperm(nTrials), 1:permut_n, 'UniformOutput', false);
    
    cfg = cell(1,64);
    for ichan = 1:64
        cfg{ichan}.ichan        = ichan;
        cfg{ichan}.permut_n     = permut_n;
        cfg{ichan}.infile       = fullfile(output_folder);
        cfg{ichan}.outfile      = fullfile(output_folder_s, ...
            'loc_difflev_all', ...
            sprintf('loc%d', loc), ...
            sprintf('%d_%d', min_th, max_th));
        cfg{ichan}.perm_indices = perm_indices;
        cfg{ichan}.trial_idx    = trial_idx;
    end
    
    % Launch jobs
    slurmfun(@circlin_corr_select_trl_lfp, cfg, ...
        'partition',   '8GB', ...
        'stopOnError', false, ...
        'useUserPath', true);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting all locations and difficulty levels

save_root = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/correlation', animalName, 'cp10_till_100', 'lfp', 'all_loc_difflev');
if ~exist(save_root, 'dir')
    mkdir(save_root);
end

cd(output_folder)
load('frequency.mat')
freq = frequency;
nCh = 64;

limit_max = nan(1, nCh);   % initialize with NaN
Corr_chan = false(nCh, numel(freq));

f1 = figure(1);
for ch = 1:nCh
    ch_folder = fullfile(output_folder,'lfp','all_loc_difflev', num2str(ch));
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

set(f1, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
set(f1, 'PaperPositionMode', 'auto');
set(f1, 'Renderer', 'opengl');
print(f1, fullfile(save_root, 'all_channels.png'),'-dpng', '-r0'); 

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
    cd(fullfile(output_folder,'lfp','all_loc_difflev', num2str(ch)));
    load correlation
    load corr_perm
    
    corr_all = [corr_all; correlation];
    corr_perm_all = cat(3, corr_perm_all, corr_perm);
end

if isempty(corr_all)
    warning('No valid channels for combined plot');
else
    corr_avg = nanmean(corr_all, 1);
    % method 1
    tmax_all = zeros(permut_n, 1);
    for p = 1:permut_n
        tmax_all(p) = nanmax(corr_perm_all(p, :, :), [], 'all'); % max over freqs + channels
%         testMax(p,:) = nanmax(squeeze(corr_perm_all(p,:,:)));
%         testMax2(p) = nanmax(testMax(p,:));
    end
    limit_avg = quantile(tmax_all, 0.95);
    
    % method 2
%     corr_perm_avg = nanmean(corr_perm_all, 3);
%     tmax_all = nanmax(corr_perm_avg, [], 2);
%     limit_avg = quantile(tmax_all, 0.95);
    
    figure;
    plot_sigfreq(freq, corr_avg, limit_avg);
    ylim([0 0.1]); title('All Channels Combined');
    xlabel('Frequency (Hz)'); ylabel('Correlation');
end

%% Plotting combinations of difficulty level and locations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(output_folder)
load('ph_all_sess.mat')
% load the max min stuff here depending on what you want to plot

% Load frequency
cd(output_folder)
load('frequency.mat')
freq = frequency;
nCh = 64;

% Folder to save all figures
save_root = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/correlation', animalName, 'ph55_till_100', 'lfp', 'loc_difflev', 'all_dlev');
if ~exist(save_root, 'dir')
    mkdir(save_root);
end

for iloc = 1:length(targ_loc)
    loc     = targ_loc(iloc);
    min_th  = min_difflev(iloc);
    max_th  = max_difflev(iloc);

    fprintf('Plotting for location %d with difficulty %d_%d\n', loc, min_th, max_th);

    % Initialize storage
    limit_max = nan(1, nCh);    
    Corr_chan = false(nCh, numel(freq));

    %% Per-channel plot
    fig1 = figure('Name', sprintf('Loc%d Dlev%d_%d', loc, min_th, max_th), ...
                  'Position', [100 100 1200 800]);

    for ch = 1:nCh
        ch_folder = fullfile(output_folder, 'lfp',...
            'loc_difflev_all', ...
            sprintf('loc%d', loc), ...
            sprintf('%d_%d', min_th, max_th), ...
            num2str(ch));

        if ~exist(ch_folder, 'dir')
            warning('Skipping channel %d (folder missing)', ch);
            continue
        end

        cd(ch_folder);

        if ~exist('correlation.mat','file') || ~exist('corr_perm.mat','file')
            warning('Skipping channel %d (missing correlation or corr_perm)', ch);
            continue
        end

        load('correlation.mat','correlation')
        load('corr_perm.mat','corr_perm')

        % Skip if correlation or permutation data has NaN
        if any(isnan(correlation)) || any(isnan(corr_perm(:)))
            warning('Skipping channel %d (NaN values present)', ch);
            continue
        end

        % Compute threshold from permutations (95% quantile)
        tmax = max(corr_perm, [], 2);
        limit_max(ch) = quantile(tmax, 0.95);

        if isnan(limit_max(ch))
            warning('Skipping channel %d (threshold NaN)', ch);
            continue
        end

        % Plot channel-specific correlation and significant frequencies
        subplot(8, 8, ch);
        plot_sigfreq(freq, correlation, limit_max(ch));
        title(['Ch ' num2str(ch)])

        % Store significant frequencies
        Corr_chan(ch,:) = correlation >= limit_max(ch);
    end

    % Save per-channel figure
    %saveas(fig1, fullfile(save_root, sprintf('Loc%d_Dlev%d_%d_channels.png', loc, min_th, max_th)));
    close(fig1);

    %% Summary image of significant frequencies
    fig2 = figure('Name', sprintf('SignificantFreqs Loc%d (%d_%d)', loc, min_th, max_th));
    imagesc(freq, 1:nCh, Corr_chan);
    set(gca, 'YDir', 'normal');
    xlabel('Frequency (Hz)'); ylabel('Channels');
    title(sprintf('Significant Frequencies Loc%d (%d_%d)', loc, min_th, max_th));
    caxis([0 1]); colorbar;

    % Save summary figure
    %saveas(fig2, fullfile(save_root, sprintf('Loc%d_Dlev%d_%d_sigfreqs.png', loc, min_th, max_th)));
    close(fig2);

    %% Combine Across Channels
    valid_idx = ~isnan(limit_max);
    corr_all = [];
    corr_perm_all = [];

    for ch = find(valid_idx)
        ch_folder = fullfile(output_folder, 'lfp',...
            'loc_difflev_all', ...
            sprintf('loc%d', loc), ...
            sprintf('%d_%d', min_th, max_th), ...
            num2str(ch));
        cd(ch_folder);
        load correlation
        load corr_perm

        corr_all = [corr_all; correlation];
        corr_perm_all = cat(3, corr_perm_all, corr_perm);
    end

    if isempty(corr_all)
        warning('No valid channels for combined plot (Loc%d (%d_%d))', loc, min_th, max_th);
        continue
    end

    corr_avg = mean(corr_all, 1);
    % method 1
    tmax_all = zeros(permut_n, 1);
    for p = 1:permut_n
        tmax_all(p) = max(corr_perm_all(p, :, :), [], 'all'); % max over freqs + channels
    end
    limit_avg = quantile(tmax_all, 0.95);
    
    % method 2
%     corr_perm_avg = mean(corr_perm_all, 3);
%     tmax_all = max(corr_perm_avg, [], 2);
%     limit_avg = quantile(tmax_all, 0.95);

    fig3 = figure('Name', sprintf('Combined_Loc%d (%d_%d)', loc, min_th, max_th));
    plot_sigfreq(freq, corr_avg, limit_avg);
    ylim([0 0.1]);
    title(sprintf('All Channels Combined Loc%d Dlev%d-%d', loc, min_th, max_th));
    xlabel('Frequency (Hz)'); ylabel('Correlation');

    % Save combined figure
    %saveas(fig3, fullfile(save_root, sprintf('Loc%d_Dlev%d_%d_combined.png', loc, min_th, max_th)));
    close(fig3);
end

%% Plotting combinations of difficulty level and locations - hits/miss %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_folder_s = fullfile(output_folder,'lfp','only_miss');
cd(output_folder)
load('ph_all_sess.mat')
% load the max min stuff here depending on what you want to plot

% Load frequency
cd(output_folder)
load('frequency.mat')
freq = frequency;
nCh = 64;

% Folder to save all figures
save_root = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/correlation', animalName, 'ph55_till_100', 'lfp', 'only_miss', 'loc_difflev', 'all_dlev');
if ~exist(save_root, 'dir')
    mkdir(save_root);
end

for iloc = 1:length(targ_loc)
    loc     = targ_loc(iloc);
    min_th  = min_difflev(iloc);
    max_th  = max_difflev(iloc);

    fprintf('Plotting for location %d with difficulty %d_%d\n', loc, min_th, max_th);

    % Initialize storage
    limit_max = nan(1, nCh);    
    Corr_chan = false(nCh, numel(freq));

    %% Per-channel plot
    fig1 = figure('Name', sprintf('Loc%d Dlev%d_%d', loc, min_th, max_th), ...
                  'Position', [100 100 1200 800]);

    for ch = 1:nCh
        ch_folder = fullfile(output_folder_s, ...
            'loc_difflev_all', ...
            sprintf('loc%d', loc), ...
            sprintf('%d_%d', min_th, max_th), ...
            num2str(ch));

        if ~exist(ch_folder, 'dir')
            warning('Skipping channel %d (folder missing)', ch);
            continue
        end

        cd(ch_folder);

        if ~exist('correlation.mat','file') || ~exist('corr_perm.mat','file')
            warning('Skipping channel %d (missing correlation or corr_perm)', ch);
            continue
        end

        load('correlation.mat','correlation')
        load('corr_perm.mat','corr_perm')

        % Skip if correlation or permutation data has NaN
        if any(isnan(correlation)) || any(isnan(corr_perm(:)))
            warning('Skipping channel %d (NaN values present)', ch);
            continue
        end

        % Compute threshold from permutations (95% quantile)
        tmax = max(corr_perm, [], 2);
        limit_max(ch) = quantile(tmax, 0.95);

        if isnan(limit_max(ch))
            warning('Skipping channel %d (threshold NaN)', ch);
            continue
        end

        % Plot channel-specific correlation and significant frequencies
        subplot(8, 8, ch);
        plot_sigfreq(freq, correlation, limit_max(ch));
        title(['Ch ' num2str(ch)])

        % Store significant frequencies
        Corr_chan(ch,:) = correlation >= limit_max(ch);
    end

    % Save per-channel figure
    saveas(fig1, fullfile(save_root, sprintf('Loc%d_Dlev%d_%d_channels.png', loc, min_th, max_th)));
    close(fig1);

    %% Summary image of significant frequencies
    fig2 = figure('Name', sprintf('SignificantFreqs Loc%d (%d_%d)', loc, min_th, max_th));
    imagesc(freq, 1:nCh, Corr_chan);
    set(gca, 'YDir', 'normal');
    xlabel('Frequency (Hz)'); ylabel('Channels');
    title(sprintf('Significant Frequencies Loc%d (%d_%d)', loc, min_th, max_th));
    caxis([0 1]); colorbar;

    % Save summary figure
    saveas(fig2, fullfile(save_root, sprintf('Loc%d_Dlev%d_%d_sigfreqs.png', loc, min_th, max_th)));
    close(fig2);

    %% Combine Across Channels
    valid_idx = ~isnan(limit_max);
    corr_all = [];
    corr_perm_all = [];

    for ch = find(valid_idx)
        ch_folder = fullfile(output_folder_s,...
            'loc_difflev_all', ...
            sprintf('loc%d', loc), ...
            sprintf('%d_%d', min_th, max_th), ...
            num2str(ch));
        cd(ch_folder);
        load correlation
        load corr_perm

        corr_all = [corr_all; correlation];
        corr_perm_all = cat(3, corr_perm_all, corr_perm);
    end

    if isempty(corr_all)
        warning('No valid channels for combined plot (Loc%d (%d_%d))', loc, min_th, max_th);
        continue
    end

    corr_avg = mean(corr_all, 1);
    % method 1
%     tmax_all = zeros(permut_n, 1);
%     for p = 1:permut_n
%         tmax_all(p) = max(corr_perm_all(p, :, :), [], 'all'); % max over freqs + channels
%     end
%     limit_avg = quantile(tmax_all, 0.95);
    
    % method 2
    corr_perm_avg = mean(corr_perm_all, 3);
    tmax_all = max(corr_perm_avg, [], 2);
    limit_avg = quantile(tmax_all, 0.95);

    fig3 = figure('Name', sprintf('Combined_Loc%d (%d_%d)', loc, min_th, max_th));
    plot_sigfreq(freq, corr_avg, limit_avg);
    ylim([0 0.1]);
    title(sprintf('All Channels Combined Loc%d Dlev%d-%d', loc, min_th, max_th));
    xlabel('Frequency (Hz)'); ylabel('Correlation');

    % Save combined figure
    saveas(fig3, fullfile(save_root, sprintf('Loc%d_Dlev%d_%d_combined.png', loc, min_th, max_th)));
    close(fig3);
end
