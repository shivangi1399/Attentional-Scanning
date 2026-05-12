clear all
close all
clc

% Code description:
% -----------------
% Checking significance of correlation per location

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
datafolder   = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks';
cd(datafolder)

animalName = 'klecks';
temp = dir;
session_names = {};
ii = 0;
for i = 1:length(temp)
    if contains(temp(i).name, animalName)
        ii = ii+1;
        session_names{ii,1} = temp(i).name;
    end
end

session_paths_files = cellfun(@(x) fullfile(datafolder,x, 'clean_lfp.mat'), ...
    session_names, 'uniform', 0);

phase_paths = cellfun(@(x) fullfile(datafolder, x,'Phase_analysis/hit_miss'), ...
    session_names, 'uniform',0);

output_folder = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks/phase_correlation/ph55_till_100';
data_folder   = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks/multi_lin_reg/cp10_till_100';

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

% Load frequency and channels
cd(data_folder)
load('frequency.mat')
freq = frequency;
nCh = 64;

% Folder to save all figures
save_root = '/mnt/hpc/projects/MWSampling/4Shivangi/Plots/correlation/klecks/ph55_till_100/lfp/loc_difflev/all_dlev';
if ~exist(save_root, 'dir')
    mkdir(save_root);
end

%% Loop through target locations

for iloc = 1:length(targ_loc)
    loc     = targ_loc(iloc);
    min_th  = min_difflev(iloc);
    max_th  = max_difflev(iloc);
    
    fprintf('Analyzing location %d with difficulty %d_%d\n', loc, min_th, max_th);
    
    limit_max = nan(1, nCh);           % threshold per channel
    Corr_chan = false(nCh, numel(freq)); % significance per freq
    
    %% Loop through channels
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
        
        % Store significant frequencies
        Corr_chan(ch,:) = correlation >= limit_max(ch);
    end
    
    %% Determine significant channels
    
    sig_channels = any(Corr_chan, 2); % logical (nCh x 1)
    nan_channels = isnan(limit_max);  % channels with NaN threshold
    
    sig_map = double(sig_channels);
    sig_map(nan_channels) = -1;
    
    sig_grid = reshape(sig_map, [8, 8])';
    
    %% Plot summary figure
    fig_summary = figure('Name', sprintf('Summary Loc%d Dlev%d_%d', loc, min_th, max_th), ...
        'Position', [100 100 600 500]);
    
    imagesc(sig_grid);
    
    % Custom colormap: white for NaN, gray for not sig, teal for sig
    cmap = [1 1 1;    % for NaN (-1)
        0.8 0.8 0.8;  % not sig (0)
        0 0.5 0.5];   % sig (1)
    colormap(cmap);
    
    % Adjust color axis so -1 maps to first colormap entry
    caxis([-1 1])
    
    colorbar('Ticks', [-1,0,1], 'TickLabels', {'Broken','Not Sig','Sig'});
    axis equal tight
    xticks(1:8); yticks(1:8);
    xlabel('Channel X'); ylabel('Channel Y');
    title(sprintf('Significant Channels - Loc %d Dlev %d-%d', loc, min_th, max_th))
    
    % Save summary figure
    saveas(fig_summary, fullfile(save_root, ...
        sprintf('Summary_Loc%d_Dlev%d_%d.png', loc, min_th, max_th)));
    
    close(fig_summary)
    
end

disp('All summary figures generated.')
