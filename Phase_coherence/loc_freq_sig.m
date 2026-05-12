clear all
close all
clc

% Code description:
% -----------------
% Coherence Analysis:Use precalculated LFPERP phase coherence per channel and visualizes significant channels on 8×8 grids.
% Average Phase:Calculates circular mean phase per channel and band, visualized with color and arrows.
% Phase Differences:Computes inter-band phase differences per channel and visualizes across the grid.
% Traveling Waves: Estimates phase gradients to show wave direction and magnitude across the electrode array.

%% Specify paths
addpath /opt/fieldtrip_github/
ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/coherence_analysis
addpath /mnt/hpc/projects/MWSampling/4Shivangi
addpath /mnt/hpc/projects/MWSampling/4Shivangi/software_folder/CircStat2012a
clc

%% Create data paths
datafolder   = '/mnt/hpc/projects/MWSampling/4Shivangi/results_hermes';
cd(datafolder)

animalName = 'hermes';
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

output_folder = '/mnt/hpc/projects/MWSampling/4Shivangi/results_hermes/phase_coherence/cp10_till_100';
data_folder   = '/mnt/hpc/projects/MWSampling/4Shivangi/results_hermes/multi_lin_reg/cp10_till_100';

%% coherence particular locations and difficulty levels - keep all Dlev

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
save_root = '/mnt/hpc/projects/MWSampling/4Shivangi/Plots/coherence/hermes/cp10_till_100/lfp/loc_difflev_all';
if ~exist(save_root, 'dir')
    mkdir(save_root);
end

%% Coherence plots per location - all frequencies %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iloc = 1:length(targ_loc)
    loc     = targ_loc(iloc);
    min_th  = min_difflev(iloc);
    max_th  = max_difflev(iloc);
    
    fprintf('Analyzing location %d with difficulty %d_%d\n', loc, min_th, max_th);
    
    limit_maxc = nan(1, nCh);   % coherence threshold
    Coh_chan = false(nCh, numel(freq)); % significance per freq
    
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
        
        if ~exist('coherence.mat','file') || ~exist('coh_perm.mat','file') || ...
                ~exist('coherence.mat','file') || ~exist('phase_spec_perm.mat','file')
            warning(['Skipping channel ' num2str(ch) ' (missing required files)']);
            continue
        end
        
        load('coherence.mat')
        load('coh_perm.mat')
        load('phase_spec_perm')
        
        % Skip if data has NaN
        if any(isnan(coh)) || any(isnan(coh_perm(:))) || ...
                any(isnan(phase_spec)) || any(isnan(phase_spec_perm(:)))
            warning(['Skipping channel ' num2str(ch) ' (NaN values present)']);
            continue
        end
        
        % Compute threshold from permutations (95% quantile)
        tmaxc = max(coh_perm, [], 2);
        limit_maxc(ch) = quantile(tmaxc, 0.95);
        
        if isnan(limit_maxc(ch))
            warning('Skipping channel %d (threshold NaN)', ch);
            continue
        end
        
        % Store significant frequencies
        Coh_chan(ch,:) = coh >= limit_maxc(ch);
    end
    
    %% Determine significant channels
    
    sig_channels = any(Coh_chan, 2); % logical (nCh x 1)
    nan_channels = isnan(limit_maxc);  % channels with NaN threshold
    
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

%% Coherence plots per location - per frequency band %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Determining frequency bands
all_sig_freqs = false(1, numel(freq));  % accumulator across all locs + chans

for iloc = 1:length(targ_loc)
    loc     = targ_loc(iloc);
    min_th  = min_difflev(iloc);
    max_th  = max_difflev(iloc);
    
    for ch = 1:nCh
        ch_folder = fullfile(output_folder, 'lfp',...
            'loc_difflev_all', ...
            sprintf('loc%d', loc), ...
            sprintf('%d_%d', min_th, max_th), ...
            num2str(ch));
        
        if ~exist(ch_folder, 'dir')
            warning('Skipping channel %d (folder missing)', ch);
            continues
        end
        
        cd(ch_folder);
        
        if ~exist('coherence.mat','file') || ~exist('coh_perm.mat','file') || ...
                ~exist('coherence.mat','file') || ~exist('phase_spec_perm.mat','file')
            warning(['Skipping channel ' num2str(ch) ' (missing required files)']);
            continue
        end
        
        load('coherence.mat')
        load('coh_perm.mat')
        load('phase_spec_perm')
        
        % Skip if data has NaN
        if any(isnan(coh)) || any(isnan(coh_perm(:))) || ...
                any(isnan(phase_spec)) || any(isnan(phase_spec_perm(:)))
            warning(['Skipping channel ' num2str(ch) ' (NaN values present)']);
            continue
        end
        
        % Compute threshold from permutations (95% quantile)
        tmaxc = max(coh_perm, [], 2);
        limit_maxc(ch) = quantile(tmaxc, 0.95);
        
        if isnan(limit_maxc(ch))
            warning('Skipping channel %d (threshold NaN)', ch);
            continue
        end
        
        % Store significant frequencies
        Coh_chan(ch,:) = coh >= limit_maxc(ch);
    end
    
    % Collapse across channels for this location
    sig_freqs_loc = any(Coh_chan, 1);
    
    % Accumulate across all locations
    all_sig_freqs = all_sig_freqs | sig_freqs_loc;
end

% Define canonical bands
canonical_bands = {
    'theta',     [4 8];
    'alpha',     [9 12];
    'beta',      [13 30];
    'gamma',     [30 80];
    'highgamma', [80 150];
    };

% Find significant bands inside canonical ranges
significant_bands = {};
for iC = 1:size(canonical_bands,1)
    band_range = canonical_bands{iC,2};
    
    % Allow a small tolerance on edges (0.5 Hz)
    in_band = freq >= (band_range(1)-0.5) & freq <= (band_range(2)+0.5);
    sig_here = all_sig_freqs & in_band;
    
    % Set minimum width rule (alpha can be narrower)
    if strcmp(canonical_bands{iC,1},'alpha')
        min_width = 2;
    else
        min_width = 5;
    end
    
    % Detect contiguous significant "islands"
    is_in = false;
    for iF = find(in_band)
        if sig_here(iF) && ~is_in
            startF = freq(iF);
            is_in = true;
        elseif ~sig_here(iF) && is_in
            endF = freq(iF-1);
            if (endF - startF) >= min_width
                significant_bands{end+1} = [startF endF];  % store range only
            end
            is_in = false;
        end
    end
    
    % Handle case if significance continues to the edge
    if is_in
        last_idx = find(in_band,1,'last');
        endF = freq(last_idx);
        if (endF - startF) >= min_width
            significant_bands{end+1} = [startF endF];
        end
    end
end

%% Define frequency bands externally

significant_bands = {
    [4, 8]
    [9, 12]; 
    [13, 30]; 
    [31, 80];
    };

%% Plot summary figure

for iloc = 1:length(targ_loc)
    loc     = targ_loc(iloc);
    min_th  = min_difflev(iloc);
    max_th  = max_difflev(iloc);
    
    fprintf('Generating grid summary for location %d with difficulty %d_%d\n', loc, min_th, max_th);
    
    % Preallocate: channels x frequency bands
    nBands = numel(significant_bands);
    Coh_chan_band = false(nCh, nBands);
    broken_chan = false(nCh,1);  % keep track of broken channels
    
    for ch = 1:nCh
        ch_folder = fullfile(output_folder, 'lfp', ...
            'loc_difflev_all', ...
            sprintf('loc%d', loc), ...
            sprintf('%d_%d', min_th, max_th), ...
            num2str(ch));
        
        % If folder missing, mark as broken
        if ~exist(ch_folder, 'dir')
            warning('Channel %d folder missing, marking as broken', ch);
            broken_chan(ch) = true;
            continue
        end
        
        cd(ch_folder);
        
        if ~exist('coherence.mat','file') || ~exist('coh_perm.mat','file')
            warning(['Channel ' num2str(ch) ' missing required files, marking as broken']);
            broken_chan(ch) = true;
            continue
        end
        
        load('coherence.mat')
        load('coh_perm.mat')
        
        % Skip if data has NaN
        if any(isnan(coh)) || any(isnan(coh_perm(:)))
            warning(['Channel ' num2str(ch) ' has NaN values, marking as broken']);
            broken_chan(ch) = true;
            continue
        end
        
        % Compute threshold from permutations (95% quantile)
        tmaxc = max(coh_perm, [], 2);
        limit_maxc = quantile(tmaxc, 0.95);
        
        if isnan(limit_maxc)
            warning(['Channel ' num2str(ch) ' threshold is NaN, marking as broken']);
            broken_chan(ch) = true;
            continue
        end
        
        % Check significance per frequency band
        for b = 1:nBands
            f_range = significant_bands{b};
            freq_idx = freq >= f_range(1) & freq <= f_range(2);
            
            if any(coh(freq_idx) >= limit_maxc)
                Coh_chan_band(ch,b) = true;
            end
        end
    end
    
    %% Plot summary grids for each frequency band
    for b = 1:nBands
        sig_channels = Coh_chan_band(:,b);
        
        % Map channels: -1 = broken, 0 = not sig, 1 = sig
        sig_map = double(sig_channels);
        sig_map(broken_chan) = -1;
        
        sig_grid = reshape(sig_map, [8, 8])';
        
        fig_band = figure('Name', sprintf('Loc%d Band%d', loc, b), ...
            'Position', [100 100 600 500]);
        
        imagesc(sig_grid);
        
        % Custom colormap: white for NaN/broken, gray for not sig, teal for sig
        cmap = [1 1 1;    % broken (-1)
            0.8 0.8 0.8;  % not sig (0)
            0 0.5 0.5];   % sig (1)
        colormap(cmap);
        
        caxis([-1 1])
        colorbar('Ticks', [-1,0,1], 'TickLabels', {'Broken','Not Sig','Sig'});
        axis equal tight
        xticks(1:8); yticks(1:8);
        xlabel('Channel X'); ylabel('Channel Y');
        title(sprintf('Loc %d Dlev %d-%d: %.0f-%.0f Hz', loc, min_th, max_th, significant_bands{b}(1), significant_bands{b}(2)));
        
        % Save summary figure
        save_folder = fullfile(save_root, ...
            sprintf('Band%d_%.0f-%.0fHz', b, significant_bands{b}(1), significant_bands{b}(2)));
        if ~exist(save_folder, 'dir')
            mkdir(save_folder);
        end
        saveas(fig_band, fullfile(save_folder, ...
            sprintf('Summary_Loc%d_Band%d_%.0f-%.0fHz.png', loc, b, significant_bands{b}(1), significant_bands{b}(2))));
        
        close(fig_band)
    end
    
end

disp('All location x frequency-band summary grids generated.');

%% Phase per location and frequency band %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

phase_avg_all = cell(length(targ_loc),1); % store per location

% Set global phase scale for all plots
phase_min = -pi;
phase_max = pi;

for iloc = 1:length(targ_loc)
    loc     = targ_loc(iloc);
    min_th  = min_difflev(iloc);
    max_th  = max_difflev(iloc);
    
    fprintf('Computing average phase for location %d\n', loc);
    
    nBands = numel(significant_bands);
    phase_avg = nan(nCh, nBands);  % channels x bands
    
    for ch = 1:nCh
        ch_folder = fullfile(output_folder, 'lfp', ...
            'loc_difflev_all', sprintf('loc%d', loc), ...
            sprintf('%d_%d', min_th, max_th), num2str(ch));
        
        if ~exist(ch_folder, 'dir')
            continue
        end
        
        cd(ch_folder);
        
        if ~exist('coherence.mat','file')
            warning(['Channel ' num2str(ch) ' missing phase_spec.mat']);
            continue
        end
        
        load('coherence.mat') % variable should be phase_spec: freq x trials/time
        load('coh_perm.mat')
        
        tmaxc = max(coh_perm, [], 2);
        limit_maxc = quantile(tmaxc, 0.95);
        %sig_mask = coh >= limit_maxc; %comment this out to plot without
        %coherence significance mask
        
        for b = 1:nBands
            f_range = significant_bands{b};
            freq_idx = freq >= f_range(1) & freq <= f_range(2);
            
            if any(freq_idx)
                if exist('sig_mask','var')
                    valid_idx = freq_idx & sig_mask;
                else
                    valid_idx = freq_idx; % fallback: use all
                end
                
                if any(valid_idx)
                    phase_data = phase_spec(valid_idx); % only significant freqs
                    phase_c = exp(1i*phase_data);
                    phase_avg(ch,b) = angle(mean(phase_c(:))); % circular mean
                else
                    phase_avg(ch,b) = NaN; % no sig freqs ? exclude channel
                end
            end
        end
    end
    
    phase_avg_all{iloc} = phase_avg;
    
    %% Visualize as 8x8 grid per band - just visualizing phase values
    for b = 1:nBands
        sig_grid_phase = reshape(phase_avg(:,b), [8,8])'; % transpose to match previous plots
        
        fig = figure;  % store figure handle
        imagesc(sig_grid_phase);
        colorbar;
        colormap(coolCircularColormap(256));
        caxis([phase_min phase_max]); % consistent scale
        title(sprintf('Avg Phase - Loc %d Band %d: %.0f-%.0f Hz', loc, b, significant_bands{b}(1), significant_bands{b}(2)));
        axis equal tight;
        hold on;
        [X,Y] = meshgrid(1:8, 1:8);
        quiver(X, Y, cos(sig_grid_phase), sin(sig_grid_phase), 0.5, 'k');
        xticks(1:8); yticks(1:8);
        xlabel('Channel X'); ylabel('Channel Y');
        
        save_folder = fullfile(save_root, ...
            sprintf('Band%d_%.0f-%.0fHz', b, significant_bands{b}(1), significant_bands{b}(2)));
        if ~exist(save_folder, 'dir')
            mkdir(save_folder);
        end
        
        % Use figure handle when saving
        saveas(fig, fullfile(save_folder, ...
            sprintf('AvgPhase_Loc_nosig%d_Band%d_%.0f-%.0fHz.png', loc, b, significant_bands{b}(1), significant_bands{b}(2))));
        
        % Close the same figure handle
        close(fig);
    end
    
end

disp('Average phase computed for all locations and bands with consistent scale.');

%% Difference betwween frequency bands per channel

num_locations = length(phase_avg_all);
num_channels = size(phase_avg_all{1},1);
num_bands = size(phase_avg_all{1},2);

% Preallocate
phase_diff_all = cell(num_locations,1);

% Compute differences
for loc = 1:num_locations
    data = phase_avg_all{loc}; % 64 x 3
    diff_mat = [data(:,2)-data(:,1), ... % Beta - Alpha
        data(:,3)-data(:,1), ... % Gamma - Alpha
        data(:,3)-data(:,2)];    % Gamma - Beta
    phase_diff_all{loc} = diff_mat; % 64 x 3
end

% Find global min/max (ignoring NaNs)
all_vals = cell2mat(cellfun(@(x) x(:), phase_diff_all, 'UniformOutput', false));
global_min = min(all_vals(~isnan(all_vals)));
global_max = max(all_vals(~isnan(all_vals)));

% Define colormap with black for NaNs
cmap = parula(256);
cmap = [0 0 0; cmap]; % prepend black

% Create single figure with 3x3 subplots
fig = figure('Position',[100 100 1200 900]);

for loc = 1:num_locations
    subplot(4,4,loc);
    data = phase_diff_all{loc};
    
    % Replace NaNs with a sentinel value
    nan_mask = isnan(data);
    data_for_plot = data;
    data_for_plot(nan_mask) = global_min - 1; % put NaNs below range
    
    % Plot
    imagesc(data_for_plot);
    colormap(cmap);
    colormap(coolCircularColormap(256));
    caxis([global_min global_max]);
    
    xlabel('Band comparison');
    ylabel('Channel');
    title(sprintf('Location %d', loc));
    xticks(1:3);
    xticklabels({'Beta-Alpha','Gamma-Alpha','Gamma-Beta'});
end

% Add one shared colorbar
h = colorbar('Position',[0.93 0.1 0.02 0.8]); % adjust placement
ylabel(h, 'Phase difference');

% Save the combined figure
saveas(fig, fullfile(save_root, 'PhaseDiff_AllLocations.png'));

%% Traveling wave analysis across locations and frequency bands %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wave_results = struct(); % Store traveling wave results

for iloc = 1:length(targ_loc)
    loc     = targ_loc(iloc);
    min_th  = min_difflev(iloc);
    max_th  = max_difflev(iloc);
    
    fprintf('Analyzing traveling wave for location %d\n', loc);
    
    nBands = numel(significant_bands);
    phase_avg = nan(nCh, nBands);  % channels x bands
    
    for ch = 1:nCh
        ch_folder = fullfile(output_folder, 'lfp', ...
            'loc_difflev_all', sprintf('loc%d', loc), ...
            sprintf('%d_%d', min_th, max_th), num2str(ch));
        
        if ~exist(ch_folder, 'dir')
            continue
        end
        
        cd(ch_folder);
        
        if ~exist('coherence.mat','file')
            warning(['Channel ' num2str(ch) ' missing coherence.mat']);
            continue
        end
        
        load('coherence.mat');
        load('coh_perm.mat')
        
        tmaxc = max(coh_perm, [], 2);
        limit_maxc = quantile(tmaxc, 0.95);
        %sig_mask = coh >= limit_maxc; % comment this out to plot without
        %coherence significance mask
        
        for b = 1:nBands
            f_range = significant_bands{b};
            freq_idx = freq >= f_range(1) & freq <= f_range(2);
            
            if any(freq_idx)
                if exist('sig_mask','var')
                    valid_idx = freq_idx & sig_mask;
                else
                    valid_idx = freq_idx; % fallback: use all
                end
                
                if any(valid_idx)
                    phase_data = phase_spec(valid_idx); % only significant freqs
                    phase_c = exp(1i*phase_data);
                    phase_avg(ch,b) = angle(mean(phase_c(:))); % circular mean
                else
                    phase_avg(ch,b) = NaN; % no sig freqs ? exclude channel
                end
            end
        end
        
    end
    
    %% plotting the gradient direction and magnitude - shows phase changes
    
    loc_results = struct();
    
    for b = 1:nBands
        sig_grid_phase = reshape(phase_avg(:,b), [8,8])'; % transpose to match grid
        
        % Wrap phase to [-pi, pi]
        sig_grid_phase = angle(exp(1i*sig_grid_phase));
        
        % Compute gradient
        [phase_dx, phase_dy] = gradient(sig_grid_phase);
        
        % Mean gradient ? wave vector
        mean_dx = mean(phase_dx(:));
        mean_dy = mean(phase_dy(:));
        
        % Magnitude and direction
        wave_magnitude = sqrt(mean_dx^2 + mean_dy^2);
        wave_direction = atan2(mean_dy, mean_dx); % radians, 0 = along X axis
        
        % Store results
        loc_results(b).phase_grid    = sig_grid_phase;
        loc_results(b).wave_magnitude = wave_magnitude;
        loc_results(b).wave_direction = wave_direction;
        loc_results(b).grad_x = phase_dx;
        loc_results(b).grad_y = phase_dy;
        
        %% Visualization
        
        % 1. Phase map with gradient vectors - Local direction of phase
        % change - Wave direction - gives both direction and magnitude
        fig1 = figure;
        imagesc(sig_grid_phase);
        colormap(coolCircularColormap(256));
        colorbar;
        caxis([-pi pi]);
        axis equal tight;
        title(sprintf('Phase & Gradient Vectors - Loc %d Band %d: %.0f-%.0f Hz', ...
            loc, b, significant_bands{b}(1), significant_bands{b}(2)));
        hold on;
        [X,Y] = meshgrid(1:8, 1:8);
        quiver(X, Y, phase_dx, phase_dy, 'k'); % <-- plot gradients
        xlabel('X'); ylabel('Y');
        
        save_folder = fullfile(save_root, ...
            sprintf('Band%d_%.0f-%.0fHz', b, significant_bands{b}(1), significant_bands{b}(2)));
        if ~exist(save_folder, 'dir')
            mkdir(save_folder);
        end
        saveas(fig1, fullfile(save_folder, ...
            sprintf('PhaseGradient_Loc%d_Band%d.png', loc, b)));
        
        % 2. Gradient magnitude heatmap - shows Local strength of phase
        % change (Wavefront steepness)
        grad_magnitude = sqrt(phase_dx.^2 + phase_dy.^2);
        
        fig2 = figure;
        imagesc(grad_magnitude);
        colorbar;
        axis equal tight;
        title(sprintf('Gradient Magnitude - Loc %d Band %d: %.0f-%.0f Hz', ...
            loc, b, significant_bands{b}(1), significant_bands{b}(2)));
        xlabel('X'); ylabel('Y');
        
        saveas(fig2, fullfile(save_folder, ...
            sprintf('GradientMag_Loc%d_Band%d.png', loc, b)));
        
        close(fig1); close(fig2);
    end
    
    wave_results(iloc).location = loc;
    wave_results(iloc).bands = loc_results;
end

disp('Traveling wave analysis completed.');












