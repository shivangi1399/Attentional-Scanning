clear all
close all
clc

% Code description:
% -----------------
% Frequency Bands: Defines frequency bands for analysis  (eg. 4-40 Hz, 41-80 Hz)
% Average Phase: Computes circular mean phase per channel and band for each location across sessions
% Phase vs Location Video: Generates a video showing per-channel phase changes across locations for each band.
% Phase vs Location Average: Plots average phase across channels for each band across all locations.

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

output_folder = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks/phase_coherence/cp10_till_100';
data_folder   = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks/multi_lin_reg/cp10_till_100';

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
save_root = '/mnt/hpc/projects/MWSampling/4Shivangi/Plots/phase_coherence/klecks/cp10_till_100/lfp/loc_difflev_all';
if ~exist(save_root, 'dir')
    mkdir(save_root);
end

%% Define frequency bands

% alpha, beta and gamma bands
% significant_bands = {
%     [3, 7]             %1: 4-8, 2:3-7 ,3:3-8
%     [8, 12];           %2: 8-12
%     [13, 30]; 
%     [31, 80];
%     };

significant_bands = {
    [4, 40]; 
    [41, 80];
    };

% 10Hz bands
% startFreq = 2;
% endFreq = 80;
% bandWidth = 10;
% 
% % Create bands
% edges = startFreq:bandWidth:endFreq;
% significant_bands = cell(length(edges)-1,1);
% 
% for i = 1:length(edges)-1
%     significant_bands{i} = [edges(i), edges(i+1)];
% end

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
                    phase_c = exp(1i*phase_data); %try complex notation?
                    phase_avg(ch,b) = angle(mean(phase_c(:))); % circular mean
                else
                    phase_avg(ch,b) = NaN; % no sig freqs ? exclude channel
                end
            end
        end
    end
    
    phase_avg_all{iloc} = phase_avg;
    
    %% Visualize as 8x8 grid per band - just visualizing phase values
%     for b = 1:nBands
%         sig_grid_phase = reshape(phase_avg(:,b), [8,8])'; % transpose to match previous plots
%         
%         fig = figure;  % store figure handle
%         imagesc(sig_grid_phase);
%         colorbar;
%         colormap(coolCircularColormap(256));
%         caxis([phase_min phase_max]); % consistent scale
%         title(sprintf('Avg Phase - Loc %d Band %d: %.0f-%.0f Hz', loc, b, significant_bands{b}(1), significant_bands{b}(2)));
%         axis equal tight;
%         hold on;
%         [X,Y] = meshgrid(1:8, 1:8);
%         quiver(X, Y, cos(sig_grid_phase), sin(sig_grid_phase), 0.5, 'k');
%         xticks(1:8); yticks(1:8);
%         xlabel('Channel X'); ylabel('Channel Y');
%         
%         save_folder = fullfile(save_root, ...
%             sprintf('Band%d_%.0f-%.0fHz', b, significant_bands{b}(1), significant_bands{b}(2)));
%         if ~exist(save_folder, 'dir')
%             mkdir(save_folder);
%         end
%         
%         % Use figure handle when saving
%         saveas(fig, fullfile(save_folder, ...
%             sprintf('AvgPhase_Loc_nosig%d_Band%d_%.0f-%.0fHz.png', loc, b, significant_bands{b}(1), significant_bands{b}(2))));
%         
%         % Close the same figure handle
%         close(fig);
%     end
    
end

disp('Average phase computed for all locations and bands with consistent scale.');

%% Phase vs location %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

videoFile = fullfile(save_root, 'Phase_vs_Location.mp4');
v = VideoWriter(videoFile, 'Motion JPEG AVI');
v.FrameRate = 2; % frames per second
open(v);

% Loop through channels and plot
nCh = size(phase_avg_all{1},1); % total channels
nBands = numel(significant_bands);
nLocs = length(targ_loc);

% Define subplot grid size (square-ish layout)
nRows = ceil(sqrt(nBands));
nCols = ceil(nBands / nRows);

legendLabels = cell(1,nBands);
for b = 1:nBands
    legendLabels{b} = sprintf('Band %d: %.0f-%.0f Hz', ...
        b, significant_bands{b}(1), significant_bands{b}(2));
end

fig = figure('Name','Phase vs Location','NumberTitle','off','Position',[100 100 1400 900]);

for ch = 1:nCh
    clf(fig); % clear figure each frame
    
    for b = 1:nBands
        ax = subplot(nRows, nCols, b);
        hold(ax,'on');
        
        % Compute phase values for this channel and band
        phase_vals = nan(nLocs,1);
        for iloc = 1:nLocs
            phase_vals(iloc) = phase_avg_all{iloc}(ch,b);
        end
        
        % Plot
        plot(ax, 1:nLocs, phase_vals, '-o', 'Color', 'b', 'LineWidth', 2);
        
        % Labels and formatting
        xlabel(ax,'Location'); 
        ylabel(ax,'Phase (rad)');
        xticks(ax,1:nLocs); 
        xticklabels(ax,targ_loc);
        title(ax, legendLabels{b});
        grid(ax,'on');
    end
    
    % Add a super title for the channel
    sgtitle(sprintf('Channel %d', ch));
    
    drawnow;
    frame = getframe(fig);
    writeVideo(v, frame);
end

close(v);
disp(['Video saved as ' videoFile]);

%% Phase vs location - avg across channels

nBands = numel(significant_bands);
nLocs = length(targ_loc);

% Compute subplot layout automatically (square-like)
nRows = ceil(sqrt(nBands));
nCols = ceil(nBands / nRows);

f = figure('Name','Phase vs Location - Average across Channels','NumberTitle','off','Position',[100 100 1400 900]);

for b = 1:nBands
    ax = subplot(nRows, nCols, b);  % auto layout
    hold(ax,'on');
    
    % Compute average across channels
    phase_vals = nan(nLocs,1);
    for iloc = 1:nLocs
        phase_data = phase_avg_all{iloc}(:,b);
        phase_vals(iloc) = angle(mean(exp(1i*phase_data), 'omitnan'));
    end
    
    % Plot
    plot(ax, 1:nLocs, phase_vals, '-o', 'Color', 'b', 'LineWidth', 2);
    
    % Labels and formatting
    xlabel(ax,'Location');
    ylabel(ax,'Phase (rad)');
    xticks(ax,1:nLocs); 
    xticklabels(ax,targ_loc);
    title(ax, sprintf('Band %d: %.0f-%.0f Hz', b, ...
        significant_bands{b}(1), significant_bands{b}(2)));
    grid(ax,'on');
end

% Save figure
saveas(f, fullfile(save_root, 'Phase_loc_allChan_40HzBands.png'));
%saveas(f, fullfile(save_root, 'Phase_loc_allChan_10HzBands.png'));


