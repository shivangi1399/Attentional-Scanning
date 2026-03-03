%% TFRs to find frequency around the ERPs (LFP data)

clear all;
close all;
clc;

%% paths

addpath /opt/fieldtrip_github/
ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf

%% Create data paths to use with slurm

datafolder = '/mnt/hpc/projects/MWSampling/4Shivangi/results';

session_names = [];
session_names = {'klecks_20170804_attentional-sampling_1';...
    'klecks_20170807_attentional-sampling_1';...
    'klecks_20170808_attentional-sampling_1';... 
    'klecks_20170810_attentional-sampling_1';...
    'klecks_20170817_attentional-sampling_1';...
    'klecks_20170818_attentional-sampling_1';...
    'klecks_20170821_attentional-sampling_1';...
    'klecks_20170822_attentional-sampling_1';...
    'klecks_20170823_attentional-sampling_1';...
    'klecks_20170824_attentional-sampling_1';...
    'klecks_20170825_attentional-sampling_1';...
    'klecks_20170828_attentional-sampling_1';...
    'klecks_20170830_attentional-sampling_1';...
    'klecks_20170831_attentional-sampling_1';...
    'klecks_20170901_attentional-sampling_1';...
    'klecks_20170904_attentional-sampling_1';...
    'klecks_20170906_attentional-sampling_1';...
    'klecks_20170908_attentional-sampling_1';...
    'klecks_20170911_attentional-sampling_1';...
    'klecks_20170913_attentional-sampling_1';...
    'klecks_20170914_attentional-sampling_1';...
    'klecks_20170915_attentional-sampling_1';...
    'klecks_20170919_attentional-sampling_1';...
    'klecks_20171020_attentional-sampling_1'};

% create cell paths to use with slurm
session_paths = [];
session_paths = cellfun(@(x) fullfile(datafolder,x), session_names, 'uniform',0);

session_paths_files = [];
session_paths_files = cellfun(@(x) fullfile(datafolder,x, 'clean_data.mat'), session_names, 'uniform',0);

% path to save (slurm) output
output_folder = '/mnt/hpc/projects/MWSampling/4Shivangi/results';
output_paths = cellfun(@(x) fullfile(output_folder, x, 'TFR'),session_names, 'uniform',0);

%% zscoring

for isess = 1:length(session_names)
    
    % load LFP data
    ESIload(session_paths_files{isess});
    disp(strcat('session- ',num2str(isess),' out of- ',num2str(length(session_names)), ', running LFP z-scoring'))
    
    zlfpTrials = fun_zscore_session(clean_data);
    
    % save zscored data
    if ~isdir(output_paths{isess})
        mkdir(output_paths{isess})
    end
    
    cd(output_paths{isess}),
    save('zlfpTrials','zlfpTrials')
end

%% TFRs for each channel across sessions

%% Find all unique channels present across sessions

all_channels = {};

for isess = 1:length(session_names)
    if exist(session_paths_files{isess}, 'file')
        load(session_paths_files{isess}, 'clean_data');
        all_channels = unique([all_channels; clean_data.label]);
    end
end

%% Compute TFR for each channel

num_channels = length(all_channels);
num_frequencies = 50;
num_timepoints = 33;
TFR_channels = nan(num_channels, num_frequencies, num_timepoints); % Storage matrix

for ichan = 1%:num_channels
    channel_name = all_channels{ichan};
    TFR_data = []; % Store session-wise data for this channel
    
    for isess = 1%:length(session_names)
        cd(output_paths{isess}),
        
        % load LFP data
        load('zlfpTrials.mat')
        dat = zlfpTrials;
        
        % Check if this session has the desired channel
        if any(strcmp(dat.label, channel_name))
            disp(['Processing ', channel_name, ' in session ', session_names{isess}]);
            
            % Extract data for the current channel
            cfg = [];
            cfg.channel = channel_name;
            data = ft_selectdata(cfg, dat);
            
            % Compute TFR
            cfg              = [];
            cfg.output       = 'pow';
            cfg.method       = 'mtmconvol';
            cfg.taper        = 'hanning';
            cfg.foi          = 2:2:100;
            cfg.t_ftimwin    = 3 ./ cfg.foi;
            cfg.toi          = -1:0.05:0.6;
            freqpow = ft_freqanalysis(cfg, data);
            
            % Store power spectrum
            TFR_data = cat(1, TFR_data, freqpow.powspctrm);
        end
    end
    
    %% Store average across available sessions
    if ~isempty(TFR_data)
        TFR_channels(ichan, :, :) = squeeze(mean(TFR_data, 1));
    end
    
end

%% plot

[num_channels, num_freqs, num_time] = size(TFR_channels);

% Compute global min and max across all channels
global_min = min(TFR_channels(:));
global_max = max(TFR_channels(:));

% Define subplot grid size
num_rows = ceil(sqrt(num_channels)); 
num_cols = ceil(num_channels / num_rows); 

figure;
tiledlayout(num_rows, num_cols, 'Padding', 'none', 'TileSpacing', 'none'); 

for ch = 1:num_channels
    row = ceil(ch / num_cols); 
    col = mod(ch - 1, num_cols) + 1; 
    
    nexttile;
    imagesc(freqpow.time, freqpow.freq, squeeze(TFR_channels(ch, :, :)));
    set(gca, 'YDir', 'normal'); 
    caxis([global_min, global_max]); % Set uniform color scale
    
    % Hide axis labels by default
    xticklabels([]); 
    yticklabels([]); 
    xticks([]); 
    yticks([]);
    
    % Show x-axis labels on the bottom row
    if row == num_rows
        xticks(linspace(min(freqpow.time), max(freqpow.time), 5)); 
        xticklabels(round(linspace(min(freqpow.time), max(freqpow.time), 5), 1)); 
        xlabel('Time');
    end
    
    % Show y-axis labels on the leftmost column
    if col == 1
        yticks(linspace(min(freqpow.freq), max(freqpow.freq), 5)); 
        yticklabels(round(linspace(min(freqpow.freq), max(freqpow.freq), 5), 1)); 
        ylabel('Frequency');
    end
end

% Add one shared colorbar
cb = colorbar;
cb.Layout.Tile = 'east'; 
colormap jet; 
sgtitle('Time-Frequency Representations for All Channels');  

%% plotting individual channels

ITC_TFR(TFR_channels, freqpow);

function ITC_TFR(TFR_channels, freqpow)

    [num_channels, num_freqs, num_time] = size(TFR_channels);

    % Compute global min and max across all channels
    %global_min = min(TFR_channels(:));
    %global_max = max(TFR_channels(:));

    % Create figure
    fig = figure;
    colormap jet; 

    % Initialize channel index
    channel_idx = 1; 

    % Display the first channel
    update_plot();

    % Assign key press function to figure
    set(fig, 'KeyPressFcn', @keyPress);

    % ======= Function to Update the Plot =======
    function update_plot()
        clf; % Clear figure before updating
        imagesc(freqpow.time, freqpow.freq, squeeze(TFR_channels(channel_idx, :, :)));
        set(gca, 'YDir', 'normal');
        %caxis([global_min, global_max]); % Set uniform color scale
        colorbar;

        % Set x and y ticks
        xticks(linspace(min(freqpow.time), max(freqpow.time), 10)); 
        xticklabels(round(linspace(min(freqpow.time), max(freqpow.time), 10), 1)); 
        xlabel('Time');

        yticks(linspace(min(freqpow.freq), max(freqpow.freq), 10)); 
        yticklabels(round(linspace(min(freqpow.freq), max(freqpow.freq), 10), 1)); 
        ylabel('Frequency');

        title(['Channel ' num2str(channel_idx) ' of ' num2str(num_channels)]);
        drawnow;
    end

    % ======= Function to Handle Key Presses =======
    function keyPress(~, event)
        if strcmp(event.Key, 'rightarrow') % Move to next channel
            if channel_idx < num_channels
                channel_idx = channel_idx + 1;
                update_plot();
            end
        elseif strcmp(event.Key, 'leftarrow') % Move to previous channel
            if channel_idx > 1
                channel_idx = channel_idx - 1;
                update_plot();
            end
        end
    end

end


