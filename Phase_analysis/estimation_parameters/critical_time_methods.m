clear all
close all
clc

%% paths

addpath /opt/fieldtrip_github/
ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf

%% Create data paths to use with slurm

datafolder = '/mnt/hpc/projects/MWSampling/4Shivangi/results';

session_names = [];
session_names = {%'klecks_20170804_attentional-sampling_1';...
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
    %'klecks_20170908_attentional-sampling_1';...
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

%% check time length

for isess = 1:length(session_names)
    cd(output_paths{isess});
    load('zlfpTrials.mat');
    
    % Display the dimensions of zlfpTrials.time
    fprintf('Session %d: zlfpTrials.time has dimensions [%d x %d]\n', ...
        isess, size(zlfpTrials.time, 1), size(zlfpTrials.time{1,1}, 2));
end

%% Find all unique channels present across sessions

all_channels = {};

for isess = 1:length(session_names)
    if exist(session_paths_files{isess}, 'file')
        load(session_paths_files{isess}, 'clean_data');
        all_channels = unique([all_channels; clean_data.label]);
    end
end

%% Collect Channel-wise LFP

collect_trials = struct();

for ichan = 1:length(all_channels)
    ichan
    channel_name = all_channels{ichan};
    
    % Ensure valid field name
    valid_channel_name = matlab.lang.makeValidName(channel_name);
    
    % Initialize an empty cell array for the current channel
    if ~isfield(collect_trials, valid_channel_name)
        collect_trials.(valid_channel_name) = {};
    end
    
    for isess = 1:length(session_names)
        cd(output_paths{isess});
        
        % Load LFP data
        load('zlfpTrials.mat');
        
        % Check if this session has the desired channel
        if any(strcmp(zlfpTrials.label, channel_name))
            disp(['Processing ', channel_name, ' in session ', session_names{isess}]);
            
            % Get unique difficulty levels for this session
            unique_DL = flip(unique(zlfpTrials.trialinfo(:,18)));
            num_DL = length(unique_DL);
            first_quarter_DL = unique_DL(1:ceil(num_DL / 4));
            
            % Select trials for the current channel
            cfg = [];
            cfg.channel = channel_name;
            cfg.trials = find(ismember(zlfpTrials.trialinfo(:,18), first_quarter_DL));
            lfpTrials = ft_selectdata(cfg, zlfpTrials);
            
            % Select hit trials
            cfg = [];
            cfg.trials = find(lfpTrials.trialinfo(:,20) == 1);
            itcTrials = ft_selectdata(cfg, lfpTrials);
            
            % Skip session if no valid trials
            if isempty(itcTrials.trial)
                disp(['Skipping session ', session_names{isess}, ' due to no valid trials.']);
                continue;
            end
            
            % Collect trials per channel
            for tr = 1:length(itcTrials.trial)
                if ~isempty(itcTrials.trial{tr})
                    collect_trials.(valid_channel_name){end+1} = itcTrials.trial{tr};
                else
                    warning(['Skipping empty trial ', num2str(tr), ' in session ', session_names{isess}]);
                end
            end
        end
    end
end

%%% Ni et al implementation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate ITC - freq wise

frequencies = 6:2:18;
sample_rate = 1000;

num_channels = length(fieldnames(collect_trials));
channel_fields = fieldnames(collect_trials);

hFig = figure;
idx = 1;

while true
    clf(hFig);
    
    % Get channel name and trial data
    channel_name = channel_fields{idx};
    trial_data = collect_trials.(channel_name);
    
    if iscell(trial_data)
        trial_data = vertcat(trial_data{:});
    end
    
    % Compute ITC
    itc = calculate_itc(trial_data, frequencies, sample_rate);
    
    % Ensure no NaNs
    itc(isnan(itc)) = 0;
    
    % Plot ITC for the current channel
    imagesc(itcTrials.time{1,1}, flip(frequencies), flipud(itc));
    set(gca, 'YDir', 'normal');
    colorbar;
    
    xlabel('Time (samples)');
    ylabel('Frequency (Hz)');
    title(['ITC - Channel: ', channel_name]);
    
    % Wait for key press
    waitforbuttonpress;
    key = get(hFig, 'CurrentCharacter');
    
    % Navigate channels
    if key == char(29) % Right arrow key
        idx = min(idx + 1, num_channels);
    elseif key == char(28) % Left arrow key
        idx = max(idx - 1, 1);
    elseif key == char(27) % Escape key to exit
        break;
    end
end

close(hFig);

%% Calculate ITC - all freq

frequencies = 6:2:18;
sample_rate = 1000;

num_channels = length(fieldnames(collect_trials));
channel_fields = fieldnames(collect_trials);

hFig = figure;
idx = 1;

while true
    clf(hFig);
    
    % Get channel name and trial data
    channel_name = channel_fields{idx};
    trial_data = collect_trials.(channel_name);
    
    if iscell(trial_data)
        trial_data = vertcat(trial_data{:});
    end
    
    % Compute ITC
    itc_f = calculate_itc(trial_data, frequencies, sample_rate);
    itc = mean(itc_f);
    
    % Ensure no NaNs
    itc(isnan(itc)) = 0;
    
    % Plot ITC vs. time for the current channel
    plot(itcTrials.time{1,1}, itc);
    xlabel('Time');
    ylabel('ITC');
    title(['ITC vs. Time - Channel: ', channel_name], 'Interpreter', 'none');
    
    % Wait for key press
    waitforbuttonpress;
    key = get(hFig, 'CurrentCharacter');
    
    % Navigate channels
    if key == char(29) % Right arrow key
        idx = min(idx + 1, num_channels);
    elseif key == char(28) % Left arrow key
        idx = max(idx - 1, 1);
    elseif key == char(27) % Escape key to exit
        break;
    end
end

close(hFig);

%% permutation test for itc (all freq)

frequencies = 6:2:18;
sample_rate = 1000;

num_channels = length(fieldnames(collect_trials));
channel_fields = fieldnames(collect_trials);
num_permutations = 1000;
max_perm_values = zeros(num_permutations, num_channels);
min_perm_values = zeros(num_permutations, num_channels);

for idx = 1:num_channels
    idx
    
    % Get channel name and trial data
    channel_name = channel_fields{idx};
    trial_data = collect_trials.(channel_name);
    
    if iscell(trial_data)
        trial_data = vertcat(trial_data{:});
    end
    
    % Permutation Test for ITC Significance
    [n_trials, trial_length] = size(trial_data);
    
    for perm = 1:num_permutations
        
        for i = 1:n_trials
            permuted_data(i, :) = trial_data(i, randperm(trial_length));
        end
        
        % Compute ITC
        permuted_itc_f = calculate_itc(permuted_data, frequencies, sample_rate);
        permuted_itc = mean(permuted_itc_f);
        permuted_itc(isnan(permuted_itc)) = 0;
        max_perm_values(perm, idx) = max(permuted_itc);
        min_perm_values(perm, idx) = min(permuted_itc);
        
    end
    
end

cd('/mnt/hpc/projects/MWSampling/4Shivangi/results/ITC')
save max_perm_values max_perm_values
save min_perm_values min_perm_values

%% plotting with statistics

cd('/mnt/hpc/projects/MWSampling/4Shivangi/results/ITC')
load('max_perm_values')
load('min_perm_values')
load('collect_trials')
load('time')

limit_max = quantile([max_perm_values], 0.975, 1);
limit_min = quantile([min_perm_values], 0.025, 1);

%limit_max = quantile([max_perm_values;min_perm_values], 0.975, 1);
%limit_min = quantile([max_perm_values;min_perm_values], 0.025, 1);

frequencies = 6:2:18;
sample_rate = 1000;

num_channels = length(fieldnames(collect_trials));
channel_fields = fieldnames(collect_trials);

hFig = figure;
idx = 1;

while true
    clf(hFig);
    
    % Get channel name and trial data
    channel_name = channel_fields{idx};
    trial_data = collect_trials.(channel_name);
    
    if iscell(trial_data)
        trial_data = vertcat(trial_data{:});
    end
    
    % Compute ITC
    itc_f = calculate_itc(trial_data, frequencies, sample_rate);
    itc = mean(itc_f);
    
    % Ensure no NaNs
    itc(isnan(itc)) = 0;
    
    % Plot ITC vs. time for the current channel
    plot(time, itc, 'k');
    hold on;
    
    % Plot limit max and min as horizontal lines
    yline(limit_max(idx), '--m', 'Max');
    yline(limit_min(idx), '--b', 'Min');
    
    xlabel('Time');
    ylabel('ITC');
    title(['ITC vs. Time - Channel: ', channel_name], 'Interpreter', 'none');
    legend('ITC', 'Max Limit', 'Min Limit');
    hold off;
    
    % Wait for key press
    waitforbuttonpress;
    key = get(hFig, 'CurrentCharacter');
    
    % Navigate channels
    if key == char(29) % Right arrow key
        idx = min(idx + 1, num_channels);
    elseif key == char(28) % Left arrow key
        idx = max(idx - 1, 1);
    elseif key == char(27) % Escape key to exit
        break;
    end
end

close(hFig);

%% ERP based method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

channel_fields = fieldnames(collect_trials);
bsltimewindow_start_ind = 1;
bsltimewindow_end_ind = 1000;

n_consec = 5;

for ichan = 1:length(all_channels)
    channel_name = channel_fields{ichan};
    trial_data = collect_trials.(channel_name);
    
    if iscell(trial_data)
        trial_data = vertcat(trial_data{:});
    end
    
    % mean subtraction and threshold
    bsldata = trial_data(:, bsltimewindow_start_ind:bsltimewindow_end_ind);
    bsln_mean = mean(bsldata,2);
    SD = std(bsldata(:));
    threshold = bsln_mean + 5 * SD;
    erp_data = trial_data - bsln_mean;
    
    % finding the critical time
    n_trials = size(erp_data, 1);
    n_timepoints = size(erp_data, 2);
    
    timeofresponse_ind = nan(n_trials, 1);
    
    for trial = 1:n_trials
        search_range = bsltimewindow_end_ind+5 : 1101; % after stim on and within 100ms after stim on
        above_thresh = abs(erp_data(trial, search_range)) > threshold(trial);
        consec_counts = conv(double(above_thresh), ones(1, n_consec), 'valid');
        ind = find(consec_counts == n_consec, 1);
        
        if ~isempty(ind)
            timeofresponse_ind(trial) = search_range(ind);  % map back to original time indices
        end
    end
    mean_timeofresponse(1,ichan) = mean(timeofresponse_ind, 'omitnan');
    
end

% critical time per channel
indices = round(mean_timeofresponse);
timeVec = time;
CriticalTime = NaN(1, length(all_channels));
validIdx = indices > 0 & indices <= length(timeVec);
CriticalTime(validIdx) = timeVec(indices(validIdx));

% Histogram
valid_critical_times = CriticalTime(~isnan(CriticalTime));
figure;
n_bins = 20;
h = histogram(valid_critical_times, n_bins, 'FaceColor', [0.2, 0.6, 0.8], 'EdgeColor', 'k');
xlabel('Time from target onset(s)');
ylabel('Number of channels');
title('Critical Time');
meanv = mean(valid_critical_times, 'omitnan');
modev = mode(valid_critical_times(~isnan(valid_critical_times)));
medianv = median(valid_critical_times(~isnan(valid_critical_times)));
xline(meanv, '--k', 'LineWidth', 1.5);
xline(modev, '--b', 'LineWidth', 1.5);
xline(medianv, '--m', 'LineWidth', 1.5);

%% Critical time for each channel - threshold on average of trials - SD of mean of trials

channel_fields = fieldnames(collect_trials);
bsltimewindow_start_ind = 700;
bsltimewindow_end_ind = 1000;

n_consec = 15;  
th = 2; 
timeofresponse_ind = nan(length(all_channels), 1);
figure;
for ichan = 1:length(all_channels)
    channel_name = channel_fields{ichan};
    trial_data = collect_trials.(channel_name);
    
    if iscell(trial_data)
        trial_data = vertcat(trial_data{:});
    end
    
    % mean subtraction and threshold
    bsldata = trial_data(:, bsltimewindow_start_ind:bsltimewindow_end_ind);
    bsln_mean = mean(bsldata(:)); 
    
    mean_erp_data = mean((trial_data - bsln_mean),1);
    sd_data = mean_erp_data(:, bsltimewindow_start_ind:bsltimewindow_end_ind);
    SD = std(sd_data(:));
    
    threshold = bsln_mean + (th * SD);
    
    % finding critical time
    search_range = bsltimewindow_end_ind+20 : 1101; % after stim on and within 100ms after stim on
    above_thresh = abs(mean_erp_data(1, search_range)) > threshold;
    consec_counts = conv(double(above_thresh), ones(1, n_consec), 'valid');
    ind = find(consec_counts == n_consec, 1);
    
    if ~isempty(ind)
        timeofresponse_ind(ichan,1) = search_range(ind);  % map back to original time indices
    end
    
    subplot(8,8,ichan)
    plot(time(1,700:1400),mean_erp_data(1,700:1400))
    hold on
    yline(threshold, '--k', 'LineWidth', 1.5);
    if ~isnan(timeofresponse_ind(ichan))
        xline(time(timeofresponse_ind(ichan)), '--k', 'LineWidth', 1.5);
    end
    
end

% critical time per channel
indices = timeofresponse_ind;
timeVec = time;
CriticalTime = NaN(1, length(all_channels));
validIdx = indices > 0 & indices <= length(timeVec);
CriticalTime(validIdx) = timeVec(indices(validIdx));

% Histogram
valid_critical_times = CriticalTime(~isnan(CriticalTime));
figure;
n_bins = 20;
h = histogram(valid_critical_times, n_bins, 'FaceColor', [0.2, 0.6, 0.8], 'EdgeColor', 'k');
xlabel('Time from target onset(s)');
ylabel('Number of channels');
title('Critical Time');
meanv = mean(valid_critical_times, 'omitnan');
modev = mode(valid_critical_times(~isnan(valid_critical_times)));
medianv = median(valid_critical_times(~isnan(valid_critical_times)));
xline(meanv, '--k', 'LineWidth', 1.5);
xline(modev, '--b', 'LineWidth', 1.5);
xline(medianv, '--m', 'LineWidth', 1.5);

%% Visualize to check - avg over trials

channel_fields = fieldnames(collect_trials);
bsltimewindow_start_ind = 1;
bsltimewindow_end_ind = 1000;
n_consec = 5;

ichan = 1;  % Start with the first channel
n_channels = length(channel_fields);  % Total number of channels
ct_select = 0.06;

% Create a figure to display the plots
hFig = figure('KeyPressFcn', @onKeyPress);

% Store the initial 'ichan' and 'n_channels' in the figure's app data
setappdata(hFig, 'ichan', ichan);
setappdata(hFig, 'n_channels', n_channels);

while ishandle(hFig)
    % Get the current channel from the figure's app data
    ichan = getappdata(hFig, 'ichan');
    n_channels = getappdata(hFig, 'n_channels');  % Get total number of channels
    
    channel_name = channel_fields{ichan};
    trial_data = collect_trials.(channel_name);
    
    if iscell(trial_data)
        trial_data = vertcat(trial_data{:});
    end
    
    % mean subtraction and threshold
    bsldata = trial_data(:, bsltimewindow_start_ind:bsltimewindow_end_ind);
    bsln_mean = mean(bsldata, 2);
    SD = std(bsldata(:));
    threshold = bsln_mean + 5 * SD;
    erp_data = trial_data - bsln_mean;
    
    % finding the critical time
    n_trials = size(erp_data, 1);
    n_timepoints = size(erp_data, 2);
    timeofresponse_ind = nan(n_trials, 1);
    
    for trial = 1:n_trials
        search_range = bsltimewindow_end_ind + 5 : 1101; % after stim on and within 100ms after stim on
        above_thresh = (erp_data(trial, search_range)) > threshold(trial);
        consec_counts = conv(double(above_thresh), ones(1, n_consec), 'valid');
        ind = find(consec_counts == n_consec, 1);
        
        if ~isempty(ind)
            timeofresponse_ind(trial) = search_range(ind);  % map back to original time indices
        end
    end
    
    avg_trial_data = mean(trial_data, 1);
    mean_timeofresponse = mean(timeofresponse_ind, 'omitnan');
    
    % Plot the current channel's ERP data
    clf(hFig);
    plot(zlfpTrials.time{1, 1}, avg_trial_data, 'LineWidth', 2)
    hold on
    
    % Plot the vertical line
    %     if ~isnan(mean_timeofresponse)
    %         xline(zlfpTrials.time{1, 1}(1, round(mean_timeofresponse)), 'k--', 'LineWidth', 2);
    %
    %         text(zlfpTrials.time{1, 1}(1, round(mean_timeofresponse)), max(avg_trial_data), sprintf('%.2f ms', zlfpTrials.time{1, 1}(1, round(mean_timeofresponse))), ...
    %             'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
    %             'FontSize', 10, 'Color', 'black');
    %     end
    
    % Plot with the selected critical time
    if ~isnan(mean_timeofresponse)
        xline(ct_select, 'k--', 'LineWidth', 2);
        
        text(ct_select, max(avg_trial_data), sprintf('%.2f ms', ct_select), ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
            'FontSize', 10, 'Color', 'black');
    end
    
    % Channel name
    title(['Channel: ' channel_name], 'Interpreter', 'none');
    xlabel('Time (ms)')
    ylabel('ERP')
    grid on;
    
    % Pause for a short time before waiting for the key press
    pause(0.2);
end

% Callback function for key press
function onKeyPress(~, event)
% Get the current channel index and number of channels from the figure's app data
ichan = getappdata(gcf, 'ichan');
n_channels = getappdata(gcf, 'n_channels');

if strcmp(event.Key, 'rightarrow')  % Right arrow key
    ichan = mod(ichan, n_channels) + 1;  % Move to next channel
elseif strcmp(event.Key, 'leftarrow')  % Left arrow key
    ichan = mod(ichan - 2, n_channels) + 1;  % Move to previous channel
elseif strcmp(event.Key, 'escape')  % Escape key to exit the loop
    close(gcf);  % Close the figure and exit the loop
end

% Update the channel index in the figure's app data
setappdata(gcf, 'ichan', ichan);
end

%%
% %% Visualize to check - per trial - too noisy - try average
%
% % --- Setup ---
% channel_fields = fieldnames(collect_trials);
% channel_name = channel_fields{4};  % Use first channel
% trial_data = collect_trials.(channel_name);
%
% if iscell(trial_data)
%     trial_data = vertcat(trial_data{:});
% end
%
% bsltimewindow_start_ind = 1;
% bsltimewindow_end_ind = 1000;
% n_consec = 3;
%
% % --- Threshold computation ---
% bsldata = trial_data(:, bsltimewindow_start_ind:bsltimewindow_end_ind);
% bsln_mean = mean(bsldata, 2);
% SD = std(bsldata(:));
% threshold = bsln_mean + 3 * SD;
%
% erp_data = trial_data - bsln_mean;
% n_trials = size(erp_data, 1);
% n_timepoints = size(erp_data, 2);
%
% % --- Detect response time per trial ---
% timeofresponse_ind = nan(n_trials, 1);
%
% for trial = 1:n_trials
%     search_range = bsltimewindow_end_ind+1 : 1131;  % Post-stim and within 100 ms
%     above_thresh = abs(erp_data(trial, search_range)) > threshold(trial);
%     consec_counts = conv(double(above_thresh), ones(1, n_consec), 'valid');
%     ind = find(consec_counts == n_consec, 1);
%
%     if ~isempty(ind)
%         timeofresponse_ind(trial) = search_range(ind);  % Map back to full index
%     end
% end
%
% % --- Initialize interactive viewer ---
% current_trial = 1;
%
% % Create plot for first trial
% hFig = figure('Name', 'ERP Trial Viewer');
% hPlot = plot(trial_data(current_trial, :), 'b');
% hold on;
% hLine = xline(timeofresponse_ind(current_trial), '--k', 'LineWidth', 1.5);
% title(['Trial ', num2str(current_trial)]);
% xlabel('Time (samples)');
% ylabel('ERP amplitude');
% adjustY(trial_data, current_trial);  % Call adjustY with necessary data
%
% % Store variables in guidata along with the figure handle and n_trials
% guidata(hFig, struct('current_trial', current_trial, 'erp_data', erp_data, ...
%     'timeofresponse_ind', timeofresponse_ind, 'hPlot', hPlot, 'hLine', hLine, ...
%     'n_trials', n_trials));
%
% % Set callback after creating everything (so variables are in scope)
% set(hFig, 'KeyPressFcn', @arrowKeyCallback);
%
% % --- Key press callback function ---
%     function arrowKeyCallback(~, event)
%         % Retrieve the current state from guidata using hFig
%         data = guidata(gcf());  % Use gcf() to get the current figure handle
%         current_trial = data.current_trial;
%         erp_data = data.erp_data;
%         timeofresponse_ind = data.timeofresponse_ind;
%         hPlot = data.hPlot;
%         hLine = data.hLine;
%         n_trials = data.n_trials;  % Get n_trials from guidata
%
%         switch event.Key
%             case 'rightarrow'
%                 current_trial = min(current_trial + 1, n_trials);
%             case 'leftarrow'
%                 current_trial = max(current_trial - 1, 1);
%             otherwise
%                 return;
%         end
%
%         % Update ERP trace
%         set(hPlot, 'YData', erp_data(current_trial, :));
%
%         % Update vertical line
%         if ~isnan(timeofresponse_ind(current_trial))
%             set(hLine, 'Value', timeofresponse_ind(current_trial), 'Visible', 'on');
%         else
%             set(hLine, 'Visible', 'off');
%         end
%
%         adjustY(erp_data, current_trial);  % Update Y limits
%         title(['Trial ', num2str(current_trial)]);
%
%         % Save the updated state back into guidata
%         data.current_trial = current_trial;
%         guidata(gcf(), data);  % Use gcf() to get the figure handle
%     end
%
% % --- Adjust Y-limits function ---
%     function adjustY(erp_data, current_trial)
%         y = erp_data(current_trial, :);
%         padding = range(y) * 0.1;
%         ylim([min(y)-padding, max(y)+padding]);
%     end

