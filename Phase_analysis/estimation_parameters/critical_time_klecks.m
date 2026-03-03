clear all
close all
clc

% Code description:
% -----------------
% Code to find critical time or time when stimulus reaches the brain
% This is done on artifact rejected and zscored data
% Critical time was calculated on linear regression subtracted erp 
% threshold was done on average of trials - SD of the mean was used

%% paths

addpath /opt/fieldtrip_github/
ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
clc

%% paths

datafolder   = '/mnt/hpc/projects/MWSampling/4Shivangi/data_Klecks';

cd(datafolder),
animalName = 'klecks';
temp = dir;
session_names = [];
ii = 0;
for i = 1:length(temp)
    if strfind(temp(i).name,animalName)
        ii = ii+1;
        session_names{ii,1} = temp(i).name;
    end
end

% paths to RS4 preprocessed data
data_path = cellfun(@(S) fullfile(datafolder, S), session_names, 'Uniform', 0);

% Path to cleaned data
output_folder = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks';
output_path = cellfun(@(x) fullfile(output_folder, x),session_names, 'uniform',0);

%% use artifact rejected data and zscore sessions

for isess = 1:length(session_names)
    
    % load LFP data
    cd(data_path{isess})
    
    % check if file exists
    if ~isfile('lfpTrials_cleanf.mat')
        fprintf('Skipping %s: lfpTrials_cleanf.mat not found\n', session_names{isess});
        continue
    end
    
    ESIload('lfpTrials_cleanf.mat');
    
    % create session folder
    if ~isdir(fullfile(output_path{isess}))
        mkdir(fullfile(output_path{isess}))
    end
    
    channels = 1:64;
    trials = 1:length(lfpTrials_cleanf.trial);
    A = cellfun(@(x) isnan(x),lfpTrials_cleanf.trial,'UniformOutput',false);
    remove_channels = cellfun(@(x) find(x(:,1)==1),A,'UniformOutput',false);
    B = cellfun(@(x) size(x,1),remove_channels,'UniformOutput',false);
    remove_trials = find(cell2mat(B)==64);
    
    % selecting clean data
    cfg = [];
    cfg.channel = setdiff(channels,remove_channels{1,1});
    cfg.trials  = setdiff(trials,remove_trials);
    clean_data  = ft_selectdata(cfg,lfpTrials_cleanf);
    cd(fullfile(output_path{isess}))
    save('clean_data','clean_data')
    
    % zscore data in order to be able to pool over sessions
    disp(strcat('session- ',num2str(isess),' out of- ',num2str(length(session_names)), ', running LFP z-scoring')) % #ok<UNRCH>
    zlfptrials = fun_zscore_session(clean_data);
    
    % save data
    cd(fullfile(output_path{isess}))
    save('zlfptrials','zlfptrials')
    
    clear lfpTrials_cleanf
    
end

%% Find all unique channels present across sessions

all_channels = {};

for isess = 1:length(session_names)
    if exist(output_path{isess}, 'file')
        cd(output_path{isess}), 
        load('clean_data');
        all_channels = unique([all_channels; clean_data.label]);
    end
end
cd('/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks/critical_time')
save all_channels all_channels

%% Collect Channel-wise LFP - all diff level

collect_trials_all = struct();

for ichan = 1:length(all_channels)
    ichan
    channel_name = all_channels{ichan};
    
    % Ensure valid field name
    valid_channel_name = matlab.lang.makeValidName(channel_name);
    
    % Initialize an empty cell array for the current channel
    if ~isfield(collect_trials_all, valid_channel_name)
        collect_trials_all.(valid_channel_name) = {};
    end
    
    for isess = 1:length(session_names)
        session_dir = output_path{isess};
        
        % Check if the session folder exists
        if ~exist(session_dir, 'dir')
            disp(['Directory not found for session ', session_names{isess}, ', skipping.']);
            continue;
        end
        
        cd(session_dir);
        
        % Load LFP data
        load('zlfptrials.mat');
        
        % Check if this session has the desired channel
        if any(strcmp(zlfptrials.label, channel_name))
            disp(['Processing ', channel_name, ' in session ', session_names{isess}]);
            
            % Select trials for the current channel
            cfg = [];
            cfg.channel = channel_name;
            lfpTrials = ft_selectdata(cfg, zlfptrials);
            
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
                    collect_trials_all.(valid_channel_name){end+1} = itcTrials.trial{tr};
                else
                    warning(['Skipping empty trial ', num2str(tr), ' in session ', session_names{isess}]);
                end
            end
        end
    end
end

cd('/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks/critical_time')
save collect_trials_all collect_trials_all

%% linear regression on erp for critical time calculation - threshold on average of trials - SD of avg

channel_fields = fieldnames(collect_trials_all);
bsltimewindow_start_ind = 1;
bsltimewindow_end_ind = 1000;

n_consec = 10;
th = 2.5;
timeofresponse_ind = nan(length(all_channels), 1);

figure;
for ichan = 1:length(all_channels)
    channel_name = channel_fields{ichan};
    trial_data = collect_trials_all.(channel_name);
    
    if iscell(trial_data)
        trial_data = vertcat(trial_data{:});
    end
    
    % Average ERP for the channel
    mean_erp_data = mean(trial_data, 1);
    
    % Linear regression on baseline using polyfit
    bsltime = bsltimewindow_start_ind:bsltimewindow_end_ind;
    t_baseline = time(bsltime);
    y_baseline = mean_erp_data(bsltime);
    p = polyfit(t_baseline, y_baseline, 1);  % Linear fit
    
    % Extrapolate baseline until 400ms post-stim
    extrap_time_inds = bsltimewindow_start_ind:find(time >= 0.4, 1);  % up to 400ms
    t_extrap = time(extrap_time_inds);
    extrap_baseline = polyval(p, t_extrap);  % Evaluate linear fit
    
    % Adjust ERP
    adj_mean_erp = mean_erp_data(extrap_time_inds) - extrap_baseline;
    
    % Compute SD and threshold from adjusted baseline segment
    bsldata_adj = adj_mean_erp(1:(bsltimewindow_end_ind - bsltimewindow_start_ind + 1));
    bsln_mean = mean(bsldata_adj);
    SD = std(bsldata_adj);
    threshold = bsln_mean + (th * SD);
    
    % Search range (within extrapolated data)
    search_start = bsltimewindow_end_ind + 20 - bsltimewindow_start_ind + 1;
    search_range = search_start : length(adj_mean_erp);
    above_thresh = abs(adj_mean_erp(search_range)) > threshold;
    consec_counts = conv(double(above_thresh), ones(1, n_consec), 'valid');
    ind = find(consec_counts == n_consec, 1);
    
    if ~isempty(ind)
        timeofresponse_ind(ichan, 1) = extrap_time_inds(search_range(ind));
    end
    
    % Plot original ERP with critical time
    subplot(8, 8, ichan)
    %figure;
    plot(time(1:end), adj_mean_erp(1:end))
    hold on
    plot(time(1:end), mean_erp_data(1:end))
    hold on
    if ~isnan(timeofresponse_ind(ichan))
        xline(time(timeofresponse_ind(ichan)), '--k', 'LineWidth', 1.5);
    end
    title(channel_name, 'Interpreter', 'none');
end

% Critical time per channel
indices = timeofresponse_ind;
timeVec = time;
CriticalTime = NaN(1, length(all_channels));
validIdx = indices > 0 & indices <= length(timeVec);
CriticalTime(validIdx) = timeVec(indices(validIdx));
cd('/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks/critical_time')
save CriticalTime CriticalTime %remove 10 sec from this to use as a critical time

% Histogram
valid_critical_times = CriticalTime(~isnan(CriticalTime));
figure;
n_bins = 20;
h = histogram(valid_critical_times, n_bins, 'FaceColor', [0.2, 0.6, 0.8], 'EdgeColor', 'k');
xlabel('Time from target onset(s)');
ylabel('Number of channels');
title('Critical Time');
meanv = mean(valid_critical_times, 'omitnan');
modev = mode(valid_critical_times);
medianv = median(valid_critical_times);
xline(meanv, '--k', 'LineWidth', 1.5);
xline(modev, '--b', 'LineWidth', 1.5);
xline(medianv, '--m', 'LineWidth', 1.5);

%% Visualize - threshold on the average of trials

% Setup
channel_fields = fieldnames(collect_trials_all);
bsltimewindow_start_ind = 1;
bsltimewindow_end_ind = 1000;
n_consec = 10;
th = 2.5;
time_limit = 0.4; % in seconds
stimulus_delay = 20; % in samples after baseline for search start

ichan = 1;
n_channels = length(channel_fields);

% Prepare figure
hFig = figure('KeyPressFcn', @onKeyPress);
setappdata(hFig, 'ichan', ichan);
setappdata(hFig, 'n_channels', n_channels);
setappdata(hFig, 'channel_fields', channel_fields);
setappdata(hFig, 'collect_trials_all', collect_trials_all);
setappdata(hFig, 'bsltimewindow_start_ind', bsltimewindow_start_ind);
setappdata(hFig, 'bsltimewindow_end_ind', bsltimewindow_end_ind);
setappdata(hFig, 'n_consec', n_consec);
setappdata(hFig, 'th', th);
setappdata(hFig, 'time_limit', time_limit);
setappdata(hFig, 'stimulus_delay', stimulus_delay);
setappdata(hFig, 'time', time);

plot_channel();

function plot_channel()
    hFig = gcf;
    ichan = getappdata(hFig, 'ichan');
    channel_fields = getappdata(hFig, 'channel_fields');
    collect_trials_all = getappdata(hFig, 'collect_trials_all');
    bsltimewindow_start_ind = getappdata(hFig, 'bsltimewindow_start_ind');
    bsltimewindow_end_ind = getappdata(hFig, 'bsltimewindow_end_ind');
    n_consec = getappdata(hFig, 'n_consec');
    th = getappdata(hFig, 'th');
    time_limit = getappdata(hFig, 'time_limit');
    stimulus_delay = getappdata(hFig, 'stimulus_delay');
    time = getappdata(hFig, 'time');

    channel_name = channel_fields{ichan};
    trial_data = collect_trials_all.(channel_name);

    if iscell(trial_data)
        trial_data = vertcat(trial_data{:});
    end

    mean_erp_data = mean(trial_data, 1);

    bsltime = bsltimewindow_start_ind:bsltimewindow_end_ind;
    t_baseline = time(bsltime);
    y_baseline = mean_erp_data(bsltime);
    p = polyfit(t_baseline, y_baseline, 1);

    extrap_time_inds = bsltimewindow_start_ind:find(time >= time_limit, 1);
    t_extrap = time(extrap_time_inds);
    extrap_baseline = polyval(p, t_extrap);

    adj_mean_erp = mean_erp_data(extrap_time_inds) - extrap_baseline;

    bsldata_adj = adj_mean_erp(1:(bsltimewindow_end_ind - bsltimewindow_start_ind + 1));
    bsln_mean = mean(bsldata_adj);
    SD = std(bsldata_adj);
    threshold = bsln_mean + (th * SD);

    search_start = bsltimewindow_end_ind + stimulus_delay - bsltimewindow_start_ind + 1;
    search_range = search_start : length(adj_mean_erp);
    above_thresh = abs(adj_mean_erp(search_range)) > threshold;
    consec_counts = conv(double(above_thresh), ones(1, n_consec), 'valid');
    ind = find(consec_counts == n_consec, 1);

    if ~isempty(ind)
        response_time_idx = extrap_time_inds(search_range(ind));
    else
        response_time_idx = NaN;
    end

    clf(hFig);
    plot(time(1:end), adj_mean_erp(1:end), 'LineWidth', 1.5)
    hold on
    plot(time(1:end), mean_erp_data(1:end), 'LineWidth', 1.5)
    if ~isnan(response_time_idx)
        %xline(time(response_time_idx), '--k', 'LineWidth', 1.5);
        xline(time(response_time_idx)-0.01, '--k', 'LineWidth', 1.5);
    end
    title(channel_name, 'Interpreter', 'none');
    xlabel('Time (s)');
    ylabel('ERP');
    grid on;
end

function onKeyPress(~, event)
    hFig = gcf;
    ichan = getappdata(hFig, 'ichan');
    n_channels = getappdata(hFig, 'n_channels');

    if strcmp(event.Key, 'rightarrow')
        ichan = mod(ichan, n_channels) + 1;
    elseif strcmp(event.Key, 'leftarrow')
        ichan = mod(ichan - 2, n_channels) + 1;
    elseif strcmp(event.Key, 'escape')
        close(hFig);
        return;
    end

    setappdata(hFig, 'ichan', ichan);
    plot_channel();
end

%% Visualize to check - threshold per trial
% 
% channel_fields = fieldnames(collect_trials_all);
% bsltimewindow_start_ind = 1;
% bsltimewindow_end_ind = 1000;
% n_consec = 5;
% 
% ichan = 1;  % Start with the first channel
% n_channels = length(channel_fields);  % Total number of channels
% ct_select = 0.072;
% 
% % Create a figure to display the plots
% hFig = figure('KeyPressFcn', @onKeyPress);
% 
% % Store the initial 'ichan' and 'n_channels' in the figure's app data
% setappdata(hFig, 'ichan', ichan);
% setappdata(hFig, 'n_channels', n_channels);
% 
% while ishandle(hFig)
%     % Get the current channel
%     ichan = getappdata(hFig, 'ichan');
%     n_channels = getappdata(hFig, 'n_channels');
%     
%     channel_name = channel_fields{ichan};
%     trial_data = collect_trials_all.(channel_name);
%     
%     if iscell(trial_data)
%         trial_data = vertcat(trial_data{:});
%     end
%     
%     % mean subtraction and threshold
%     bsldata = trial_data(:, bsltimewindow_start_ind:bsltimewindow_end_ind);
%     bsln_mean = mean(bsldata, 2);
%     SD = std(bsldata(:));
%     threshold = bsln_mean + 5 * SD;
%     erp_data = trial_data - bsln_mean;
%     
%     % finding the critical time
%     n_trials = size(erp_data, 1);
%     n_timepoints = size(erp_data, 2);
%     timeofresponse_ind = nan(n_trials, 1);
%     
%     for trial = 1:n_trials
%         search_range = bsltimewindow_end_ind + 5 : 1101; % after stim on and within 100ms after stim on
%         above_thresh = (erp_data(trial, search_range)) > threshold(trial);
%         consec_counts = conv(double(above_thresh), ones(1, n_consec), 'valid');
%         ind = find(consec_counts == n_consec, 1);
%         
%         if ~isempty(ind)
%             timeofresponse_ind(trial) = search_range(ind);  % map back to original time indices
%         end
%     end
%     
%     avg_trial_data = mean(trial_data, 1);
%     mean_timeofresponse = mean(timeofresponse_ind, 'omitnan');
%     
%     % Plot the current channel's ERP data
%     clf(hFig);
%     plot(time, avg_trial_data, 'LineWidth', 2)
%     hold on
%     
%     % Plot the vertical line
%     %     if ~isnan(mean_timeofresponse)
%     %         xline(zlfptrials.time{1, 1}(1, round(mean_timeofresponse)), 'k--', 'LineWidth', 2);
%     %
%     %         text(zlfptrials.time{1, 1}(1, round(mean_timeofresponse)), max(avg_trial_data), sprintf('%.2f ms', zlfptrials.time{1, 1}(1, round(mean_timeofresponse))), ...
%     %             'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
%     %             'FontSize', 10, 'Color', 'black');
%     %     end
%     
%     % Plot with the selected critical time
%     if ~isnan(mean_timeofresponse)
%         xline(ct_select, 'k--', 'LineWidth', 2);
%         
%         text(ct_select, max(avg_trial_data), sprintf('%.3f ms', ct_select), ...
%             'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
%             'FontSize', 10, 'Color', 'black');
%     end
%     
%     % Channel name
%     title(['Channel: ' channel_name], 'Interpreter', 'none');
%     xlabel('Time (ms)')
%     ylabel('ERP')
%     grid on;
%     
%     pause(0.2);
% end
% 
% % Callback function for key press
% function onKeyPress(~, event)
% % Get the current channel index and number of channels from the figure's app data
% ichan = getappdata(gcf, 'ichan');
% n_channels = getappdata(gcf, 'n_channels');
% 
% if strcmp(event.Key, 'rightarrow')  % Right arrow key
%     ichan = mod(ichan, n_channels) + 1;  % Move to next channel
% elseif strcmp(event.Key, 'leftarrow')  % Left arrow key
%     ichan = mod(ichan - 2, n_channels) + 1;  % Move to previous channel
% elseif strcmp(event.Key, 'escape')  % Escape key to exit the loop
%     close(gcf);  % Close the figure and exit the loop
% end
% 
% % Update the channel index in the figure's app data
% setappdata(gcf, 'ichan', ichan);
% end
% 
% 
