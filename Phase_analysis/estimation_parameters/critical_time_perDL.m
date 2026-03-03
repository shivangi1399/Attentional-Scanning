clear all
close all
clc

%% paths

addpath /opt/fieldtrip_github/
ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
clc

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Collect Channel-wise LFP - all diff level

collect_trials_all = struct();
dfi = 1;

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
        cd(output_paths{isess});
        
        % Load LFP data
        load('zlfpTrials.mat');
        
        % Check if this session has the desired channel
        if any(strcmp(zlfpTrials.label, channel_name))
            disp(['Processing ', channel_name, ' in session ', session_names{isess}]);
            
            % Get unique difficulty levels for this session
            unique_DL = flip(unique(zlfpTrials.trialinfo(:,18)));
            num_DL = length(unique_DL);
            %first_quarter_DL = unique_DL(1:ceil(num_DL / 4));
            DL = unique_DL(dfi);
            
            % Select trials for the current channel
            cfg = [];
            cfg.channel = channel_name;
            %cfg.trials = find(ismember(zlfpTrials.trialinfo(:,18), DL));
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
                    collect_trials_all.(valid_channel_name){end+1} = itcTrials.trial{tr};
                else
                    warning(['Skipping empty trial ', num2str(tr), ' in session ', session_names{isess}]);
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Critical time per trial and channel

channel_fields = fieldnames(collect_trials_all);
bsltimewindow_start_ind = 1;
bsltimewindow_end_ind = 1000;

n_consec = 5;
n_channels = length(channel_fields);
timeofresponse_all = struct();

% Step 1: Compute and store timeofresponse_ind for each channel
for ichan = 1:n_channels
    channel_name = channel_fields{ichan};
    trial_data = collect_trials_all.(channel_name);
    
    if iscell(trial_data)
        trial_data = vertcat(trial_data{:});
    end
    
    bsldata = trial_data(:, bsltimewindow_start_ind:bsltimewindow_end_ind);
    bsln_mean = mean(bsldata,2);
    SD = std(bsldata(:));
    threshold = bsln_mean + 5 * SD;
    erp_data = trial_data - bsln_mean;
    
    n_trials = size(erp_data, 1);
    timeofresponse_ind = nan(n_trials, 1);
    
    for trial = 1:n_trials
        search_range = bsltimewindow_end_ind+5 : 1101;
        above_thresh = abs(erp_data(trial, search_range)) > threshold(trial);
        consec_counts = conv(double(above_thresh), ones(1, n_consec), 'valid');
        ind = find(consec_counts == n_consec, 1);
        
        if ~isempty(ind)
            timeofresponse_ind(trial) = search_range(ind);
        end
    end
    
    timeofresponse_all.(channel_name) = timeofresponse_ind;
    
    time_vals = nan(size(timeofresponse_ind));
    valid_idx = ~isnan(timeofresponse_ind);
    timeofresponse_all.(channel_name) = time(round(timeofresponse_ind(valid_idx)));
end

% Step 2: Plot histograms

channel_fields = fieldnames(timeofresponse_all);
n_channels = length(channel_fields);

% global min and max values for x and y axes
x_min = inf;
x_max = -inf;
y_max = -inf;

for ichan = 1:n_channels
    channel_name = channel_fields{ichan};
    timeofresponse_ind = timeofresponse_all.(channel_name);
    
    % Find the minimum and maximum x values (time)
    x_min = min(x_min, min(timeofresponse_ind(~isnan(timeofresponse_ind))));
    x_max = max(x_max, max(timeofresponse_ind(~isnan(timeofresponse_ind))));
    
    % Calculate histogram counts to find the maximum y value
    [counts, edges] = histcounts(timeofresponse_ind, 'BinEdges', linspace(x_min, x_max, 20)); % Create 20 bins
    y_max = max(y_max, max(counts));
end

figure;
for ichan = 1:n_channels
    channel_name = channel_fields{ichan};
    timeofresponse_ind = timeofresponse_all.(channel_name);
    
    % Plot current channel
    subplot(8, 8, ichan);
    histogram(timeofresponse_ind, 'BinEdges', linspace(x_min, x_max, 20), 'EdgeColor', 'k', 'FaceColor', 'b', 'FaceAlpha', 0.5);
    hold on;
    
    % mean line
    mean_time = nanmean(timeofresponse_ind);
    y_limits = [0,y_max];
    plot([mean_time, mean_time], y_limits, 'm--', 'LineWidth', 1);
    
    % axis limits
    xlim([x_min, x_max]);
    ylim([0, y_max]);
    hold off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Critical time for each channel - threshold per trial

channel_fields = fieldnames(collect_trials_all);
bsltimewindow_start_ind = 1;
bsltimewindow_end_ind = 1000;

n_consec = 5;

for ichan = 1:length(all_channels)
    channel_name = channel_fields{ichan};
    trial_data = collect_trials_all.(channel_name);
    
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
mean_val = mean(valid_critical_times, 'omitnan');
xline(mean_val, '--k', 'LineWidth', 1.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Critical time for each channel - threshold on average of trials - SD after mean subtraction

channel_fields = fieldnames(collect_trials_all);
bsltimewindow_start_ind = 800;
bsltimewindow_end_ind = 1000;

n_consec = 3;
th = 0.5;
timeofresponse_ind = nan(length(all_channels), 1);

for ichan = 1:length(all_channels)
    channel_name = channel_fields{ichan};
    trial_data = collect_trials_all.(channel_name);
    
    if iscell(trial_data)
        trial_data = vertcat(trial_data{:});
    end
    
    % mean subtraction and threshold
    bsldata = trial_data(:, bsltimewindow_start_ind:bsltimewindow_end_ind);
    bsln_mean = mean(bsldata(:)); 
    erp_data = (trial_data - bsln_mean);
    sd_data = erp_data(:, bsltimewindow_start_ind:bsltimewindow_end_ind);
    SD = std(sd_data(:));
    threshold = bsln_mean + (th * SD);
    
    mean_erp_data = mean((trial_data - bsln_mean),1);
    
    % finding critical time
    search_range = bsltimewindow_end_ind+20 : 1101; % after stim on and within 100ms after stim on
    above_thresh = abs(mean_erp_data(1, search_range)) > threshold;
    consec_counts = conv(double(above_thresh), ones(1, n_consec), 'valid');
    ind = find(consec_counts == n_consec, 1);
    
    if ~isempty(ind)
        timeofresponse_ind(ichan,1) = search_range(ind);  % map back to original time indices
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

%% Critical time for each channel - threshold on average of trials - SD of mean of trials

channel_fields = fieldnames(collect_trials_all);
bsltimewindow_start_ind = 600;
bsltimewindow_end_ind = 1000;

n_consec = 5;  
th = 2; 
timeofresponse_ind = nan(length(all_channels), 1);

figure;
for ichan = 1:length(all_channels)
    channel_name = channel_fields{ichan};
    trial_data = collect_trials_all.(channel_name);
    
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
    %upper_threshold = bsln_mean + (th * SD);
    %lower_threshold = bsln_mean - (th * SD);
    
    % finding critical time
    search_range = bsltimewindow_end_ind+20 : 1101; % after stim on and within 100ms after stim on
    above_thresh = abs(mean_erp_data(1, search_range)) > threshold;
    %above_thresh = (mean_erp_data(1, search_range) > upper_threshold) | (mean_erp_data(1, search_range) < lower_threshold);
    consec_counts = conv(double(above_thresh), ones(1, n_consec), 'valid');
    ind = find(consec_counts == n_consec, 1);
    
    if ~isempty(ind)
        timeofresponse_ind(ichan,1) = search_range(ind);  % map back to original time indices
    end
    
    subplot(8,8,ichan)
    plot(time(1,500:end),mean_erp_data(1,500:end))
    hold on
    %yline(threshold, '--k', 'LineWidth', 1.5);
    %yline(upper_threshold, '--k', 'LineWidth', 1.5);
    %yline(lower_threshold, '--k', 'LineWidth', 1.5);
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

%% parameter search for critial time of each channel - threshold on average of trials - SD of mean of trials

channel_fields = fieldnames(collect_trials_all);
bsltimewindow_end_ind = 1000;
timeVec = time;

best_params = struct('start', 0, 'n_consec', 0, 'th', 0, 'early_hits', Inf);

for bsltimewindow_start_ind = 600:50:950
    for n_consec = 5:20
        for th = 2.0:0.25:3.5

            timeofresponse_ind = nan(length(all_channels), 1);

            for ichan = 1:length(all_channels)
                channel_name = channel_fields{ichan};
                trial_data = collect_trials_all.(channel_name);

                if iscell(trial_data)
                    trial_data = vertcat(trial_data{:});
                end

                bsldata = trial_data(:, bsltimewindow_start_ind:bsltimewindow_end_ind);
                bsln_mean = mean(bsldata(:));

                mean_erp_data = mean((trial_data - bsln_mean),1);
                sd_data = mean_erp_data(:, bsltimewindow_start_ind:bsltimewindow_end_ind);
                SD = std(sd_data(:));

                threshold = bsln_mean + (th * SD);

                search_range = bsltimewindow_end_ind+20 : 1101;
                above_thresh = abs(mean_erp_data(1, search_range)) > threshold;
                consec_counts = conv(double(above_thresh), ones(1, n_consec), 'valid');
                ind = find(consec_counts == n_consec, 1);

                if ~isempty(ind)
                    timeofresponse_ind(ichan,1) = search_range(ind);
                end
            end

            % evaluate results
            indices = timeofresponse_ind;
            CriticalTime = NaN(1, length(all_channels));
            validIdx = indices > 0 & indices <= length(timeVec);
            CriticalTime(validIdx) = timeVec(indices(validIdx));
            n_early_hits = sum(CriticalTime == 0.02);

            if n_early_hits < best_params.early_hits
                best_params.start = bsltimewindow_start_ind;
                best_params.n_consec = n_consec;
                best_params.th = th;
                best_params.early_hits = n_early_hits;
            end
        end
    end
end

fprintf('Best Parameters:\nBaseline Start: %d\nConsec: %d\nThreshold: %.2f\n20ms Hits: %d\n', ...
    best_params.start, best_params.n_consec, best_params.th, best_params.early_hits);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% linear regression on erp for critical time calculation - threshold on average of trials - SD of avg

channel_fields = fieldnames(collect_trials_all);
bsltimewindow_start_ind = 1;
bsltimewindow_end_ind = 1000;

n_consec = 10;
th = 2.5;
timeofresponse_ind = nan(length(all_channels), 1);

%figure;
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
    %subplot(8, 8, ichan)
    figure;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Collect Channel-wise LFP - separate diff level

cd('/mnt/hpc/projects/MWSampling/4Shivangi/results/critical_time')
load('all_channels.mat')

collect_erp_perDL = struct();

% Step 1: Gather all unique difficulty levels across sessions
all_difficulty_levels = [];

for isess = 1:length(session_names)
    cd(output_paths{isess});
    load('zlfpTrials.mat');
    
    all_difficulty_levels = [all_difficulty_levels; unique(zlfpTrials.trialinfo(:,18))];
end

% Unique and sorted difficulty levels
all_difficulty_levels = unique(all_difficulty_levels);

% Step 2: Collect ERP trials per channel and difficulty level
for ichan = 1:length(all_channels)
    channel_name = all_channels{ichan};
    valid_channel_name = matlab.lang.makeValidName(channel_name);
    
    % substructure for the channel
    if ~isfield(collect_erp_perDL, valid_channel_name)
        collect_erp_perDL.(valid_channel_name) = struct();
    end

    % cell arrays for each difficulty level
    for iDL = 1:length(all_difficulty_levels)
        DL = all_difficulty_levels(iDL);
        collect_erp_perDL.(valid_channel_name).(['DL' num2str(DL)]) = {};
    end

    for isess = 1:length(session_names)
        cd(output_paths{isess});
        load('zlfpTrials.mat');

        if any(strcmp(zlfpTrials.label, channel_name))
            disp(['Processing ', channel_name, ' in session ', session_names{isess}]);

            for iDL = 1:length(all_difficulty_levels)
                DL = all_difficulty_levels(iDL);

                % Select trials with this difficulty level
                trial_idx_DL = find(zlfpTrials.trialinfo(:,18) == DL);
                if isempty(trial_idx_DL)
                    continue;
                end

                cfg = [];
                cfg.channel = channel_name;
                cfg.trials = trial_idx_DL;
                lfpTrials = ft_selectdata(cfg, zlfpTrials);

                % Select only hit trials
                cfg = [];
                cfg.trials = find(lfpTrials.trialinfo(:,20) == 1);
                itcTrials = ft_selectdata(cfg, lfpTrials);

                if isempty(itcTrials.trial)
                    disp(['No valid hit trials for ', channel_name, ' at DL=', num2str(DL), ' in ', session_names{isess}]);
                    continue;
                end

                % Collect ERP trials for the current difficulty level
                for tr = 1:length(itcTrials.trial)
                    if ~isempty(itcTrials.trial{tr})
                        collect_erp_perDL.(valid_channel_name).(['DL' num2str(DL)]){end+1} = itcTrials.trial{tr};
                    else
                        warning(['Empty trial ', num2str(tr), ' at DL=', num2str(DL)]);
                    end
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Average for each channel - threshold per trial

% Parameters
bsltimewindow_start_ind = 1;
bsltimewindow_end_ind = 1000;
n_consec = 5;
n_bins = 20;

% Structure to store critical time per DL
CriticalTime_perDL = struct();
all_times_combined = [];

% Loop through all difficulty levels
for iDL = 1:length(all_difficulty_levels)
    DL = all_difficulty_levels(iDL);
    DL_field = ['DL' num2str(DL)];
    CriticalTime_perDL.(DL_field) = [];

    % Loop through all channels
    for ichan = 1:length(all_channels)
        channel_name = all_channels{ichan};
        valid_channel_name = matlab.lang.makeValidName(channel_name);

        % Check if data exists for this channel and DL
        if ~isfield(collect_erp_perDL, valid_channel_name) || ...
           ~isfield(collect_erp_perDL.(valid_channel_name), DL_field)
            continue;
        end

        trial_data = collect_erp_perDL.(valid_channel_name).(DL_field);

        if isempty(trial_data)
            continue;
        end

        % Convert to matrix if possible
        if iscell(trial_data)
            try
                trial_data = vertcat(trial_data{:});
            catch
                warning(['Skipping ', channel_name, ' DL=', num2str(DL), ' due to inconsistent trial sizes.']);
                continue;
            end
        end

        % Baseline and threshold
        bsldata = trial_data(:, bsltimewindow_start_ind:bsltimewindow_end_ind);
        bsln_mean = mean(bsldata, 2);
        SD = std(bsldata(:));
        threshold = bsln_mean + 5 * SD;
        erp_data = trial_data - bsln_mean;

        % Critical time detection
        n_trials = size(erp_data, 1);
        timeofresponse_ind = nan(n_trials, 1);

        for trial = 1:n_trials
            search_range = bsltimewindow_end_ind + 5 : min(size(erp_data, 2), 1101);
            above_thresh = abs(erp_data(trial, search_range)) > threshold(trial);
            consec_counts = conv(double(above_thresh), ones(1, n_consec), 'valid');
            ind = find(consec_counts == n_consec, 1);
            if ~isempty(ind)
                timeofresponse_ind(trial) = search_range(ind);
            end
        end

        % Convert to time, if valid
        mean_idx = round(mean(timeofresponse_ind, 'omitnan'));
        if ~isnan(mean_idx) && mean_idx > 0 && mean_idx <= length(time)
            time_val = time(mean_idx);
            CriticalTime_perDL.(DL_field)(end+1) = time_val;
            all_times_combined(end+1) = time_val;
        end
    end
end

% Determine global x-axis limits
xmin = min(all_times_combined);
xmax = max(all_times_combined);

% Plot histograms per DL with consistent x-axis
DL_fields = fieldnames(CriticalTime_perDL);
figure;
tiledlayout('flow');

for i = 1:length(DL_fields)
    thisDL = DL_fields{i};
    times = CriticalTime_perDL.(thisDL);

    if isempty(times)
        continue;
    end

    nexttile;
    histogram(times, n_bins, ...
        'FaceColor', [0.2, 0.6, 0.8], ...
        'EdgeColor', 'k');
    xlim([xmin, xmax]);

    % Add mean line
    mean_val = mean(times, 'omitnan');
    xline(mean_val, '--k', 'LineWidth', 1.2);

    title(['DL = ', extractAfter(thisDL, 'DL')]);
    xlabel('Time from target onset (s)');
    ylabel('Number of channels');
end

sgtitle('Critical Time per Difficulty Level');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Critical time for sections of difficulty levels
%% for first quarter of DL
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
    %yline(threshold, '--k', 'LineWidth', 1.5);
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

%% for last quarter of DL
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
            sorted_DL = sort(unique_DL);             
            num_DL = length(sorted_DL);             
            lowest_quartile_DL = sorted_DL(1:ceil(num_DL / 4)); 
            
            % Select trials for the current channel
            cfg = [];
            cfg.channel = channel_name;
            cfg.trials = find(ismember(zlfpTrials.trialinfo(:,18), lowest_quartile_DL));
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

%% Critical time for each channel - threshold on average of trials - SD of mean of trials

channel_fields = fieldnames(collect_trials);
bsltimewindow_start_ind = 800;
bsltimewindow_end_ind = 1000;

n_consec = 15;  
th = 4; 
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
    %yline(threshold, '--k', 'LineWidth', 1.5);
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualize to check - threshold per trial

channel_fields = fieldnames(collect_trials_all);
bsltimewindow_start_ind = 1;
bsltimewindow_end_ind = 1000;
n_consec = 5;

ichan = 1;  % Start with the first channel
n_channels = length(channel_fields);  % Total number of channels
ct_select = 0.068; %0.056;

% Create a figure to display the plots
hFig = figure('KeyPressFcn', @onKeyPress);

% Store the initial 'ichan' and 'n_channels' in the figure's app data
setappdata(hFig, 'ichan', ichan);
setappdata(hFig, 'n_channels', n_channels);

while ishandle(hFig)
    % Get the current channel
    ichan = getappdata(hFig, 'ichan');
    n_channels = getappdata(hFig, 'n_channels');  
    
    channel_name = channel_fields{ichan};
    trial_data = collect_trials_all.(channel_name);
    
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
    plot(time, avg_trial_data, 'LineWidth', 2)
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
        
        text(ct_select, max(avg_trial_data), sprintf('%.3f ms', ct_select), ...
            'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
            'FontSize', 10, 'Color', 'black');
    end
    
    % Channel name
    title(['Channel: ' channel_name], 'Interpreter', 'none');
    xlabel('Time (ms)')
    ylabel('ERP')
    grid on;
    
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

