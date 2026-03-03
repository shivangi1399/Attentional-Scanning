clear all
close all
clc

%% Setup paths

addpath /opt/fieldtrip_github/
ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/coherence_analysis
addpath /mnt/hpc/projects/MWSampling/4Shivangi
clc

%% Paths and parameters

datafolder    = '/mnt/hpc/projects/MWSampling/4Shivangi/results_hermes';
output_folder = '/mnt/hpc/projects/MWSampling/4Shivangi/results_hermes/coherence/per_location';
animalName    = 'hermes';

% Get session names
temp = dir(datafolder);
session_names = {};
for i = 1:length(temp)
    if contains(temp(i).name, animalName)
        session_names{end+1,1} = temp(i).name;
    end
end
session_paths = cellfun(@(x) fullfile(datafolder, x), session_names, 'UniformOutput', 0);

%% Compute coherence per session and per location

all_coh = struct();
all_freq = [];

for isess = 1:length(session_names)
    fprintf('\n===== Session %d / %d: %s =====\n', isess, length(session_names), session_names{isess});
    cd(session_paths{isess});

    % Load LFP and MUA data
    load('clean_lfp.mat'); load('clean_mua.mat');
    lfpTrials = clean_data; clear clean_data
    muaTrials = clean_mua;  clear clean_mua

    % Add prefixes
    lfpTrials.label = cellfun(@(x) ['lfp_' x], lfpTrials.label, 'UniformOutput', false);
    muaTrials.label = cellfun(@(x) ['mua_' x], muaTrials.label, 'UniformOutput', false);

    % Collect unique locations from trialinfo (column 16)
    locations = unique(lfpTrials.trialinfo(:,16));
    locations(isnan(locations)) = [];

    %% Loop through each location
    for iloc = 1:length(locations)
        loc = locations(iloc);
        fprintf('--- Location %.0f (%d/%d) ---\n', loc, iloc, length(locations));

        % Select trials for this location
        trials_idx = find(lfpTrials.trialinfo(:,16) == loc);
        if isempty(trials_idx)
            fprintf('No trials for location %.0f\n', loc);
            continue;
        end

        cfg = [];
        cfg.trials = trials_idx;
        lfpLoc = ft_selectdata(cfg, lfpTrials);
        muaLoc = ft_selectdata(cfg, muaTrials);

        % Trim trials
        t1 = 1001; % start index
        for t = 1:min(numel(lfpLoc.trial), numel(muaLoc.trial))
            nMin = 1200; % fixed window length
            lfpLoc.trial{t} = lfpLoc.trial{t}(:, t1:nMin);
            muaLoc.trial{t} = muaLoc.trial{t}(:, t1:nMin);
            lfpLoc.time{t}  = lfpLoc.time{t}(t1:nMin);
            muaLoc.time{t}  = muaLoc.time{t}(t1:nMin);
        end

        % Combine LFP + MUA
        cfg = []; cfg.keepsampleinfo = 'no';
        dataLoc = ft_appenddata(cfg, lfpLoc, muaLoc);

        %% Frequency analysis
        cfg_freq = [];
        cfg_freq.output     = 'fourier';
        cfg_freq.method     = 'mtmfft';
        cfg_freq.foilim     = [2 100];
        cfg_freq.tapsmofrq  = 5;
        cfg_freq.keeptrials = 'yes';
        cfg_freq.channel    = {'lfp_*' 'mua_*'};
        freqLoc = ft_freqanalysis(cfg_freq, dataLoc);

        %% Coherence
        cfg_coh = [];
        cfg_coh.method     = 'coh';
        cfg_coh.channelcmb = {'lfp_*' 'mua_*'};
        fdLoc = ft_connectivityanalysis(cfg_coh, freqLoc);

        %% Identify matching pairs (lfp_1 - mua_1, etc.)
        nPairs = min(length(lfpLoc.label), length(muaLoc.label));
        same_pairs_idx = zeros(nPairs,1);
        for p = 1:nPairs
            idx = find(strcmp(fdLoc.labelcmb(:,1), lfpLoc.label{p}) & ...
                       strcmp(fdLoc.labelcmb(:,2), muaLoc.label{p}));
            if ~isempty(idx)
                same_pairs_idx(p) = idx;
            end
        end
        same_pairs_idx = same_pairs_idx(same_pairs_idx>0);

        %% Store results
        all_coh(isess).location(iloc).locID          = loc;
        all_coh(isess).location(iloc).cohspctrm      = fdLoc.cohspctrm;
        all_coh(isess).location(iloc).labelcmb       = fdLoc.labelcmb;
        all_coh(isess).location(iloc).freq           = fdLoc.freq;
        all_coh(isess).location(iloc).same_pairs_idx = same_pairs_idx;
        all_freq = fdLoc.freq;
    end
end

fprintf('\n Location-wise coherence analysis completed.\n');

if ~isfolder(output_folder)
    mkdir(output_folder)
end
cd(output_folder)
save all_coh_locationwise all_coh -v7.3

%% Average across sessions per location and channel

save_folder = '/mnt/hpc/projects/MWSampling/4Shivangi/Plots/coherence/hermes/per_location/';
if ~isfolder(save_folder)
    mkdir(save_folder);
end

nSessions = length(all_coh);
nFreqs = length(all_freq);

% Collect all unique location IDs across all sessions
all_locIDs = [];
for s = 1:nSessions
    all_locIDs = [all_locIDs [all_coh(s).location.locID]]; %#ok<AGROW>
end
all_locIDs = unique(all_locIDs);

% Loop through each unique location ID
for ilocID = all_locIDs
    fprintf('\nPlotting Location ID %d\n', ilocID);

    % Collect all matching pairs across sessions for this location ID
    matching_pairs = {};
    for s = 1:nSessions
        loc_idx = find([all_coh(s).location.locID] == ilocID, 1);
        if isempty(loc_idx)
            continue;
        end
        fd = all_coh(s).location(loc_idx);
        if isempty(fd.same_pairs_idx)
            continue;
        end
        matching_pairs = [matching_pairs; fd.labelcmb(fd.same_pairs_idx,:)]; %#ok<AGROW>
    end

    % Skip if no pairs found for this location
    if isempty(matching_pairs)
        fprintf('No matching pairs found for Location ID %d  skipping.\n', ilocID);
        continue;
    end

    % Unique matching LFP-MUA pairs across sessions
    pair_strings = strcat(matching_pairs(:,1), '_', matching_pairs(:,2));
    [~, ia] = unique(pair_strings, 'stable');
    matching_pairs = matching_pairs(ia,:);
    nPairs = size(matching_pairs,1);

    % Combine coherence across sessions for each pair
    combined_coh = nan(nPairs, nFreqs, nSessions);
    for s = 1:nSessions
        loc_idx = find([all_coh(s).location.locID] == ilocID, 1);
        if isempty(loc_idx)
            continue;
        end
        fd = all_coh(s).location(loc_idx);
        if isempty(fd.same_pairs_idx)
            continue;
        end
        session_pairs = fd.labelcmb(fd.same_pairs_idx,:);
        for p = 1:nPairs
            idx = find(strcmp(session_pairs(:,1), matching_pairs{p,1}) & ...
                       strcmp(session_pairs(:,2), matching_pairs{p,2}));
            if ~isempty(idx)
                combined_coh(p,:,s) = fd.cohspctrm(fd.same_pairs_idx(idx),:);
            end
        end
    end

    % Average coherence across sessions
    mean_coh = nanmean(combined_coh, 3);

    % Plotting
    fig = figure('Color','w','Name',sprintf('Location %d', ilocID),'Position',[100 100 1200 900]);
    nCols = ceil(sqrt(nPairs)); 
    nRows = ceil(nPairs / nCols);
    freq = all_freq;

    for p = 1:nPairs
        subplot(nRows, nCols, p);
        plot(freq, mean_coh(p,:), 'k', 'LineWidth', 1.5);
        ylim([0 1]);
        xlim([2 100]);
        grid on;

        xlabel('Frequency (Hz)', 'FontSize', 8);
        ylabel('Coherence', 'FontSize', 8);

        lfp_name = erase(matching_pairs{p,1}, 'lfp_');
        title(lfp_name, 'Interpreter', 'none', 'FontSize', 8);
    end

    sgtitle(sprintf('Location %d  Mean LFPMUA Coherence', ilocID), ...
        'FontWeight','bold','FontSize',12);

    % Save figure
    save_name = fullfile(save_folder, sprintf('Location_%d_mean_coherence.png', ilocID));
    saveas(fig, save_name);
    close(fig);
end

%% Rank channels by mean coherence

disp('Ranking channels across all locations...');

nSessions = length(all_coh);
nFreqs = length(all_freq);
ranked_channels = struct();

% Collect all unique location IDs across sessions 
all_locIDs = [];
for s = 1:nSessions
    all_locIDs = [all_locIDs [all_coh(s).location.locID]]; %#ok<AGROW>
end
all_locIDs = unique(all_locIDs);

% Loop through each location ID
count = 0;  % struct index counter
for ilocID = all_locIDs
    fprintf('\nRanking Channels for Location %d\n', ilocID);

    % Collect all matching pairs across sessions for this location
    matching_pairs = {};
    for s = 1:nSessions
        loc_idx = find([all_coh(s).location.locID] == ilocID, 1);
        if isempty(loc_idx) || isempty(all_coh(s).location(loc_idx).same_pairs_idx)
            continue;
        end
        fd = all_coh(s).location(loc_idx);
        matching_pairs = [matching_pairs; fd.labelcmb(fd.same_pairs_idx,:)]; %#ok<AGROW>
    end

    % Skip if no pairs found
    if isempty(matching_pairs)
        fprintf('No matching pairs found for Location %d  skipping.\n', ilocID);
        continue;
    end

    % Unique matching pairs across sessions
    pair_strings = strcat(matching_pairs(:,1), '_', matching_pairs(:,2));
    [~, ia] = unique(pair_strings, 'stable');
    matching_pairs = matching_pairs(ia,:);
    nPairs = size(matching_pairs,1);

    % Combine coherence across sessions for each pair
    combined_coh = nan(nPairs, nFreqs, nSessions);
    for s = 1:nSessions
        loc_idx = find([all_coh(s).location.locID] == ilocID, 1);
        if isempty(loc_idx)
            continue;
        end
        fd = all_coh(s).location(loc_idx);
        if isempty(fd.same_pairs_idx)
            continue;
        end
        session_pairs = fd.labelcmb(fd.same_pairs_idx,:);
        for p = 1:nPairs
            idx = find(strcmp(session_pairs(:,1), matching_pairs{p,1}) & ...
                       strcmp(session_pairs(:,2), matching_pairs{p,2}));
            if ~isempty(idx)
                combined_coh(p,:,s) = fd.cohspctrm(fd.same_pairs_idx(idx),:);
            end
        end
    end

    % Compute mean coherence across sessions and frequencies
    mean_coh_freq = squeeze(nanmean(combined_coh, 2));  % [nPairs x nSessions]
    mean_coh_all = nanmean(mean_coh_freq, 2);           % average across sessions

    % Rank pairs
    [sorted_coh, sort_idx] = sort(mean_coh_all, 'descend');

    % Increment counter and store in struct
    count = count + 1;
    ranked_channels(count).locID = ilocID;
    ranked_channels(count).channels = matching_pairs(sort_idx, :);
    ranked_channels(count).mean_coherence = sorted_coh;

    % Display top 5
    fprintf('Top 5 pairs for Location %d:\n', ilocID);
    for k = 1:min(5, nPairs)
        fprintf('  %s-%s : %.3f\n', ranked_channels(count).channels{k,1}, ...
                                     ranked_channels(count).channels{k,2}, ...
                                     ranked_channels(count).mean_coherence(k));
    end
end

cd(output_folder)
save ranked_channels ranked_channels -v7.3

