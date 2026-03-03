%% Hit ERPs for combinations of locations and channels (0-100 ms analysis)

clear all
close all
clc

% Code description:
% -----------------
% Aggregates and averages hit ERPs across sessions, computes ERP metrics (peak, latency, AUC, SNR), 
% ranks channels using a weighted combination of metrics and MUALFP coherence, selects top 
% channels per location, and visualizes results both interactively and as spatial clusters on the 
% 8×8 electrode grid

%% paths and config
datafolder = '/mnt/hpc/projects/MWSampling/4Shivangi/results_hermes';
animalName = 'hermes';

% gather sessions
cd(datafolder)
temp = dir;
session_names = {};
ii = 0;
for i = 1:length(temp)
    if contains(temp(i).name, animalName)
        ii = ii+1;
        session_names{ii,1} = temp(i).name;
    end
end

output_paths = cellfun(@(x) fullfile(datafolder, x, 'ERP_LFP', 'per_loc'), session_names, 'uniform', 0);
output_plot_dir = '/mnt/hpc/projects/MWSampling/4Shivangi/Plots/phase_coherence';
if ~isfolder(output_plot_dir), mkdir(output_plot_dir); end

%% get all unique locations
all_locs = [];
for isess = 1:length(output_paths)
    trialfile = fullfile(output_paths{isess}, 'trial_perm_ind.mat');
    if exist(trialfile,'file')
        load(trialfile,'trial_perm_ind');
        if isfield(trial_perm_ind,'locations')
            all_locs = union(all_locs, trial_perm_ind.locations(:));
        end
    end
end

%% aggregate hit ERPs per channel per location --------------------------------------------------------------------

disp('Collecting hit ERPs for all locations...')
all_hitERPs = struct();
for isess = 1:length(output_paths)
    fprintf('Session %d/%d\n', isess, length(output_paths))
    cd(output_paths{isess})
    if ~exist('trial_perm_ind.mat','file'), continue; end
    load('trial_perm_ind.mat')
    if ~isfield(trial_perm_ind,'loc'), continue; end
    for iloc = 1:length(trial_perm_ind.loc)
        locid = trial_perm_ind.loc(iloc).locationID;
        locERPdir = fullfile(output_paths{isess}, sprintf('loc%d', locid), 'ERP_real');
        hitfile = fullfile(locERPdir, 'norm_hit_timelock.mat');
        if ~exist(hitfile,'file'), continue; end
        data = load(hitfile);
        hit = data.norm_hit_timelock;
        labels = hit.label;
        time = hit.time;
        avg = hit.avg;  % channels x time
        for ichan = 1:length(labels)
            ch = labels{ichan};
            safe_ch = matlab.lang.makeValidName(ch);
            if ~isfield(all_hitERPs,safe_ch), all_hitERPs.(safe_ch) = struct(); end
            loc_field = sprintf('loc%d',locid);
            if ~isfield(all_hitERPs.(safe_ch),loc_field), all_hitERPs.(safe_ch).(loc_field) = []; end
            all_hitERPs.(safe_ch).(loc_field) = [all_hitERPs.(safe_ch).(loc_field); avg(ichan,:)];
        end
    end
end

%% average across sessions

disp('Averaging across sessions...')
chan_names = fieldnames(all_hitERPs);
for ic = 1:length(chan_names)
    ch = chan_names{ic};
    locfields = fieldnames(all_hitERPs.(ch));
    for il = 1:length(locfields)
        mat = all_hitERPs.(ch).(locfields{il});
        if isempty(mat)
            all_hitERPs.(ch).(locfields{il}) = [];
        else
            all_hitERPs.(ch).(locfields{il}) = nanmean(mat,1);
        end
    end
end

%% Interactive ERP Channel Browser (0-100 ms)

disp('Launching interactive ERP channel browser...');

% Number of locations
nLocs = length(all_locs);

% Teal shades for plotting (light to dark)
colors = [linspace(0.5,0,nLocs)', linspace(0.9,0.2,nLocs)', linspace(0.9,0.2,nLocs)'];

% Map safe channel names back to original labels if available
orig_chan_names = {};
safe_chan_names = fieldnames(all_hitERPs);
for i = 1:length(safe_chan_names)
    safe = safe_chan_names{i};
    orig_chan_names{i} = safe;
    for j = 1:length(labels)
        if strcmp(matlab.lang.makeValidName(labels{j}), safe)
            orig_chan_names{i} = labels{j};
            break;
        end
    end
end

% Initialize figure
fig = figure('Name','ERP Channel Browser','Units','normalized','Position',[0.1 0.1 0.6 0.5]);
movegui(fig,'center');
chan_idx = 1; 
nChans = length(safe_chan_names); 
running = true;

while running && isvalid(fig)
    clf(fig,'reset'); 
    safe_ch = safe_chan_names{chan_idx}; 
    ch_display = orig_chan_names{chan_idx};
    hold on; 
    legend_entries = {}; 
    y_all = [];

    for iloc = 1:nLocs
        locid = all_locs(iloc); 
        locfield = sprintf('loc%d', locid);
        if isfield(all_hitERPs.(safe_ch), locfield)
            y = all_hitERPs.(safe_ch).(locfield);
            if isempty(y), continue; end
            tmask = time >= 0 & time <= 0.1; % 0-100 ms
            plot(time(tmask), y(tmask), 'LineWidth', 1.5, 'Color', colors(iloc,:));
            legend_entries{end+1} = sprintf('Loc %d', locid); %#ok<AGROW>
            y_all = [y_all, y(tmask)]; %#ok<AGROW>
        end
    end

    xlabel('Time (s)'); 
    ylabel('Normalized ERP');
    title(sprintf('Hit ERPs - %s (%d/%d)', ch_display, chan_idx, nChans));

    if ~isempty(legend_entries)
        legend(legend_entries, 'Location','northeastoutside');
    end

    grid on; 
    xlim([0 0.1]);
    if ~isempty(y_all)
        padding = 0.1 * (max(y_all)-min(y_all)); 
        ylim([min(y_all)-padding, max(y_all)+padding]);
    end
    xline(0,'--k','LineWidth',1);

    % Instructions
    annotation('textbox',[0.35 0.01 0.4 0.05],'String','Left/Right or P/N: change channel | Esc: exit','EdgeColor','none','HorizontalAlignment','center','FontSize',10);

    % Wait for key press
    waitforbuttonpress; 
    key = get(fig,'CurrentKey');
    switch key
        case {'rightarrow','n'}, chan_idx = min(chan_idx+1, nChans);
        case {'leftarrow','p'}, chan_idx = max(chan_idx-1, 1);
        case 'escape', running = false;
    end
end

if isvalid(fig), close(fig); end
disp('Interactive browsing finished.');

%% ERP metric analysis (0-100 ms) ---------------------------------------------------------------------------------------------------------------------------

analysis_window = [0 0.1];
t_analysis = time >= analysis_window(1) & time <= analysis_window(2);
baseline_mask = time < 0; % pre-stimulus baseline

% Channel grid (8x8)
[X,Y] = meshgrid(1:8,1:8); chan_coords = [X(:) Y(:)];

safe_chan_names = fieldnames(all_hitERPs);
nChannels = length(safe_chan_names); 
nLocs = length(all_locs);

display_names = safe_chan_names;
if exist('labels','var')
    for i=1:nChannels
        for j=1:length(labels)
            if strcmp(matlab.lang.makeValidName(labels{j}),safe_chan_names{i})
                display_names{i}=labels{j}; break;
            end
        end
    end
end

PeakAmp=nan(nChannels,nLocs); 
PeakLat=nan(nChannels,nLocs);
AUC=nan(nChannels,nLocs); 
SNR=nan(nChannels,nLocs);
HasData=false(nChannels,nLocs);

distMat = squareform(pdist(chan_coords));

for ic=1:nChannels
    safe_ch = safe_chan_names{ic};
    for iloc=1:nLocs
        locfield = sprintf('loc%d',all_locs(iloc));
        if isfield(all_hitERPs.(safe_ch),locfield)
            y = all_hitERPs.(safe_ch).(locfield);
            if isempty(y) || all(isnan(y)), continue; end
            HasData(ic,iloc)=true;

            % Peak amplitude & latency (time at which abs amp is  the largest)
            y_peak = y(t_analysis);
            [~,idx_rel] = max(abs(y_peak));
            idx_all = find(t_analysis); peak_idx = idx_all(idx_rel);
            PeakAmp(ic,iloc) = y(peak_idx);
            PeakLat(ic,iloc) = time(peak_idx);

            % AUC
            y_auc = y(t_analysis);
            AUC(ic,iloc) = trapz(time(t_analysis), abs(y_auc));

            % SNR
            basevals = y(baseline_mask);
            baseSD = nanstd(basevals(:));
            if baseSD ~=0, SNR(ic,iloc)=PeakAmp(ic,iloc)/baseSD; end
        end
    end
end

cd(fullfile(datafolder, 'coherence', 'per_location'));
load('ranked_channels.mat')
disp('ERP metrics computed.')

%% Select channels based on amplitude, latency and SNR criteria per location

% User-defined weights
w_amp = 0.7;     % weight for amplitude
w_lat = 0;     % weight for latency
w_snr = 0.1;     % weight for SNR
w_coh = 0.2;     % weight for coherence ranking

pct_threshold = 60;  % top X% of channels to select
selected_ch_per_loc = cell(1, nLocs);

for iloc = 1:nLocs
    ch_candidates = find(HasData(:, iloc));
    if isempty(ch_candidates), continue; end

    amps = abs(PeakAmp(ch_candidates, iloc));
    lats = PeakLat(ch_candidates, iloc);
    snr_vals = SNR(ch_candidates, iloc);

    % Coherence ranking info
    mean_coh = ranked_channels(iloc).mean_coherence;

    % Align coherence values with the same channel indices )
    coh_vals = mean_coh(ch_candidates);

    % Normalize all metrics to [0,1]
    norm_amp = (amps - min(amps)) / (max(amps) - min(amps) + eps);
    norm_lat = (max(lats) - lats) / (max(lats) - min(lats) + eps); % faster latency = higher score
    norm_snr = (snr_vals - min(snr_vals)) / (max(snr_vals) - min(snr_vals) + eps);
    norm_coh = (coh_vals - min(coh_vals)) / (max(coh_vals) - min(coh_vals) + eps);

    % Weighted score including coherence
    score = w_amp*norm_amp + w_lat*norm_lat + w_snr*norm_snr + w_coh*norm_coh;

    % Select top pct_threshold channels based on weighted score
    thresh = prctile(score, 100 - pct_threshold);
    sel_idx = ch_candidates(score >= thresh);

    % Optional: sort by score
    [~, sort_idx] = sort(score(sel_idx), 'descend');
    sel_idx = sel_idx(sort_idx);

    selected_ch_per_loc{iloc} = sel_idx;
end

%% Display selected channels

for iloc = 1:nLocs
    sel_idx = selected_ch_per_loc{iloc};
    if isempty(sel_idx)
        fprintf('Location %d: no channels meet criteria\n', all_locs(iloc));
    else
        fprintf('Location %d: ', all_locs(iloc));
        fprintf('%s ', display_names{sel_idx});
        fprintf('\n');
    end
end

%% Interactive cluster inspection per location

dist_thresh = 2;  % distance threshold for clustering (in channel grid units)

% Custom cluster colors
custom_colors = [0 0.5 0.5;      % dark teal 
                 0.5 0 0.5;      % dark magenta 
                 0 0 0.75;       % blue
                 0.5 0.5 0.5];   % gray
custom_colors = repmat(custom_colors, ceil(20/size(custom_colors,1)), 1);

% Main loop
running = true;
while running
    %% Ask user for location
    fprintf('Available locations: %s\n', strjoin(arrayfun(@num2str, all_locs, 'UniformOutput', false), ', '));
    loc_choice = input('Enter location number to inspect (or 0 to exit): ');

    if loc_choice == 0
        running = false;
        break;
    end

    %% Validate input
    loc_idx = find(all_locs == loc_choice, 1);
    if isempty(loc_idx)
        fprintf('Invalid location number. Try again.\n');
        continue;
    end

    sel_idx = selected_ch_per_loc{loc_idx};

    %% Keep only channels that exist in all_hitERPs and have data for this location
    valid_sel_idx = sel_idx(arrayfun(@(k) isfield(all_hitERPs.(safe_chan_names{k}), sprintf('loc%d', loc_choice)), sel_idx));

    if isempty(valid_sel_idx)
        fprintf('No selected channels with data for location %d.\n', loc_choice);
        continue;
    end

    %% Build adjacency matrix and cluster thresholded channels
    adjMat = distMat(valid_sel_idx, valid_sel_idx) < dist_thresh;
    G = graph(adjMat);
    bins = conncomp(G);
    nClusters = max(bins);

    %% ERP PLOT
    figure('Name', sprintf('Location %d - ERPs', loc_choice), ...
           'Units', 'normalized', 'Position', [0.1 0.1 0.7 0.7]);
    hold on;

    remaining_idx = setdiff(1:nChannels, valid_sel_idx);

    % Plot gray ERPs for non-thresholded channels
    for k = 1:length(remaining_idx)
        ch = safe_chan_names{remaining_idx(k)};
        loc_field = sprintf('loc%d', loc_choice);
        if isfield(all_hitERPs, ch) && isfield(all_hitERPs.(ch), loc_field)
            y = all_hitERPs.(ch).(loc_field);
            tmask = time >= 0 & time <= 0.1;
            plot(time(tmask), y(tmask), 'Color', [0.8 0.8 0.8], 'LineWidth', 1, 'HandleVisibility', 'off');
        end
    end

    % Plot colored ERPs for thresholded clusters
    legend_entries = {};
    for icl = 1:nClusters
        cluster_idx = valid_sel_idx(bins == icl);
        color = custom_colors(icl,:);
        for k = 1:length(cluster_idx)
            ch = safe_chan_names{cluster_idx(k)};
            loc_field = sprintf('loc%d', loc_choice);
            if isfield(all_hitERPs.(ch), loc_field)
                y = all_hitERPs.(ch).(loc_field);
                tmask = time >= 0 & time <= 0.1;
                plot(time(tmask), y(tmask), 'LineWidth', 1.5, 'Color', color);
                legend_entries{end+1} = ch; %#ok<AGROW>
            end
        end
    end

    xlabel('Time (s)');
    ylabel('Normalized ERP');
    title(sprintf('Location %d - ERPs', loc_choice));
    if ~isempty(legend_entries)
        legend(legend_entries, 'Location', 'northeastoutside', 'Interpreter', 'none');
    end
    xline(0, '--k');
    grid on;
    hold off;

    %% CLUSTER GRID PLOT 
    figure('Name', sprintf('Location %d - Cluster Grid', loc_choice), ...
           'Units', 'normalized', 'Position', [0.2 0.2 0.6 0.6]);
    hold on;

    % Helper: extract numeric channel ID from names like 'V4_37'
    get_ch_number = @(name) sscanf(name, '%*[^_]_%d');

    % All channel numbers
    all_ch_numbers = arrayfun(@(name) get_ch_number(safe_chan_names{name}), 1:nChannels);

    % Grid layout: 8 rows x 8 columns
    nRows = 8;
    nCols = 8;

    % Column-wise bottom-to-top
    col_idx = ceil(all_ch_numbers / nRows);           % column number
    row_idx = nRows - mod(all_ch_numbers-1, nRows);  % row number

    % Shade specific columns
    light_pink = [1 0.9 0.95]; % RGB light pink
    highlight_cols = [1 4 5 8]; % columns to shade
    for c = highlight_cols
        rectangle('Position',[c-0.5,0.5,1,nRows],'FaceColor',light_pink,'EdgeColor','none');
    end

    % Plot all channels
    for k = 1:nChannels
        if ismember(k, valid_sel_idx)
            % Thresholded: color by cluster
            cluster_idx = find(valid_sel_idx == k,1);
            cluster_id = bins(cluster_idx);
            color = custom_colors(cluster_id,:);
        else
            % Non-thresholded: gray
            color = [0.8 0.8 0.8];
        end
        scatter(col_idx(k), row_idx(k), 300, 'o', 'MarkerFaceColor', color, ...
                'MarkerEdgeColor','k','LineWidth',1.5);
        text(col_idx(k)+0.2, row_idx(k), num2str(all_ch_numbers(k)), ...
             'FontSize', 10, 'FontWeight', 'bold');
    end

    xlabel('Column');
    ylabel('Row');
    title(sprintf('Location %d - Cluster Grid ', loc_choice));
    axis([0 nCols+1 0 nRows+1]);
    axis equal;
    grid on;
    hold off;
end

