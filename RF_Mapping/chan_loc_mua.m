clearvars
close all
clc

% Code description:
% -----------------
% Organizes V4 receptive field (RF) mapping data for two monkeys
% Aligns the RFs with behavioral target locations, and saves them.
% It then visualizes RF centers using the max activity method - but its not reliable
% Then we use Gaussian-fitted RFs, and overlap the center and radius on the location map

%% load RF mapping data

datafolder = '/mnt/hpc/projects/MWSampling/RF_mappings';
cd(datafolder);

items = dir(datafolder);
isSession = [items.isdir] & ~startsWith({items.name}, '.');
sessionList = {items(isSession).name};

RF_data = struct();

for s = 1:numel(sessionList)
    sessionName = sessionList{s};
    sessionPath = fullfile(datafolder, sessionName);
    
    % Find .RF file(s)
    rfFiles = dir(fullfile(sessionPath, '*.RF'));
    
    if isempty(rfFiles)
        fprintf('No .RF file found in %s\n', sessionName);
        continue;
    end
    
    % Load the first .RF file
    RFfile = fullfile(sessionPath, rfFiles(1).name);
    safeName = matlab.lang.makeValidName(sessionName);
    
    try
        RF_data.(safeName) = load(RFfile, '-mat');
    catch ME
        fprintf('Error loading %s: %s\n', RFfile, ME.message);
        continue;
    end
end

%% Map the channels according to the scanning data

sessionNames = fieldnames(RF_data);
RF_data_V4 = struct();  % initialize new struct

for isess = 1:numel(sessionNames)
    D = RF_data.(sessionNames{isess});
    
    if ~isfield(D, 'sessInfo') || ~isfield(D.sessInfo, 'monkey')
        error('Session %s does not have D.sessInfo.monkey field.', sessionNames{isess});
    end
    
    monkey = D.sessInfo.monkey;
    
    % ----- Identify V4 labels and indices -----
    isV4 = startsWith(D.label, 'V4-');
    V4_labels = D.label(isV4);
    V4_indices = find(isV4);
    
    % ----- Convert label numbers -----
    parts = erase(V4_labels, 'V4-');
    V4_nums = str2double(parts);
    
    switch lower(monkey)
        case 'hermes'
            V4_nums(isnan(V4_nums)) = 64;
            V4_labels(isnan(V4_nums)) = {'V4-64'};
        case 'klecks'
            V4_nums(isnan(V4_nums)) = 127;
            V4_labels(isnan(V4_nums)) = {'V4-127'};
        otherwise
            error('Unknown monkey name: %s', monkey);
    end
    
    [sorted_nums, sort_order] = sort(V4_nums);
    V4_sorted_indices = V4_indices(sort_order);
    
    % Increment numbers for Klecks
    if strcmpi(monkey, 'klecks')
        incremented_nums = sorted_nums + 1;
        incremented_nums(incremented_nums == 128) = 128;   % preserve special label
    else
        incremented_nums = sorted_nums;
    end
    
    V4_sorted_labels = cellstr(strcat('V4-', string(incremented_nums)));
    
    % ----- Reorder D.RF according to V4 -----
    if isfield(D, 'RF')
        RF_V4_sorted = D.RF(:, V4_indices);          % select V4 columns
        RF_V4_sorted = RF_V4_sorted(:, sort_order);  % reorder according to sorted labels
    else
        RF_V4_sorted = [];  % if no RF field exists
    end
    
    % ----- Update -----
    RF_data_V4.(sessionNames{isess}) = struct( ...
        'sessInfo', D.sessInfo, ...
        'V4_labels', table(V4_sorted_labels), ...
        'RF_V4_sorted', RF_V4_sorted ...
        );
end

% Save 

RF_hermes = struct();
RF_klecks = struct();
fields = fieldnames(RF_data_V4);

for i = 1:numel(fields)
    fname = fields{i};
    
    if startsWith(fname,'hermes')
        RF_hermes.(fname) = RF_data_V4.(fname);
    elseif startsWith(fname,'klecks')
        RF_klecks.(fname) = RF_data_V4.(fname);
    end
end

outdir_hermes = '/mnt/hpc/projects/MWSampling/4Shivangi/results_hermes/RF_mapping';
outdir_klecks = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks/RF_mapping';

if ~exist(outdir_hermes, 'dir'); mkdir(outdir_hermes); end
if ~exist(outdir_klecks, 'dir'); mkdir(outdir_klecks); end

% save
save(fullfile(outdir_hermes, 'RF_hermes.mat'), 'RF_hermes');
save(fullfile(outdir_klecks, 'RF_klecks.mat'), 'RF_klecks');

%% Mapping target on RF %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot max activity RF centers and target locations in RF coordiantes - example session

% Get session and data
fn = fieldnames(RF_data_V4);
session_name = fn{13}; %select one session
S = RF_data_V4.(session_name);

% Target locations (already RF-centered pixels)
x_target_pix = clean_data.trialinfo(:,16);
y_target_pix = clean_data.trialinfo(:,17);

% Get all channels
nChannels = length(S.RF_V4_sorted);

% Pixels per degree (make sure this is defined per session)
pixPerDeg = S.sessInfo.ppd;

% Prepare figure
figure('Color','w');
scatter(x_target_pix, y_target_pix, 50, [0.8 0.8 0.8], 'filled'); % targets in light gray
hold on;

colors = lines(nChannels);  % unique colors for each channel

for ichan = 1:nChannels
    % RF map
    RF_map_struct = S.RF_V4_sorted(ichan);
    RF_map_matrix = RF_map_struct.map;
    [nY, nX] = size(RF_map_matrix);
    
    % Pixel axes centered at RF map center
    x_pix_axis = (1:nX) - (nX+1)/2;      % origin at map center
    y_pix_axis = (nY+1)/2 - (1:nY);      % +Y is up
    
    % Convert RF center from degrees ? pixels
    xRF_center_pix = RF_map_struct.xMax * pixPerDeg;
    yRF_center_pix = RF_map_struct.yMax * pixPerDeg;
    
    % Find closest indices in the RF map
    [~, x_idx] = min(abs(x_pix_axis - xRF_center_pix));
    [~, y_idx] = min(abs(y_pix_axis - yRF_center_pix));
    
    % Plot RF center
    scatter(x_pix_axis(x_idx), y_pix_axis(y_idx), 100, colors(ichan,:), 'filled');
end

% Set axes equal and fixed limits
axis equal;
xlim([-300 300]);
ylim([-200 200]);

xlabel('X (pixels, RF-centered)');
ylabel('Y (pixels, RF-centered)');
title(['Target positions and RF centers (pixel space) - All channels - ' session_name]);
grid on;

%% plot max activity RF centers and target locations in screen coordinates - example session - stick to this

% Screen info
fn = fieldnames(RF_data_V4);
screenXpix = 1680;
screenYpix = 1050;
centerX = screenXpix / 2;  % 840
centerY = screenYpix / 2;  % 525

% Get session and data
session_name = fn{5}; %select one session
S = RF_data_V4.(session_name);

% Target locations (RF-centered pixels)
x_target_pix = clean_data.trialinfo(:,16);
y_target_pix = clean_data.trialinfo(:,17);

% Convert target locations to screen coordinates (origin at center)
x_target_screen = x_target_pix + centerX;
y_target_screen = y_target_pix + centerY;

% Get all channels
nChannels = length(S.RF_V4_sorted);

% Pixels per degree (defined per session)
pixPerDeg = S.sessInfo.ppd;

% Prepare figure
figure('Color','w');
scatter(x_target_screen, y_target_screen, 50, [0.8 0.8 0.8], 'filled'); % targets in light gray
hold on;

colors = lines(nChannels);  % unique colors for each channel

for ichan = 1:nChannels
    % RF map
    RF_map_struct = S.RF_V4_sorted(ichan);
    RF_map_matrix = RF_map_struct.map;
    [nY, nX] = size(RF_map_matrix);
    
    % Pixel axes centered at RF map center
    x_pix_axis = (1:nX) - (nX+1)/2;      % origin at map center
    y_pix_axis = (nY+1)/2 - (1:nY);      % +Y is up
    
    % Convert RF center from degrees to pixels
    xRF_center_pix = RF_map_struct.xMax * pixPerDeg;
    yRF_center_pix = RF_map_struct.yMax * pixPerDeg;
    
    % Find closest indices in the RF map
    [~, x_idx] = min(abs(x_pix_axis - xRF_center_pix));
    [~, y_idx] = min(abs(y_pix_axis - yRF_center_pix));
    
    % Convert RF center to screen coordinates (origin at center)
    xRF_screen = x_pix_axis(x_idx) + centerX;
    yRF_screen = y_pix_axis(y_idx) + centerY;
    
    % Plot RF center
    scatter(xRF_screen, yRF_screen, 100, colors(ichan,:), 'filled');
end

% Set axes equal and limits to full screen
axis equal;
xlim([0 screenXpix]);
ylim([0 screenYpix]);

xlabel('X (screen pixels, origin at center)');
ylabel('Y (screen pixels, origin at center)');
title(['Target positions and RF centers (screen coordinates) - All channels - ' session_name]);
grid on;

%% Max Activity RF map in pixel space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('/mnt/hpc/projects/MWSampling/4Shivangi/results_hermes/RF_mapping/RF_hermes.mat');  % loads RF_hermes
load('/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks/RF_mapping/RF_klecks.mat');  % loads RF_klecks

% Screen info
screenXpix = 1680;
screenYpix = 1050;
centerX = screenXpix / 2;
centerY = screenYpix / 2;

% Base results folder
results_base = '/mnt/hpc/projects/MWSampling/4Shivangi/results_';

% Monkey data structure
monkeyData = struct('hermes', RF_hermes, 'klecks', RF_klecks);

% Loop over monkeys 
monkeys = fieldnames(monkeyData);  
for iMonkey = 1:numel(monkeys)
    monkey = monkeys{iMonkey};          
    RF_struct = monkeyData.(monkey);

    % Create output folder for PDFs
    outdir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/RF_Mapping', monkey, 'loc_RF_map', 'max_activity');
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    % List all session folders for this monkey
    session_folders = dir(fullfile([results_base monkey], [monkey '_*']));
    folder_names = {session_folders([session_folders(:).isdir]).name};

    % Loop over all RF sessions
    sessions = fieldnames(RF_struct);
    for s = 1:numel(sessions)
        session_name = sessions{s};
        S = RF_struct.(session_name);

        % --- Extract RF mapping date ---
        rf_date_str = regexp(session_name, '\d{8}', 'match');
        if isempty(rf_date_str)
            error('No date found in session name %s', session_name);
        end
        rf_date = datenum(rf_date_str{1}, 'yyyymmdd');

        % --- Find closest behavioral session folder ---
        folder_dates = nan(1, numel(folder_names));
        for f = 1:numel(folder_names)
            date_match = regexp(folder_names{f}, '\d{8}', 'match');
            if ~isempty(date_match)
                folder_dates(f) = datenum(date_match{1}, 'yyyymmdd');
            end
        end

        [~, closest_idx] = min(abs(folder_dates - rf_date));
        closest_folder = folder_names{closest_idx};
        fprintf('RF session %s matched to behavioral folder %s\n', session_name, closest_folder);

        % --- Load clean_lfp.mat ---
        lfp_path = fullfile([results_base monkey], closest_folder, 'clean_lfp.mat');
        if exist(lfp_path, 'file')
            load(lfp_path, 'clean_data');  % loads clean_data
        else
            error('clean_lfp.mat not found in folder %s', closest_folder);
        end

        % --- Pixels per degree ---
        pixPerDeg = S.sessInfo.ppd;

        % --- Target locations ---
        x_target_pix = clean_data.trialinfo(:,16);
        y_target_pix = clean_data.trialinfo(:,17);

        x_target_screen = x_target_pix + centerX;
        y_target_screen = y_target_pix + centerY;

        % --- Channels & layout ---
        nChannels = length(S.RF_V4_sorted);
        nRows = 4; nCols = 4; maxPerFig = nRows * nCols;
        nFigs = min(ceil(nChannels / maxPerFig), 4);

        % --- Output PDF path ---
        outPDF = fullfile(outdir, ['RF_targets_' session_name '.pdf']);

        % --- Plot RFs ---
        for figIdx = 1:nFigs
            startIdx = (figIdx-1)*maxPerFig + 1;
            endIdx   = min(figIdx*maxPerFig, nChannels);
            chans = startIdx:endIdx;

            figure('Color','w','Units','normalized','Position',[0 0 1 1]);
            t = tiledlayout(nRows, nCols, 'TileSpacing','compact', 'Padding','compact');

            for ii = 1:length(chans)
                ichan = chans(ii);
                nexttile; hold on;

                % Plot targets
                scatter(x_target_screen, y_target_screen, 8, 'k', 'filled');

                % RF
                RF = S.RF_V4_sorted(ichan);
                RFmap = RF.map;
                [nY, nX] = size(RFmap);

                x_rf_pix = (1:nX) - (nX + 1)/2;
                y_rf_pix = (nY + 1)/2 - (1:nY);

                xRF_center_pix = RF.xMax * pixPerDeg;
                yRF_center_pix = RF.yMax * pixPerDeg;

                x_screen = x_rf_pix + xRF_center_pix + centerX;
                y_screen = y_rf_pix + yRF_center_pix + centerY;

                h = imagesc(x_screen, y_screen, RFmap);
                set(h, 'AlphaData', 0.9 * mat2gray(RFmap));

                axis equal;
                xlim([0 screenXpix]);
                ylim([0 screenYpix]);
                set(gca, 'YDir','normal', 'XTick',[], 'YTick',[]);

                title(sprintf('Ch %d', ichan), 'FontSize', 8);
            end

            sgtitle(sprintf('RF maps (%d%d)  %s', chans(1), chans(end), strrep(session_name,'_','\_')), ...
                'Interpreter','tex', 'FontSize', 12);

            colormap jet;
            exportgraphics(gcf, outPDF, 'ContentType','vector', 'Append', true);
            close(gcf);
        end
    end
end

%% Gaussian fit RF map in pixel space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load RF data
load('/mnt/hpc/projects/MWSampling/4Shivangi/results_hermes/RF_mapping/RF_hermes.mat');  % loads RF_hermes
load('/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks/RF_mapping/RF_klecks.mat');  % loads RF_klecks

% Screen info
screenXpix = 1680;
screenYpix = 1050;
centerX = screenXpix / 2;
centerY = screenYpix / 2;

% Monkey data structure
monkeyData = struct('hermes', RF_hermes, 'klecks', RF_klecks);

% Threshold for Gaussian fit loading
ThresZ = 3;

% Loop over monkeys 
monkeys = fieldnames(monkeyData);  
for iMonkey = 1:numel(monkeys)
    monkey = monkeys{iMonkey};          
    RF_struct = monkeyData.(monkey);
    
    % Create output folder for PDFs
    outdir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/RF_Mapping', monkey, 'loc_RF_map', 'gaussian_fit');
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end
    
    % Base directory for Gaussian fit results
    gauss_base_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' monkey], 'RF_mapping');
    thres_dir = fullfile(gauss_base_dir, sprintf('ThresZ_%d', ThresZ));
    
    % List all session folders for this monkey
    behavioral_base = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' monkey]);
    session_folders = dir(fullfile(behavioral_base, [monkey '_*']));
    folder_names = {session_folders([session_folders(:).isdir]).name};
    
    % Loop over all RF sessions
    sessions = fieldnames(RF_struct);
    for s = 1:numel(sessions)
        session_name = sessions{s};
        S = RF_struct.(session_name);
        
        fprintf('Processing session: %s\n', session_name);
        
        % --- Extract RF mapping date ---
        rf_date_str = regexp(session_name, '\d{8}', 'match');
        if isempty(rf_date_str)
            warning('No date found in session name %s - skipping', session_name);
            continue;
        end
        rf_date = datenum(rf_date_str{1}, 'yyyymmdd');
        date_token = rf_date_str{1};
        
        % --- Load Gaussian fit results ---
        gauss_results_file = fullfile(thres_dir, sprintf('%s_centerResults.mat', date_token));
        gauss_fit_file = fullfile(thres_dir, sprintf('%s.mat', date_token));
        gauss_extrap_file = fullfile(thres_dir, sprintf('%s_centerExtrap.mat', date_token));
        
        if ~exist(gauss_results_file, 'file') || ~exist(gauss_fit_file, 'file') || ~exist(gauss_extrap_file, 'file')
            warning('Gaussian fit files not found for session %s (date: %s) - skipping', session_name, date_token);
            fprintf('  Looking for files in: %s\n', thres_dir);
            continue;
        end
        
        fprintf('  Loading Gaussian fit files...\n');
        load(gauss_results_file, 'CenterResults');
        load(gauss_fit_file, 'GFit');
        load(gauss_extrap_file, 'CenterExtrap');
        
        % --- Find closest behavioral session folder ---
        folder_dates = nan(1, numel(folder_names));
        for f = 1:numel(folder_names)
            date_match = regexp(folder_names{f}, '\d{8}', 'match');
            if ~isempty(date_match)
                folder_dates(f) = datenum(date_match{1}, 'yyyymmdd');
            end
        end
        [~, closest_idx] = min(abs(folder_dates - rf_date));
        closest_folder = folder_names{closest_idx};
        fprintf('  Matched to behavioral folder: %s\n', closest_folder);
        
        % --- Load clean_lfp.mat ---
        lfp_path = fullfile(behavioral_base, closest_folder, 'clean_lfp.mat');
        if ~exist(lfp_path, 'file')
            warning('clean_lfp.mat not found in folder %s - skipping', closest_folder);
            continue;
        end
        load(lfp_path, 'clean_data');  % loads clean_data
        
        % --- Pixels per degree ---
        pixPerDeg = S.sessInfo.ppd;
        
        % --- Get degree axes ---
        xDeg = S.sessInfo.xDeg;  % corresponds to columns of RF map
        yDeg = S.sessInfo.yDeg;  % corresponds to rows of RF map
        
        % --- Get target locations (RF-centered pixels) ---
        x_target_pix = clean_data.trialinfo(:,16);
        y_target_pix = clean_data.trialinfo(:,17);
        
        % Convert target locations to screen coordinates (origin at center)
        x_target_screen = x_target_pix + centerX;
        y_target_screen = y_target_pix + centerY;
        
        % --- Get number of channels ---
        nChannels = length(S.RF_V4_sorted);
        
        % --- Create figure ---
        fig = figure('Color','w', 'Position', [100 100 1200 900]);
        scatter(x_target_screen, y_target_screen, 50, [0.8 0.8 0.8], 'filled', 'DisplayName', 'Targets');
        hold on;
        colors = lines(nChannels);
        
        % Track which types of centers we plotted (for legend)
        has_gaussian = false;
        has_extrap = false;
        
        % --- Plot Gaussian-fitted RF centers ---
        for ichan = 1:nChannels
            % RF map
            RF_map_struct = S.RF_V4_sorted(ichan);
            RF_map_matrix = RF_map_struct.map;
            
            if isempty(RF_map_matrix)
                continue;
            end
            
            [nY, nX] = size(RF_map_matrix);
            
            % Pixel axes centered at RF map center (SAME AS MAX ACTIVITY)
            x_pix_axis = (1:nX) - (nX+1)/2;      % origin at map center
            y_pix_axis = (nY+1)/2 - (1:nY);      % +Y is up
            
            % Get Gaussian-fit center in RF map indices
            cx_gauss = GFit.FitResults(ichan, 4);  % x0 - column index
            cy_gauss = GFit.FitResults(ichan, 5);  % y0 - row index
            
            % Check if this is a valid center
            hasValidCenter = CenterExtrap.validGaussian(ichan) || CenterExtrap.isExtrapolated(ichan);
            
            if hasValidCenter && ~isnan(cx_gauss) && ~isnan(cy_gauss)
                % Convert from indices to degrees (SAME AS MAX ACTIVITY BUT FROM GAUSSIAN)
                xRF_deg = interp1(1:length(xDeg), xDeg, cx_gauss, 'linear', 'extrap');
                yRF_deg = interp1(1:length(yDeg), yDeg, cy_gauss, 'linear', 'extrap');
                
                % Convert RF center from degrees to pixels (SAME AS MAX ACTIVITY)
                xRF_center_pix = xRF_deg * pixPerDeg;
                yRF_center_pix = yRF_deg * pixPerDeg;
                
                % Find closest indices in the RF map (SAME AS MAX ACTIVITY)
                [~, x_idx] = min(abs(x_pix_axis - xRF_center_pix));
                [~, y_idx] = min(abs(y_pix_axis - yRF_center_pix));
                
                % Convert RF center to screen coordinates (SAME AS MAX ACTIVITY)
                xRF_screen = x_pix_axis(x_idx) + centerX;
                yRF_screen = y_pix_axis(y_idx) + centerY;
                
                % Plot RF center
                if CenterExtrap.validGaussian(ichan)
                    % Valid Gaussian fit - filled circle
                    scatter(xRF_screen, yRF_screen, 100, colors(ichan,:), 'filled', 'HandleVisibility', 'off');
                    has_gaussian = true;
                else
                    % Extrapolated - X marker
                    scatter(xRF_screen, yRF_screen, 100, colors(ichan,:), 'x', 'LineWidth', 2, 'HandleVisibility', 'off');
                    has_extrap = true;
                end
            end
        end
        
        % Add dummy plots for legend
        if has_gaussian
            scatter(nan, nan, 100, 'k', 'filled', 'DisplayName', 'Gaussian centers');
        end
        if has_extrap
            scatter(nan, nan, 100, 'k', 'x', 'LineWidth', 2, 'DisplayName', 'Extrapolated centers');
        end
        
        % --- Finalize plot ---
        axis equal;
        xlim([0 screenXpix]);
        ylim([0 screenYpix]);
        xlabel('X (screen pixels, origin at center)');
        ylabel('Y (screen pixels, origin at center)');
        title(sprintf('%s - Target positions and RF centers (Gaussian fit)', session_name), 'Interpreter', 'none');
        grid on;
        legend('Location', 'best');
        
        % --- Save figure ---
        out_pdf = fullfile(outdir, sprintf('%s_RFcenters_gaussian.pdf', session_name));
        exportgraphics(fig, out_pdf, 'ContentType', 'vector');
        fprintf('  Saved: %s\n', out_pdf);
        
        close(fig);
    end
end

%% Plot Gaussian RF FWHM with target overlap 

% Load RF data
load('/mnt/hpc/projects/MWSampling/4Shivangi/results_hermes/RF_mapping/RF_hermes.mat');
load('/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks/RF_mapping/RF_klecks.mat');

% Screen info
screenXpix = 1680;
screenYpix = 1050;
centerX = screenXpix / 2;
centerY = screenYpix / 2;

% Monkey data structure
monkeyData = struct('hermes', RF_hermes, 'klecks', RF_klecks);

% Threshold for Gaussian fit loading
ThresZ = 3;

% Loop over monkeys 
monkeys = fieldnames(monkeyData);  
for iMonkey = 1%1:numel(monkeys)
    monkey = monkeys{iMonkey};          
    RF_struct = monkeyData.(monkey);
    
    % Create output folder for PDFs
    outdir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/RF_Mapping', monkey, 'loc_RF_map', 'gaussian_overlap');
    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end
    
    % Base directory for Gaussian fit results
    gauss_base_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' monkey], 'RF_mapping');
    thres_dir = fullfile(gauss_base_dir, sprintf('ThresZ_%d', ThresZ));
    
    % List all session folders for this monkey
    behavioral_base = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' monkey]);
    session_folders = dir(fullfile(behavioral_base, [monkey '_*']));
    folder_names = {session_folders([session_folders(:).isdir]).name};
    
    % Loop over all RF sessions
    sessions = fieldnames(RF_struct);
    for s = 3:numel(sessions)
        session_name = sessions{s};
        S = RF_struct.(session_name);
        
        fprintf('Processing session: %s\n', session_name);
        
        % --- Extract RF mapping date ---
        rf_date_str = regexp(session_name, '\d{8}', 'match');
        if isempty(rf_date_str)
            warning('No date found in session name %s - skipping', session_name);
            continue;
        end
        rf_date = datenum(rf_date_str{1}, 'yyyymmdd');
        date_token = rf_date_str{1};
        
        % --- Load Gaussian fit results ---
        gauss_results_file = fullfile(thres_dir, sprintf('%s_centerResults.mat', date_token));
        gauss_fit_file = fullfile(thres_dir, sprintf('%s.mat', date_token));
        gauss_extrap_file = fullfile(thres_dir, sprintf('%s_centerExtrap.mat', date_token));
        
        if ~exist(gauss_results_file, 'file') || ~exist(gauss_fit_file, 'file') || ~exist(gauss_extrap_file, 'file')
            warning('Gaussian fit files not found for session %s (date: %s) - skipping', session_name, date_token);
            continue;
        end
        
        fprintf('  Loading Gaussian fit files...\n');
        load(gauss_results_file, 'CenterResults');
        load(gauss_fit_file, 'GFit');
        load(gauss_extrap_file, 'CenterExtrap');
        
        % --- Find closest behavioral session folder ---
        folder_dates = nan(1, numel(folder_names));
        for f = 1:numel(folder_names)
            date_match = regexp(folder_names{f}, '\d{8}', 'match');
            if ~isempty(date_match)
                folder_dates(f) = datenum(date_match{1}, 'yyyymmdd');
            end
        end
        [~, closest_idx] = min(abs(folder_dates - rf_date));
        closest_folder = folder_names{closest_idx};
        fprintf('  Matched to behavioral folder: %s\n', closest_folder);
        
        % --- Load clean_lfp.mat ---
        lfp_path = fullfile(behavioral_base, closest_folder, 'clean_lfp.mat');
        if ~exist(lfp_path, 'file')
            warning('clean_lfp.mat not found in folder %s - skipping', closest_folder);
            continue;
        end
        load(lfp_path, 'clean_data');
        
        % --- Pixels per degree ---
        pixPerDeg = S.sessInfo.ppd;
        
        % --- Get degree axes ---
        xDeg = S.sessInfo.xDeg;
        yDeg = S.sessInfo.yDeg;
        
        % --- Get target locations ---
        x_target_pix = clean_data.trialinfo(:,16);
        y_target_pix = clean_data.trialinfo(:,17);
        x_target_screen = x_target_pix + centerX;
        y_target_screen = y_target_pix + centerY;
        
        % --- GRID-BASED LOCATION NUMBERING (UPPER-LEFT = 1) ---
        target_coords = [x_target_pix, y_target_pix];
        [unique_targets, ~, orig_ids] = unique(target_coords, 'rows', 'stable');
        
        % Sort by Y descending (top to bottom), then X ascending (left to right)
        [~, sort_idx] = sortrows(unique_targets, [-2 1]);
        unique_targets = unique_targets(sort_idx);
        
        % Remap target IDs to new spatial order
        target_ids = zeros(size(orig_ids));
        for k = 1:numel(sort_idx)
            target_ids(orig_ids == sort_idx(k)) = k;
        end
        
        n_unique_targets = size(unique_targets, 1);
        target_colors = lines(n_unique_targets);
        
        % --- Compute median FWHM from valid Gaussians for extrapolated channels ---
        valid_fwhm = CenterResults.FWHM(CenterExtrap.validGaussian);
        median_fwhm = median(valid_fwhm(~isnan(valid_fwhm)));
        
        % --- Get number of channels ---
        nChannels = length(S.RF_V4_sorted);
        
        % --- Create PDF ---
        out_pdf = fullfile(outdir, sprintf('%s_RF_FWHM_targets.pdf', session_name));
        
        nRows = 4;
        nCols = 4;
        channels_per_page = nRows * nCols;
        n_pages = ceil(nChannels / channels_per_page);
        
        for page = 1:n_pages
            fig = figure('Color','w', 'Position', [50 50 1400 1400]);
            sgtitle(sprintf('%s - Page %d/%d - FWHM and Targets', session_name, page, n_pages), ...
                'Interpreter', 'none', 'FontSize', 14);
            
            start_ch = (page - 1) * channels_per_page + 1;
            end_ch = min(page * channels_per_page, nChannels);
            
            for idx = 1:channels_per_page
                ichan = start_ch + idx - 1;
                if ichan > nChannels
                    break;
                end
                
                subplot(nRows, nCols, idx);
                hold on;
                
                % Plot all targets with colors based on location ID
                for t = 1:length(x_target_screen)
                    loc_id = target_ids(t);
                    scatter(x_target_screen(t), y_target_screen(t), 40, ...
                        target_colors(loc_id, :), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
                end
                
                % Get RF map
                RF_map_struct = S.RF_V4_sorted(ichan);
                RF_map_matrix = RF_map_struct.map;
                
                if ~isempty(RF_map_matrix)
                    [nY, nX] = size(RF_map_matrix);
                    
                    % Check if has any center (valid Gaussian OR extrapolated)
                    hasValidCenter = CenterExtrap.validGaussian(ichan) || CenterExtrap.isExtrapolated(ichan);
                    
                    if hasValidCenter
                        cx_gauss = GFit.FitResults(ichan, 4);
                        cy_gauss = GFit.FitResults(ichan, 5);
                        
                        % Convert center to screen coordinates
                        xRF_deg = interp1(1:length(xDeg), xDeg, cx_gauss, 'linear', 'extrap');
                        yRF_deg = interp1(1:length(yDeg), yDeg, cy_gauss, 'linear', 'extrap');
                        xRF_center_pix = xRF_deg * pixPerDeg;
                        yRF_center_pix = yRF_deg * pixPerDeg;
                        
                        % Pixel axes
                        x_pix_axis = (1:nX) - (nX+1)/2;
                        y_pix_axis = (nY+1)/2 - (1:nY);
                        [~, x_idx] = min(abs(x_pix_axis - xRF_center_pix));
                        [~, y_idx] = min(abs(y_pix_axis - yRF_center_pix));
                        xRF_screen = x_pix_axis(x_idx) + centerX;
                        yRF_screen = y_pix_axis(y_idx) + centerY;
                        
                        % Determine FWHM radius
                        if CenterExtrap.validGaussian(ichan)
                            % Use actual FWHM
                            fwhm_radius_pix = CenterResults.FWHM(ichan) * mean(diff(xDeg)) * pixPerDeg;
                            circle_color = 'b';
                            cross_color = 'b';
                            status_text = 'Valid';
                        else
                            % Use median FWHM for extrapolated
                            fwhm_radius_pix = median_fwhm * mean(diff(xDeg)) * pixPerDeg;
                            circle_color = 'r';
                            cross_color = 'r';
                            status_text = 'Extrap (median FWHM)';
                        end
                        
                        % Draw FWHM circle
                        th = linspace(0, 2*pi, 100);
                        xcirc = xRF_screen + fwhm_radius_pix * cos(th);
                        ycirc = yRF_screen + fwhm_radius_pix * sin(th);
                        
                        if CenterExtrap.validGaussian(ichan)
                            plot(xcirc, ycirc, 'Color', circle_color, 'LineWidth', 2);
                        else
                            plot(xcirc, ycirc, 'Color', circle_color, 'LineWidth', 2, 'LineStyle', '--');
                        end
                        
                        % Mark center
                        plot(xRF_screen, yRF_screen, '+', 'Color', cross_color, 'MarkerSize', 10, 'LineWidth', 2);
                        
                        % Find targets inside FWHM
                        dist_to_center = sqrt((x_target_screen - xRF_screen).^2 + ...
                                             (y_target_screen - yRF_screen).^2);
                        inside_fwhm = dist_to_center <= fwhm_radius_pix;
                        n_inside = sum(inside_fwhm);
                        
                        % Get unique locations inside FWHM
                        if n_inside > 0
                            inside_loc_ids = unique(target_ids(inside_fwhm));
                            loc_str = sprintf('Loc:%s', sprintf('%d,', inside_loc_ids));
                            loc_str = loc_str(1:end-1); % remove trailing comma
                        else
                            loc_str = 'Loc: none';
                        end
                        
                        % Title with channel number and info
                        if CenterExtrap.validGaussian(ichan)
                            title(sprintf('Ch %d | %s\nCenter:(%.0f,%.0f) FWHM:%.0f', ...
                                ichan, loc_str, xRF_screen, yRF_screen, fwhm_radius_pix), ...
                                'FontSize', 9, 'Color', 'k');
                        else
                            title(sprintf('Ch %d [EXTRAP] | %s\nCenter:(%.0f,%.0f) FWHM:%.0f*', ...
                                ichan, loc_str, xRF_screen, yRF_screen, fwhm_radius_pix), ...
                                'FontSize', 9, 'Color', 'r');
                        end
                    else
                        title(sprintf('Ch %d - No fit', ichan), ...
                            'FontSize', 9, 'Color', [0.5 0.5 0.5]);
                    end
                else
                    title(sprintf('Ch %d - No data', ichan), ...
                        'FontSize', 9, 'Color', [0.5 0.5 0.5]);
                end
                
                xlim([0 screenXpix]);
                ylim([0 screenYpix]);
                axis equal;
                grid on;
                set(gca, 'FontSize', 8);
                xlabel('X (pixels)', 'FontSize', 8);
                ylabel('Y (pixels)', 'FontSize', 8);
                hold off;
            end
            
            % Save page
            if page == 1
                exportgraphics(fig, out_pdf, 'ContentType', 'vector');
            else
                exportgraphics(fig, out_pdf, 'ContentType', 'vector', 'Append', true);
            end
            
            close(fig);
        end
        
        fprintf('  Saved: %s\n', out_pdf);
        
        % --- Also create a summary table ---
        summary_file = fullfile(outdir, sprintf('%s_channel_target_summary.txt', session_name));
        fid = fopen(summary_file, 'w');
        fprintf(fid, 'Channel\tStatus\tRF_Center_X\tRF_Center_Y\tFWHM_radius\tLocations_Inside\tN_Targets_Inside\n');
        
        for ichan = 1:nChannels
            RF_map_matrix = S.RF_V4_sorted(ichan).map;
            
            if isempty(RF_map_matrix)
                fprintf(fid, '%d\tNo_Data\t-\t-\t-\t-\t-\n', ichan);
                continue;
            end
            
            [nY, nX] = size(RF_map_matrix);
            
            hasValidCenter = CenterExtrap.validGaussian(ichan) || CenterExtrap.isExtrapolated(ichan);
            
            if hasValidCenter
                cx_gauss = GFit.FitResults(ichan, 4);
                cy_gauss = GFit.FitResults(ichan, 5);
                
                xRF_deg = interp1(1:length(xDeg), xDeg, cx_gauss, 'linear', 'extrap');
                yRF_deg = interp1(1:length(yDeg), yDeg, cy_gauss, 'linear', 'extrap');
                xRF_center_pix = xRF_deg * pixPerDeg;
                yRF_center_pix = yRF_deg * pixPerDeg;
                
                x_pix_axis = (1:nX) - (nX+1)/2;
                y_pix_axis = (nY+1)/2 - (1:nY);
                [~, x_idx] = min(abs(x_pix_axis - xRF_center_pix));
                [~, y_idx] = min(abs(y_pix_axis - yRF_center_pix));
                xRF_screen = x_pix_axis(x_idx) + centerX;
                yRF_screen = y_pix_axis(y_idx) + centerY;
                
                if CenterExtrap.validGaussian(ichan)
                    fwhm_radius_pix = CenterResults.FWHM(ichan) * mean(diff(xDeg)) * pixPerDeg;
                    status = 'Valid_Gaussian';
                else
                    fwhm_radius_pix = median_fwhm * mean(diff(xDeg)) * pixPerDeg;
                    status = 'Extrapolated';
                end
                
                dist_to_center = sqrt((x_target_screen - xRF_screen).^2 + ...
                                     (y_target_screen - yRF_screen).^2);
                inside_fwhm = dist_to_center <= fwhm_radius_pix;
                n_inside = sum(inside_fwhm);
                
                if n_inside > 0
                    inside_loc_ids = unique(target_ids(inside_fwhm));
                    loc_str = sprintf('%d,', inside_loc_ids);
                    loc_str = loc_str(1:end-1);
                else
                    loc_str = 'none';
                end
                
                fprintf(fid, '%d\t%s\t%.1f\t%.1f\t%.1f\t%s\t%d\n', ...
                    ichan, status, xRF_screen, yRF_screen, fwhm_radius_pix, loc_str, n_inside);
            else
                fprintf(fid, '%d\tNo_Fit\t-\t-\t-\t-\t-\n', ichan);
            end
        end
        
        fclose(fid);
        fprintf('  Summary saved: %s\n', summary_file);
    end
end
