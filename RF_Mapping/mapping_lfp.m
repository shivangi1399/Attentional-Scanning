clearvars
close all
clc

% Code description:
% -----------------
% Loads session RF mapping data, extracts and sorts V4 channels, normalizes and reorders RF maps, 
% fits Gaussian surfaces to each channel's RF to determine centers and FWHM, extrapolates missing 
% or poor fits using neighbors or plane fits, and visualizes RFs with Gaussian centers and FWHM 
% masks across sessions, saving results and plots separately for Hermes and Klecks.

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
    rfFiles = dir(fullfile(sessionPath, '*.stimOnlfpx.RF'));
    
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

    sessStr = D.sessInfo.session;
    tok = regexp(sessStr, '^(hermes|klecks)_', 'tokens', 'once');
    monkey = tok{1};

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

%% Save and plot RFs

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

fn = fieldnames(RF_data_V4);

for i = 1:numel(fn)

    session_name = fn{i};
    S = RF_data_V4.(session_name);
    probe_fields = S.V4_labels{:, 1};

    % ----- Determine if Hermes or Klecks -----

    if contains(session_name, 'hermes', 'IgnoreCase', true)
        isHermes = true;
        outdir = outdir_hermes;
        ch_offset = 0;
    elseif contains(session_name, 'klecks', 'IgnoreCase', true)
        isHermes = false;
        outdir = outdir_klecks;
        ch_offset = 64; % Klecks (65-128) maps to 1-64 positions
    else
        warning('Skipping unknown session type: %s', session_name);
        continue
    end

    % ----- Plotting in a loop -----
    fig = figure('Position',[100 100 1200 900]); % Slightly taller for 8 rows

    % Define grid dimensions
    nRows = 8;
    nCols = 8;

    for k = 1:numel(probe_fields)
        probeName = probe_fields{k};
        probeData = S.RF_V4_sorted(k);

        % 1. Extract the channel number
        tmp = regexp(probeName, '\d+$', 'match');
        if isempty(tmp)
            warning('Could not extract number from %s', probeName);
            continue;
        end
        chNum = str2double(tmp{1});

        % 2. Normalize channel number for the grid (1-64)
        relCh = chNum - ch_offset;

        % Safety check to ensure we stay within grid bounds
        if relCh < 1 || relCh > 64
            warning('Channel %d is out of range for the 8x8 grid logic.', chNum);
            continue;
        end

        % 3. Subplot Index
        col_idx = ceil(relCh / nRows);
        row_idx = mod(relCh - 1, nRows) + 1;
        subplot_idx = (row_idx - 1) * nCols + col_idx;

        % 4. Plot
        subplot(nRows, nCols, subplot_idx);

        if isfield(probeData, 'map')
            imagesc(S.sessInfo.xDeg , S.sessInfo.yDeg, probeData.map);% find the xDeg and yDeg for lfp map
            axis image;
            colormap jet;

            set(gca, 'XTick', [], 'YTick', []);
            title(probeName, 'Interpreter', 'none', 'FontSize', 8);
        end
    end

    sgtitle(sprintf('%s RF Mapping', session_name), 'Interpreter','none');

    pdfname = fullfile(outdir, sprintf('%s_RFmapping_lfp.pdf', session_name));
    exportgraphics(fig, pdfname, 'ContentType','vector');

    close(fig);
end

%% Mapping the max activity centers on screen

fn = fieldnames(RF_data_V4);

for isess = 1%:numel(fn)

    % 1. Extract Labels for this Session
    session_name = fn{isess};
    S = RF_data_V4.(session_name);
    probe_fields = S.V4_labels{:, 1};

    % 2. Setup Parameters - check these, might change depending on session!
    screenXpix = 1680;
    screenYpix = 1050;
    pixPerDeg  = S.sessInfo.ppd;

    centerX = screenXpix / 2;
    centerY = screenYpix / 2;

    % 3. Visualization Setup
    figure('Color', 'w', 'Name', ['RF Map: ' session_name]);
    rectangle('Position', [0, 0, screenXpix, screenYpix], 'EdgeColor', 'k', 'LineWidth', 2);
    hold on;
    set(gca, 'YDir', 'reverse'); % 0,0 is top-left
    axis equal;
    grid on;
    xlim([-100, screenXpix + 100]);
    ylim([-100, screenYpix + 100]);
    xlabel('Screen X (pixels)');
    ylabel('Screen Y (pixels)');
    title(['RF Positions - ' session_name], 'Interpreter', 'none');

    % 4. Calculate Positions and Plot with Extracted Labels
    numRFs = size((S.V4_labels),1);

    for ind = 1:numRFs

        % --- Calculate Position ---
        degX = S.RF_V4_sorted(ind).xMax;
        degY = S.RF_V4_sorted(ind).yMax;

        x_pix = centerX + (degX * pixPerDeg);
        y_pix = centerY - (degY * pixPerDeg);

        % --- Determine Label ---
        % Default to index
        labelStr = num2str(ind);

        % Try to grab the specific label from the list we extracted
        if ~isempty(probe_fields) && ind <= length(probe_fields)
            rawLabel = probe_fields(ind);

            % Handle if the label inside the cell is a number or a string
            if iscell(rawLabel)
                labelStr = char(rawLabel);
            elseif isnumeric(rawLabel)
                labelStr = num2str(rawLabel);
            else
                labelStr = string(rawLabel);
            end
        end

        % --- Plotting ---
        tealColor = [0, 0.7, 0.7];
        plot(x_pix, y_pix, '.', 'Color', tealColor, 'MarkerSize', 25);

        % Add text label
        text(x_pix + 15, y_pix, labelStr, ...
            'FontSize', 9, ...
            'Color', 'k', ...
            'Interpreter', 'none');
    end

    hold off;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fitting the gaussian to find center of RFs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fn = fieldnames(RF_data_V4);
outdir_hermes = '/mnt/hpc/projects/MWSampling/4Shivangi/results_hermes/RF_mapping';
outdir_klecks = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks/RF_mapping';

for isess = 1:numel(fn)

    session_name = fn{isess};
    S = RF_data_V4.(session_name);

    % Determine if Hermes or Klecks
    if contains(session_name,'hermes','IgnoreCase',true)
        outdir = outdir_hermes;
    elseif contains(session_name,'klecks','IgnoreCase',true)
        outdir = outdir_klecks;
    else
        continue
    end

    %% Fitting the Gaussian

    numChan = numel(S.RF_V4_sorted);   % should be 64

    FitResults = nan(numChan, 6);      % [xIdx, yIdx, sigma, x0, y0, r2]
    GaussChan  = cell(numChan, 1);     % store Gaussian surface for each channel

    ThresZ = 3; % Z-threshold for noise

    for ichan = 1:numChan

        % 1) Extract RF map (nX × nY)-----
        Zraw = S.RF_V4_sorted(ichan).map;

        if isempty(Zraw) || all(isnan(Zraw), 'all')
            continue;
        end

        % Dimensions
        [nX, nY] = size(Zraw);

        % XY grid for fitting
        [xi, yi] = meshgrid(1:nY, 1:nX);

        % 2) Z-score RF -----
        mu = mean(Zraw(:), 'omitnan');
        sd = std(Zraw(:), 0, 'omitnan');

        if sd == 0 || isnan(sd)
            continue;
        end

        Z = (Zraw - mu) ./ sd;

        % 3) Threshold noise -----
        Z(Z < ThresZ) = 0;
        Z(isnan(Z))   = 0;

        if ~any(Z(:))
            continue;
        end

        % 4) Fit Gaussian -----
        opts = struct();
        opts.iso = true;
        opts.positive = true;

        results = autoGaussianSurf(xi, yi, Z, opts);

        % 5) Convert Gaussian center to pixel indices -----
        [~, xIdx] = min(abs(xi(1,:) - results.x0));
        [~, yIdx] = min(abs(yi(:,1) - results.y0));

        % 6) Save results -----
        FitResults(ichan, :) = [
            xIdx, ...        % 1
            yIdx, ...        % 2
            results.sigma, ...  % 3
            results.x0, ...     % 4
            results.y0, ...     % 5
            results.r2 ...      % 6
            ];

        GaussChan{ichan} = results.G;   % the fitted Gaussian surface

    end

    % Pack output
    GFit = struct();
    GFit.FitResults = FitResults;
    GFit.GaussChan  = GaussChan;
    GFit.nX = nX;
    GFit.nY = nY;

    % --- Save ---
    date_token = regexp(session_name, '\d{8}', 'match', 'once');
    thres_dir = fullfile(outdir, sprintf('ThresZ_%d', ThresZ));
    if ~exist(thres_dir,'dir'); mkdir(thres_dir); end

    outfile = sprintf('%s.mat', date_token);
    save(fullfile(thres_dir, outfile), 'GFit');
    clear  FitResults GaussChan

    %% Finding the center

    % Preallocate
    numChan = size(GFit.FitResults, 1);
    masked = nan(numChan, GFit.nX, GFit.nY);
    FWHM   = nan(numChan, 1);

    % Grid
    [Y, X] = meshgrid(1:GFit.nX, 1:GFit.nY);

    for ichan = 1:numChan

        r2 = GFit.FitResults(ichan, 6);
        if r2 <= 0
            continue;
        end

        % Extract parameters
        sigma   = GFit.FitResults(ichan, 3);
        centerX = GFit.FitResults(ichan, 1);   % xIdx (pixel)
        centerY = GFit.FitResults(ichan, 2);   % yIdx (pixel)

        % Convert Gaussian surface from cell
        gaussianSurface = GFit.GaussChan{ichan};
        if isempty(gaussianSurface)
            continue;
        end

        % ---- Compute FWHM ----
        % FWHM = 2*sqrt(2*ln(2))*sigma
        FWHM(ichan) = 2 * sqrt(2 * log(2)) * sigma;

        % Use FWHM radius as threshold
        threshold = FWHM(ichan);

        % ---- Create circular mask ----
        circularMask = ((X - centerX).^2 + (Y - centerY).^2) <= threshold^2;

        % ---- Apply mask ----
        tmp = gaussianSurface;
        tmp(~circularMask) = nan;

        masked(ichan,:,:) = tmp;
    end

    % --- Save ---
    CenterResults.masked = masked;
    CenterResults.FWHM   = FWHM;

    save(fullfile(thres_dir, sprintf('%s_centerResults.mat', date_token)), ...
        'CenterResults');

    %% Extrapolate RF centers for channels with failed Gaussian fits

    numChan = size(GFit.FitResults,1);   % should be 64
    nRows = 8;
    nCols = 8;

    % Extract existing centers and fit quality
    x0 = GFit.FitResults(:,4);
    y0 = GFit.FitResults(:,5);
    r2 = GFit.FitResults(:,6);

    validFit = r2 > 0 & ~isnan(x0) & ~isnan(y0);

    % --- Manual override: refit of selected channels ---

    ForceExtrapBySession = struct();

    % Hermes
    ForceExtrapBySession.hermes_20170418 = [2 5 6 7 19];
    ForceExtrapBySession.hermes_20170711 = [5 9 64];
    ForceExtrapBySession.hermes_20170808 = [27 32 37 41 49 64];
    ForceExtrapBySession.hermes_20170823 = [5 11];
    ForceExtrapBySession.hermes_20170829 = [3 5 11 64];

    % Klecks
    ForceExtrapBySession.klecks_20170829 = [7 15 37 55];
    ForceExtrapBySession.klecks_20171020 = ...
        [4 11 12 14 15 28 31 32 36 37 40 43 46 47 50 53 55 64];

    % Identify monkey and date
    monkey = '';
    if contains(session_name,'hermes','IgnoreCase',true), monkey='hermes'; end
    if contains(session_name,'klecks','IgnoreCase',true), monkey='klecks'; end
    date_token = regexp(session_name, '\d{8}', 'match', 'once');

    forcedChannels = [];

    if ~isempty(date_token)
        key = sprintf('%s_%s', monkey, date_token);
        if isfield(ForceExtrapBySession, key)
            forcedChannels = ForceExtrapBySession.(key);
        end
    end

    % Safety + formatting
    forcedChannels = forcedChannels(:);
    forcedChannels = forcedChannels(forcedChannels >= 1 & forcedChannels <= numChan);

    % Force these channels to be treated as bad Gaussian fits
    if ~isempty(forcedChannels)
        validFit(forcedChannels) = false;
        fprintf('[%s] Forced extrapolation: %s\n', ...
            key, mat2str(forcedChannels'));
    end


    % --- Preallocate extrapolated centers ---
    x0_ext = x0;
    y0_ext = y0;
    isExtrapolated = false(numChan,1);

    % --- Helper to convert channel index to grid position ---
    chan2grid = @(ch) deal( ...
        mod(ch-1, nRows) + 1, ...         % row
        ceil(ch / nRows) );               % col

    % --- Helper to convert grid position to channel index ---
    grid2chan = @(r,c) (c-1)*nRows + r;

    % --- Loop over channels with bad fits ---
    for ch = 1:numChan

        if validFit(ch)
            continue
        end

        [r, c] = chan2grid(ch);

        neighborCenters = [];

        % 8-connected neighbors
        for dr = -1:1
            for dc = -1:1
                if dr == 0 && dc == 0
                    continue
                end

                rr = r + dr;
                cc = c + dc;

                if rr < 1 || rr > nRows || cc < 1 || cc > nCols
                    continue
                end

                nch = grid2chan(rr,cc);

                if validFit(nch)
                    neighborCenters(end+1,:) = [x0(nch), y0(nch)]; %#ok<SAGROW>
                end
            end
        end

        % --- If neighbors exist, use median ---
        if ~isempty(neighborCenters)
            x0_ext(ch) = median(neighborCenters(:,1));
            y0_ext(ch) = median(neighborCenters(:,2));
            isExtrapolated(ch) = true;
        end
    end

    % Global fallback for any remaining failures (plane fit)
    stillBad = ~validFit & ~isExtrapolated;

    if any(stillBad)

        % Grid coordinates
        [rowGrid, colGrid] = arrayfun(chan2grid, (1:numChan)');
        rowGrid = rowGrid(:);
        colGrid = colGrid(:);

        % Fit planes x0 = a*r + b*c + d
        px = fit([rowGrid(validFit), colGrid(validFit)], x0(validFit), 'poly11');
        py = fit([rowGrid(validFit), colGrid(validFit)], y0(validFit), 'poly11');

        for ch = find(stillBad)'
            x0_ext(ch) = px(rowGrid(ch), colGrid(ch));
            y0_ext(ch) = py(rowGrid(ch), colGrid(ch));
            isExtrapolated(ch) = true;
        end
    end

    % Store results in a familiar structure
    CenterExtrap = struct();
    CenterExtrap.x0 = x0_ext;
    CenterExtrap.y0 = y0_ext;
    CenterExtrap.isExtrapolated = isExtrapolated;
    CenterExtrap.validGaussian = validFit;

    % Update GFit
    GFit.FitResults(:,4) = x0_ext;
    GFit.FitResults(:,5) = y0_ext;

    save(fullfile(thres_dir, sprintf('%s.mat', date_token)), 'GFit');
    save(fullfile(thres_dir, sprintf('%s_centerExtrap.mat', date_token)), ...
        'CenterExtrap');

end

%% Plot RFs with Gaussian-fit centers and threshold masks

for isess = 1:numel(fn)

    session_name = fn{isess};
    S = RF_data_V4.(session_name);

    ThresZ = 3;

    % ----- Determine the correct directory for this session -----
    if contains(session_name, 'hermes', 'IgnoreCase', true)
        base_dir = '/mnt/hpc/projects/MWSampling/4Shivangi/results_hermes/RF_mapping';
        plot_dir = '/mnt/hpc/projects/MWSampling/4Shivangi/Plots/RF_Mapping/hermes';
    elseif contains(session_name, 'klecks', 'IgnoreCase', true)
        base_dir = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks/RF_mapping';
        plot_dir = '/mnt/hpc/projects/MWSampling/4Shivangi/Plots/RF_Mapping/klecks';
    else
        warning('Unknown session type: %s  skipping.', session_name);
        continue
    end

    thres_dir = fullfile(base_dir, sprintf('ThresZ_%d', ThresZ));

    % ----- Load previously saved Gaussian fit results -----
    date_token = regexp(session_name, '\d{8}', 'match', 'once');

    load(fullfile(thres_dir, sprintf('%s_centerResults.mat', date_token)), ...
        'CenterResults');
    load(fullfile(thres_dir, sprintf('%s.mat', date_token)), 'GFit');
    load(fullfile(thres_dir, sprintf('%s_centerExtrap.mat', date_token)), ...
        'CenterExtrap');


    masked = CenterResults.masked;
    FWHM   = CenterResults.FWHM;

    % ----- Plot all channels -----
    fig = figure('Position',[200 100 1400 900]);
    sgtitle(sprintf('%s  RFs with Gaussian Centers + FWHM Mask', session_name), ...
        'Interpreter','none');

    nRows = 8;
    nCols = 8;

    for ch = 1:size(GFit.FitResults,1)

        subplot(nRows, nCols, ch)

        % Raw RF map
        Zraw = S.RF_V4_sorted(ch).map;
        if isempty(Zraw)
            continue
        end

        % Gaussian Fit Info
        cx    = GFit.FitResults(ch,4);   % x0
        cy    = GFit.FitResults(ch,5);   % y0
        sigma = GFit.FitResults(ch,3);
        r2    = GFit.FitResults(ch,6);

        imagesc(Zraw);
        hold on
        colormap parula
        axis image off

        % --- CASE 1: Valid Gaussian fit ---
        if CenterExtrap.validGaussian(ch)

            % plot Gaussian-fit center
            plot(cx, cy, 'wo', 'MarkerSize', 2, 'LineWidth', 1.5);

            % plot FWHM mask
            fwhm_radius = FWHM(ch);
            th = linspace(0, 2*pi, 180);
            xcirc = cx + fwhm_radius * cos(th);
            ycirc = cy + fwhm_radius * sin(th);
            plot(xcirc, ycirc, 'w-', 'LineWidth', 1.0);

            title(sprintf('Ch %d', ch), 'FontSize', 8)

            % --- CASE 2: Extrapolated channel ---
        elseif CenterExtrap.isExtrapolated(ch)

            % plot extrapolated center ONLY
            plot(cx, cy, 'rx', 'MarkerSize', 6, 'LineWidth', 1.5);
            title(sprintf('Ch %d (extrap)', ch), 'FontSize', 8)

            % --- CASE 3: No center at all ---
        else
            title(sprintf('Ch %d  no center', ch), 'FontSize', 8)
        end

    end

    % ----- Save -----
    plotdirt = fullfile(plot_dir, sprintf('extr_ThresZ_%d', ThresZ));
    if ~exist(plotdirt,'dir'); mkdir(plotdirt); end
    out_pdf = fullfile(plotdirt, sprintf('%s_RFGaussian.pdf', session_name));
    exportgraphics(fig, out_pdf, 'ContentType', 'vector');
    close(fig)

end


