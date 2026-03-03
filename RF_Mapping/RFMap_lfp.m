clearvars
close all
clc

% Code description:
% -----------------
% This code extracts high-gamma LFP power, aligns it to bar stimulus geometry, and back-projects it 
% to estimate spatial receptive fields for V4 electrodes, automatically finding the best response 
% latency

%% paths
rootDir = '/mnt/hpc/projects/MWSampling/RF_mappings';
addpath /mnt/hpc/opt/fieldtrip_esi
ft_defaults

%% discover sessions
items = dir(rootDir);
isSession = [items.isdir] & ~startsWith({items.name}, '.');
sessionList = {items(isSession).name};

%% loop sessions

for s = 1:numel(sessionList)

    sessionName = sessionList{s};
    sessionPath = fullfile(rootDir, sessionName);
    fprintf('\n=== %s ===\n', sessionName);

    %% find stimOn lfp
    lfpFile = dir(fullfile(sessionPath, '*.stimOn.lfp'));

    if isempty(lfpFile)
        fprintf('  No stimOn.lfp found, skipping\n');
        continue
    end

    lfpPath = fullfile(sessionPath, lfpFile(1).name);
    load(lfpPath, '-mat');

    %% trial selection
    trial_endtime = 3;

    cfg = [];
    cfg.latency = [0 trial_endtime];
    data2 = ft_selectdata(cfg, data);

    %% frequency analysis
    cfg = [];
    cfg.method     = 'mtmconvol';
    cfg.tapsmofrq  = 20;
    cfg.output     = 'pow';
    cfg.toi        = 0:0.01:trial_endtime;
    cfg.foi        = 40:1:100;
    cfg.t_ftimwin  = ones(size(cfg.foi)) * 0.1;
    cfg.keeptrials = 'yes';

    pow = ft_freqanalysis(cfg, data2);

    %% average over frequency
    cfg = [];
    cfg.avgoverfreq = 'yes';
    trials = ft_selectdata(cfg, pow);
    trials.fsample = data.fsample;

    %% RF backprojection
    latencyRng      = 0:0.01:0.5;
    directionColumn = 3;

    [map, bestLatency] = ft_barmap_backproject(trials, directionColumn, latencyRng);

    %% build RF struct
    nChan = size(map, 3);

    RF = struct( ...
        'map',         cell(1, nChan), ...
        'latency',     cell(1, nChan), ...
        'xMax',        cell(1, nChan), ...
        'yMax',        cell(1, nChan) ...
    );

    for ch = 1:nChan
        m = squeeze(map(:,:,ch));
        [~, idx] = max(m(:));
        [yMax, xMax] = ind2sub(size(m), idx);

        RF(ch).map         = m;
        RF(ch).latency     = bestLatency(ch);
        RF(ch).xMax        = xMax;
        RF(ch).yMax        = yMax;
    end

    %% sessInfo
    sessInfo = struct();
    sessInfo.session   = sessionName;
    sessInfo.fsample   = trials.fsample;
    sessInfo.nChannels = nChan;

    %% assemble + save
    RF_data = struct();
    RF_data.sessInfo = sessInfo;
    RF_data.RF       = RF;
    RF_data.label    = trials.label(:);

    outFile = fullfile(sessionPath, [sessionName '.stimOnlfpx.RF']);
    save(outFile, '-struct', 'RF_data', '-mat');

    fprintf('  Saved %s\n', outFile);
end

