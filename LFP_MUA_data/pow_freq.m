clear all
close all
clc

% Code description:
% -----------------
% Power spectrum analysis across all channels and sessions using FieldTrip
% use fooof and baseline trials to remove 1/f

%% Specify paths

addpath /opt/fieldtrip_github/
ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/eye_data
clc

%% Paths for data

datafolder = '/mnt/hpc/projects/MWSampling/4Shivangi/data_Klecks';
outputfolder   = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks';

cd(outputfolder),
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

data_paths = [];
data_paths = cellfun(@(x) fullfile(datafolder,x), session_names, 'uniform',0);

session_paths = [];
session_paths = cellfun(@(x) fullfile(outputfolder,x), session_names, 'uniform',0);

session_paths_files = [];
session_paths_files = cellfun(@(x) fullfile(outputfolder,x, 'clean_data.mat'), session_names, 'uniform',0);

%% Using fooof to remove 1/f

all_pow = containers.Map;
all_labels = {};

for isess = 1:length(session_paths_files)

    ESIload(session_paths_files{isess});  
    load(fullfile(session_paths{isess}, 'zlfptrials.mat')); 

    % Trial selection
    hitIdx  = find(clean_data.trialinfo(:,20) == 1);
    missIdx = find(clean_data.trialinfo(:,20) == 5);
    
    cfg=[];
    cfg.toilim    = [-1 0]; % change this for bp (baseline period) or rp (response period)
    dat = ft_redefinetrial(cfg, zlfptrials);

    % Frequency analysis (aperiodic)
    cfg = [];
    cfg.method     = 'mtmfft';
    cfg.taper      = 'hann';
    cfg.keeptrials = 'no';
    cfg.foi        = 2:2:80;
    cfg.pad        = 'nextpow2';
    cfg.trials     = [hitIdx; missIdx];
    cfg.output     = 'fooof_aperiodic';  % 1/f only
    fractal        = ft_freqanalysis(cfg, dat);

    % Frequency analysis
    cfg.output     = 'pow';
    original       = ft_freqanalysis(cfg, dat);

    % Flattened spectra = log10(original) - log10(fractal)
    flatten = log10(original.powspctrm) - log10(fractal.powspctrm);

    % Store per channel
    for ch = 1:length(original.label)
        chName = original.label{ch};
        if ~isKey(all_pow, chName)
            all_pow(chName) = [];
            all_labels{end+1} = chName;
        end
        tmp = all_pow(chName);
        tmp = cat(1, tmp, flatten(ch,:));
        all_pow(chName) = tmp;
    end
end

% Average across sessions
nChans = numel(all_labels);
freqs  = original.freq;
mean_flatten = zeros(nChans, length(freqs));
se_flatten   = zeros(nChans, length(freqs));

for i = 1:nChans
    chName = all_labels{i};
    data = all_pow(chName);          % sessions x frequencies
    mean_flatten(i,:) = mean(data,1);
    se_flatten(i,:)   = std(data,0,1) ./ sqrt(size(data,1)); % SE
end

% plot
peakFreqs = zeros(nChans, 1);
figure;
for ch = 1:nChans
    subplot(8,8,ch);
    hold on;
    
    xVals = freqs(1:end);
    meanVals = mean_flatten(ch,1:end);
    seVals   = se_flatten(ch,1:end);
    
    % Shaded error patch
    fill([xVals fliplr(xVals)], ...
         [meanVals+seVals fliplr(meanVals-seVals)], ...
         [0.8 0.8 0.8], 'EdgeColor','blue', 'FaceAlpha', 0.4);
    
    % Mean line
    plot(xVals, meanVals, 'k', 'LineWidth', 1.5);
    
    % Peak frequency
    [peakPower, peakIdx] = max(meanVals);
    peakFreq = xVals(peakIdx);
    peakFreqs(ch) = peakFreq;
    plot(peakFreq, peakPower, 'ro', 'MarkerSize', 6, 'LineWidth', 1.5);
    
    title(all_labels{ch}, 'Interpreter','none', 'FontSize', 6);
    xlim([xVals(1) xVals(end)]);
    
    if mod(ch,8)==1
        ylabel('Power (log10)');
    end
    if ch > nChans-8
        xlabel('Frequency (Hz)');
    end
    set(gca,'FontSize',6);
    hold off;
end
sgtitle('1/f Removed Power Spectra with Peak Frequencies');

% Average across all channels
overall_mean = mean(mean_flatten,1);   % mean across channels
overall_se   = std(mean_flatten,0,1) ./ sqrt(size(mean_flatten,1)); % SE across channels
figure;
hold on;
fill([freqs fliplr(freqs)], ...
     [overall_mean+overall_se fliplr(overall_mean-overall_se)], ...
     [0.8 0.8 0.8], 'EdgeColor','blue', 'FaceAlpha', 0.4);
plot(freqs, overall_mean, 'k', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Power (log10)');
title('1/f Removed Spectrum (Mean ± SE across sessions and channels)');
set(gca,'FontSize',10);
box on;

%% Using baseline trials to remove 1/f

all_pow = containers.Map;
all_labels = {};

for isess = 1:length(session_paths_files)
    
    %% baseline trials
    cd(data_paths{isess})
    ESIload('V4_lfp_data_noTarg.mat')
    
    % zscoring
    disp(strcat('session- ',num2str(isess),' out of- ',num2str(length(session_names)), ', running z-scoring for baseline trials')) %#ok<UNRCH>
    zlfpTrials_b = fun_zscore_session(lfpTrials);
    
    % keep only correct trials and baseline condition (cond 2: catch trials, cond 3: baseline)
    hits = [];
    hits = find(zlfpTrials_b.trialinfo(:,20)==1 & zlfpTrials_b.trialinfo(:,15)==3);
    misses = [];
    misses = find(zlfpTrials_b.trialinfo(:,20)==5 & zlfpTrials_b.trialinfo(:,15)==3);
    
    % redefine trials
    cfg=[];
    cfg.trials = [hits' misses']; 
    cfg.toilim    = [0.3 1.3];
    datb = ft_redefinetrial(cfg, zlfpTrials_b);
    
    % power spectra
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmfft';
    cfg.taper        = 'hann';
    cfg.keeptrials = 'no';
    cfg.foi          = 2:2:80;
    freqpow_b = ft_freqanalysis(cfg,datb);
    
    %% Baseline period of target trials
    
    ESIload(session_paths_files{isess});  
    load(fullfile(session_paths{isess}, 'zlfptrials.mat')); 

    % Trial selection
    hitIdx  = find(clean_data.trialinfo(:,20) == 1);
    missIdx = find(clean_data.trialinfo(:,20) == 5);
    
    cfg=[];
    cfg.toilim    = [-1 0]; % change this for bp (baseline period) or rp (response period)
    dat = ft_redefinetrial(cfg, zlfptrials);

    % Frequency analysis
    cfg = [];
    cfg.method     = 'mtmfft';
    cfg.taper      = 'hann';
    cfg.keeptrials = 'no';
    cfg.output     = 'pow';
    cfg.foi        = 2:2:80;
    cfg.trials     = [hitIdx; missIdx];
    freqpow_tbp       = ft_freqanalysis(cfg, dat);

    % Restrict to common channels
    commonChans = intersect(freqpow_b.label, freqpow_tbp.label, 'stable');
    [~, idxB]   = ismember(commonChans, freqpow_b.label);
    [~, idxTBP] = ismember(commonChans, freqpow_tbp.label);

    pow_b_common   = freqpow_b.powspctrm(idxB, :);
    pow_tbp_common = freqpow_tbp.powspctrm(idxTBP, :);

    %% Remove 1/f component (difference)
    PF = pow_tbp_common-pow_b_common;

    % Store results per channel
    for ch = 1:length(commonChans)
        chName = commonChans{ch};
        if ~isKey(all_pow, chName)
            all_pow(chName) = [];
            all_labels{end+1} = chName;
        end
        tmp = all_pow(chName);
        tmp = cat(1, tmp, PF(ch,:)); % append new session
        all_pow(chName) = tmp;
    end
end

%% Aggregate results across sessions 
freqs = freqpow_b.freq; % use baseline structure for frequencies
nChans = numel(all_labels);
mean_pf = zeros(nChans, length(freqs));
se_pf   = zeros(nChans, length(freqs));

for i = 1:nChans
    chName = all_labels{i};
    data = all_pow(chName); % sessions x frequencies
    mean_pf(i,:) = mean(data, 1);
    se_pf(i,:)   = std(data, 0, 1) ./ sqrt(size(data,1));
end

% Plot per channel
peakFreqs = zeros(nChans, 1);
figure;
for ch = 1:nChans
    subplot(8,8,ch); hold on;
    xVals = freqs(1:16);
    meanVals = mean_pf(ch,1:16);
    seVals   = se_pf(ch,1:16);

    % Shaded error patch
    fill([xVals fliplr(xVals)], ...
         [meanVals+seVals fliplr(meanVals-seVals)], ...
         [0.8 0.8 0.8], 'EdgeColor','blue', 'FaceAlpha', 0.4);

    % Mean line
    plot(xVals, meanVals, 'k', 'LineWidth', 1.5);

    % Peak frequency
    [peakPower, peakIdx] = max(meanVals);
    peakFreqs(ch) = xVals(peakIdx);
    plot(xVals(peakIdx), peakPower, 'ro', 'MarkerSize', 6, 'LineWidth', 1.5);

    title(all_labels{ch}, 'Interpreter','none', 'FontSize', 6);
    xlim([xVals(1) xVals(end)]);
    if mod(ch,8)==1
        ylabel('Power');
    end
    if ch > nChans-8
        xlabel('Frequency (Hz)');
    end
    set(gca,'FontSize',6); hold off;
end
sgtitle('1/f Removed Power Spectra with Peak Frequencies');

% Average across all channels
mean_flatten = mean_pf; % rename for clarity
overall_mean = mean(mean_flatten, 1);
overall_se   = std(mean_flatten, 0, 1) ./ sqrt(size(mean_flatten,1));

figure; hold on;
fill([freqs fliplr(freqs)], ...
     [overall_mean+overall_se fliplr(overall_mean-overall_se)], ...
     [0.8 0.8 0.8], 'EdgeColor','blue', 'FaceAlpha', 0.4);
plot(freqs, overall_mean, 'k', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Powe');
title('1/f Removed Spectrum (Mean ± SE across sessions and channels)');
set(gca,'FontSize',10); box on;

