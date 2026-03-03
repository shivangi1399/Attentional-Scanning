clear all
close all
clc

% Code description:
% -----------------
% Power spectrum analysis across all channels and sessions using FieldTrip
% Checking significance or correlation per location based on the different frequency windows

%% Specify paths

addpath /opt/fieldtrip_github/
ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/eye_data
clc

%% Paths for data

datafolder   = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks';

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

session_paths = [];
session_paths = cellfun(@(x) fullfile(datafolder,x), session_names, 'uniform',0);

session_paths_files = [];
session_paths_files = cellfun(@(x) fullfile(datafolder,x, 'clean_data.mat'), session_names, 'uniform',0);

output_paths = cellfun(@(x) fullfile(datafolder, x,'ERP'),session_names, 'uniform',0);

%% Defining frequency bands - change once the method is defined

all_pow = containers.Map;
all_labels = {};

for isess = 1:length(session_paths_files)

    ESIload(session_paths_files{isess});  
    load(fullfile(session_paths{isess}, 'zlfptrials.mat')); 

    % Trial selection
    hitIdx  = find(clean_data.trialinfo(:,20) == 1);
    missIdx = find(clean_data.trialinfo(:,20) == 5);
    
    cfg=[];
    cfg.toilim    = [-1 0];
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
for i = 1:nChans
    chName = all_labels{i};
    mean_flatten(i,:) = mean(all_pow(chName),1);
end

% Plot
peakFreqs = zeros(nChans, 1);  
figure;
for ch = 1:nChans
    subplot(8,8,ch);
    plot(freqs(12:end), mean_flatten(ch,12:end), 'k', 'LineWidth', 1.5);
    hold on;
    
    % Find peak frequency and power for this channel
    [peakPower, peakIdx] = max(mean_flatten(ch,12:end));
    peakFreq = freqs(peakIdx);
    
    % Store peak frequency
    peakFreqs(ch) = peakFreq;
    
    % Mark peak frequency with a red circle
    plot(peakFreq, peakPower, 'ro', 'MarkerSize', 6, 'LineWidth', 1.5);
    
    title(all_labels{ch}, 'Interpreter','none', 'FontSize', 6);
    xlim([freqs(12) freqs(end)]);
    
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






