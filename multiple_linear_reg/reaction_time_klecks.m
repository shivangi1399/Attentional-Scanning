% Description:
% Clean the eye data like lfp data and compute reaction time for each of
% the sessions. The reaction time calculated here is the saccade initiation
% time and done in Zhang2025
clearvars
close all
clc

%% Dependencies
addpath /opt/fieldtrip_github/
ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
clc

%% Paths
datafolder = '/mnt/hpc/projects/MWSampling/4Shivangi/data_Klecks';
resultfolder   = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks';

cd(resultfolder)
animalName = 'klecks';
temp = dir;

session_names = {};
ii = 0;
for i = 1:length(temp)
    if contains(temp(i).name, animalName)
        ii = ii + 1;
        session_names{ii,1} = temp(i).name; 
    end
end

eye_paths = cellfun(@(x) ...
    fullfile(datafolder, x, 'eye_data.mat'), ...
    session_names, 'uniformoutput', false);
data_path = cellfun(@(x) fullfile(datafolder, x), session_names, 'Uniform', 0);
session_paths_RT = cellfun(@(x) fullfile(resultfolder,x), session_names, 'uniform',0);

%% Saccade detection config

cfg_detectSac.method = 'Engbert2003';
cfg_detectSac.params = {6};  % lambda scaling factor for threshold
cfg_selSac = struct();
cfg_selSac.starttime_range  = [0.1 0.5];    % seconds after target onset
cfg_selSac.magnitude_range  = [1 10];       % degrees

%% Selecting trials and channels corresponding to neural data

for isess = 1:length(session_names)
    
    % load LFP data
    cd(data_path{isess})
    
    % check if file exists
    if ~isfile('lfpTrials_cleanf.mat')
        fprintf('Skipping %s: lfpTrials_cleanf.mat not found\n', session_names{isess});
        continue
    end
    
    ESIload('lfpTrials_cleanf.mat');
    
    trials = 1:length(lfpTrials_cleanf.trial);
    A = cellfun(@(x) isnan(x),lfpTrials_cleanf.trial,'UniformOutput',false);
    remove_channels = cellfun(@(x) find(x(:,1)==1),A,'UniformOutput',false);
    B = cellfun(@(x) size(x,1),remove_channels,'UniformOutput',false);
    remove_trials = find(cell2mat(B)==64);

    ESIload(eye_paths{isess})
    cfg = [];
    cfg.trials  = setdiff(trials,remove_trials);
    eyeData  = ft_selectdata(cfg,eyeTrials);

    eye_file = fullfile(session_paths_RT{isess}, 'eyeData.mat');
    save(eye_file, 'eyeData'); %cleaned eye data

%% Compute RTs

    cfg = [];
    cfg.eye_file      = eye_file;
    cfg.cfg_detectSac = cfg_detectSac;
    cfg.cfg_selSac    = cfg_selSac;

    RT_sess = compute_RT_eye(cfg);
    cd(fullfile(session_paths_RT{isess}));
    save('RT_sess','RT_sess')

end