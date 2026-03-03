% Description:
% Clean pupil data like lfp and mua data

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

pup_paths = cellfun(@(x) ...
    fullfile(datafolder, x, 'pup_data.mat'), ...
    session_names, 'uniformoutput', false);
data_path = cellfun(@(x) fullfile(datafolder, x), session_names, 'Uniform', 0);
session_paths_RT = cellfun(@(x) fullfile(resultfolder,x), session_names, 'uniform',0);

%% Selecting trials and channels corresponding to neural data

for isess = 1:length(session_names)
    isess
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

    ESIload(pup_paths{isess})
    cfg = [];
    cfg.trials  = setdiff(trials,remove_trials);
    pupData  = ft_selectdata(cfg,pupTrials);

    eye_file = fullfile(session_paths_RT{isess}, 'pupData.mat');
    save(eye_file, 'pupData'); %clean pupil data

end