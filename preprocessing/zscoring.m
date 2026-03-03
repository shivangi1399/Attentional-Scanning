clear all
close all
clc

% Code description:
% -----------------
% zscoring artifect rejected data

%% paths

addpath /opt/fieldtrip_github/
ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
clc

%% paths

datafolder   = '/mnt/hpc/projects/MWSampling/4Shivangi/data_Hermes';

cd(datafolder),
animalName = 'hermes';
temp = dir;
session_names = [];
ii = 0;
for i = 1:length(temp)
    if strfind(temp(i).name,animalName)
        ii = ii+1;
        session_names{ii,1} = temp(i).name;
    end
end

% paths to RS4 preprocessed data
data_path = cellfun(@(S) fullfile(datafolder, S), session_names, 'Uniform', 0);

% Path to cleaned data
output_folder = '/mnt/hpc/projects/MWSampling/4Shivangi/results_hermes';
output_path = cellfun(@(x) fullfile(output_folder, x),session_names, 'uniform',0);

%% use artifact rejected data and zscore sessions

for isess = 1:length(session_names)
    
    % load LFP data
    cd(data_path{isess})
    
    % check if file exists
    if ~isfile('lfpTrials_clean.mat')
        fprintf('Skipping %s: lfpTrials_clean.mat not found\n', session_names{isess});
        continue
    end
    
    ESIload('lfpTrials_clean.mat');
    
    % create session folder
    if ~isdir(fullfile(output_path{isess}))
        mkdir(fullfile(output_path{isess}))
    end
    
    channels = 1:64;
    trials = 1:length(lfpTrials_clean.trial);
    A = cellfun(@(x) isnan(x),lfpTrials_clean.trial,'UniformOutput',false);
    remove_channels = cellfun(@(x) find(x(:,1)==1),A,'UniformOutput',false);
    B = cellfun(@(x) size(x,1),remove_channels,'UniformOutput',false);
    remove_trials = find(cell2mat(B)==64);
    
    % selecting clean data
    cfg = [];
    cfg.channel = setdiff(channels,remove_channels{1,1});
    cfg.trials  = setdiff(trials,remove_trials);
    clean_data  = ft_selectdata(cfg,lfpTrials_clean);
    cd(fullfile(output_path{isess}))
    save('clean_data','clean_data')
    
    % zscore data in order to be able to pool over sessions
    disp(strcat('session- ',num2str(isess),' out of- ',num2str(length(session_names)), ', running LFP z-scoring')) % #ok<UNRCH>
    zlfptrials = fun_zscore_session(clean_data);
    
    % save data
    cd(fullfile(output_path{isess}))
    save('zlfptrials','zlfptrials')
    
    clear lfpTrials_clean
    
end