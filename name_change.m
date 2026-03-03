clear all
close all
clc

%% Add necessary paths
addpath /opt/fieldtrip_github/
ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/eye_data

%% Define root data directory
datafolder = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks';
animalName = 'klecks';

% Get session folders
cd(datafolder)
temp = dir;
session_names = {};
ii = 0;
for i = 1:length(temp)
    if contains(temp(i).name, animalName)
        ii = ii + 1;
        session_names{ii,1} = temp(i).name;
    end
end

% Build full paths
data_paths = cellfun(@(x) fullfile(datafolder, x), session_names, 'UniformOutput', false);

%% Rename folder inside each session
for i = 1:length(data_paths)
    phase_path = fullfile(data_paths{i}) %, 'Phase_analysis', 'hit_miss');
    target_folder = fullfile(phase_path, 'ERP');
    new_folder = fullfile(phase_path, 'ERP_LFP');

    if exist(target_folder, 'dir')
        fprintf('Renaming in session %s...\n', session_names{i});
        try
            movefile(target_folder, new_folder);
        catch ME
            warning('Failed to rename folder in %s: %s', session_names{i}, ME.message);
        end
    else
        fprintf('No target folder in session %s\n', session_names{i});
    end
end

disp('Folder renaming complete.');

%% Rename files inside each session

for i = 1:length(data_paths)
    phase_path = fullfile(data_paths{i}); %, 'Phase_analysis', 'hit_miss');
    old_mat_file = fullfile(phase_path, 'clean_data.mat');
    new_mat_file = fullfile(phase_path, 'clean_lfp.mat');

    if exist(old_mat_file, 'file')
        fprintf('Renaming .mat file in session %s...\n', session_names{i});
        try
            movefile(old_mat_file, new_mat_file);
        catch ME
            warning('Failed to rename .mat file in %s: %s', session_names{i}, ME.message);
        end
    else
        fprintf('No clean_data.mat file in session %s\n', session_names{i});
    end
end

disp('File renaming complete.');


