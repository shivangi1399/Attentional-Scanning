% code to clean up data from artifacts, bad channels and eye-blinks.
% This is for Klecks data
% Elena's script here /mnt/hpc/projects/MWSampling/4Shivangi/code/Elena/run_cleanUp_data.m

% Description of clean-up procedure:
% ------------------------------
% Run this code after the data have been cut into trials

% Step A:
% plot and save all trials per channel (gives a good overview of problematic
% channels and trials with artifacts)

% Step B:
% Run semi-automated, visual artifact rejection
% Save ID of channels and trials that were rejected
% Data to clean:
% 1. LFP
% 2. MUA
% 3. Eye data (look for blinks)

% Step C:
% RE-plot the summary plots for all trials per channel (gives a good overview of
% problematic channels and trials with artifacts) -> find the remaining bad
% trials from summary plot and remove them for the respective sessions

% Step D:
% Load clean data Detect and remove eye blinks using ft_databrowser

clear all
close all
clc

%% create paths

animalName = 'klecks';

addpath /opt/fieldtrip_github/
ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf

% path to data
datafolder   = '/mnt/hpc/projects/MWSampling/4Shivangi/data_Klecks';

% path to save plots
figfolder = '/mnt/hpc/projects/MWSampling/4Shivangi/Plots/preprocessing/artifact_summary_klecks';

%% important variables

saveOutput = 0;
plottingInput  = 1;
plottingoutput = 1;

stepA = 1;
stepB = 0;
stepC = 0;
stepD = 0;

LFP = 0; % Clean up LFP data
MUA = 1; % Clean up MUA data
EYE = 0;

%% get file names and create paths to use

cd(datafolder),

temp = dir;
sesslist = [];
ii = 0;
for i = 1:length(temp)
    if strfind(temp(i).name,animalName)
        ii = ii+1;
        sesslist{ii,1} = temp(i).name;
    end
end

% paths to online preprocessed data
% sesspaths = cellfun(@(S) fullfile(datafolder, S,'RZ2'), sesslist, 'Uniform', 0);

% paths to RS4 preprocessed data
sesspaths = cellfun(@(S) fullfile(datafolder, S), sesslist, 'Uniform', 0);

%% Step A:
% Remove the same trials and channels in the mua data, that was removed in
% already cleaned lfp data. Plot and save all trials per channel
% This gives a good overview of problematic channels and trials with artifacts)

if stepA
    for isess = 1:length(sesspaths)
        
        disp(['reviewing session: ', sesslist{isess}])
        
        % create folder to save figures
        if ~isdir(fullfile(figfolder, sesslist{isess}))
            mkdir(fullfile(figfolder, sesslist{isess}))
        end
        
        try
            % Load data
            ESIload(fullfile(sesspaths{isess}, 'V4_mua_data.mat')); % mua data
            cd(sesspaths{isess});
            ESIload('lfpTrials_cleanf'); % cleaned lfp data
            
            %% plot original data
%             close all
%             FigH = figure('Position', get(0, 'Screensize'));
%             for ichan = 1:length(muaTrials.label)
%                 subplot(8,8,ichan)
%                 
%                 for itrial = 1:size(muaTrials.trialinfo,1)
%                     plot(muaTrials.time{1,1},muaTrials.trial{1,itrial}(ichan,:)),
%                     hold on
%                 end
%                 % title(['Chan ',lfpTrials.label{ichan}])
%                 title(num2str(ichan))
%                 xlim([min(muaTrials.time{1,1}),max(muaTrials.time{1,1})])
%             end
%             
%             % save figure
%             cd(fullfile(figfolder,sesslist{isess}))
%             set(FigH,'Units','Inches');
%             pos = get(FigH,'Position');
%             set(FigH,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%             print(FigH,[sesslist{isess},'_orgn_mua'],'-dpdf','-r0')
            
            %% Remove channels/trials in muaTrials based on NaNs in lfpTrials_cleanf
            
            nTrials = length(lfpTrials_cleanf.trial);
            for t = 1:nTrials
                nanCh = any(isnan(lfpTrials_cleanf.trial{t}), 2); % channels with NaNs in this trial
                
                if any(nanCh)
                    % Set corresponding channels to NaN or remove from muaTrials
                    muaTrials.trial{t}(nanCh, :) = NaN;
                end
            end
            muaTrials_clean = [];
            muaTrials_clean = muaTrials;
            
            cd(sesspaths{isess}),
            ESIsave('muaTrials_clean','muaTrials_clean')
                
            %% plot cleaned data
%             close all
%             FigH = figure('Position', get(0, 'Screensize'));
%             for ichan = 1:length(muaTrials.label)
%                 subplot(8,8,ichan)
%                 
%                 for itrial = 1:size(muaTrials.trialinfo,1)
%                     plot(muaTrials.time{1,1},muaTrials.trial{1,itrial}(ichan,:)),
%                     hold on
%                 end
%                 % title(['Chan ',lfpTrials.label{ichan}])
%                 title(num2str(ichan))
%                 xlim([min(muaTrials.time{1,1}),max(muaTrials.time{1,1})])
%             end
%             
%             % save figure
%             cd(fullfile(figfolder,sesslist{isess}))
%             set(FigH,'Units','Inches');
%             pos = get(FigH,'Position');
%             set(FigH,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%             print(FigH,[sesslist{isess},'_1cn_mua'],'-dpdf','-r0')
            
        catch
            continue
        end
    end
end

%% Step B:
% Run semi-automated, visual artifact rejection
% Save ID of channels and trials that were rejected

if stepB
    
    %% Clean up LFP data
    if LFP
        for isess = 1:length(sesspaths)
            
            disp(['reviewing session: ',sesslist{isess}])
            
            % load data into trials
            
            try
                %% LFP data
                cd(sesspaths{isess}),
                ESIload(fullfile(sesspaths{isess},'V4_lfp_data.mat'));
                
                cfg = [];
                cfg.method = 'summary';
                cfg.keepchannel = 'nan';
                cfg.keeptrial   = 'nan';
                cfg.latency     = [-1 0.4];
                lfpTrials_cleanf = ft_rejectvisual(cfg,lfpTrials);
                
                cd(sesspaths{isess}),
                ESIsave('lfpTrials_cleanf','lfpTrials_cleanf')
                
                clear lfpTrials_cleanf lfpTrials
                close all
                
            catch
                continue
            end
            
        end
        
        %% plot LFP summaries after artifact rejection
        for isess = [6,23] %1:length(sesspaths)
            
            disp(['reviewing session: ',sesslist{isess}])
            
            % load data into trials
            
            try
                
                % create folder to save figures
                if ~isdir(fullfile(figfolder,sesslist{isess}))
                    mkdir(fullfile(figfolder,sesslist{isess}))
                end
                
                %% LFP data
                ESIload(fullfile(sesspaths{isess},'lfpTrials_cleanf.mat'));
                
                close all
                FigH = figure('Position', get(0, 'Screensize'));
                for ichan = 1:length(lfpTrials_cleanf.label)
                    subplot(8,8,ichan)
                    
                    for itrial = 1:size(lfpTrials_cleanf.trialinfo,1)
                        plot(lfpTrials_cleanf.time{1,1},lfpTrials_cleanf.trial{1,itrial}(ichan,:)),
                        hold on
                    end
                    % title(['Chan ',lfpTrials.label{ichan}])
                    title(num2str(ichan))
                    xlim([min(lfpTrials_cleanf.time{1,1}),max(lfpTrials_cleanf.time{1,1})])
                end
                
                
                % save figure
                cd(fullfile(figfolder,sesslist{isess}))
                set(FigH,'Units','Inches');
                pos = get(FigH,'Position');
                set(FigH,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
                print(FigH,[sesslist{isess},'_cleaned_lfp'],'-dpdf','-r0')
                
                clear lfpTrials_cleanf
                
                
            end
            
        end
    end
end

%% Step C:
% plot LFP summaries after artifact rejection - so we can remove the remaining bad trials

if stepC
    if LFP
        for isess = [4,6,10,16,23] %1:length(sesspaths)
            
            disp(['reviewing session: ',sesslist{isess}])
            
            % load data into trials
            try
                
                % create folder to save figures
                if ~isdir(fullfile(figfolder,sesslist{isess}))
                    mkdir(fullfile(figfolder,sesslist{isess}))
                end
                
                %% LFP data
                ESIload(fullfile(sesspaths{isess},'lfpTrials_cleanf.mat'));
                
                close all
                FigH = figure('Position', get(0, 'Screensize'));
                for ichan = 60:63%length(lfpTrials_cleanf.label)
                    subplot(2,2,ichan)
                    
                    for itrial = 1:size(lfpTrials_cleanf.trialinfo,1)
                        plot(lfpTrials_cleanf.time{1,1},lfpTrials_cleanf.trial{1,itrial}(ichan,:)),
                        hold on
                        end_time = lfpTrials_cleanf.time{1,1}(end);
                        end_val = lfpTrials_cleanf.trial{1, itrial}(ichan,end);
                        text(end_time, end_val, num2str(itrial), 'FontSize', 6)
                    end
                    % title(['Chan ',lfpTrials.label{ichan}])
                    title(num2str(ichan))
                    xlim([min(lfpTrials_cleanf.time{1,1}),max(lfpTrials_cleanf.time{1,1})])
                end
                
                
                % save figure
                cd(fullfile(figfolder,sesslist{isess}))
                set(FigH,'Units','Inches');
                pos = get(FigH,'Position');
                set(FigH,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
                print(FigH,[sesslist{isess},'_trial_cleaning'],'-dpdf','-r0')
                
                clear lfpTrials_cleanf
                
                
            end
        end
        
        %% removing the remaining bad trials
        
        for isess = [6,23]
            
            disp(['reviewing session: ',sesslist{isess}])
            
            try
                %% LFP data
                cd(sesspaths{isess}),
                ESIload(fullfile(sesspaths{isess},'lfpTrials_cleanf'));
                lfpt = lfpTrials_cleanf;
                
                cfg = [];
                cfg.method = 'summary';
                cfg.keepchannel = 'nan';
                cfg.keeptrial   = 'nan';
                cfg.latency     = [-1 0.4];
                lfpTrials_cleanf = ft_rejectvisual(cfg,lfpt);
                
                cd(sesspaths{isess}),
                ESIsave('lfpTrials_cleanf','lfpTrials_cleanf')
                
                clear lfpTrials_cleanf lfpt
                close all
                
            catch
                continue
            end
        end
    end
end

% rerun the summaries to check artifact rejection

%% Step D:
% Run semi-automated, visual artifact rejection
% using ft_databrowser

if stepD
    
    %% Clean up LFP data
    if LFP
        for isess = [6,23] %1:length(sesspaths)
            
            disp(['reviewing session: ',sesslist{isess}])
            
            % load data into trials
            
            try
                %% LFP data
                cd(sesspaths{isess}),
                ESIload(fullfile(sesspaths{isess},'lfpTrials_cleanf.mat'));
                lfpTrials = lfpTrials_cleanf;
                cfg = [];
                %cfg.method = 'summary';
                %cfg.keepchannel = 'nan';
                %cfg.keeptrial   = 'nan';
                %cfg.latency     = [-1 0.4];
                lfpTrials_cleanf = ft_databrowser(cfg,lfpTrials);
                
                close all
                
            catch
                continue
            end
            
        end
    end
end
