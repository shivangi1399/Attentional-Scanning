clear all
close all
clc

%% Specify paths

addpath /opt/fieldtrip_github/
ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/
clc

%% Define important variables

zscore_run   = 1;     % if first_run, z-score session (over session) and save this z-scored 
                      % data should be later on used for averaging across sessions
run_iter     = 0;
find_order   = 1;
iter_n       = 100;   % number of iterations to run
t            = -0.5:0.25:0; % the critical time is the new zero
f            = logspace(log10(2), log10(80), 40);
exms         = 800;
arord        = 51;

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
session_paths_files = cellfun(@(x) fullfile(datafolder,x, 'clean_lfp.mat'), session_names, 'uniform',0);

output_paths = cellfun(@(x) fullfile(datafolder, x,'Phase_analysis/hit_miss'),session_names, 'uniform',0);

L = linspace(-1,0.4,1401);

%% load data per session and redefine it

% Steps:
%--------
% 1) load session
% 2) redefine data such that the critical point in the new zero (this is done channel wise)
% 2) zscore  and detrend data
% 3) extrapolate data using AR for iter_n
% 4) phase analysis using fourier

if run_iter
    
    if zscore_run
        for isess = 1:length(session_names)
            
            % load LFP data
            ESIload(session_paths_files{isess});
            
            % create session folder to save output if it does not exist already
            if ~isdir(fullfile(output_paths{isess},'100iter_cut@30'))
                mkdir(fullfile(output_paths{isess},'100iter_cut@30'))
            end
            
            % redefine trials using channel specific critical time
            cd(fullfile(datafolder, 'critical_time'))
            load CriticalTime.mat
            load all_channels.mat
            
            CriticalTime_n = 0.040*ones(1,length(CriticalTime)); %min(CriticalTime)*ones(1,length(CriticalTime));
            
            ch_cp_dat = []; % data structure for data aligned to critical time
            ch_cp_dat.time_actual = [];
            ch_cp_dat.sampleinfo_actual = [];
            ch_cp_dat.sampleinfo_actual = cell(1, length(clean_data.label));
            ch_cp_dat.label = clean_data.label;
            nTrials = sum(clean_data.trialinfo(:,20) == 1 | clean_data.trialinfo(:,20) == 5);
            ch_cp_dat.trial = cell(1, nTrials);

            for ichan = 1:length(clean_data.label)
                label = clean_data.label{ichan}; % current channel
                idx = find(strcmp(all_channels, label)); % match in all_channels
                    
                if ~isempty(idx)
                    cp = CriticalTime_n(idx)-0.01; % get corresponding critical time
                    
                    cfg = []; % select channel specific data
                    cfg.channel = label;
                    ch_dat = ft_selectdata(cfg, clean_data);
                    
                    cfg = []; % cut the trials at cp such that there is one sec of data before the cut
                    cfg.toilim = [(-1+cp)  cp];
                    cfg.trials = find(clean_data.trialinfo(:,20) == 1 | clean_data.trialinfo(:,20) == 5);
                    dat = ft_redefinetrial(cfg, ch_dat);
                else
                    warning('Channel %s not found in all_channels.', label);
                end
                ch_cp_dat.sampleinfo_actual{ichan} = dat.sampleinfo;
                ch_cp_dat.time_actual(ichan,:) = dat.time{1, 1};
                
                for itrial = 1:length(dat.trial)
                    if ichan == 1
                        ch_cp_dat.trial{itrial} = zeros(length(clean_data.label), length(dat.time{itrial}));
                    end
                    ch_cp_dat.trial{itrial}(ichan, :) = dat.trial{itrial};
                end
                
            end
            
            ch_cp_dat.fsample = dat.fsample;
            ch_cp_dat.hdr = dat.hdr;
            ch_cp_dat.sampleinfo = dat.sampleinfo;
            ch_cp_dat.cfg = dat.cfg;
            ch_cp_dat.trialinfo = dat.trialinfo;
            ch_cp_dat.time = repmat({linspace(-1, 0, 1001)}, 1, numel(ch_cp_dat.trial)); % 0 is the new critical time or point
            
            % save the data with the actual times and sampleinfo
            cd(fullfile(output_paths{isess},'100iter_cut@30'))
            save('ch_cp_dat','ch_cp_dat')
            
            % zscore data in order to be able to pool over sessions
            fieldsToRemove = {'time_actual', 'sampleinfo_actual'};
            existingFields = isfield(ch_cp_dat, fieldsToRemove);
            
            if any(existingFields)
                ch_cp_datz = rmfield(ch_cp_dat, fieldsToRemove(existingFields));
            end
            
            disp(strcat('session- ',num2str(isess),' out of- ',num2str(length(session_names)), ', running LFP z-scoring')) % #ok<UNRCH>
            zdat = fun_zscore_session(ch_cp_datz);

            % detrend data
            cfg = [];
            cfg.detrend = 'yes';
            zdat = ft_preprocessing(cfg, zdat);
            
            % save data
            cd(fullfile(output_paths{isess},'100iter_cut@30'))
            save('zdat','zdat')
            
            % Find correct model order
            if find_order
                arordmin = 1;
                arordmax = 300;
                
                for ichan = 1:length(zdat.label)
                    
                    data = cellfun(@(x) x(ichan,:),zdat.trial,'UniformOutput',false);
                    
                    [Aest,Cest,SBC,FPE] = subfunc_ARord(data,arordmin,arordmax);
                    
                    min_sbc = min(SBC);
                    order_sbc(isess,ichan) = find(min_sbc==SBC);
                    order_chan(1,ichan)=find(min_sbc==SBC);
                end
                cd(fullfile(output_paths{isess},'100iter_cut@30'))
                save order_chan order_chan
            end
        end
    end
     
    %% Find model order for the sessions
    
    if find_order
        corder_all = [];  
        for isess = 1:length(session_names)
            cd(fullfile(output_paths{isess}, '100iter_cut@30'))
            load order_chan  
            corder_all = [corder_all order_chan(:)'];
        end
        arord = mean(corder_all,'all');
    end
    
    %% AR extrapolation and fourier
    
    for isess = 1%:length(session_names)
        isess
        cfg = [];
        cfg = cell(1,iter_n);
        
        for iter = 1:iter_n
            
            if ~isdir(fullfile(output_paths{isess},'100iter_cut@30',num2str(iter)))
                mkdir(fullfile(output_paths{isess},'100iter_cut@30',num2str(iter)))
            end

            cfg{iter}.inputfile   = fullfile(output_paths{isess},'100iter_cut@30','zdat.mat');
            cfg{iter}.exms        = exms;
            cfg{iter}.arord       = arord;
            cfg{iter}.cycles      = 3;
            cfg{iter}.toi         = t;
            cfg{iter}.foi         = f;
            cfg{iter}.fsample     = 1000; % Sampling rate of original data
            cfg{iter}.nfsample    = 1000; % Sampling rate after downsampling
            cfg{iter}.outputfile  = fullfile(output_paths{isess},'100iter_cut@30',num2str(iter));

        end


        slurmfun(@AR_fourier, cfg, ...
            'partition',     '8GB', ...
            'stopOnError',   false,  ...
            'useUserPath',   true    );
        %         cellfun(@(x) AR_fourier(x),cfg,'UniformOutput',false);

    end
    pause(400)
end

%% Phase estimation using transform

for isess= 1:length(session_names)
    isess
    ESIload(session_paths_files{isess});
    
    cfg = [];
    cfg = cell(1,length(clean_data.label));
    
    for ichan = 1:length(clean_data.label)
        
        
        if ~isdir(fullfile(output_paths{isess},'100iter_cut@30','phase',num2str(ichan)))
            mkdir(fullfile(output_paths{isess},'100iter_cut@30','phase',num2str(ichan)))
        end
        
        cfg{ichan}.inputfile   = fullfile(output_paths{isess},'100iter_cut@30');
        cfg{ichan}.ichan       = ichan;
        cfg{ichan}.iter_n      = 100;
        cfg{ichan}.toi         = 3;
        cfg{ichan}.outputfile  = fullfile(output_paths{isess},'100iter_cut@30','phase',num2str(ichan));
        
    end
    
    slurmfun(@phase_transform, cfg, ...
        'partition',     '8GB', ...
        'stopOnError',   false,  ...
        'useUserPath',   true    );
    
    %         cellfun(@(x) phase_transform(x), cfg,'UniformOutput',false);

end

%% check upto which session has the estimation run

iter = 100;
for isess = 1:length(session_names)
    isess
    cd(fullfile(output_paths{isess},'100iter_cut@30',num2str(iter)))
    load freqpow.mat 
    disp(['channel Num = ', num2str(length(freqpow.label))])
end




