% Description:
% -------------
% Goal: combine phases from all sessions that have been estimated before
% ----
% if concatenate_phases:
% Load the pre-computed spectra per trial and session and save in a common
% structure that also includes the trialinfo
% -------------
% if concatenate_erp_amp:
% estimate erp amplitude and save in a common structure with the
% phase
% -------------
% if concatenate_RT:
% load pre estimated response time and save in a common structure with the
% phase

clear all
close all
clc

%% Dependecies

addpath /opt/fieldtrip_github/
ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/
clc

%% logicals

concatenate_phases = 1;
concatenate_erp_amp = 1;
concatenate_RT     = 0;

%% Create data paths
animalName = 'hermes';  % Change this to switch animal (e.g. 'klecks')

datafolder   = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName]);

cd(datafolder),
temp = dir;
session_names = [];
ii = 0;
for i = 1:length(temp)
    if strfind(temp(i).name,animalName)
        ii = ii+1;
        session_names{ii,1} = temp(i).name;
    end
end

session_paths_files = [];
session_paths_files = cellfun(@(x) fullfile(datafolder,x, 'clean_data.mat'), session_names, 'uniform',0);

phase_paths = cellfun(@(x) fullfile(datafolder, x,'Phase_analysis/hit_miss'),session_names, 'uniform',0);
RT_paths = cellfun(@(x) fullfile(datafolder, x,'RT'),session_names, 'uniform',0);
phase_folder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName], 'phase_coherence', 'cp_plus_100');
% cp_till/plus_100 shows the erp range used, ph55_till_100 means that
% phase est was done at 55 and the 10 or 25 after cp shows how many ms after
% the defelction point the signal was cut
%% load channel per session matrix

cd(datafolder)
load('V4_lfp_data.mat')

all_channels = lfpTrials.label(:);
num_channels = length(all_channels);
num_sessions = length(phase_paths);

channels_sessions = NaN(num_channels, num_sessions);

for isess = 1:num_sessions
    ESIload(session_paths_files{isess});
    
    for ichan = 1:num_channels
        if any(strcmp(all_channels{ichan}, clean_data.label))
            channels_sessions(ichan, isess) = 1;
        end
    end
end

cd(datafolder)
save channels_sessions channels_sessions

%% load critical time per channel matrix

cd(fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animalName], 'critical_time'))
load('CriticalTime.mat')
load('all_channels.mat')

cd(datafolder)
load('V4_lfp_data.mat')
full_channels = lfpTrials.label(:);

CriticalTime_mat = NaN(length(full_channels), 1);

for i = 1:length(all_channels)
    idx = find(strcmp(full_channels, all_channels{i}));
    if ~isempty(idx)
        CriticalTime_mat(idx) = CriticalTime(i);
    end
end

%% combine phases across session along with erp amplitude


cd(datafolder)
load('channels_sessions.mat')
chan_orig = 1:64;

if concatenate_phases
    
    ph_all = []; trlinfo_all = [];
    session_name_all = {};
    ERP_ampl_all = []; ERP_trialinfo = [];
    RT_time  = [];  RT_sample = [];
    
    % loop through sessions and concatenate the phases at cp over all
    % sessions
    
    for isess = 1:length(phase_paths)
        isess
        
        ESIload(session_paths_files{isess}) % load broadband data for ERP-amplitude estimation
        
        %% load phase estimations
        
        cd(fullfile(phase_paths{isess}, '100iter_cut@cp/phase'));
        
        % list of available channels
        dr = dir;
        dr = dr(3:end);
        chan_nums = sort(str2double({dr.name}));
        
        % Preallocate for all channels
        cd(fullfile(phase_paths{isess}, '100iter_cut@cp/phase', num2str(chan_nums(1))));
        load('phase.mat'); 
        ph = nan([size(phase.ar_phase), length(chan_orig)]); % size: (samples × trials × 64)
        
        % Loop only over active channels
        active_ch = find(channels_sessions(:, isess) == 1);
        for idx = 1:length(active_ch)
            ichan = active_ch(idx);
            cd(fullfile(phase_paths{isess}, '100iter_cut@cp/phase', num2str(chan_nums(idx))));
            load('phase.mat');
            ph(:, :, ichan) = phase.ar_phase;
        end
        
        ph_all = [ph_all; ph];
        trlinfo_all = [trlinfo_all; phase.trialinfo];
        session_name_all = [session_name_all; repmat({session_names{isess}}, size(phase.trialinfo, 1), 1)];

        %% ERP amplitude
        
        if concatenate_erp_amp
            [~, ia_phase, ib_clean] = intersect(phase.trialinfo(:,14), clean_data.trialinfo(:,14));
            nTrials = length(ia_phase);
            
            ERP_ampl = nan(length(chan_orig), nTrials);
           
            active_ch = find(channels_sessions(:, isess) == 1);
            for ichan = 1:length(active_ch')
                ct = CriticalTime_mat(active_ch(ichan))-0.01;
                if isnan(ct), continue; end
                
                cfg = [];
                cfg.channel = ichan;
                cfg.latency = [ct ct+0.1];
                data_chan = ft_selectdata(cfg, clean_data);
                
                % ERP amplitude (RMS)
                ampl = cellfun(@(x) sqrt(mean((x(:) - mean(x(:))).^2)), data_chan.trial);
                
                % Align trials
                [~, ia_chan, ib_chan] = intersect(phase.trialinfo(:,14), data_chan.trialinfo(:,14));
                
                % Store
                ch_data = nan(1, nTrials);
                ch_data(ia_chan) = ampl(ib_chan);
                ERP_ampl(active_ch(ichan),:) = ch_data;
            end
            
            ERP_ampl_all  = [ERP_ampl_all; ERP_ampl'];
            ERP_trialinfo = [ERP_trialinfo; clean_data.trialinfo(ib_clean,:)];
        end
         
        %% RT
        if concatenate_RT
            load(fullfile(RT_paths{isess},'responsetime.mat'))
            
            % find trials in common
            idxA = []; idxB = [];
            [~,idxA,idxB] = intersect(transf.trialinfo(:,14),responsetime.trialinfo(:,14),'stable');
            RT_time   = [RT_time;responsetime.time(idxB,1)];
            RT_sample = [RT_sample;responsetime.sample(idxB,1)];
        end
        
    end
    
    ph_comb =[];
    ph_comb.phase_all = ph_all;
    ph_comb.trialinfo = trlinfo_all;
    ph_comb.session_names = session_name_all;
    ph_comb.dimord    = 'trlx freq x chan';
    
    if concatenate_erp_amp
        ph_comb.ERP_ampl_all  = ERP_ampl_all;
        ph_comb.ERP_trialinfo = ERP_trialinfo;
    end
    
    if concatenate_RT
        ph_comb.RT_time   = RT_time;
        ph_comb.RT_sample = RT_sample;
    end
    
    if ~isdir(phase_folder)
        mkdir(fullfile(phase_folder))
    end
    
    cd(phase_folder)
    save('ph_all_sess','ph_comb')
    clear ph_all_sess ph_comb
    
end

