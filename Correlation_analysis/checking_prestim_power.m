%% To check power spectra of pre-stimulus data that was used to estimate phase

% Since pre-stimulus data is used to calculate the phase, we need to check if the
% data has enough power in frequencies which later show us the significant correlations

% Find freq power spectra and normalize by baseline trial data otherwise we will have
% very high 1/f noise and that will make detection of power difficult

clear all
close all
clc
%% paths

addpath /opt/fieldtrip_github/
ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf

% path to temporary data
datafolder = '/mnt/hpc/projects/MWSampling/4Shivangi/data'; % in windows:\\cs\projects etc

session_names = [];
session_names = {'klecks_20170804_attentional-sampling_1';...
    'klecks_20170807_attentional-sampling_1';...
    'klecks_20170808_attentional-sampling_1';...
    'klecks_20170810_attentional-sampling_1';...
    'klecks_20170817_attentional-sampling_1';...
    'klecks_20170818_attentional-sampling_1';...
    'klecks_20170821_attentional-sampling_1';...
    'klecks_20170822_attentional-sampling_1';...
    'klecks_20170823_attentional-sampling_1';...
    'klecks_20170824_attentional-sampling_1';...
    'klecks_20170825_attentional-sampling_1';...
    'klecks_20170828_attentional-sampling_1';...
    'klecks_20170830_attentional-sampling_1';...
    'klecks_20170831_attentional-sampling_1';...
    'klecks_20170901_attentional-sampling_1';...
    'klecks_20170904_attentional-sampling_1';...
    'klecks_20170906_attentional-sampling_1';...
    'klecks_20170908_attentional-sampling_1';...
    'klecks_20170911_attentional-sampling_1';...
    'klecks_20170913_attentional-sampling_1';...
    'klecks_20170914_attentional-sampling_1';...
    'klecks_20170915_attentional-sampling_1';...
    'klecks_20170919_attentional-sampling_1';...
    'klecks_20171020_attentional-sampling_1'};

% create cell paths to use with slurm
session_paths = [];
session_paths = cellfun(@(x) fullfile(datafolder,x), session_names, 'uniform',0);

%%  

PF_h = nan(length(session_names),64,30);
PF_m = nan(length(session_names),64,30);

for isess = 1:length(session_names)
    cd(session_paths{isess})
    
    %% baseline trials
    
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
    cfg.toilim    = [1  2];
    datb = ft_redefinetrial(cfg, zlfpTrials_b);
    
    %plot(dat.time{1, 1},dat.trial{1, 1}(1,:))
    
    % power spectra
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmfft';
    cfg.taper        = 'hanning';
    cfg.foi          = 2:2:60;
    freqpow_b = ft_freqanalysis(cfg,datb);
    
    % for ichan = 1:64
    %     subplot(8,8,ichan)
    %     plot(freqpow_Hit_b.freq,freqpow_Hit_b.powspctrm(ichan,:),'k')
    % end
    
    %%  power in the baseline period of the target trials
    
    % load LFP data
    ESIload('V4_lfp_data.mat')
    
    % zscoring
    disp(strcat('session- ',num2str(isess),' out of- ',num2str(length(session_names)), ', running z-scoring for target trials')) %#ok<UNRCH>
    zlfpTrials = fun_zscore_session(lfpTrials);
    
    % redefine trial
    cfg=[];
    cfg.toilim    = [-1 0];
    dat = ft_redefinetrial(cfg, zlfpTrials);
    
    % detrending
    cfg=[];
    cfg.detrend = 'yes';
    cfg.demean = 'yes';
    dat = ft_preprocessing(cfg, dat);
    
    hits = []; misses = [];
    hits= find(zlfpTrials.trialinfo(:,20)==1);
    misses = find(zlfpTrials.trialinfo(:,20)==5);
    
    trial_no = min(length(hits),length(misses));
    
    % hits
    hit = randperm(max(length(hits),length(misses)),trial_no);
    cfg = [];
    cfg.trials = hits(hit)';
    H_dat = ft_selectdata(cfg,dat);
    
    % misses
    cfg = [];
    cfg.trials = misses';
    M_dat = ft_selectdata(cfg,dat);
    
    
    % spectral analysis
    
    %%% hits
    
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmfft';
    cfg.taper        = 'hanning';
    cfg.foi          = 2:2:60;      % analysis 2 to 60 Hz in steps of 2 Hz
    freqpow_Hit = ft_freqanalysis(cfg,H_dat);
    
    
    %%% misses
    
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmfft';
    cfg.taper        = 'hanning';
    cfg.foi          = 2:2:60;      % analysis 2 to 60 Hz in steps of 2 Hz
    freqpow_Miss = ft_freqanalysis(cfg, M_dat);
    
    % plot
    
    %     for ichan = 1:64
    %         subplot(8,8,ichan)
    %         plot(freqpow_Hit.freq,freqpow_Hit.powspctrm(ichan,:),'k')
    %         plot(freqpow_Hit.freq,freqpow_Miss.powspctrm(ichan,:),'r')
    %     end
    %
    
    PF_h(isess,:,:) = (freqpow_Hit.powspctrm-freqpow_b.powspctrm);
    PF_m(isess,:,:) = (freqpow_Miss.powspctrm-freqpow_b.powspctrm);
    
end

%% plot - corrected by baseline trials

PF_avg_h = reshape(mean(PF_h,1),64,30);
PF_avg_m = reshape(mean(PF_m,1),64,30);

for ichan = 1:64
    subplot(8,8,ichan)
    plot(freqpow_Hit.freq,PF_avg_h(ichan,:),'k')
    hold on
    plot(freqpow_Hit.freq,PF_avg_m(ichan,:),'r')
end

PF_avgc_h = mean(PF_avg_h,1);
PF_avgc_m = mean(PF_avg_m,1);

figure;
plot(freqpow_Hit.freq,PF_avgc_h,'k')
hold on
plot(freqpow_Hit.freq,PF_avgc_m,'r')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% eye data - to confirm the time window in baseline trial

ESIload('eye_data_noTarg.mat')

% keep only correct trials and baseline condition (cond 2: catch trials, cond 3: baseline)
hits_e = [];
hits_e = find(eyeTrials.trialinfo(:,20)==1 & eyeTrials.trialinfo(:,15)==3);

% redefine trials
cfg=[];
cfg.trials = hits_e';
cfg.toilim    = [1  2];
dat_eye = ft_redefinetrial(cfg, eyeTrials);

figure;
for tr = 1:48
    subplot(8,6,tr)
    plot(dat_eye.trial{1, tr}(2,:),dat_eye.trial{1, tr}(1,:))
end


