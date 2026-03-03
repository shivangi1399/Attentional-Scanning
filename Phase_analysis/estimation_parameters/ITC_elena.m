
% Description:
% -------------
% last update: 20220510
% Goal: compute intertial coherence and phase opposition metrics (like POS)
% ----
% Load the pre-computed spectra per trial and session and compute
% The ITC-R coherence and other phase opposition
% metrics are computed between hits and misses
% if MS_control: the eye velocity per trial is loaded and only trials with eye
% ------------- velocity < eye_thres are used in the ITC analysis.
%               Otherwise, all hits and misses are included
% if loc_control: include only trials with specific target location
% -------------- Otherwise, all hits and misses included

% Still TO DO:
% make code time resolved! pos(channel x freq x toi) - wih SLURM?
% imagesc(pos,squeeze(ichan,toi,freq)
% imagesc(coh,squeeze(ichan,toi,freq)

% fix error with no RT_control!
% run with MS control, next with RT control, next with loc_control and save

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

%% Create data paths to use with slurm
% path to temporary data
datafolder = '/mnt/hpc/projects/MWSampling/4Shivangi/results';

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
    %'klecks_20170825_attentional-sampling_1';... %no trials taken
    'klecks_20170828_attentional-sampling_1';...
    'klecks_20170830_attentional-sampling_1';...
    'klecks_20170831_attentional-sampling_1';...
    'klecks_20170901_attentional-sampling_1';...
    'klecks_20170904_attentional-sampling_1';...
    'klecks_20170906_attentional-sampling_1';...
    'klecks_20170908_attentional-sampling_1';...
    %'klecks_20170911_attentional-sampling_1';... %problem in phase calculation
    'klecks_20170913_attentional-sampling_1';...
    'klecks_20170914_attentional-sampling_1';...
    'klecks_20170915_attentional-sampling_1';...
    'klecks_20170919_attentional-sampling_1';...
    'klecks_20171020_attentional-sampling_1'};


% path to save (slurm) output
output_folder = '/mnt/hpc/projects/MWSampling/4Shivangi/results';
phase_paths = cellfun(@(x) fullfile(output_folder, x,'Phase_analysis/hit_miss/50'),session_names, 'uniform',0);
erp_ph_paths = cellfun(@(x) fullfile(output_folder, x,'ERP_Phase_corr'),session_names, 'uniform',0);

phase_folder = '/mnt/hpc/projects/MWSampling/4Shivangi/results/phase_coherence/phase';

%%
% load channel list per session
cd('/mnt/hpc/projects/MWSampling/4Shivangi/results')
load('channel_sessions_ph.mat'); % var name: channels_sessions
chan_orig = 1:64;

timepoints = 0; % vector with timepoint of phase estimation

% logicals
concatenate_phsaes = 0;
compute_grandgrand = 0;
compute_within_sess =1;

%%
if concatenate_phsaes
    for it = 1
        
        toi = timepoints(it);
        phase_outfolder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/results/phase_coherence/phase/',num2str(toi));
        
        ph_all =[]; trlinfo_all = [];
        
        % loop through sessions and concatenate the phases at toi over all
        % sessions
        
        for isess = 1:length(phase_paths)
            
            cd(fullfile(phase_paths{isess},'final_100iter_cut@0.03/phase_100/stim_onset')); % CHANGE when more tois are used
            
            % Get names of existing channels in each session folder and sort them in assending order
            dr = [];
            dr = dir;
            dr = dr(3:end);
            k = [];
            for ichan = 1:length(dr), k(ichan) = str2num(dr(ichan).name); end
            k = sort(k);
            
            % compare existing channels with original vector with 64 entries
            j = 0; ph =[];
            for ichan = 1:length(chan_orig)
                
                lg = channels_sessions(ichan,isess)==1;
                
                if lg
                    j = j+1
                    cd(fullfile(phase_paths{isess},'final_100iter_cut@0.03/phase_100/stim_onset',num2str(k(j)))) % CHANGE when more tois are used
                    load('transf.mat')
                    
                    ph(:,:,ichan) = transf.ar_transform; % transf.ar_transform: trials x frequency
                else
                    ph(:,:,ichan) = nan(size(transf.ar_transform));
                end
            end
            
            ph_all = [ph; ph_all];
            
            trlinfo_all = [trlinfo_all; transf.trialinfo];
            clear transf
            
        end
        
        ph =[];
        ph.phase_all = ph_all;
        ph.trialinfo = trlinfo_all;
        ph.dimord    = 'trlx freq x chan';
        ph.toi       = 0; % timepoint of phase estimation
        
        if ~isdir(phase_outfolder)
            mkdir(fullfile(phase_outfolder))
        end
        
        cd(phase_outfolder)
        save('ph_all_sess','ph')
        
        clear ph_all_sess ph
    end
    
end


%% calculate phase opposition metrics within each session

if compute_within_sess
    for it = 1
        
        toi = timepoints(it);
        phase_outfolder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/results/phase_coherence/phase/',num2str(toi));
        
        % load the matrix containg data from all sessions
        cd(fullfile(phase_folder,num2str(timepoints(it))))
        load('ph_all_sess.mat')
        
        % define transisions from one session to another using the the
        % recording sample info (increases within session and abruptly
        % decreases afterwards)
        [~, sess_ind] = findpeaks(ph.trialinfo(:,1));
        sess_ind = [1; sess_ind+1; length(ph.phase_all)];
        
        target_loc = unique(ph.trialinfo(:,16));
        
        % combine neighboring target locations to increase number of trials
        % per location
        target_ind = 1:9;
        target_groups = [];
        for igroup = 1:max(target_ind)-2
            
            target_groups(igroup,:) = [target_ind(igroup) target_ind(igroup)+1 target_ind(igroup)+2]
            
        end
        
        ITCr_loc = []; POS_loc = [];
        
        for isess = 1:length(phase_paths)
                
                str =  sess_ind(isess);% beginning of session
                stp =  sess_ind(isess+1)-1; % end of session
                
                %%  compute intertrial phase opposition metrics per target location group
                for igroup = 1:max(target_ind)-2
                    
                    hits = []; misses = [];
                    hits   = find(ph.trialinfo(str:stp,20)==1 & (ph.trialinfo(str:stp,16))==target_loc(target_groups(igroup,1))|ph.trialinfo(str:stp,16)==target_loc(target_groups(igroup,2))...
                        | ph.trialinfo(str:stp,16)==target_loc(target_groups(igroup,3)));
                    misses = find(ph.trialinfo(str:stp,20)==5 & (ph.trialinfo(str:stp,16))==target_loc(target_groups(igroup,1))|ph.trialinfo(str:stp,16)==target_loc(target_groups(igroup,2))...
                        | ph.trialinfo(str:stp,16)==target_loc(target_groups(igroup,3)));
                    
                    % create subselection to equalize number of trials per condition
                    min_trlnum = min([length(hits),length(misses)]);
                    hits2use = []; miss2use =[]; % maybe use the same across trials
                    
                    for isel = 1:1000
                        hits2use(:,isel)= randsample(hits, min_trlnum);
                        miss2use(:,isel)= randsample(misses, min_trlnum);
                    end
                    
                    %% compute ITC
                    for ichan = 1:64
                        ichan
                        ITC_subsel = [];
                        for isel = 1:1000
                        
                        ph_hits = []; ph_hits = ph.phase_all(hits2use(:,isel),:,ichan);
                        ph_miss = []; ph_miss = ph.phase_all(miss2use(:,isel),:,ichan);
                        
                        temp(isel,:) = fun_itc_coherence([ph_hits; -1.*ph_miss]);
                        
                        
                        % POS
                        coh_hits = fun_itc_coherence([ph_hits;]);
                        coh_miss = fun_itc_coherence([ph_miss]);
                        coh_all  = fun_itc_coherence([ph_hits; ph_miss]);
                        
                        temp_pos(isel,:)= coh_hits+coh_miss -(2.*(coh_all));
                        coh_all = []; coh_miss = []; coh_hits = [];
                    end
                    
                    ITCr_loc.coh{isess}(ichan,:,igroup)= nanmean(temp,1); temp = [];
                    POS_loc.pos{isess}(ichan,:) = nanmean(temp_pos,1); temp_pos = [];
                    
                end
                
                %         ITCr_loc.de(:,1)= de;
                %         ITCr_loc.de_thres  = [min_thres;max_thres];
                ITCr_loc.trialinfo{isess,igroup}(:,:,igroup) =  ph.trialinfo([hits;misses],:);
                
                %         POS_loc.de(ichan,:,1) = de;
                %         POS_loc.de_thres  = [min_thres;max_thres];
                POS_loc.trialinfo{isess,igroup}(:,:,igroup) = ph.trialinfo([hits;misses],:);
                

            end
            
            
        end
    end
end

%%

figure; 

% Process ITC data
subplot(1,2,1); 
coh_data = cellfun(@(x) squeeze(nanmean(x, 3)), ITCr_loc.coh, 'UniformOutput', false);
coh_matrix = cell2mat(coh_data(:)); 
plot(2:2:60, nanmean(coh_matrix, 1)); 
title('ITC Coherence'); %x axis freq bins

% Process POS data
subplot(1,2,2); 
pos_data = cellfun(@(x) squeeze(nanmean(x, 1)), POS_loc.pos, 'UniformOutput', false);
pos_matrix = cell2mat(pos_data(:));
plot(2:2:60, nanmean(pos_matrix, 1)); 
title('POS Coherence');

%% calculate intertrial coherence and phase opposition metrics
coh_outfolder = '/mnt/hpc/projects/MWSampling/4Shivangi/results/phase_coherence/coherence_results';


for it = 1
    
    toi = timepoints(it);
    phase_outfolder   = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/results/phase_coherence/phase/',num2str(toi));
    coh_outfolder_toi = fullfile(coh_outfolder,num2str(toi));
    if ~isdir(coh_outfolder_toi)
        mkdir(coh_outfolder_toi)
    end
    
    cd(phase_outfolder)
    load('ph_all_sess.mat')
    
    targ_loc = unique(ph.trialinfo(:,16));
    
    if compute_grandgrand % CHANGE
        %%
        % phase opposition across all target locations and DE
        hits = []; misses = [];
        hits   = find(ph.trialinfo(:,20)==1);
        misses = find(ph.trialinfo(:,20)==5); % TO DO! IMPLEMENT TRIAL SUBSELECTION!
        
        % 1) compute ITCr
        ITCr_all = [];
        
        for ichan = 1:64
            ph_hits = []; ph_hits = ph_all(hits,:,ichan);
            ph_miss = []; ph_miss = ph_all(misses,:,ichan);
            
            ITCr_all.coh(ichan,:) = fun_itc_coherence([ph_hits; -1.*ph_miss]);
        end
        
        ITCr_all.trialinfo = ph.trialinfo;
        cd(coh_outfolder_toi)
        save('ITCr_all','ITCr_all')
        
        % 1) compute POS, POS = ITCa + ITCb - 2ITCall
        for ichan = 1:64
            ph_hits = []; ph_hits = ph_all(hits,:,ichan);
            ph_miss = []; ph_miss = ph_all(misses,:,ichan);
            
            POS.coh_hits(ichan,:)= fun_itc_coherence([ph_hits;]);
            POS.coh_miss(ichan,:)= fun_itc_coherence([ph_miss]);
            POS.coh_all(ichan,:) = fun_itc_coherence([ph_hits; ph_miss]);
            POS.pos(ichan,:) =  (POS.coh_hits(ichan,:)+POS.coh_miss(ichan,:))-2.*(POS.coh_all(ichan,:));
        end
        
        POS.trialinfo = ph.trialinfo;
        cd(coh_outfolder_toi)
        save('POS_all','POS')
        
    end
    
    POS_loc = []; ITCr_loc = [];
    
    for iloc = 1:length(targ_loc)
        
        iloc
        disp('--------------------')
        hits = []; misses = [];
        hits   = find(ph.trialinfo(:,20)==1 & ph.trialinfo(:,16)==targ_loc(iloc));
        misses = find(ph.trialinfo(:,20)==5 & ph.trialinfo(:,16)==targ_loc(iloc));
        
        % constrain difficulty levels
        de = [];
        de = ph.trialinfo([hits;misses],18);
        min_thres = quantile(de,0.05);
        max_thres = quantile(de,1-0.05);
        
        hits = []; misses = [];
        hits   = find(ph.trialinfo(:,20)==1 & ph.trialinfo(:,16)==targ_loc(iloc)& ph.trialinfo(:,18)>min_thres & ph.trialinfo(:,18)<max_thres);
        misses = find(ph.trialinfo(:,20)==5 & ph.trialinfo(:,16)==targ_loc(iloc)& ph.trialinfo(:,18)>min_thres & ph.trialinfo(:,18)<max_thres);
        
        %     figure; hist(de, unique(de))
        %     hold on; plot([min_thres min_thres],ylim,'r')
        %     plot([max_thres max_thres],ylim,'r')
        
        % create subselection to equalize number of trials per condition
        min_trlnum = min([length(hits),length(misses)]);
        hits2use = []; miss2use =[]; % maybe use the same across trials
        
        for isel = 1:1000
            hits2use(:,isel)= randsample(hits, min_trlnum);
            miss2use(:,isel)= randsample(misses, min_trlnum);
        end
        
        %% compute ITC
        for ichan = 1:64
            ichan
            ITC_subsel = [];
            for isel = 1:1000
                
                ph_hits = []; ph_hits = ph.phase_all(hits2use(:,isel),:,ichan);
                ph_miss = []; ph_miss = ph.phase_all(miss2use(:,isel),:,ichan);
                
                temp(isel,:) = fun_itc_coherence([ph_hits; -1.*ph_miss]);
                
                
                % POS
                coh_hits = fun_itc_coherence([ph_hits;]);
                coh_miss = fun_itc_coherence([ph_miss]);
                coh_all  = fun_itc_coherence([ph_hits; ph_miss]);
                
                temp_pos(isel,:)= coh_hits+coh_miss -(2.*(coh_all));
                coh_all = []; coh_miss = []; coh_hits = [];
            end
            
            ITCr_loc.coh(ichan,:,1)= nanmean(temp,1); temp = [];
            POS_loc.pos(ichan,:,1) = nanmean(temp_pos,1); temp_pos = [];
                      
        end
        
        ITCr_loc.de(:,1)= de;
        ITCr_loc.de_thres  = [min_thres;max_thres];
        ITCr_loc.trialinfo =  ph.trialinfo([hits;misses],:);
        
        POS_loc.de(ichan,:,1) = de;
        POS_loc.de_thres  = [min_thres;max_thres];
        POS_loc.trialinfo = ph.trialinfo([hits;misses],:);
        
          
        cd(coh_outfolder_toi)
        save(['ITCr_loc','_',num2str(iloc)],'ITCr_loc')
        clear ITCr_loc
        
        cd(coh_outfolder_toi)
        save(['POS_loc','_',num2str(iloc)],'POS_loc')
        clear POS_loc
        
    end

    
end

%% Plot phase distributions
figure; 
ichan = 10
for ifreq =1:30
    subplot(10,3,ifreq)
    tt =[]; tm =[];
    tt = squeeze(angle(ph.phase_all(hits,ifreq,ichan)));
    tm = squeeze(angle(ph.phase_all(misses,ifreq,ichan)));

    histogram(tt,'FaceAlpha',0.4,'NumBins',5,'Normalization','probability'); 
    hold on,histogram(tm,'FaceAlpha',0.4,'NumBins',5,'Normalization','probability')
    
end


%% Compute Intertrial phase clustering per target location
i = 0;
temp =[];
 for iloc = 1:length(targ_loc)
        
        iloc
        disp('--------------------')
        hits = []; misses = [];
        hits   = find(ph.trialinfo(:,20)==1 & ph.trialinfo(:,16)==targ_loc(iloc));
        misses = find(ph.trialinfo(:,20)==5 & ph.trialinfo(:,16)==targ_loc(iloc));
        
        % constrain difficulty levels
        de = [];
        de = ph.trialinfo([hits;misses],18);
        min_thres = quantile(de,0.01);
        max_thres = quantile(de,1-0.01);
        
        hits = []; misses = [];
        hits   = find(ph.trialinfo(:,20)==1 & ph.trialinfo(:,16)==targ_loc(iloc)& ph.trialinfo(:,18)>min_thres & ph.trialinfo(:,18)<max_thres);
        misses = find(ph.trialinfo(:,20)==5 & ph.trialinfo(:,16)==targ_loc(iloc)& ph.trialinfo(:,18)>min_thres & ph.trialinfo(:,18)<max_thres);
        
        % compare distributions before and after removing extreme de
        % figure; histogram(de,unique(de),'FaceAlpha',0.4), hold on, histogram(ph.trialinfo([hits;misses],18), unique(de))
        
        
        % create subselection to equalize number of trials per condition
        min_trlnum = min([length(hits),length(misses)]);
        hits2use = []; miss2use =[]; % maybe use the same across trials
        
        for isel = 1:1000
            hits2use(:,isel)= randsample(hits, min_trlnum);
            miss2use(:,isel)= randsample(misses, min_trlnum);
        end
        
        
        
        
        %% compute ITPC (weighted by DE)
        for ichan = 1:64
            ichan
            ITC_subsel = [];
            for isel = 1:1000
            
                % weight based on DE
                w_hits = [];
                w_hits = 1./ph.trialinfo(hits2use(:,isel),18);
        
                ph_hits = []; ph_hits = ph.phase_all(hits2use(:,isel),:,ichan);
                
%                 ITPC_h_sel(isel,:) = abs(mean(exp(1i*angle(ph_hits))));
                ITPC_h_sel(isel,:) = abs(nanmean(w_hits.*exp(1i*angle(ph_hits))));
            end
            
            ITPC_hits_loc(ichan,:,iloc) = nanmean(ITPC_h_sel,1);
            ITPC_h_sel =[];
            
            w_miss = [];
            w_miss = 1./ph.trialinfo(miss2use(:,1),18);

            ph_miss = []; ph_miss = ph.phase_all(miss2use(:,1),:,ichan);
            ITPC_miss_loc(ichan,:,iloc)= abs(nanmean(w_miss.*exp(1i*angle(ph_miss))));
        end
        
 end
 
% % % figure; 
% % % for ichan = 1:64
% % %     subplot(8,8,ichan)
% % %     plot(ITPC_hits(ichan,:)),
% % %     hold on,
% % %     plot(ITPC_miss(ichan,:))
% % % end

figure
for iloc = 1:9
    
    %     cd(coh_outfolder_toi)
    %     load(['ITCr_loc','_',num2str(iloc),'.mat'])
    %
    %     subplot(3,3,iloc)
    %     plot(nanmean(ITCr_loc.coh(:,:),1))
    
    cd(coh_outfolder_toi)
    load(['POS_loc','_',num2str(iloc),'.mat'])
    
    subplot(3,3,iloc)
    plot(nanmean(POS_loc.pos(:,:),1))
    
end


close all
for iloc = 1:9
   figure
    
    cd(coh_outfolder_toi)
    load(['POS_loc','_',num2str(iloc),'.mat'])
    for ichan = 1:64
        subplot(8,8,ichan)
        plot(POS_loc.pos(ichan,:))
    end
end


close all
for iloc = 1:9
    
        cd(coh_outfolder_toi)
        load(['ITCr_loc','_',num2str(iloc),'.mat'])
    %
    %     subplot(3,3,iloc)
    %     plot(nanmean(ITCr_loc.coh(:,:),1))
    
    cd(coh_outfolder_toi)
%     load(['POS_loc','_',num2str(iloc),'.mat'])

    figure
    for ichan = 1:64
        subplot(8,8,ichan)
        %         plot(POS_loc.pos(ichan,:))
        plot(ITCr_loc.coh(ichan,:))
        
    end
end





figure
for iloc = 1:9
    subplot(3,3,iloc)
    plot(freqpow.freq, nanmean(ITC_r_hits(:,:,iloc)))
    hold on
    plot(freqpow.freq, nanmean(ITC_r(:,:,iloc)))
    
    plot(freqpow.freq, nanmean(ITC_r_miss(:,:,iloc)))
    
    
    title(num2str(iloc))
    
    %     hold on
    %     pause
end
c = parula(9);

figure
for ichan = 1:64
    subplot(8,8,ichan)
    %     plot(freqpow.freq, nanmean(ITC_r_hits(:,:,iloc)))
    hold on
    
    for iloc = 1:8
        plot(freqpow.freq, squeeze(ITC_r(ichan,:,iloc)),'color', c(iloc,:))
        %         plot(freqpow.freq, squeeze(ITC_r_miss(ichan,:,iloc)),'color', c(iloc,:))
        
        hold on
    end
    
    title(num2str(iloc))
    
    %     hold on
    %     pause
end

c =[]; c= parula(3);

figure
for ichan = 1:64
    subplot(8,8,ichan)
    %     plot(freqpow.freq, nanmean(ITC_r_hits(:,:,iloc)))
    hold on
    
    plot(freqpow.freq, squeeze(nanmean(ITC_r(ichan,:,1:3),3)),'color', c(1,:))
    plot(freqpow.freq, squeeze(nanmean(ITC_r(ichan,:,4:6),3)),'color', c(2,:))
    plot(freqpow.freq, squeeze(nanmean(ITC_r(ichan,:,7:9),3)),'color', c(3,:))
    
    title(num2str(iloc))
    
    %     hold on
    %     pause
end


% average over channels
group1 = nanmean(nanmean(ITC_r(:,:,1:3),1),3);
group2 = nanmean(nanmean(ITC_r(:,:,4:6),1),3);
group3 = nanmean(nanmean(ITC_r(:,:,7:9),1),3);

group1_h = nanmean(nanmean(ITC_r_hits(:,:,1:3),1),3);
group2_h = nanmean(nanmean(ITC_r_hits(:,:,4:6),1),3);
group3_h = nanmean(nanmean(ITC_r_hits(:,:,7:9),1),3);


figure;
plot(freqpow.freq,group1,'color', c(1,:))
hold on
plot(freqpow.freq,group2,'color', c(2,:))
plot(freqpow.freq,group3,'color', c(3,:))

figure;
plot(freqpow.freq,group1_h,'color', c(1,:))
hold on
plot(freqpow.freq,group2_h,'color', c(2,:))
plot(freqpow.freq,group3_h,'color', c(3,:))
plot(freqpow.freq,mean([group1_h;group2_h;group3_h]),'color', 'r')



for isess = 1:length(phase_paths)
    
    cd(fullfile(phase_paths{isess},'final_100iter_cut@0.03/phase_100/stim_onset'))
    
    % look for available channels (note that some channels might have been
    % excluded due to artifacts
    dr = dir;
    ITC_r = [];
    
    for iloc = 1:9
        
        cd(fullfile(phase_paths{isess},'final_100iter_cut@0.03/phase_100/stim_onset','1'))
        load('transf.mat') % -> transf.ar_transform % trials x frequency
        
        loc = unique(transf.trialinfo(:,16));
        
        %% trial selection
        hits = []; misses = [];
        hits   = find(transf.trialinfo(:,20)==1 & transf.trialinfo(:,16)==loc(iloc));
        misses = find(transf.trialinfo(:,20)==5 & transf.trialinfo(:,16)==loc(iloc));
        min_trlnum = min([length(hits),length(misses)]);
        
        % create subselection
        hits2use = []; miss2use =[]; % maybe use the same across trials
        for isel = 1:1000
            hits2use(:,isel)= randsample(hits, min_trlnum);
            miss2use(:,isel)= randsample(misses, min_trlnum);
        end
        
        
        for ichan = 1:length(dr)-2
            i = 2+ichan;
            chan_name = dr(i).name
            
            cd(fullfile(phase_paths{isess},'final_100iter_cut@0.03/phase_100/stim_onset',chan_name))
            load('transf.mat') % -> transf.ar_transform % trials x frequency
            
            
            %        coh_hits = fun_itc_coherence(transf.ar_transform(hits2use));
            %        coh_misses = fun_itc_coherence(transf.ar_transform(miss2use));
            
            
            ITC_subsel = [];
            for isel = 1:1000
                
                ITC_subsel(isel,:) = fun_itc_coherence([transf.ar_transform(hits2use(:,isel),:); -1*(transf.ar_transform(miss2use(:,isel),:))]);
                ITC_subsel_test(isel,:) = fun_itc_coherence([transf.ar_transform(hits2use(:,isel),:)]);
                
            end
            
            ITC_r(ichan,:,iloc) = nanmean(ITC_subsel);
            ITC_r_test(ichan,:,iloc) = nanmean(ITC_subsel_test);
            
        end
        
        %         figure
        %         for ichan = 1:size(ITC_r)
        %             subplot(8,8,ichan)
        %             plot(ITC_r(ichan,:))
        %         end
        % %
        %         figure; plot( freqpow.freq,nanmean(ITC_r)), hold on
        %         plot( freqpow.freq,nanmean(ITC_r_test))
        %
        
    end
end

