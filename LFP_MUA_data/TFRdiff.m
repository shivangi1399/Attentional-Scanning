clear all
close all
clc
%% paths linux

addpath /opt/fieldtrip_github/
ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf

% path to temporary data
datafolder = '/mnt/hpc/projects/MWSampling/4Shivangi/data'; % in windows:\\cs\projects etc
outputfolder = '/mnt/hpc/projects/MWSampling/4Shivangi/results';

s1 = 'klecks_20170804_attentional-sampling_1';
s2 = 'klecks_20170807_attentional-sampling_1';
s3 = 'klecks_20170808_attentional-sampling_1';
s4 = 'klecks_20170810_attentional-sampling_1';
s5 = 'klecks_20170825_attentional-sampling_1';
s6 = 'klecks_20170828_attentional-sampling_1';
s7 = 'klecks_20170901_attentional-sampling_1';

%% variables

run_perm_real=0;
per_session=0; %isess=6;
pooled_data=1; remove_sess=1;

%% save permutations and real TFR

if run_perm_real
    
    %% TFR Perm
 
    for i = 1:7
        i
        v = ['cd(fullfile(datafolder,s' int2str(i) ')),'];
        eval(v);
        
        % load LFP data
        load('zlfpTrials.mat')
        lfpTrials = zlfpTrials;
        dat = lfpTrials;
        %load('freqpow.mat'); %var window
        
        cfg              = [];
        cfg.output       = 'pow';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.foi          = 2:2:100;
        cfg.t_ftimwin    = 3./cfg.foi;
        cfg.toi          = -1:0.05:0.6;
        cfg.keeptrials  = 'yes';
        freqpow = ft_freqanalysis(cfg,dat);
        v = ['cd(fullfile(outputfolder,s' int2str(i) ')),'];
        eval(v);
        cd('TFR')
        save freqpow.mat
        
        % get triial info
        trials_ind = [];
        trials_ind = find(lfpTrials.trialinfo(:,20)==1);
        trials_ind = [trials_ind ; find(lfpTrials.trialinfo(:,20)==5)];
        
        hit_n = sum(lfpTrials.trialinfo(:,20)==1);
        miss_n  = sum(lfpTrials.trialinfo(:,20)==5);
        
        %shuffle
        rand_matrix = nan(1000, length(trials_ind));
        for iter = 1:1000
            rand_matrix(iter,:)= permutate(trials_ind');
        end
        
        % TFR diff
        
        TFR_diff = nan(1000,64,50,33);
        
        
        for perm=1:1000
            perm
            
            
            hits = []; misses = [];
            
            hits = rand_matrix(perm,1:hit_n);
            misses = rand_matrix(perm,hit_n+1:end);
            trial_no = min(length(hits),length(misses));
            hits = hits(randperm(max(length(hits),length(misses)),trial_no))';
            
            % hits
            
            freqpow_H = freqpow.powspctrm(hits,:,:,:);
            freqpowavg_H = mean(freqpow_H,1); %avg over all trials
            
            % misses
            freqpow_M = freqpow.powspctrm(misses,:,:,:);
            freqpowavg_M = mean(freqpow_M,1); %avg over all trials
            
            
            % diff
            TFR_diff(perm,:,:,:) = log10((freqpowavg_H)./(freqpowavg_M));
            % imagesc(freqpow.time,freqpow.freq,reshape(TFR_diff(perm,ichan,:,:),30,33))
            % set(gca,'YDir','normal')
            
        end
        v = ['cd(fullfile(outputfolder,s' int2str(i) ')),'];
        eval(v);
        cd('TFR')
        save TFR_diff.mat
    end
    
    
    %% TFR real
    
    TFRdiff_r = nan(7,64,50,33);
    
    for i=1:7
        
        v = ['cd(fullfile(datafolder,s' int2str(i) ')),'];
        eval(v);
        
        % load LFP data
        load('zlfpTrials.mat')
        lfpTrials = zlfpTrials;
        dat = lfpTrials;
        
        % average evoked response of hit trials
        hits = find(lfpTrials.trialinfo(:,20)==1);
        misses = find(lfpTrials.trialinfo(:,20)==5);
        trial_no = min(length(hits),length(misses));
        
        dat = lfpTrials;
        
        % avg after freq analysis
        % hits
        hit = randperm(max(length(hits),length(misses)),trial_no);
        cfg = [];
        cfg.trials = hits(hit)';
        timelock_H = ft_selectdata(cfg,dat);
        
        % misses
        cfg = [];
        cfg.trials = misses';
        timelock_M = ft_selectdata(cfg,dat);
        
        
        %%% hits
        
        cfg              = [];
        cfg.output       = 'pow';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.foi          = 2:2:100;
        cfg.t_ftimwin    = 3./cfg.foi;
        cfg.toi          = -1:0.05:0.6;
        freqpow_H = ft_freqanalysis(cfg, timelock_H);
        
        
        %%% misses
        
        cfg              = [];
        cfg.output       = 'pow';
        cfg.method       = 'mtmconvol';
        cfg.taper        = 'hanning';
        cfg.foi          = 2:2:100;
        cfg.t_ftimwin    = 3./cfg.foi;
        cfg.toi          = -1:0.05:0.6;
        freqpow_M = ft_freqanalysis(cfg, timelock_M);
        
        %TFRdiff_r = log10((freqpow_H.powspctrm)./(freqpow_M.powspctrm));
        TFRdiff_r(i,:,:,:) = log10((freqpow_H.powspctrm)./(freqpow_M.powspctrm));
        
    end
    
    cd('/mnt/hpx/projects/MWSampling/4Shivangi/Plots/freqanalysis/TFRdiff')
    save TFRdiff_r.mat
    
end


%% Multiple Comparison

if per_session
    
    
    v = ['cd(fullfile(outputfolder,s' int2str(isess) ')),'];
    eval(v);
    cd('TFR')
    load('TFR_diff.mat')
    
    cd('/mnt/hpx/projects/MWSampling/4Shivangi/Plots/freqanalysis/TFRdiff')
    load('TFRdiff_r.mat')
    
    limit= zeros(64,2);
    
    h1 = figure('name','TFR diff','Unit','centimeter', 'Position', [0 0 30 20]);
    
    for ichan=1:64
        
        max_r = max((reshape(TFR_diff(:,ichan,:,:),1000,50,33)),[],[2,3]);
        min_r = min((reshape(TFR_diff(:,ichan,:,:),1000,50,33)),[],[2,3]);
        
        t_dist = [max_r;min_r];
        limit(ichan,1) = quantile(t_dist(:,1), 0.975);
        limit(ichan,2) = quantile(t_dist(:,1), 0.025);
        
        % plotting
        plotm = zeros(50,33);
        
        A = reshape(TFRdiff_r(isess,ichan,:,:),50,33);
        
        for i=1:50
            for j=1:33
                if A(i,j)<limit(ichan,2)|A(i,j)>limit(ichan,1)
                    plotm(i,j)=1;
                end
            end
        end
        
        subplot(8,8,ichan)
        imagesc(freqpow.time,freqpow.freq,reshape(TFRdiff_r(isess,ichan,:,:),50,33))
        set(gca,'YDir','normal')
        hold on
        imagesc(freqpow.time,freqpow.freq,plotm)
        set(gca,'YDir','normal')
        xlim([-0.5 0.15])
        clims = [limit(ichan,1) limit(ichan,2)];
        
    end
    cd('/mnt/hpx/projects/MWSampling/4Shivangi/Plots/freqanalysis/TFRdiff')
    v=['savefig("session_' int2str(isess) '.fig")'];
    eval(v);
end



%% Multiple comparisons for pooled data

if pooled_data
    cd('/mnt/hpx/projects/MWSampling/4Shivangi/Plots/freqanalysis/TFRdiff')
    load('TFRdiff_r.mat')
    if remove_sess
        remove_sessions = [2,5,7];
        TFRdiff_r(2,:) = [];
        TFRdiff_r(4,:) = [];
        TFRdiff_r(5,:) = [];
    end
    TFRdiff_ravg = reshape(mean(TFRdiff_r,1),64,50,33);
    
    
    maxmin_matrix = zeros(64,2000,7);
    
    for isess=[1 3 4 5 6]
        v = ['cd(fullfile(outputfolder,s' int2str(isess) ')),'];
        eval(v);
        cd('TFR')
        load('TFR_diff.mat')
        
        for ichan=1:64
            
            max_r = max((reshape(TFR_diff(:,ichan,:,:),1000,50,33)),[],[2,3]);
            min_r = min((reshape(TFR_diff(:,ichan,:,:),1000,50,33)),[],[2,3]);
            
            maxmin_matrix(ichan,:,isess) = [max_r;min_r];
            
        end
    end
    
    if remove_sess
        remove_sessions = [2,5,7];
        maxmin_matrix(:,:,2) = [];
        maxmin_matrix(:,:,4) = [];
        maxmin_matrix(:,:,5) = [];
    end
    
    t_dist =[];
    t_dist=[maxmin_matrix(:,:,1) maxmin_matrix(:,:,2) maxmin_matrix(:,:,3) maxmin_matrix(:,:,4)];
    limit= zeros(64,2);
    
        %% plotting method 1
    
        h1 = figure('name','TFR diff','Unit','centimeter', 'Position', [0 0 30 20]);
        for ichan=1:64
            
            %             limit(ichan,1) = quantile(t_dist(ichan,:), 0.975);
            %             limit(ichan,2) = quantile(t_dist(ichan,:), 0.025);
            %
            %             plotm = zeros(50,33);
            %             A = reshape(TFRdiff_ravg(ichan,:,:),50,33);
            %
            %             for i=1:50
            %                 for j=1:33
            %                     if A(i,j)<limit(ichan,2)|A(i,j)>limit(ichan,1)
            %                         plotm(i,j)=1;
            %                     end
            %                 end
            %             end
            
            subplot(8,8,ichan)
            imagesc(freqpow.time,freqpow.freq,reshape(TFRdiff_ravg(ichan,:,:),50,33))
            set(gca,'YDir','normal')
            %             hold on
            %             imagesc(freqpow.time,freqpow.freq,plotm)
            %             set(gca,'YDir','normal')
            hold on
            a = gca;
            CLr(ichan,1)= a.CLim(1,1);
            CLr(ichan,2)= a.CLim(1,2);
            xlim([-0.1 0.15])
            caxis([-0.376 0.573])
        end
        colorbar('south')
        %     cd('/mnt/hpx/projects/MWSampling/4Shivangi/Plots/freqanalysis/TFRdiff')
        %     savefig('TFR_diff_pooled')
    
    
    %% plotting method 2
    
    CLr = zeros(64,2);
    h1 = figure('name','TFR diff','Unit','centimeter', 'Position', [0 0 30 20]);
    for ichan=1%:64
        
        limit(ichan,1) = quantile(t_dist(ichan,:), 0.975);
        limit(ichan,2) = quantile(t_dist(ichan,:), 0.025);
        
        plotm = reshape(TFRdiff_ravg(ichan,:,:),50,33);
        A = reshape(TFRdiff_ravg(ichan,:,:),50,33);
        
        for i=1:50
            for j=1:33
                if (limit(ichan,1)> A(i,j)) && (A(i,j) >limit(ichan,2))
                    plotm(i,j)= NaN;
                end
            end
        end
        
        %subplot(8,8,ichan)
        %         imagesc(freqpow.time,freqpow.freq,reshape(TFRdiff_ravg(ichan,:,:),50,33))
        %         set(gca,'YDir','normal')
        %         hold on
        imagesc(freqpow.time,freqpow.freq,plotm)
        set(gca,'YDir','normal')
        %         hold on
        %         a = gca;
        %         CLr(ichan,1)= a.CLim(1,1);
        %         CLr(ichan,2)= a.CLim(1,2);
        xlim([-0.1 0.15])
        caxis([-0.37 1])
        xline(0.1,'w')
    end
    colorbar('south')
%     cd('/mnt/hpx/projects/MWSampling/4Shivangi/Plots/freqanalysis/TFRdiff')
%     savefig('TFR_diff_pooled')
    
end



















