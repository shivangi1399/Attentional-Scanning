% Multiple linear regression: phase predicts RT / MUA / LFP / hit-miss
% H3 (abs_per_pos_diff): group trials by (stimulus position x difficulty bin)
% cells; difficulty (trialinfo col 18) is binned into nDiffBins quantile bins
% within each position. Run regression within each cell, average R^2 across
% cells, then across channels and animals.
clearvars; close all; clc

%% Dependencies
addpath /opt/fieldtrip_github/; ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/multiple_linear_reg/functions
addpath /mnt/hpc/projects/MWSampling/4Shivangi/
clc

%% Settings
nPerm     = 1000;
alpha     = 0.05;
nDiffBins = 4;          % within-position quantile bins of difficulty (col 18)
Y_vars    = {'RT','MUA_ERP_ampl_all','LFP_ERP_ampl_all','hit_miss'};
animals   = {'hermes', 'klecks'};

%% Loop over animals

for a = 1:numel(animals)
    animalName = animals{a};
    fprintf('\n=== Processing %s ===\n', animalName);

    data_folder    = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ...
        ['results_' animalName], 'multi_lin_reg', 'cp10_till_100');
    results_folder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ...
        ['results_' animalName], 'multi_lin_reg', 'abs_per_pos_diff', 'cp10_till_100');
    if ~exist(results_folder,'dir'), mkdir(results_folder); end
    cd(data_folder); load('ph_all_sess.mat')

    ph_all  = ph_comb.phase_all;
    amp_all = ph_comb.amp_all;
    [numTrials, numFreq, numCh] = size(ph_all);

    reg_results = struct();
    reg_results.nDiffBins = nDiffBins;

    for d = 1:length(Y_vars)
        depVarName = Y_vars{d};
        fprintf('Processing %s — %s...\n', animalName, depVarName);

        isLogistic = strcmp(depVarName,'hit_miss');

        if isLogistic
            Y_all_orig  = ph_comb.trialinfo(:,20);
            Y_all_orig(Y_all_orig==5) = 0;
            trlInfo_all = ph_comb.trialinfo;
        else
            Y_all_orig  = ph_comb.(depVarName);
            switch depVarName
                case 'RT',               trlInfo_all = ph_comb.RT_trialinfo;
                case 'MUA_ERP_ampl_all', trlInfo_all = ph_comb.MUA_ERP_trialinfo;
                case 'LFP_ERP_ampl_all', trlInfo_all = ph_comb.LFP_ERP_trialinfo;
            end
        end

        positions = unique(trlInfo_all(:,16));
        nPos      = numel(positions);

        % Within-position quantile bins of difficulty (dE00, col 18)
        diff_bin = bin_difficulty_per_pos(trlInfo_all(:,18), trlInfo_all(:,16), positions, nDiffBins);

        % Enumerate (position x difficulty) cells
        cell_pos = repmat(positions(:), 1, nDiffBins);
        cell_dif = repmat(1:nDiffBins, nPos, 1);
        cell_pos = cell_pos(:);
        cell_dif = cell_dif(:);
        nCell    = numel(cell_pos);

        reg_results.(depVarName).R2_phase    = NaN(numCh, numFreq);
        reg_results.(depVarName).R2_MUA      = NaN(numCh, numFreq);
        reg_results.(depVarName).R2_Amp      = NaN(numCh, numFreq);
        reg_results.(depVarName).R2_AmpPhase = NaN(numCh, numFreq);
        reg_results.(depVarName).R2_any      = NaN(numCh, numFreq);
        reg_results.(depVarName).R_phase     = NaN(numCh, numFreq);
        reg_results.(depVarName).phi_pref    = NaN(numCh, numFreq);
        reg_results.(depVarName).cell_pos    = cell_pos;
        reg_results.(depVarName).cell_dif    = cell_dif;

        for ch = 1:numCh
            fprintf('Channel %d/%d\n', ch, numCh);

            if ~isLogistic
                Y_all = Y_all_orig(:,ch);
                if all(isnan(Y_all)), continue; end
            else
                Y_all = Y_all_orig;
            end

            %% Pre-compute per-cell, per-frequency X and Y for real data and perms

            X_all_freq_cell = cell(nCell, numFreq);
            Y_all_freq_cell = cell(nCell, numFreq);

            Robs_phase    = zeros(numFreq,1);
            Robs_MUA      = zeros(numFreq,1);
            Robs_Amp      = zeros(numFreq,1);
            Robs_AmpPhase = zeros(numFreq,1);
            Robs_any      = zeros(numFreq,1);
            R_phase       = zeros(numFreq,1);
            phi_pref      = zeros(numFreq,1);

            for f = 1:numFreq
                X = [ph_comb.pup_baseline(:,ch), ph_comb.MUA_baseline(:,ch), ...
                     amp_all(:,f,ch), sin(ph_all(:,f,ch)), cos(ph_all(:,f,ch))];

                if strcmp(depVarName,'RT')
                    keepIdx = trlInfo_all(:,20) == 1;
                else
                    keepIdx = true(size(Y_all));
                end

                R2_cell_phase    = NaN(nCell,1);
                R2_cell_MUA      = NaN(nCell,1);
                R2_cell_Amp      = NaN(nCell,1);
                R2_cell_AmpPhase = NaN(nCell,1);
                R2_cell_any      = NaN(nCell,1);
                beta_sin_cell    = NaN(nCell,1);
                beta_cos_cell    = NaN(nCell,1);

                for c = 1:nCell
                    cell_mask = (trlInfo_all(:,16) == cell_pos(c)) & ...
                                (diff_bin == cell_dif(c));
                    nanIdx    = any(isnan(X),2) | isnan(Y_all) | ~keepIdx | ~cell_mask;
                    X_clean   = X(~nanIdx,:);
                    Y_clean   = Y_all(~nanIdx);

                    X_all_freq_cell{c,f} = X_clean;
                    Y_all_freq_cell{c,f} = Y_clean;

                    if length(Y_clean) < size(X_clean,2)+1, continue; end

                    if isLogistic
                        b_full  = glmfit(X_clean,Y_clean,'binomial','link','logit');
                        p_full  = glmval(b_full,X_clean,'logit');
                        LL_full = sum(Y_clean.*log(p_full+eps)+(1-Y_clean).*log(1-p_full+eps));

                        b_null  = glmfit(ones(size(Y_clean,1),1),Y_clean,'binomial','link','logit');
                        p_null  = glmval(b_null,ones(size(Y_clean,1),1),'logit');
                        LL_null = sum(Y_clean.*log(p_null+eps)+(1-Y_clean).*log(1-p_null+eps));
                        R2_cell_any(c) = 1-(LL_full/LL_null);

                        b_red = glmfit(X_clean(:,1:3),Y_clean,'binomial','link','logit');
                        p_red = glmval(b_red,X_clean(:,1:3),'logit');
                        LL_red = sum(Y_clean.*log(p_red+eps)+(1-Y_clean).*log(1-p_red+eps));
                        R2_cell_phase(c) = 1-(LL_full/LL_red);

                        b_red = glmfit(X_clean(:,[1 3 4 5]),Y_clean,'binomial','link','logit');
                        p_red = glmval(b_red,X_clean(:,[1 3 4 5]),'logit');
                        LL_red = sum(Y_clean.*log(p_red+eps)+(1-Y_clean).*log(1-p_red+eps));
                        R2_cell_MUA(c) = 1-(LL_full/LL_red);

                        b_red = glmfit(X_clean(:,[1 2 4 5]),Y_clean,'binomial','link','logit');
                        p_red = glmval(b_red,X_clean(:,[1 2 4 5]),'logit');
                        LL_red = sum(Y_clean.*log(p_red+eps)+(1-Y_clean).*log(1-p_red+eps));
                        R2_cell_Amp(c) = 1-(LL_full/LL_red);

                        b_red = glmfit(X_clean(:,1:2),Y_clean,'binomial','link','logit');
                        p_red = glmval(b_red,X_clean(:,1:2),'logit');
                        LL_red = sum(Y_clean.*log(p_red+eps)+(1-Y_clean).*log(1-p_red+eps));
                        R2_cell_AmpPhase(c) = 1-(LL_full/LL_red);

                    else
                        X_full   = [ones(size(X_clean,1),1), X_clean];
                        b_full   = regress(Y_clean,X_full);
                        RSS_full = sum((Y_clean - X_full*b_full).^2);
                        RSS_null = sum((Y_clean - mean(Y_clean)).^2);

                        R2_cell_any(c) = max(0,(RSS_null-RSS_full)/RSS_null);

                        X_red = [ones(size(X_clean,1),1), X_clean(:,1:3)];
                        b_red = regress(Y_clean,X_red); RSS_red = sum((Y_clean-X_red*b_red).^2);
                        R2_cell_phase(c) = max(0,(RSS_red-RSS_full)/RSS_null);

                        X_red = [ones(size(X_clean,1),1), X_clean(:,[1 3 4 5])];
                        b_red = regress(Y_clean,X_red); RSS_red = sum((Y_clean-X_red*b_red).^2);
                        R2_cell_MUA(c) = max(0,(RSS_red-RSS_full)/RSS_null);

                        X_red = [ones(size(X_clean,1),1), X_clean(:,[1 2 4 5])];
                        b_red = regress(Y_clean,X_red); RSS_red = sum((Y_clean-X_red*b_red).^2);
                        R2_cell_Amp(c) = max(0,(RSS_red-RSS_full)/RSS_null);

                        X_red = [ones(size(X_clean,1),1), X_clean(:,1:2)];
                        b_red = regress(Y_clean,X_red); RSS_red = sum((Y_clean-X_red*b_red).^2);
                        R2_cell_AmpPhase(c) = max(0,(RSS_red-RSS_full)/RSS_null);
                    end

                    beta_sin_cell(c) = b_full(end-1);
                    beta_cos_cell(c) = b_full(end);
                end

                % Average R^2 across (position x difficulty) cells
                Robs_phase(f)    = mean(R2_cell_phase,    'omitnan');
                Robs_MUA(f)      = mean(R2_cell_MUA,      'omitnan');
                Robs_Amp(f)      = mean(R2_cell_Amp,      'omitnan');
                Robs_AmpPhase(f) = mean(R2_cell_AmpPhase, 'omitnan');
                Robs_any(f)      = mean(R2_cell_any,       'omitnan');

                % Preferred phase: circular mean of per-cell beta vectors
                mean_sin = mean(beta_sin_cell, 'omitnan');
                mean_cos = mean(beta_cos_cell, 'omitnan');
                R_phase(f)  = sqrt(mean_sin^2 + mean_cos^2);
                phi_pref(f) = atan2(mean_sin, mean_cos);
            end

            %% Permutation (SLURM)

            nJobs = ceil(nPerm/10);
            cfg_array = cell(nJobs,1);
            output_dir = fullfile(results_folder,'perm_R_pos_diff',depVarName,num2str(ch));
            if ~exist(output_dir,'dir'), mkdir(output_dir); end
            for j = 1:nJobs
                perm_start = (j-1)*10 + 1;
                perm_end   = min(j*10, nPerm);
                cfg_array{j}.X_all_freq_cell = X_all_freq_cell;
                cfg_array{j}.Y_all_freq_cell = Y_all_freq_cell;
                cfg_array{j}.numFreq         = numFreq;
                cfg_array{j}.nCell           = nCell;
                cfg_array{j}.isLogistic      = isLogistic;
                cfg_array{j}.perm_idx        = perm_start:perm_end;
                cfg_array{j}.output_dir      = output_dir;
            end

            slurmfun(@regress_perm_R_per_pos_diff, cfg_array, ...
                'partition','8GB','stopOnError',false,'useUserPath',true, ...
                'waitForToolboxes',{'statistics_toolbox'});

            %% Collect permutation nulls

            null_max_phase    = zeros(nPerm,1);
            null_max_MUA      = zeros(nPerm,1);
            null_max_Amp      = zeros(nPerm,1);
            null_max_AmpPhase = zeros(nPerm,1);
            null_max_any      = zeros(nPerm,1);

            for p = 1:nPerm
                perm_file = fullfile(output_dir, sprintf('perm_%04d.mat',p));
                load(perm_file,'results');
                null_max_phase(p)    = results.null_max_phase;
                null_max_MUA(p)      = results.null_max_MUA;
                null_max_Amp(p)      = results.null_max_Amp;
                null_max_AmpPhase(p) = results.null_max_AmpPhase;
                null_max_any(p)      = results.null_max_any;
            end

            thresholds.thresh_phase    = prctile(null_max_phase,    100*(1-alpha));
            thresholds.thresh_MUA      = prctile(null_max_MUA,      100*(1-alpha));
            thresholds.thresh_Amp      = prctile(null_max_Amp,      100*(1-alpha));
            thresholds.thresh_AmpPhase = prctile(null_max_AmpPhase, 100*(1-alpha));
            thresholds.thresh_any      = prctile(null_max_any,      100*(1-alpha));

            stats.p_phase    = mean(null_max_phase    >= max(Robs_phase));
            stats.p_MUA      = mean(null_max_MUA      >= max(Robs_MUA));
            stats.p_Amp      = mean(null_max_Amp      >= max(Robs_Amp));
            stats.p_AmpPhase = mean(null_max_AmpPhase >= max(Robs_AmpPhase));
            stats.p_any      = mean(null_max_any      >= max(Robs_any));

            reg_results.(depVarName).R2_phase(ch,:)    = Robs_phase;
            reg_results.(depVarName).R2_MUA(ch,:)      = Robs_MUA;
            reg_results.(depVarName).R2_Amp(ch,:)      = Robs_Amp;
            reg_results.(depVarName).R2_AmpPhase(ch,:) = Robs_AmpPhase;
            reg_results.(depVarName).R2_any(ch,:)      = Robs_any;
            reg_results.(depVarName).stats(ch)         = stats;
            reg_results.(depVarName).thresholds(ch)    = thresholds;
            reg_results.(depVarName).R_phase(ch,:)     = R_phase;
            reg_results.(depVarName).phi_pref(ch,:)    = phi_pref;
        end

        %% Channel-average permutation null

        fprintf('Computing channel-average null for %s...\n', animalName);

        null_avg_R_phase_freq    = zeros(nPerm, numFreq);
        null_avg_R_MUA_freq      = zeros(nPerm, numFreq);
        null_avg_R_Amp_freq      = zeros(nPerm, numFreq);
        null_avg_R_AmpPhase_freq = zeros(nPerm, numFreq);
        null_avg_R_any_freq      = zeros(nPerm, numFreq);

        for perm = 1:nPerm
            R_phase_ch    = NaN(numCh, numFreq);
            R_MUA_ch      = NaN(numCh, numFreq);
            R_Amp_ch      = NaN(numCh, numFreq);
            R_AmpPhase_ch = NaN(numCh, numFreq);
            R_any_ch      = NaN(numCh, numFreq);

            for ch = 1:numCh
                perm_file = fullfile(results_folder,'perm_R_pos_diff',depVarName,num2str(ch), ...
                    sprintf('perm_%04d.mat',perm));
                if ~isfile(perm_file), continue; end
                load(perm_file,'results');
                R_phase_ch(ch,:)    = results.null_R_phase';
                R_MUA_ch(ch,:)      = results.null_R_MUA';
                R_Amp_ch(ch,:)      = results.null_R_Amp';
                R_AmpPhase_ch(ch,:) = results.null_R_AmpPhase';
                R_any_ch(ch,:)      = results.null_R_any';
            end

            null_avg_R_phase_freq(perm,:)    = mean(R_phase_ch,    1,'omitnan');
            null_avg_R_MUA_freq(perm,:)      = mean(R_MUA_ch,      1,'omitnan');
            null_avg_R_Amp_freq(perm,:)      = mean(R_Amp_ch,      1,'omitnan');
            null_avg_R_AmpPhase_freq(perm,:) = mean(R_AmpPhase_ch, 1,'omitnan');
            null_avg_R_any_freq(perm,:)      = mean(R_any_ch,      1,'omitnan');
        end

        reg_results.(depVarName).channel_avg_thresh.phase    = prctile(max(null_avg_R_phase_freq,   [],2), 100*(1-alpha));
        reg_results.(depVarName).channel_avg_thresh.MUA      = prctile(max(null_avg_R_MUA_freq,     [],2), 100*(1-alpha));
        reg_results.(depVarName).channel_avg_thresh.Amp      = prctile(max(null_avg_R_Amp_freq,     [],2), 100*(1-alpha));
        reg_results.(depVarName).channel_avg_thresh.AmpPhase = prctile(max(null_avg_R_AmpPhase_freq,[],2), 100*(1-alpha));
        reg_results.(depVarName).channel_avg_thresh.any      = prctile(max(null_avg_R_any_freq,     [],2), 100*(1-alpha));

        reg_results.(depVarName).channel_avg_R.phase    = mean(reg_results.(depVarName).R2_phase,    1,'omitnan');
        reg_results.(depVarName).channel_avg_R.MUA      = mean(reg_results.(depVarName).R2_MUA,      1,'omitnan');
        reg_results.(depVarName).channel_avg_R.Amp      = mean(reg_results.(depVarName).R2_Amp,      1,'omitnan');
        reg_results.(depVarName).channel_avg_R.AmpPhase = mean(reg_results.(depVarName).R2_AmpPhase, 1,'omitnan');
        reg_results.(depVarName).channel_avg_R.any      = mean(reg_results.(depVarName).R2_any,      1,'omitnan');

        chan_avg_save_dir = fullfile(results_folder,'perm_R_pos_diff',depVarName);
        if ~exist(chan_avg_save_dir,'dir'), mkdir(chan_avg_save_dir); end
        obs_avg = reg_results.(depVarName).channel_avg_R;
        save(fullfile(chan_avg_save_dir,'channel_avg_results.mat'), ...
            'null_avg_R_phase_freq','null_avg_R_MUA_freq','null_avg_R_Amp_freq', ...
            'null_avg_R_AmpPhase_freq','null_avg_R_any_freq','obs_avg','-v7.3');
        fprintf('Saved channel-average results for %s cross-animal analysis.\n', animalName);
    end

    cd(results_folder)
    save('multi_regression_channelwise_R2_abs_per_pos_diff.mat','reg_results','-v7.3');
    fprintf('Saved per-animal results for %s.\n', animalName);
end

%% Monkey-average

nAnimals = numel(animals);
R2_types = {'phase','MUA','Amp','AmpPhase','any'};

for d = 1:length(Y_vars)
    depVarName = Y_vars{d};
    fprintf('\n=== Monkey-average for %s ===\n', depVarName);

    obs_monkey  = struct(); perm_monkey = struct();
    for t = 1:numel(R2_types)
        obs_monkey.(R2_types{t})  = [];
        perm_monkey.(R2_types{t}) = [];
    end

    for a = 1:nAnimals
        avg_file = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ...
            ['results_' animals{a}],'multi_lin_reg','abs_per_pos_diff','cp10_till_100', ...
            'perm_R_pos_diff',depVarName,'channel_avg_results.mat');
        if ~isfile(avg_file)
            warning('Channel-average results not found for %s (%s).', depVarName, animals{a}); continue
        end
        tmp = load(avg_file);
        for t = 1:numel(R2_types)
            rtype = R2_types{t};
            obs_monkey.(rtype)  = cat(1,obs_monkey.(rtype),tmp.obs_avg.(rtype));
            perm_monkey.(rtype) = cat(3,perm_monkey.(rtype),tmp.(['null_avg_R_' rtype '_freq']));
        end
    end

    if size(obs_monkey.phase,1) < nAnimals
        warning('Not all animals have results for %s. Skipping.', depVarName); continue
    end

    monkey_avg_obs = struct(); tmax_monkey = struct();
    thresh_monkey  = struct(); p_monkey    = struct();

    for t = 1:numel(R2_types)
        rtype = R2_types{t};
        monkey_avg_obs.(rtype) = mean(obs_monkey.(rtype),1);
        perm_avg = mean(perm_monkey.(rtype),3);
        tmax_monkey.(rtype)   = max(perm_avg,[],2);
        thresh_monkey.(rtype) = quantile(tmax_monkey.(rtype),1-alpha);
        p_monkey.(rtype)      = mean(tmax_monkey.(rtype) >= max(monkey_avg_obs.(rtype)));
    end

    monkey_save_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/results_combined', ...
        'multi_lin_reg','abs_per_pos_diff','cp10_till_100',depVarName);
    if ~exist(monkey_save_dir,'dir'), mkdir(monkey_save_dir); end
    save(fullfile(monkey_save_dir,'monkey_avg_results.mat'), ...
        'monkey_avg_obs','tmax_monkey','thresh_monkey','p_monkey','obs_monkey','animals','-v7.3');

    fprintf('Saved monkey-average for %s\n', depVarName);
    fprintf('  Thresholds: phase=%.4f  MUA=%.4f  Amp=%.4f  AmpPhase=%.4f  any=%.4f\n', ...
        thresh_monkey.phase,thresh_monkey.MUA,thresh_monkey.Amp,thresh_monkey.AmpPhase,thresh_monkey.any);
    fprintf('  p-values:   phase=%.4f  MUA=%.4f  Amp=%.4f  AmpPhase=%.4f  any=%.4f\n', ...
        p_monkey.phase,p_monkey.MUA,p_monkey.Amp,p_monkey.AmpPhase,p_monkey.any);
end
