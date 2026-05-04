% =====================================================================
% Multiple linear regression: phase predicts DV
% Hypothesis H2 (abs_per_pos/)
%
% Claim: each stimulus position has its own phase-DV relationship;
% positions are NOT required to share a preferred phase.
%
% Recipe: group trials by stimulus position (trialinfo col 16); fit
%   DV ~ pupil + MUA_baseline + amp + sin(φ) + cos(φ) within each
%   position per channel; R_phase = |β_cos + i·β_sin| and R² averaged
%   across positions (Way 2); arithmetic mean across channels/animals.
%   Permutation null: shuffle DV within position; max-stat correction.
%
% See sampling_compare/README.md for the Way-1 / Way-2 framing.
% =====================================================================
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
nPerm  = 1000;
alpha  = 0.05;
Y_vars  = {'RT','MUA_ERP_ampl_all','LFP_ERP_ampl_all','hit_miss'};
animals = {'hermes', 'klecks'};

%% Loop over animals

for a = 1:numel(animals)
    animalName = animals{a};
    fprintf('\n=== Processing %s ===\n', animalName);

    data_folder    = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ...
        ['results_' animalName], 'multi_lin_reg', 'cp10_till_100');
    results_folder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ...
        ['results_' animalName], 'multi_lin_reg', 'abs_per_pos', 'cp10_till_100');
    if ~exist(results_folder,'dir'), mkdir(results_folder); end
    cd(data_folder); load('ph_all_sess.mat')

    ph_all  = ph_comb.phase_all;
    amp_all = ph_comb.amp_all;
    [numTrials, numFreq, numCh] = size(ph_all);

    reg_results = struct();

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

        reg_results.(depVarName).R2_phase    = NaN(numCh, numFreq);
        reg_results.(depVarName).R2_MUA      = NaN(numCh, numFreq);
        reg_results.(depVarName).R2_Amp      = NaN(numCh, numFreq);
        reg_results.(depVarName).R2_AmpPhase = NaN(numCh, numFreq);
        reg_results.(depVarName).R2_any      = NaN(numCh, numFreq);
        reg_results.(depVarName).R_phase     = NaN(numCh, numFreq);
        reg_results.(depVarName).phi_pref    = NaN(numCh, numFreq);

        for ch = 1:numCh
            fprintf('Channel %d/%d\n', ch, numCh);

            if ~isLogistic
                Y_all = Y_all_orig(:,ch);
                if all(isnan(Y_all)), continue; end
            else
                Y_all = Y_all_orig;
            end

            %% Pre-compute per-position, per-frequency X and Y for real data and perms

            X_all_freq_pos = cell(nPos, numFreq);
            Y_all_freq_pos = cell(nPos, numFreq);

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

                R2_pos_phase    = NaN(nPos,1);
                R2_pos_MUA      = NaN(nPos,1);
                R2_pos_Amp      = NaN(nPos,1);
                R2_pos_AmpPhase = NaN(nPos,1);
                R2_pos_any      = NaN(nPos,1);
                beta_sin_pos    = NaN(nPos,1);
                beta_cos_pos    = NaN(nPos,1);

                for p = 1:nPos
                    pos_mask = trlInfo_all(:,16) == positions(p);
                    nanIdx   = any(isnan(X),2) | isnan(Y_all) | ~keepIdx | ~pos_mask;
                    X_clean  = X(~nanIdx,:);
                    Y_clean  = Y_all(~nanIdx);

                    X_all_freq_pos{p,f} = X_clean;
                    Y_all_freq_pos{p,f} = Y_clean;

                    if length(Y_clean) < size(X_clean,2)+1, continue; end

                    if isLogistic
                        b_full  = glmfit(X_clean,Y_clean,'binomial','link','logit');
                        p_full  = glmval(b_full,X_clean,'logit');
                        LL_full = sum(Y_clean.*log(p_full+eps)+(1-Y_clean).*log(1-p_full+eps));

                        b_null  = glmfit(ones(size(Y_clean,1),1),Y_clean,'binomial','link','logit');
                        p_null  = glmval(b_null,ones(size(Y_clean,1),1),'logit');
                        LL_null = sum(Y_clean.*log(p_null+eps)+(1-Y_clean).*log(1-p_null+eps));
                        R2_pos_any(p) = 1-(LL_full/LL_null);

                        b_red = glmfit(X_clean(:,1:3),Y_clean,'binomial','link','logit');
                        p_red = glmval(b_red,X_clean(:,1:3),'logit');
                        LL_red = sum(Y_clean.*log(p_red+eps)+(1-Y_clean).*log(1-p_red+eps));
                        R2_pos_phase(p) = 1-(LL_full/LL_red);

                        b_red = glmfit(X_clean(:,[1 3 4 5]),Y_clean,'binomial','link','logit');
                        p_red = glmval(b_red,X_clean(:,[1 3 4 5]),'logit');
                        LL_red = sum(Y_clean.*log(p_red+eps)+(1-Y_clean).*log(1-p_red+eps));
                        R2_pos_MUA(p) = 1-(LL_full/LL_red);

                        b_red = glmfit(X_clean(:,[1 2 4 5]),Y_clean,'binomial','link','logit');
                        p_red = glmval(b_red,X_clean(:,[1 2 4 5]),'logit');
                        LL_red = sum(Y_clean.*log(p_red+eps)+(1-Y_clean).*log(1-p_red+eps));
                        R2_pos_Amp(p) = 1-(LL_full/LL_red);

                        b_red = glmfit(X_clean(:,1:2),Y_clean,'binomial','link','logit');
                        p_red = glmval(b_red,X_clean(:,1:2),'logit');
                        LL_red = sum(Y_clean.*log(p_red+eps)+(1-Y_clean).*log(1-p_red+eps));
                        R2_pos_AmpPhase(p) = 1-(LL_full/LL_red);

                    else
                        X_full   = [ones(size(X_clean,1),1), X_clean];
                        b_full   = regress(Y_clean,X_full);
                        RSS_full = sum((Y_clean - X_full*b_full).^2);
                        RSS_null = sum((Y_clean - mean(Y_clean)).^2);

                        R2_pos_any(p) = max(0,(RSS_null-RSS_full)/RSS_null);

                        X_red = [ones(size(X_clean,1),1), X_clean(:,1:3)];
                        b_red = regress(Y_clean,X_red); RSS_red = sum((Y_clean-X_red*b_red).^2);
                        R2_pos_phase(p) = max(0,(RSS_red-RSS_full)/RSS_null);

                        X_red = [ones(size(X_clean,1),1), X_clean(:,[1 3 4 5])];
                        b_red = regress(Y_clean,X_red); RSS_red = sum((Y_clean-X_red*b_red).^2);
                        R2_pos_MUA(p) = max(0,(RSS_red-RSS_full)/RSS_null);

                        X_red = [ones(size(X_clean,1),1), X_clean(:,[1 2 4 5])];
                        b_red = regress(Y_clean,X_red); RSS_red = sum((Y_clean-X_red*b_red).^2);
                        R2_pos_Amp(p) = max(0,(RSS_red-RSS_full)/RSS_null);

                        X_red = [ones(size(X_clean,1),1), X_clean(:,1:2)];
                        b_red = regress(Y_clean,X_red); RSS_red = sum((Y_clean-X_red*b_red).^2);
                        R2_pos_AmpPhase(p) = max(0,(RSS_red-RSS_full)/RSS_null);
                    end

                    beta_sin_pos(p) = b_full(end-1);
                    beta_cos_pos(p) = b_full(end);
                end

                % Average R² across positions (Way 2 across positions —
                % each position's strength contributes regardless of direction).
                Robs_phase(f)    = mean(R2_pos_phase,    'omitnan');
                Robs_MUA(f)      = mean(R2_pos_MUA,      'omitnan');
                Robs_Amp(f)      = mean(R2_pos_Amp,      'omitnan');
                Robs_AmpPhase(f) = mean(R2_pos_AmpPhase, 'omitnan');
                Robs_any(f)      = mean(R2_pos_any,       'omitnan');

                % R_phase / phi_pref: Way 2 across positions — magnitude per
                % position, then arithmetic mean of magnitudes; angle per
                % position, then circular mean of angles. Positions are
                % allowed to disagree on preferred phase (H2 framing).
                R_phase_pos  = sqrt(beta_sin_pos.^2 + beta_cos_pos.^2);
                R_phase(f)   = mean(R_phase_pos, 'omitnan');

                phi_pref_pos = atan2(beta_sin_pos, beta_cos_pos);
                ok_pos       = ~isnan(phi_pref_pos);
                if any(ok_pos)
                    phi_pref(f) = angle(mean(exp(1i * phi_pref_pos(ok_pos))));
                else
                    phi_pref(f) = NaN;
                end
            end

            %% Permutation (SLURM)

            nJobs = ceil(nPerm/10);
            cfg_array = cell(nJobs,1);
            output_dir = fullfile(results_folder,'perm_R_pos',depVarName,num2str(ch));
            if ~exist(output_dir,'dir'), mkdir(output_dir); end
            for j = 1:nJobs
                perm_start = (j-1)*10 + 1;
                perm_end   = min(j*10, nPerm);
                cfg_array{j}.X_all_freq_pos = X_all_freq_pos;
                cfg_array{j}.Y_all_freq_pos = Y_all_freq_pos;
                cfg_array{j}.numFreq        = numFreq;
                cfg_array{j}.nPos           = nPos;
                cfg_array{j}.isLogistic     = isLogistic;
                cfg_array{j}.perm_idx       = perm_start:perm_end;
                cfg_array{j}.output_dir     = output_dir;
            end

%             slurmfun(@regress_perm_R_pos, cfg_array, ...
%                 'partition','8GB','stopOnError',false,'useUserPath',true, ...
%                 'waitForToolboxes',{'statistics_toolbox'});

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

        % H2 paired-test channel-average for R_phase = |complex β|.
        % Way 2 across positions: |β_pos| then mean across positions per
        % channel; then mean across channels. Stays NaN if perm files
        % lack betas (older runs, before update).
        null_avg_R_phase_mag_freq = nan(nPerm, numFreq);

        for perm = 1:nPerm
            R_phase_ch    = NaN(numCh, numFreq);
            R_MUA_ch      = NaN(numCh, numFreq);
            R_Amp_ch      = NaN(numCh, numFreq);
            R_AmpPhase_ch = NaN(numCh, numFreq);
            R_any_ch      = NaN(numCh, numFreq);
            R_phase_mag_ch = NaN(numCh, numFreq);

            for ch = 1:numCh
                perm_file = fullfile(results_folder,'perm_R_pos',depVarName,num2str(ch), ...
                    sprintf('perm_%04d.mat',perm));
                if ~isfile(perm_file), continue; end
                load(perm_file,'results');
                R_phase_ch(ch,:)    = results.null_R_phase';
                R_MUA_ch(ch,:)      = results.null_R_MUA';
                R_Amp_ch(ch,:)      = results.null_R_Amp';
                R_AmpPhase_ch(ch,:) = results.null_R_AmpPhase';
                R_any_ch(ch,:)      = results.null_R_any';
                if isfield(results,'null_b_sin_pos')
                    % Way 2 across positions: |β_pos| per pos, mean over pos.
                    R_pos_mag = sqrt(results.null_b_sin_pos.^2 + results.null_b_cos_pos.^2);
                    R_phase_mag_ch(ch,:) = mean(R_pos_mag, 1, 'omitnan');
                end
            end

            null_avg_R_phase_freq(perm,:)    = mean(R_phase_ch,    1,'omitnan');
            null_avg_R_MUA_freq(perm,:)      = mean(R_MUA_ch,      1,'omitnan');
            null_avg_R_Amp_freq(perm,:)      = mean(R_Amp_ch,      1,'omitnan');
            null_avg_R_AmpPhase_freq(perm,:) = mean(R_AmpPhase_ch, 1,'omitnan');
            null_avg_R_any_freq(perm,:)      = mean(R_any_ch,      1,'omitnan');
            null_avg_R_phase_mag_freq(perm,:) = mean(R_phase_mag_ch, 1, 'omitnan');
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

        % R_phase / phi_pref channel-average: Way 2 across channels (channels
        % may disagree on preferred phase — H4 framing layered on H2).
        R_phase_ch  = reg_results.(depVarName).R_phase;
        phi_pref_ch = reg_results.(depVarName).phi_pref;
        reg_results.(depVarName).channel_avg_R.R_phase  = mean(R_phase_ch, 1, 'omitnan');
        reg_results.(depVarName).channel_avg_R.phi_pref = angle(mean(exp(1i * phi_pref_ch), 1, 'omitnan'));

        chan_avg_save_dir = fullfile(results_folder,'perm_R_pos',depVarName);
        if ~exist(chan_avg_save_dir,'dir'), mkdir(chan_avg_save_dir); end
        obs_avg = reg_results.(depVarName).channel_avg_R;
        save(fullfile(chan_avg_save_dir,'channel_avg_results.mat'), ...
            'null_avg_R_phase_freq','null_avg_R_MUA_freq','null_avg_R_Amp_freq', ...
            'null_avg_R_AmpPhase_freq','null_avg_R_any_freq', ...
            'null_avg_R_phase_mag_freq','obs_avg','-v7.3');
        fprintf('Saved channel-average results for %s cross-animal analysis.\n', animalName);
    end

    cd(results_folder)
    save('multi_regression_channelwise_R2_abs_per_pos.mat','reg_results','-v7.3');
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
    obs_monkey.R_phase  = [];
    obs_monkey.phi_pref = [];

    perm_monkey_R_phase = [];   % [nPerm × nFreq × nAnimals] for paired test

    for a = 1:nAnimals
        avg_file = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ...
            ['results_' animals{a}],'multi_lin_reg','abs_per_pos','cp10_till_100','perm_R_pos',depVarName,'channel_avg_results.mat');
        if ~isfile(avg_file)
            warning('Channel-average results not found for %s (%s).', depVarName, animals{a}); continue
        end
        tmp = load(avg_file);
        for t = 1:numel(R2_types)
            rtype = R2_types{t};
            obs_monkey.(rtype)  = cat(1,obs_monkey.(rtype),tmp.obs_avg.(rtype));
            perm_monkey.(rtype) = cat(3,perm_monkey.(rtype),tmp.(['null_avg_R_' rtype '_freq']));
        end
        if isfield(tmp.obs_avg, 'R_phase')
            obs_monkey.R_phase  = cat(1, obs_monkey.R_phase,  tmp.obs_avg.R_phase);
            obs_monkey.phi_pref = cat(1, obs_monkey.phi_pref, tmp.obs_avg.phi_pref);
        end
        if isfield(tmp,'null_avg_R_phase_mag_freq')
            perm_monkey_R_phase = cat(3, perm_monkey_R_phase, tmp.null_avg_R_phase_mag_freq);
        end
    end

    if size(obs_monkey.phase,1) < nAnimals
        warning('Not all animals have results for %s. Skipping.', depVarName); continue
    end

    monkey_avg_obs   = struct(); tmax_monkey = struct();
    thresh_monkey    = struct(); p_monkey    = struct();
    perm_monkey_avg  = struct();

    for t = 1:numel(R2_types)
        rtype = R2_types{t};
        monkey_avg_obs.(rtype) = mean(obs_monkey.(rtype),1);
        perm_avg = mean(perm_monkey.(rtype),3);
        tmax_monkey.(rtype)   = max(perm_avg,[],2);
        thresh_monkey.(rtype) = quantile(tmax_monkey.(rtype),1-alpha);
        p_monkey.(rtype)      = mean(tmax_monkey.(rtype) >= max(monkey_avg_obs.(rtype)));

        % Persist per-perm null curve for paired test.
        perm_monkey_avg.(rtype) = perm_avg;
    end

    % R_phase / phi_pref monkey-average: Way 2 across animals (H2+H4 — animals
    % may disagree on preferred phase). Mean of magnitudes; circular mean of
    % angles.
    if size(obs_monkey.R_phase, 1) >= 1
        monkey_avg_obs.R_phase  = mean(obs_monkey.R_phase, 1, 'omitnan');
        monkey_avg_obs.phi_pref = angle(mean(exp(1i * obs_monkey.phi_pref), 1, 'omitnan'));
    end

    % R_phase paired null monkey-average: Way 2 across animals (mean of
    % per-animal channel-avg magnitude curves). Stays absent if perm
    % files lacked betas.
    if ~isempty(perm_monkey_R_phase)
        perm_monkey_avg.R_phase = mean(perm_monkey_R_phase, 3, 'omitnan');
    end

    monkey_save_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/results_combined', ...
        'multi_lin_reg','abs_per_pos','cp10_till_100',depVarName);
    if ~exist(monkey_save_dir,'dir'), mkdir(monkey_save_dir); end
    save(fullfile(monkey_save_dir,'monkey_avg_results.mat'), ...
        'monkey_avg_obs','tmax_monkey','thresh_monkey','p_monkey', ...
        'perm_monkey_avg','obs_monkey','animals','-v7.3');

    fprintf('Saved monkey-average for %s\n', depVarName);
    fprintf('  Thresholds: phase=%.4f  MUA=%.4f  Amp=%.4f  AmpPhase=%.4f  any=%.4f\n', ...
        thresh_monkey.phase,thresh_monkey.MUA,thresh_monkey.Amp,thresh_monkey.AmpPhase,thresh_monkey.any);
    fprintf('  p-values:   phase=%.4f  MUA=%.4f  Amp=%.4f  AmpPhase=%.4f  any=%.4f\n', ...
        p_monkey.phase,p_monkey.MUA,p_monkey.Amp,p_monkey.AmpPhase,p_monkey.any);
end
