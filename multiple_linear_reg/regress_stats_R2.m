clearvars;
close all;
clc;

% Description:
% -------------

% Goal: Conduct channel-wise regression analysis relating independent variables
% to each dependent variable across the frequency spectrum. Statistical
% significance is evaluated using non-parametric permutation testing to control
% for family-wise error rate.
% Computes both per-channel and channel-average max-stat thresholds.
% F normalizes by residual variance and df and R normalizes by total variance (RSS0)

%% Dependencies
addpath /opt/fieldtrip_github/; ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/
clc

%% Settings
nPerm  = 1000;       % number of permutations
alpha  = 0.05;
Y_vars = {'RT','MUA_ERP_ampl_all','LFP_ERP_ampl_all','hit_miss'};
animals = {'hermes', 'klecks'};

%% Loop over animals

for a = 1:numel(animals)
    animalName = animals{a};
    fprintf('\n=== Processing %s ===\n', animalName);

    %% Load data
    info_folder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ...
        ['results_' animalName], 'multi_lin_reg', 'cp10_till_100');
    cd(info_folder)
    load('ph_all_sess.mat') % loads ph_comb

    ph_all  = ph_comb.phase_all; % trials x freqs x channels
    amp_all = ph_comb.amp_all;   % trials x freqs x channels
    [numTrials, numFreq, numCh] = size(ph_all);

    %% Initialize results
    reg_results = struct();

    %% Loop over dependent variables

    for d = 1:length(Y_vars)
        depVarName = Y_vars{d};
        fprintf('Processing %s — %s...\n', animalName, depVarName);

        isLogistic = strcmp(depVarName,'hit_miss');

        % Load dependent variable
        if isLogistic
            Y_all_orig = ph_comb.trialinfo(:,20);
            Y_all_orig(Y_all_orig==5) = 0;
            trlInfo_all = ph_comb.trialinfo;
        else
            Y_all_orig = ph_comb.(depVarName);
            switch depVarName
                case 'RT', trlInfo_all = ph_comb.RT_trialinfo;
                case 'MUA_ERP_ampl_all', trlInfo_all = ph_comb.MUA_ERP_trialinfo;
                case 'LFP_ERP_ampl_all', trlInfo_all = ph_comb.LFP_ERP_trialinfo;
            end
        end

        %% Pre-allocate with NaN so skipped channels are excluded by 'omitnan'
        reg_results.(depVarName).R2_phase    = NaN(numCh, numFreq);
        reg_results.(depVarName).R2_MUA      = NaN(numCh, numFreq);
        reg_results.(depVarName).R2_Amp      = NaN(numCh, numFreq);
        reg_results.(depVarName).R2_AmpPhase = NaN(numCh, numFreq);
        reg_results.(depVarName).R2_any      = NaN(numCh, numFreq);
        reg_results.(depVarName).R_phase     = NaN(numCh, numFreq);
        reg_results.(depVarName).phi_pref    = NaN(numCh, numFreq);

        %% Loop over channels
        for ch = 1:numCh
            fprintf('Channel %d/%d\n', ch, numCh);

            % Channel-specific DV
            if ~isLogistic
                Y_all = Y_all_orig(:,ch);
                if all(isnan(Y_all)), continue; end
            else
                Y_all = Y_all_orig;
            end

            %% Pre-compute clean data per frequency
            X_all_freq = cell(numFreq,1);
            Y_all_freq = cell(numFreq,1);
            Robs_phase    = zeros(numFreq,1);
            Robs_MUA      = zeros(numFreq,1);
            Robs_Amp      = zeros(numFreq,1);
            Robs_AmpPhase = zeros(numFreq,1);
            Robs_any      = zeros(numFreq,1);
            R_phase       = zeros(numFreq,1);
            phi_pref      = zeros(numFreq,1);

            for f = 1:numFreq
                % Design matrix: pupil, MUA baseline, amplitude, sin/cos phase
                X = [ph_comb.pup_baseline(:,ch), ph_comb.MUA_baseline(:,ch), ...
                     amp_all(:,f,ch), sin(ph_all(:,f,ch)), cos(ph_all(:,f,ch))];

                % RT: keep only hit trials
                if strcmp(depVarName,'RT')
                    keepIdx = trlInfo_all(:,20) == 1;
                else
                    keepIdx = true(size(Y_all));
                end

                nanIdx = any(isnan(X),2) | isnan(Y_all) | ~keepIdx;
                X_clean = X(~nanIdx,:);
                Y_clean = Y_all(~nanIdx);

                X_all_freq{f} = X_clean;
                Y_all_freq{f} = Y_clean;

                if isempty(Y_clean), continue; end

                %% Full model
                if isLogistic
                    b_full = glmfit(X_clean,Y_clean,'binomial','link','logit');
                else
                    X_full = [ones(size(X_clean,1),1), X_clean];
                    b_full = regress(Y_clean,X_full);
                    RSS_full = sum((Y_clean - X_full*b_full).^2);
                end

                %% Compute observed R² / pseudo-R²
                if ~isLogistic
                    % Null model
                    X_null = ones(size(X_clean,1),1);
                    b_null = regress(Y_clean,X_null);
                    RSS_null = sum((Y_clean - X_null*b_null).^2);

                    % Any predictor
                    Robs_any(f) = max(0,(RSS_null - RSS_full)/RSS_null);

                    % Phase
                    X_red = [ones(size(X_clean,1),1), X_clean(:,1:3)];
                    b_red = regress(Y_clean,X_red);
                    RSS_red = sum((Y_clean - X_red*b_red).^2);
                    Robs_phase(f) = max(0,(RSS_red - RSS_full)/RSS_null);

                    % MUA
                    X_red = [ones(size(X_clean,1),1), X_clean(:,[1 3 4 5])];
                    b_red = regress(Y_clean,X_red);
                    RSS_red = sum((Y_clean - X_red*b_red).^2);
                    Robs_MUA(f) = max(0,(RSS_red - RSS_full)/RSS_null);

                    % Amp
                    X_red = [ones(size(X_clean,1),1), X_clean(:,[1 2 4 5])];
                    b_red = regress(Y_clean,X_red);
                    RSS_red = sum((Y_clean - X_red*b_red).^2);
                    Robs_Amp(f) = max(0,(RSS_red - RSS_full)/RSS_null);

                    % Amp+Phase
                    X_red = [ones(size(X_clean,1),1), X_clean(:,1:2)];
                    b_red = regress(Y_clean,X_red);
                    RSS_red = sum((Y_clean - X_red*b_red).^2);
                    Robs_AmpPhase(f) = max(0,(RSS_red - RSS_full)/RSS_null);

                else
                    % Logistic: pseudo-R²
                    p_full = glmval(b_full,X_clean,'logit');
                    LL_full = sum(Y_clean.*log(p_full+eps) + (1-Y_clean).*log(1-p_full+eps));

                    % Any
                    b_null = glmfit(ones(size(Y_clean,1),1),Y_clean,'binomial','link','logit');
                    p_null = glmval(b_null,ones(size(Y_clean,1),1),'logit');
                    LL_null = sum(Y_clean.*log(p_null+eps) + (1-Y_clean).*log(1-p_null+eps));
                    Robs_any(f) = 1 - (LL_full/LL_null);

                    % Phase
                    b_red = glmfit(X_clean(:,1:3),Y_clean,'binomial','link','logit');
                    p_red = glmval(b_red,X_clean(:,1:3),'logit');
                    LL_red = sum(Y_clean.*log(p_red+eps) + (1-Y_clean).*log(1-p_red+eps));
                    Robs_phase(f) = 1 - (LL_full/LL_red);

                    % MUA
                    b_red = glmfit(X_clean(:,[1 3:5]),Y_clean,'binomial','link','logit');
                    p_red = glmval(b_red,X_clean(:,[1 3:5]),'logit');
                    LL_red = sum(Y_clean.*log(p_red+eps) + (1-Y_clean).*log(1-p_red+eps));
                    Robs_MUA(f) = 1 - (LL_full/LL_red);

                    % Amp
                    b_red = glmfit(X_clean(:,[1 2 4 5]),Y_clean,'binomial','link','logit');
                    p_red = glmval(b_red,X_clean(:,[1 2 4 5]),'logit');
                    LL_red = sum(Y_clean.*log(p_red+eps) + (1-Y_clean).*log(1-p_red+eps));
                    Robs_Amp(f) = 1 - (LL_full/LL_red);

                    % Amp+Phase
                    b_red = glmfit(X_clean(:,1:2),Y_clean,'binomial','link','logit');
                    p_red = glmval(b_red,X_clean(:,1:2),'logit');
                    LL_red = sum(Y_clean.*log(p_red+eps) + (1-Y_clean).*log(1-p_red+eps));
                    Robs_AmpPhase(f) = 1 - (LL_full/LL_red);
                end

                % Phase tuning amplitude & preferred phase (valid for both linear and logistic:
                % both regress() and glmfit() return intercept first, so end-1=beta_sin, end=beta_cos)
                beta_sin = b_full(end-1);
                beta_cos = b_full(end);
                R_phase(f)  = sqrt(beta_sin^2 + beta_cos^2);
                phi_pref(f) = atan2(beta_sin, beta_cos);
            end

            %% Permutation testing: SLURM call
            nJobs = ceil(nPerm/10);
            cfg_array = cell(nJobs,1);
            output_dir = fullfile(info_folder,'perm_R',depVarName,num2str(ch));

            for j = 1:nJobs
                perm_start = (j-1)*10 + 1;
                perm_end   = min(j*10, nPerm);
                cfg_array{j}.X_all_freq = X_all_freq;
                cfg_array{j}.Y_all_freq = Y_all_freq;
                cfg_array{j}.numFreq    = numFreq;
                cfg_array{j}.isLogistic = isLogistic;
                cfg_array{j}.perm_idx   = perm_start:perm_end;
                cfg_array{j}.output_dir = output_dir;
            end

%             slurmfun(@regress_perm_R, cfg_array, ...
%                 'partition','8GB','stopOnError',false,'useUserPath',true, ...
%                 'waitForToolboxes',{'statistics_toolbox'});

            %% Collect permutation nulls
            null_max_phase    = zeros(nPerm,1);
            null_max_MUA      = zeros(nPerm,1);
            null_max_Amp      = zeros(nPerm,1);
            null_max_AmpPhase = zeros(nPerm,1);
            null_max_any      = zeros(nPerm,1);

            for p = 1:nPerm
                perm_file = fullfile(output_dir, sprintf('perm_%04d.mat', p));
                load(perm_file,'results');
                null_max_phase(p)    = results.null_max_phase;
                null_max_MUA(p)      = results.null_max_MUA;
                null_max_Amp(p)      = results.null_max_Amp;
                null_max_AmpPhase(p) = results.null_max_AmpPhase;
                null_max_any(p)      = results.null_max_any;
            end

            %% Compute thresholds and p-values (channel-wise)
            thresholds.thresh_phase    = prctile(null_max_phase,100*(1-alpha));
            thresholds.thresh_MUA      = prctile(null_max_MUA,100*(1-alpha));
            thresholds.thresh_Amp      = prctile(null_max_Amp,100*(1-alpha));
            thresholds.thresh_AmpPhase = prctile(null_max_AmpPhase,100*(1-alpha));
            thresholds.thresh_any      = prctile(null_max_any,100*(1-alpha));

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

        %% Channel-average permutation null (max-stat across frequencies)
        fprintf('Computing channel-average permutation null for %s...\n', animalName);

        null_avg_R_phase    = zeros(nPerm,1);
        null_avg_R_MUA      = zeros(nPerm,1);
        null_avg_R_Amp      = zeros(nPerm,1);
        null_avg_R_AmpPhase = zeros(nPerm,1);
        null_avg_R_any      = zeros(nPerm,1);

        null_avg_R_phase_freq    = zeros(nPerm, numFreq);
        null_avg_R_MUA_freq      = zeros(nPerm, numFreq);
        null_avg_R_Amp_freq      = zeros(nPerm, numFreq);
        null_avg_R_AmpPhase_freq = zeros(nPerm, numFreq);
        null_avg_R_any_freq      = zeros(nPerm, numFreq);

        for perm = 1:nPerm
            % Initialize container for per-channel R
            R_phase_ch    = NaN(numCh, numFreq);
            R_MUA_ch      = NaN(numCh, numFreq);
            R_Amp_ch      = NaN(numCh, numFreq);
            R_AmpPhase_ch = NaN(numCh, numFreq);
            R_any_ch      = NaN(numCh, numFreq);

            % Load each channel's permutation result
            for ch = 1:numCh
                perm_file = fullfile(info_folder,'perm_R',depVarName,num2str(ch),...
                    sprintf('perm_%04d.mat',perm));
                if ~isfile(perm_file)
                    continue;
                end
                load(perm_file,'results');
                R_phase_ch(ch,:)    = results.null_R_phase';
                R_MUA_ch(ch,:)      = results.null_R_MUA';
                R_Amp_ch(ch,:)      = results.null_R_Amp';
                R_AmpPhase_ch(ch,:) = results.null_R_AmpPhase';
                R_any_ch(ch,:)      = results.null_R_any';
            end

            % Average over channels (dim 1), result is 1 x numFreq
            avg_R_phase    = mean(R_phase_ch,    1, 'omitnan');
            avg_R_MUA      = mean(R_MUA_ch,      1, 'omitnan');
            avg_R_Amp      = mean(R_Amp_ch,      1, 'omitnan');
            avg_R_AmpPhase = mean(R_AmpPhase_ch, 1, 'omitnan');
            avg_R_any      = mean(R_any_ch,      1, 'omitnan');

            % Max-stat across frequencies
            null_avg_R_phase(perm)    = max(avg_R_phase);
            null_avg_R_MUA(perm)      = max(avg_R_MUA);
            null_avg_R_Amp(perm)      = max(avg_R_Amp);
            null_avg_R_AmpPhase(perm) = max(avg_R_AmpPhase);
            null_avg_R_any(perm)      = max(avg_R_any);

            % Store frequency-resolved channel-average nulls (for cross-animal analysis)
            null_avg_R_phase_freq(perm,:)    = avg_R_phase;
            null_avg_R_MUA_freq(perm,:)      = avg_R_MUA;
            null_avg_R_Amp_freq(perm,:)      = avg_R_Amp;
            null_avg_R_AmpPhase_freq(perm,:) = avg_R_AmpPhase;
            null_avg_R_any_freq(perm,:)      = avg_R_any;
        end

        % Compute channel-average thresholds
        reg_results.(depVarName).channel_avg_thresh.phase    = prctile(null_avg_R_phase, 100*(1-alpha));
        reg_results.(depVarName).channel_avg_thresh.MUA      = prctile(null_avg_R_MUA,   100*(1-alpha));
        reg_results.(depVarName).channel_avg_thresh.Amp      = prctile(null_avg_R_Amp,   100*(1-alpha));
        reg_results.(depVarName).channel_avg_thresh.AmpPhase = prctile(null_avg_R_AmpPhase, 100*(1-alpha));
        reg_results.(depVarName).channel_avg_thresh.any      = prctile(null_avg_R_any,   100*(1-alpha));

        % Compute observed channel-average R²
        reg_results.(depVarName).channel_avg_R.phase    = mean(reg_results.(depVarName).R2_phase,1,'omitnan');
        reg_results.(depVarName).channel_avg_R.MUA      = mean(reg_results.(depVarName).R2_MUA,1,'omitnan');
        reg_results.(depVarName).channel_avg_R.Amp      = mean(reg_results.(depVarName).R2_Amp,1,'omitnan');
        reg_results.(depVarName).channel_avg_R.AmpPhase = mean(reg_results.(depVarName).R2_AmpPhase,1,'omitnan');
        reg_results.(depVarName).channel_avg_R.any      = mean(reg_results.(depVarName).R2_any,1,'omitnan');

        % Save channel-average results for cross-animal analysis
        chan_avg_save_dir = fullfile(info_folder, 'perm_R', depVarName);
        if ~exist(chan_avg_save_dir, 'dir'), mkdir(chan_avg_save_dir); end
        obs_avg = reg_results.(depVarName).channel_avg_R;
        save(fullfile(chan_avg_save_dir, 'channel_avg_results.mat'), ...
            'null_avg_R_phase_freq', 'null_avg_R_MUA_freq', 'null_avg_R_Amp_freq', ...
            'null_avg_R_AmpPhase_freq', 'null_avg_R_any_freq', 'obs_avg', '-v7.3');
        fprintf('Saved channel-average results for %s cross-animal analysis.\n', animalName);

    end

    % Save per-animal results
    cd(info_folder)
    save('multi_regression_channelwise_R2.mat','reg_results','-v7.3');
    fprintf('Saved per-animal results for %s.\n', animalName);

end  % animal loop

%% Monkey-average: combine across animals

nAnimals = numel(animals);
R2_types = {'phase', 'MUA', 'Amp', 'AmpPhase', 'any'};

for d = 1:length(Y_vars)
    depVarName = Y_vars{d};
    fprintf('\n=== Monkey-average for %s ===\n', depVarName);

    obs_monkey  = struct();
    perm_monkey = struct();
    for t = 1:numel(R2_types)
        obs_monkey.(R2_types{t})  = [];
        perm_monkey.(R2_types{t}) = [];
    end

    for a = 1:nAnimals
        animal_folder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ...
            ['results_' animals{a}], 'multi_lin_reg', 'cp10_till_100');
        avg_file = fullfile(animal_folder, 'perm_R', depVarName, 'channel_avg_results.mat');

        if ~isfile(avg_file)
            warning('Channel-average results not found for %s (%s). Run per-animal section first.', ...
                depVarName, animals{a});
            continue
        end

        tmp = load(avg_file);

        for t = 1:numel(R2_types)
            rtype = R2_types{t};
            obs_monkey.(rtype)  = cat(1, obs_monkey.(rtype),  tmp.obs_avg.(rtype));
            perm_monkey.(rtype) = cat(3, perm_monkey.(rtype), tmp.(['null_avg_R_' rtype '_freq']));
        end
    end

    % Check all animals loaded
    if size(obs_monkey.phase, 1) < nAnimals
        warning('Not all animals have results for %s. Skipping monkey-average.', depVarName);
        continue
    end

    monkey_avg_obs = struct();
    tmax_monkey    = struct();
    thresh_monkey  = struct();
    p_monkey       = struct();

    for t = 1:numel(R2_types)
        rtype = R2_types{t};

        % Average observed R² across animals
        monkey_avg_obs.(rtype) = mean(obs_monkey.(rtype), 1);

        % Average permutation null across animals, then max-stat
        perm_avg = mean(perm_monkey.(rtype), 3);
        tmax_monkey.(rtype)   = max(perm_avg, [], 2);
        thresh_monkey.(rtype) = quantile(tmax_monkey.(rtype), 1 - alpha);

        % p-value
        p_monkey.(rtype) = mean(tmax_monkey.(rtype) >= max(monkey_avg_obs.(rtype)));
    end

    % Save
    monkey_save_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/results_combined/multi_lin_reg/cp10_till_100', depVarName);
    if ~exist(monkey_save_dir, 'dir'), mkdir(monkey_save_dir); end
    save(fullfile(monkey_save_dir, 'monkey_avg_results.mat'), ...
        'monkey_avg_obs', 'tmax_monkey', 'thresh_monkey', 'p_monkey', ...
        'obs_monkey', 'animals', '-v7.3');

    fprintf('Saved monkey-average results for %s\n', depVarName);
    fprintf('  Thresholds: phase=%.4f  MUA=%.4f  Amp=%.4f  AmpPhase=%.4f  any=%.4f\n', ...
        thresh_monkey.phase, thresh_monkey.MUA, thresh_monkey.Amp, ...
        thresh_monkey.AmpPhase, thresh_monkey.any);
    fprintf('  p-values:   phase=%.4f  MUA=%.4f  Amp=%.4f  AmpPhase=%.4f  any=%.4f\n', ...
        p_monkey.phase, p_monkey.MUA, p_monkey.Amp, ...
        p_monkey.AmpPhase, p_monkey.any);
end
