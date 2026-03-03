clearvars;
close all;
clc

% Description:
% -------------
% Goal: Conduct channel-wise regression analysis relating independent variables
% to each dependent variable across the frequency spectrum. Statistical
% significance is evaluated using non-parametric permutation testing to control
% for family-wise error rate.
% F normalizes by residual variance and df and R normalizes by total variance (RSS0)

%% Dependencies

addpath /opt/fieldtrip_github/; ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/
clc

%% Settings

nPerm = 1000; % number of permutations
alpha  = 0.05;
Y_vars = {'RT','MUA_ERP_ampl_all','LFP_ERP_ampl_all','hit_miss'};

%% Load data

info_folder = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks/multi_lin_reg/cp10_till_100';
cd(info_folder)
load('ph_all_sess.mat') % loads ph_comb

ph_all  = ph_comb.phase_all; % trials x freqs x channels
amp_all = ph_comb.amp_all;   % trials x freqs x channels
[numTrials, numFreq, numCh] = size(ph_all);

%% Loop over dependent variables

reg_results = struct();

for d = 1:length(Y_vars)
    
    depVarName = Y_vars{d};
    fprintf('Processing %s...\n', depVarName);
    
    isLogistic = strcmp(depVarName,'hit_miss');
    
    % Load DV
    if isLogistic
        Y_all_orig = ph_comb.trialinfo(:,20);
        Y_all_orig(Y_all_orig==1) = 1;
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
    
    %% Loop over channels
    
    for ch = 1:numCh
        fprintf('Channel %d/%d\n', ch, numCh);
        
        % Select channel-specific DV
        if ~isLogistic
            Y_all = Y_all_orig(:,ch);
            if all(isnan(Y_all)), continue; end
        else
            Y_all = Y_all_orig;
        end
        
        %% Pre-compute clean data for all frequencies
        
        X_all_freq = cell(numFreq, 1);
        Y_all_freq = cell(numFreq, 1);
        
        obs_F_phase    = zeros(numFreq,1);
        obs_F_MUA      = zeros(numFreq,1);
        obs_F_Amp      = zeros(numFreq,1);
        obs_F_AmpPhase = zeros(numFreq,1);
        obs_F_any      = zeros(numFreq,1);
        R_phase        = zeros(numFreq,1);
        phi_pref       = zeros(numFreq,1);
        
        for f = 1:numFreq
            % ---------------- Construct design matrix ----------------
            X = [ph_comb.pup_baseline(:,ch), ...
                ph_comb.MUA_baseline(:,ch), ...
                amp_all(:,f,ch), ...
                sin(ph_all(:,f,ch)), cos(ph_all(:,f,ch))];
            
            % RT: keep only hit trials
            if strcmp(depVarName,'RT')
                keepIdx = trlInfo_all(:,20) == 1;
            else
                keepIdx = true(size(Y_all));
            end
            
            % Remove NaNs
            nanIdx = any(isnan(X),2) | isnan(Y_all) | ~keepIdx;
            X_clean = X(~nanIdx,:);
            Y_clean = Y_all(~nanIdx);
            
            % Store clean data
            X_all_freq{f} = X_clean;
            Y_all_freq{f} = Y_clean;
            
            if isempty(Y_clean), continue; end
            
            % ---------------- Full model regression ----------------
            if isLogistic
                b_full = glmfit(X_clean, Y_clean, 'binomial', 'link', 'logit');
            else
                X_full = [ones(size(X_clean,1),1), X_clean];
                b_full = regress(Y_clean, X_full);
            end
            
            % ---------------- Compute observed stats ----------------
            if ~isLogistic
                X_full = [ones(size(X_clean,1),1), X_clean];
                RSS_full = sum((Y_clean - X_full*b_full).^2); %residual sum of squares
                df2 = size(X_clean,1) - size(X_clean,2) - 1; %degree of freedom
                
                % Phase (sin+cos)
                X_red = [ones(size(X_clean,1),1), X_clean(:,1:3)];
                b_red = regress(Y_clean, X_red);
                RSS_red = sum((Y_clean - X_red*b_red).^2);
                df1 = 2;
                obs_F_phase(f) = ((RSS_red - RSS_full)/df1) / (RSS_full/df2);
                
                % MUA_baseline
                X_red = [ones(size(X_clean,1),1), X_clean(:,[1 3 4 5])];
                b_red = regress(Y_clean, X_red);
                RSS_red = sum((Y_clean - X_red*b_red).^2);
                df1 = 1;
                obs_F_MUA(f) = ((RSS_red - RSS_full)/df1) / (RSS_full/df2);
                
                % Amp
                X_red = [ones(size(X_clean,1),1), X_clean(:,3)];  % only amplitude removed
                b_red = regress(Y_clean, X_red);
                RSS_red = sum((Y_clean - X_red*b_red).^2);
                df1 = 1;
                obs_F_Amp(f) = ((RSS_red - RSS_full)/df1) / (RSS_full/df2);
                
                % Amp+Phase together
                X_red = [ones(size(X_clean,1),1), X_clean(:,1:2)];
                b_red = regress(Y_clean, X_red);
                RSS_red = sum((Y_clean - X_red*b_red).^2);
                df1 = 3;
                obs_F_AmpPhase(f) = ((RSS_red - RSS_full)/df1) / (RSS_full/df2);
                
                % Any predictor
                X_null = ones(size(X_clean,1),1);
                b_null = regress(Y_clean, X_null);
                RSS_null = sum((Y_clean - X_null*b_null).^2);
                df1 = size(X_clean,2);
                obs_F_any(f) = ((RSS_null - RSS_full)/df1) / (RSS_full/df2);
            else
                % Logistic: deviance differences
                y_hat_full = glmval(b_full, X_clean, 'logit');
                dev_full = -2 * sum(Y_clean.*log(y_hat_full + eps) + ...
                    (1-Y_clean).*log(1 - y_hat_full + eps));
                
                % Phase
                X_red = X_clean(:,1:3); % remove phase
                b_red = glmfit(X_red, Y_clean, 'binomial', 'link', 'logit');
                y_hat_red = glmval(b_red, X_red, 'logit');
                dev_red = -2 * sum(Y_clean.*log(y_hat_red + eps) + ...
                    (1-Y_clean).*log(1 - y_hat_red + eps));
                obs_F_phase(f) = dev_red - dev_full;
                
                % MUA_baseline
                X_red = X_clean(:,[1 3:5]);
                b_red = glmfit(X_red, Y_clean, 'binomial', 'link', 'logit');
                y_hat_red = glmval(b_red, X_clean(:,[1 3:5]), 'logit');
                dev_red = -2 * sum(Y_clean.*log(y_hat_red + eps) + ...
                    (1-Y_clean).*log(1 - y_hat_red + eps));
                obs_F_MUA(f) = dev_red - dev_full;
                
                % Amp
                X_red = X_clean(:,3);
                b_red = glmfit(X_red, Y_clean, 'binomial', 'link', 'logit');
                y_hat_red = glmval(b_red, X_red, 'logit');
                dev_red = -2 * sum(Y_clean.*log(y_hat_red + eps) + ...
                    (1-Y_clean).*log(1 - y_hat_red + eps));
                obs_F_Amp(f) = dev_red - dev_full;
                
                % Amp+Phase together
                X_red = X_clean(:,1:2);
                b_red = glmfit(X_red, Y_clean, 'binomial', 'link', 'logit');
                y_hat_red = glmval(b_red, X_red, 'logit');
                dev_red = -2 * sum(Y_clean.*log(y_hat_red + eps) + ...
                    (1-Y_clean).*log(1 - y_hat_red + eps));
                obs_F_AmpPhase(f) = dev_red - dev_full;
                
                % Any predictor
                b_null = glmfit(ones(size(Y_clean,1),1), Y_clean, 'binomial', 'link', 'logit');
                y_hat_null = glmval(b_null, ones(size(Y_clean)), 'logit');
                dev_null = -2 * sum(Y_clean.*log(y_hat_null + eps) + ...
                    (1-Y_clean).*log(1 - y_hat_null + eps));
                obs_F_any(f) = dev_null - dev_full;
            end
            
            % Phase amplitude & preferred phase
            beta_sin = b_full(end-1);
            beta_cos = b_full(end);
            R_phase(f) = sqrt(beta_sin^2 + beta_cos^2);
            phi_pref(f) = atan2(beta_sin, beta_cos);
        end
        
        %% Permutation testing (max-stat)
        
        % Prepare configuration for all permutations
        
        nPerms = 1000;
        permsPerJob = 10;
        nJobs = ceil(nPerms / permsPerJob);
        
        cfg_array = cell(nJobs, 1);
        
        for j = 1:nJobs
            perm_start = (j-1)*permsPerJob + 1;
            perm_end   = min(j*permsPerJob, nPerms);
            
            cfg_array{j}.X_all_freq = X_all_freq;
            cfg_array{j}.Y_all_freq = Y_all_freq;
            cfg_array{j}.numFreq = numFreq;
            cfg_array{j}.isLogistic = isLogistic;
            
            % pass a vector of permutation indices
            cfg_array{j}.perm_idx = perm_start:perm_end;
            
            cfg_array{j}.output_dir = fullfile(info_folder, 'perm', depVarName, num2str(ch));
        end
        
        slurmfun(@regress_perm, cfg_array, ...
            'partition', '8GB', ...
            'stopOnError', false, ...
            'useUserPath', true, ...
            'waitForToolboxes', {'statistics_toolbox'});
        
    
    %% Collect permutation results from SLURM output
    
    null_max_phase    = zeros(nPerms, 1);
    null_max_MUA      = zeros(nPerms, 1);
    null_max_Amp      = zeros(nPerms, 1);
    null_max_AmpPhase = zeros(nPerms, 1);
    null_max_any      = zeros(nPerms, 1);
    
    % Load results from each permutation
    output_dir = fullfile(info_folder, 'perm', depVarName, num2str(ch));
    for p = 1:nPerms
        perm_file = fullfile(output_dir, sprintf('perm_%04d.mat', p));
        if exist(perm_file, 'file')
            load(perm_file, 'results');
            null_max_phase(p)    = results.null_max_phase;
            null_max_MUA(p)      = results.null_max_MUA;
            null_max_Amp(p)      = results.null_max_Amp;
            null_max_AmpPhase(p) = results.null_max_AmpPhase;
            null_max_any(p)      = results.null_max_any;
        else
            warning('Missing permutation file: %s', perm_file);
        end
    end
    
    %% Compute corrected p-values and thresholds
    
    alpha = 0.05; % desired significance level
    
    % FWER-corrected thresholds (95th percentile of null max-stat distribution)
    thresh_phase    = prctile(null_max_phase,    100 * (1 - alpha));
    thresh_MUA      = prctile(null_max_MUA,      100 * (1 - alpha));
    thresh_Amp      = prctile(null_max_Amp,      100 * (1 - alpha));
    thresh_AmpPhase = prctile(null_max_AmpPhase, 100 * (1 - alpha));
    thresh_any      = prctile(null_max_any,      100 * (1 - alpha));
    
    % FWER-corrected p-values (proportion of null max-stats >= observed max-stat)
    p_phase    = mean(null_max_phase    >= max(obs_F_phase));
    p_MUA      = mean(null_max_MUA      >= max(obs_F_MUA));
    p_Amp      = mean(null_max_Amp      >= max(obs_F_Amp));
    p_AmpPhase = mean(null_max_AmpPhase >= max(obs_F_AmpPhase));
    p_any      = mean(null_max_any      >= max(obs_F_any));
    
    %% Store results
    
    FWER_thresholds = struct('thresh_phase',    thresh_phase, ...
        'thresh_MUA',      thresh_MUA, ...
        'thresh_Amp',      thresh_Amp, ...
        'thresh_AmpPhase', thresh_AmpPhase, ...
        'thresh_any',      thresh_any);
    
    FWER_pvals = struct('p_phase',    p_phase, ...
        'p_MUA',      p_MUA, ...
        'p_Amp',      p_Amp, ...
        'p_AmpPhase', p_AmpPhase, ...
        'p_any',      p_any);
    
    % Store in results structure
    reg_results.(depVarName).R_phase(ch, :)        = R_phase;
    reg_results.(depVarName).phi_pref(ch, :)       = phi_pref;
    reg_results.(depVarName).obs_F_phase(ch, :)    = obs_F_phase;
    reg_results.(depVarName).obs_F_MUA(ch, :)      = obs_F_MUA;
    reg_results.(depVarName).obs_F_Amp(ch, :)      = obs_F_Amp;
    reg_results.(depVarName).obs_F_AmpPhase(ch, :) = obs_F_AmpPhase;
    reg_results.(depVarName).obs_F_any(ch, :)      = obs_F_any;
    reg_results.(depVarName).pvals(ch)             = FWER_pvals;
    reg_results.(depVarName).FWER_thresholds(ch)   = FWER_thresholds;
    
    % Store full null distributions
    reg_results.(depVarName).null_distributions.null_max_phase    = null_max_phase;
    reg_results.(depVarName).null_distributions.null_max_MUA      = null_max_MUA;
    reg_results.(depVarName).null_distributions.null_max_Amp      = null_max_Amp;
    reg_results.(depVarName).null_distributions.null_max_AmpPhase = null_max_AmpPhase;
    reg_results.(depVarName).null_distributions.null_max_any      = null_max_any;
    
    end
end

% Save results
cd(info_folder)
save('multi_regression_perm_maxstat.mat','reg_results','-v7.3')
