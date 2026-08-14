clear all
close all
clc

% Code description:
% -----------------
% code to compute the average LFP ERP per condition and the condition
% differences, for BOTH animals (klecks and hermes) in one script.
% Original per-animal versions:
% differences between them
% if zscore_run =1,  data z-scoring takes place (LFP)
% Statistics:
% -----------------
% Permutation of Condition labels per session
% Each permutation -> 1 slurm job
% After running all permutations, it gathers their results, computes
% condition difference and applies pixel based multiple comparison

%% Specify paths

addpath /opt/fieldtrip_github/
ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/eye_data
clc

%% Define important variables

zscore_run = 0;     % z-score session (over session) and save this z-scored data that should
                    % be later on used for averaging across sessions -
                    % zscoring done already while finding the critical time
first_runperm = 1;  % first time to run this permutation, if not, loads the matrix
                    % containing permuted label indices
permut_n   = 1000;  % number of iterations to run
plotting  = 1;      % plots the condition difference of session with significance limits
plotERPs = 1;       % if plotERPs, it plots the ERPs per condition

%% Animal configuration
% Both animals run through the same code below; everything that differs
% between them lives here.

animals(1).name            = 'klecks';
animals(1).resultsfolder   = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks';
animals(1).remove_sessions = [2, 15, 16, 18];
animals(1).plotfolder      = '/mnt/hpc/projects/MWSampling/4Shivangi/Plots/Hitvsmiss/erpdiff_klecks/lfp';

animals(2).name            = 'hermes';
animals(2).resultsfolder   = '/mnt/hpc/projects/MWSampling/4Shivangi/results_hermes';
animals(2).remove_sessions = [1, 3, 21];
animals(2).plotfolder      = '/mnt/hpc/projects/MWSampling/4Shivangi/Plots/Hitvsmiss/erpdiff_hermes/lfp';

for ianimal = 1:length(animals)

animalName      = animals(ianimal).name;
datafolder      = animals(ianimal).resultsfolder;
remove_sessions = animals(ianimal).remove_sessions;
plotfolder      = animals(ianimal).plotfolder;
fprintf('\n=============== %s ===============\n', animalName)

% animals have different channel counts, so these must not carry over
clear perm_cond_max perm_cond_min

%% Paths for data

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

session_paths = [];
session_paths = cellfun(@(x) fullfile(datafolder,x), session_names, 'uniform',0);

session_paths_files = [];
session_paths_files = cellfun(@(x) fullfile(datafolder,x, 'clean_lfp.mat'), session_names, 'uniform',0);

output_paths = cellfun(@(x) fullfile(datafolder, x,'ERP_LFP'),session_names, 'uniform',0);

L = linspace(-1,0.4,1401);

%% load data per session, create permutation matrices of trial indices
% Steps:
% 1) load session
% 2) compute the real average ERP per condition
% 3) compute real condtion differences between hits and misses
% 4) permute condition labels and create a matrix containing the permuted
%    label indeces, save this matrix that can be then used with slurm

if first_runperm
    
    for isess = 1:length(session_names)
        
        % load LFP data
        ESIload(session_paths_files{isess});
        lfpTrials = clean_data;
        
        % Create structure contating trialinfo and  permutation matrix,
        % save this matrix that can then be used with slurm
        
        trial_perm_ind = [];
        trial_perm_ind.sessname = session_names{isess};
        trial_perm_ind.trialinfo = lfpTrials.trialinfo;
        
        trial_perm_ind.hit  = find(lfpTrials.trialinfo(:,20)==1);
        trial_perm_ind.miss = find(lfpTrials.trialinfo(:,20)==5);
        
        
        % permute condition labels
        trial_ind_vec = [];
        trial_ind_vec = [trial_perm_ind.hit ; trial_perm_ind.miss];
        
        % shuffle
        rand_matrix = nan(permut_n, length(trial_ind_vec));
        for iter = 1:permut_n
            rand_matrix(iter,:)= permutate(trial_ind_vec');
        end
        
        trial_perm_ind.rand_matrix = rand_matrix;
        
        % create session folder to save output if it does not exist
        % already
        if ~isdir(output_paths{isess})
            mkdir(output_paths{isess})
        end
        
        % save the structure containing all trialinfo - an permutation
        % matrix of this session
        cd(output_paths{isess}),
        save('trial_perm_ind','trial_perm_ind')
        
        % zscore data in order to be able to pool over sessions
        if zscore_run
            disp(strcat('session- ',num2str(isess),' out of- ',num2str(length(session_names)), ', running LFP z-scoring')) %#ok<UNRCH>
            zlfptrials = fun_zscore_session(lfpTrials);
            % save zscored data
            cd(session_paths{isess}),
            save('zlfptrials','zlfptrials')
            
        else
            % load z-scored data
            cd(session_paths{isess}),
            load('zlfptrials.mat')
            
        end
        
        clear lfpTrials
        
        %% Compute average response or real data per condition
        
        if ~isdir(fullfile(output_paths{isess},'ERP_real'))
            mkdir(fullfile(output_paths{isess},'ERP_real'))
        end
        
        % activation period ----
        % average hit response
        cfg = [];
        cfg.trials     = trial_perm_ind.hit';
        cfg.outputfile = fullfile(output_paths{isess},'ERP_real','hit_lfp_avg.mat');
        
        ft_timelockanalysis(cfg,zlfptrials)
        
        % average miss response
        cfg = [];
        cfg.trials     = trial_perm_ind.miss';
        cfg.outputfile = fullfile(output_paths{isess},'ERP_real','miss_lfp_avg.mat');
        
        ft_timelockanalysis(cfg,zlfptrials)
        
        
        % baseline period ----
        cfg = [];
        cfg.latency = [-0.4 0];
        cfg.trials     =  [trial_perm_ind.hit' trial_perm_ind.miss']; % include all hits and misses
        cfg.outputfile = fullfile(output_paths{isess},'ERP_real','bsl_lfp_avg.mat');
        
        ft_timelockanalysis(cfg,zlfptrials)
        
        clear zlfptrials trial_perm_ind
        
        
    end
    
    
    %% LFP - compute average PERMUTED responses with slurm
    
    for isess = 1:length(session_names)

        % load pre-computed matrix containing permuted condition labels of
        % session
        cd(output_paths{isess}),
        load('trial_perm_ind.mat')
        
        cfg = [];
        cfg = cell(1,permut_n);
        
        %% average permuted hits

        for iperm = 1:permut_n
            
            if ~isdir(fullfile(output_paths{isess},num2str(iperm)))
                mkdir(fullfile(output_paths{isess},num2str(iperm)))
            end
            
            cfg{iperm}.trials      = trial_perm_ind.rand_matrix(iperm,1:length(trial_perm_ind.hit));
            cfg{iperm}.inputfile   =  fullfile(session_paths{isess},'zlfptrials.mat');
            cfg{iperm}.outputfile  = fullfile(output_paths{isess},num2str(iperm),'hit_lfp_avg.mat');
            
        end
        
        %  ft_timelockanalysis(cfg{iperm})
        
        slurmfun(@ft_timelockanalysis, cfg, ...
            'partition',     '8GB', ...
            'stopOnError',   false,  ...
            'useUserPath',   true    );
        
        
        %% average permuted miss
        
        cfg = [];
        cfg = cell(1,permut_n);
        
        for iperm = 1:permut_n
            
            cfg{iperm}.trials      = trial_perm_ind.rand_matrix(iperm,length(trial_perm_ind.hit)+1:end);
            cfg{iperm}.inputfile   =  fullfile(session_paths{isess},'zlfptrials.mat');
            cfg{iperm}.outputfile  = fullfile(output_paths{isess},num2str(iperm),'miss_lfp_avg.mat');
            
        end
        
        
        %  ft_timelockanalysis(cfg{iperm})
        
        slurmfun(@ft_timelockanalysis, cfg, ...
            'partition',     '8GB', ...
            'stopOnError',   false,  ...
            'useUserPath',   true    );
        
    end
    
end

%% baseline normalization of real data
% Steps:
% 1) load average stim and average baseline responses
% 2) normalize real activity with average baseline
% 3) save average baseline per channel for each session

for isess = 1:length(session_names)
    isess
    hit_timelock  = [];
    miss_timelock = [];
    bsl_timelock  = [];
    
    % load average of real conditions
    cd(fullfile(output_paths{isess},'ERP_real')),
    hit_timelock  = load('hit_lfp_avg.mat');
    miss_timelock = load('miss_lfp_avg.mat');
    
    % load baseline
    bsl_timelock  = load('bsl_lfp_avg.mat');
    bsl_avg = mean(bsl_timelock.timelock.avg,2); % avg bsl value per channel
    
    cd(fullfile(output_paths{isess},'ERP_real')),
    save('bsl_avg','bsl_avg')
    
    % normalize hits
    norm_hit_timelock = [];
    norm_hit_timelock = hit_timelock.timelock;
    norm_hit_timelock.avg = hit_timelock.timelock.avg-bsl_avg;
    
    % save normalized hits
    cd(fullfile(output_paths{isess},'ERP_real')),
    save('norm_hit_timelock','norm_hit_timelock')
    
    % normalize misses
    norm_miss_timelock = [];
    norm_miss_timelock = miss_timelock.timelock;
    norm_miss_timelock.avg = miss_timelock.timelock.avg-bsl_avg;
    
    % save normalized misses
    cd(fullfile(output_paths{isess},'ERP_real')),
    save('norm_miss_timelock','norm_miss_timelock')
    
    clear norm_miss_timelock norm_hit_timelock
    
end

%% Condition Difference of Permuted trials
% Steps:
% 1) normalize average condition responses
% 2) compute condition difference
% 3) get max and min values per permutation and per channel, save these
%    distributions

for isess = 1:length(session_names)
    isess
    
    cd(fullfile(output_paths{isess},'ERP_real')),
    load('bsl_avg')
    
    for iperm = 1:permut_n
        
        cd(fullfile(output_paths{isess},num2str(iperm))),
        
        hit_perm = [];
        hit_perm = load('hit_lfp_avg.mat');
        norm_hit = hit_perm.timelock.avg-bsl_avg;
        
        clear hit_perm
        
        miss_perm = [];
        miss_perm = load('miss_lfp_avg.mat');
        norm_miss = miss_perm.timelock.avg-bsl_avg;
        
        clear miss_perm
        
        perm_cond_diff = [];
        perm_cond_diff = norm_hit-norm_miss;
        
        perm_cond_max(:,iperm) = max(perm_cond_diff,[],2);
        perm_cond_min(:,iperm) = min(perm_cond_diff,[],2);
        
        clear perm_cond_diff 
    end
    
    cd(fullfile(output_paths{isess})),
    save('perm_cond_max','perm_cond_max')
    save('perm_cond_min','perm_cond_min')
    
    clear perm_cond_min perm_cond_max
end

%% Plot real data difference
% Steps:
% 1) load average normalized responses of real data
% 2) compute the condition difference of real data and save it
% 3) load the permutation max-min distributions
% 4) define percentiles
% 5) if plotting = 1, plot difference and percentiles and save plot per
%    session

for isess = 1:length(session_names)

    cd(fullfile(output_paths{isess})),
    load('perm_cond_max.mat');
    load('perm_cond_min.mat');
    
    % load real data average
    cd(fullfile(output_paths{isess},'ERP_real')),
    load('norm_hit_timelock.mat');
    load('norm_miss_timelock.mat');
    
    real_cond_diff = norm_hit_timelock.avg-norm_miss_timelock.avg;
    
    % significance limits
    limit_max = quantile([perm_cond_max perm_cond_min], 0.975,2);
    limit_min = quantile([perm_cond_max perm_cond_min], 0.025,2);
    
    if plotting   
        if plotERPs
            fig = figure('Units', 'normalized', 'Position', [0, 0, 1, 1]);
            for ichan =1:length(norm_hit_timelock.label)
                sig_time = (real_cond_diff(ichan,:)>=limit_max(ichan,:) | real_cond_diff(ichan,:)<=limit_min(ichan,:));
                sigonset  = find(conv(sig_time,[1 -1])==1);
                sigoffset = find(conv(sig_time,[1 -1])==-1)-1;
                durplot=[L(sigonset); L(sigoffset)];
                
                subplot(8,8,ichan)
                plot(norm_hit_timelock.time,norm_hit_timelock.avg(ichan,:),'k')
                hold on
                plot(norm_miss_timelock.time,norm_miss_timelock.avg(ichan,:),'r')
                for i=1:size(durplot,2)
                    v = [durplot(1,i) -1.5; durplot(2,i) -1.5; durplot(2,i) 1; durplot(1,i) 1];
                    f = [1 2 3 4];
                    patch('Faces',f,'Vertices',v,'FaceColor','blue','FaceAlpha',.2,'EdgeColor','blue','EdgeAlpha',.2)
                end
                xlim([-0.1 0.15])
                ylim([-1.5 1])
                title(sprintf('Ch %s', char(norm_hit_timelock.label(ichan))))
            end
        end
    end
    
    % figure
    cd(plotfolder)
    savefig_filename = [session_names{isess} '.pdf'];
    set(fig, 'PaperPositionMode', 'auto')
    print(fig, savefig_filename, '-dpdf', '-fillpage');
end

%% pooling in all sessions
% Steps:
% 1) Remove the bad sessions and redefine output paths
% 2) load the max min values for permuations channel wise for each session
%    accounting for NaNs
% 3) compute the condition difference of real normalized data across for
%    all sessions channel wise and accounting for NaNs
% 4) Plot the pooled erps

% Remove unwanted sessions
keep_sessions = setdiff(1:length(session_names), remove_sessions);
session_namesp = session_names(keep_sessions);
output_pathsp = output_paths(keep_sessions);

% Load channels
cd(fullfile(datafolder,'critical_time'))
load('all_channels')
num_channels = length(all_channels);
num_permutations = permut_n;
num_sessions = length(session_namesp);

% Significance limits
% -------------------
% The observed statistic is the mean of the condition difference over
% sessions, so the null must be built the same way: average the PERMUTED
% session differences over sessions FIRST, then take the max and min over
% time. Pooling per-session maxima (the earlier version here) compares a
% session-averaged value against a single-session null and is far too
% conservative. This matches the construction already used in the
% correlation, coherence, regression and scanning pipelines.
[limit_max, limit_min, limit_max_all, limit_min_all] = pooled_perm_limits( ...
    output_pathsp, all_channels, num_permutations, length(L), 'hit_lfp_avg.mat', 'miss_lfp_avg.mat');

% Initialize ERP matrices
time_len = length(L);
ERPdiff_r = nan(num_sessions, num_channels, time_len);
ERP_hit = nan(num_channels, time_len, num_sessions);
ERP_miss = nan(num_channels, time_len, num_sessions);

% ERP data
for isess = 1:num_sessions
    cd(fullfile(output_pathsp{isess}, 'ERP_real'))
    load('norm_hit_timelock.mat');
    load('norm_miss_timelock.mat');

    hit_avg_mapped = nan(num_channels, time_len);
    miss_avg_mapped = nan(num_channels, time_len);

    for ichan = 1:length(norm_hit_timelock.label)
        chan_name = norm_hit_timelock.label{ichan};
        idx = find(strcmp(all_channels, chan_name));
        if ~isempty(idx)
            hit_avg_mapped(idx, :) = norm_hit_timelock.avg(ichan, :);
            miss_avg_mapped(idx, :) = norm_miss_timelock.avg(ichan, :);
        end
    end

    % Compute ERP difference, handling NaNs
    diff_mapped = nan(num_channels, time_len);
    for c = 1:num_channels
        h = hit_avg_mapped(c, :);
        m = miss_avg_mapped(c, :);
        valid = ~isnan(h) & ~isnan(m);
        diff = nan(1, time_len);
        diff(valid) = h(valid) - m(valid);
        diff_mapped(c, :) = diff;
    end

    ERPdiff_r(isess, :, :) = diff_mapped;

    if plotERPs
        ERP_hit(:, :, isess) = hit_avg_mapped;
        ERP_miss(:, :, isess) = miss_avg_mapped;
    end
end

% Average difference
ERPdiff_avg = reshape(nanmean(ERPdiff_r, 1), num_channels, time_len);

% Plotting
if plotting
    if plotERPs
        fig = figure('Units', 'normalized', 'Position', [0, 0, 1, 1]);
        for ichan = 1:num_channels
            sig_time = (ERPdiff_avg(ichan, :) >= limit_max(ichan)) | (ERPdiff_avg(ichan, :) <= limit_min(ichan));
            sigonset  = find(conv(sig_time, [1 -1]) == 1);
            sigoffset = find(conv(sig_time, [1 -1]) == -1) - 1;
            durplot = [L(sigonset); L(sigoffset)];

            subplot(8, 8, ichan)
            plot(norm_hit_timelock.time, nanmean(ERP_hit(ichan, :, :), 3), 'k')
            hold on
            plot(norm_miss_timelock.time, nanmean(ERP_miss(ichan, :, :), 3), 'r')
            for i = 1:size(durplot, 2)
                v = [durplot(1, i) -1; durplot(2, i) -1; durplot(2, i) 1; durplot(1, i) 1];
                patch('Faces', [1 2 3 4], 'Vertices', v, 'FaceColor', 'blue', 'FaceAlpha', .2, 'EdgeColor', 'blue', 'EdgeAlpha', .2)
            end
            xlim([-0.1 0.15])
            ylim([-1 1])
            title(sprintf('Ch %d', char(ichan)))
        end
    end

    % Save figure
    cd(plotfolder)
    savefig('ERPpooled_LFP')
    set(fig, 'PaperPositionMode', 'auto')
    print(fig, 'ERPpooled_LFP', '-dpdf', '-fillpage');
end

%% Pooling across all sessions and channels
% Steps:
% 1) Remove the bad sessions and redefine output paths
% 2) load the max min values for permuations channel wise for each session
%    accounting for NaNs
% 3) compute the condition difference of real normalized data across for
%    all sessions and channels ccounting for NaNs
% 4) Plot the pooled erps

% Remove unwanted sessions
keep_sessions = setdiff(1:length(session_names), remove_sessions);
session_namesp = session_names(keep_sessions);
output_pathsp = output_paths(keep_sessions);

% Load channels
cd(fullfile(datafolder,'critical_time'))
load('all_channels')
num_channels = length(all_channels);
num_permutations = permut_n;
num_sessions = length(session_namesp);
time_len = length(L);

% Initialize ERP matrices
ERPdiff_r = nan(num_sessions, num_channels, time_len);
ERP_hit = nan(num_channels, time_len, num_sessions);
ERP_miss = nan(num_channels, time_len, num_sessions);

% ERP data
for isess = 1:num_sessions
    cd(fullfile(output_pathsp{isess}, 'ERP_real'))
    load('norm_hit_timelock.mat');
    load('norm_miss_timelock.mat');

    hit_avg_mapped = nan(num_channels, time_len);
    miss_avg_mapped = nan(num_channels, time_len);

    for ichan = 1:length(norm_hit_timelock.label)
        chan_name = norm_hit_timelock.label{ichan};
        idx = find(strcmp(all_channels, chan_name));
        if ~isempty(idx)
            hit_avg_mapped(idx, :) = norm_hit_timelock.avg(ichan, :);
            miss_avg_mapped(idx, :) = norm_miss_timelock.avg(ichan, :);
        end
    end

    % Compute ERP difference, handling NaNs
    diff_mapped = nan(num_channels, time_len);
    for c = 1:num_channels
        h = hit_avg_mapped(c, :);
        m = miss_avg_mapped(c, :);
        valid = ~isnan(h) & ~isnan(m);
        diff = nan(1, time_len);
        diff(valid) = h(valid) - m(valid);
        diff_mapped(c, :) = diff;
    end

    ERPdiff_r(isess, :, :) = diff_mapped;

    if plotERPs
        ERP_hit(:, :, isess) = hit_avg_mapped;
        ERP_miss(:, :, isess) = miss_avg_mapped;
    end
end

L = norm_hit_timelock.time;  
time_len = length(L);

if plotting && plotERPs
    % Grand average ERPs
    ERP_hit_all = squeeze(nanmean(reshape(ERP_hit, num_channels, time_len, []), 1));  % time x 1
    ERP_miss_all = squeeze(nanmean(reshape(ERP_miss, num_channels, time_len, []), 1));
    
    % Difference
    ERPdiff_all_avg = nanmean(ERP_hit_all,2) - nanmean(ERP_miss_all,2);
    
    % Pooled permutation thresholds
    % Averaged over channels AND sessions, so the null is averaged the same
    % way before the max/min over time is taken.
    [~, ~, limit_max_all, limit_min_all] = pooled_perm_limits( ...
        output_pathsp, all_channels, num_permutations, time_len, 'hit_lfp_avg.mat', 'miss_lfp_avg.mat');
    
    % Identify significant timepoints
    sig_time = (ERPdiff_all_avg >= limit_max_all) | (ERPdiff_all_avg <= limit_min_all);
    sig_time = double(sig_time(:)');  
    
    sigonset  = find(conv(sig_time, [1 -1]) == 1);
    sigoffset = find(conv(sig_time, [1 -1]) == -1) - 1;
    durplot = [L(sigonset); L(sigoffset)];
    
    % Plot
    fig = figure('Name', 'Pooled ERP - All Channels & Sessions', 'Units', 'normalized', 'Position', [0.2, 0.2, 0.6, 0.5]);
    plot(norm_hit_timelock.time, nanmean(ERP_hit_all,2), 'k', 'LineWidth', 1.7); hold on;
    plot(norm_miss_timelock.time, nanmean(ERP_miss_all,2), 'r', 'LineWidth', 1.7);
    
    % Significance shading
    for i = 1:size(durplot, 2)
        v = [durplot(1, i) -1; durplot(2, i) -1; durplot(2, i) 1; durplot(1, i) 1];
        patch('Faces', [1 2 3 4], 'Vertices', v, 'FaceColor', 'blue', 'FaceAlpha', .2, 'EdgeColor', 'none')
    end
    
    xlabel('Time (s)');
    ylabel('Normalized ERP');
    title('Pooled ERP Across All Channels and Sessions');
    legend('Hit', 'Miss');
    xlim([-0.1 0.15]);
    ylim([-1 1]);
    grid on;
    cd(plotfolder)
    savefig('ERPpooled_GrandAvg');
    print('ERPpooled_GrandAvg', '-dpdf', '-fillpage');
end


end  % animal loop

%% ======================= pooled permutation null ======================= %%

function [limit_max, limit_min, limit_max_all, limit_min_all] = ...
    pooled_perm_limits(output_pathsp, all_channels, num_permutations, time_len, hitfile, missfile)
% Null distribution for the SESSION-AVERAGED condition difference.
%
% For every permutation the permuted session differences are averaged over
% sessions first, and only then reduced to a max and a min over time. That
% puts the null on the same scale as the observed session average.
%
% The common baseline cancels in the difference, (hit-bsl)-(miss-bsl) =
% hit-miss, so bsl_avg is not needed here.
%
% Reads the permuted averages written by the slurm jobs, because
% perm_cond_max / perm_cond_min have already collapsed over time and cannot
% be used to rebuild a session-averaged null. Results are cached so that
% calling this twice with the same session list costs one pass.

persistent cache_key cache_out

key = [strjoin(output_pathsp(:)', '|') '|' hitfile '|' num2str(num_permutations)];
if ~isempty(cache_key) && strcmp(key, cache_key)
    [limit_max, limit_min, limit_max_all, limit_min_all] = deal(cache_out{:});
    return
end

num_channels = length(all_channels);
num_sessions = length(output_pathsp);

perm_sum = zeros(num_permutations, num_channels, time_len, 'single');
perm_cnt = zeros(num_channels, 1);

for isess = 1:num_sessions

    fprintf('  pooled null: session %d/%d\n', isess, num_sessions)

    ref = load(fullfile(output_pathsp{isess}, 'ERP_real', hitfile));
    idx = nan(length(ref.timelock.label),1);
    for ichan = 1:length(ref.timelock.label)
        j = find(strcmp(all_channels, ref.timelock.label{ichan}));
        if ~isempty(j), idx(ichan) = j; end
    end
    ok = ~isnan(idx);
    if ~any(ok)
        warning('session %d: no channel matches all_channels', isess), continue
    end
    perm_cnt(idx(ok)) = perm_cnt(idx(ok)) + 1;

    for iperm = 1:num_permutations
        h = load(fullfile(output_pathsp{isess}, num2str(iperm), hitfile));
        m = load(fullfile(output_pathsp{isess}, num2str(iperm), missfile));
        d = h.timelock.avg - m.timelock.avg;
        perm_sum(iperm, idx(ok), :) = perm_sum(iperm, idx(ok), :) + ...
            reshape(single(d(ok,:)), 1, sum(ok), time_len);
    end
end

cnt = perm_cnt;
cnt(cnt == 0) = NaN;
perm_avg = perm_sum ./ reshape(single(cnt), [1 num_channels 1]);
clear perm_sum

% per channel
limit_max = nan(num_channels,1);
limit_min = nan(num_channels,1);
for ichan = 1:num_channels
    x = squeeze(perm_avg(:,ichan,:));
    if all(isnan(x(:))), continue, end
    dist = double([max(x,[],2); min(x,[],2)]);
    limit_max(ichan) = quantile(dist, 0.975);
    limit_min(ichan) = quantile(dist, 0.025);
end

% grand average over channels
pa   = squeeze(nanmean(perm_avg, 2));
dist = double([max(pa,[],2); min(pa,[],2)]);
limit_max_all = quantile(dist, 0.975);
limit_min_all = quantile(dist, 0.025);

cache_key = key;
cache_out = {limit_max, limit_min, limit_max_all, limit_min_all};
end
