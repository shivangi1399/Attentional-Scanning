%% per-location ERP difference pipeline
clear all
close all
clc

%% Paths

addpath /opt/fieldtrip_github/
ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/eye_data
clc

%% Data paths

datafolder   = '/mnt/hpc/projects/MWSampling/4Shivangi/results_hermes';
cd(datafolder)

animalName = 'hermes';
temp = dir;
session_names = {};
ii = 0;
for i = 1:length(temp)
    if ~isempty(strfind(temp(i).name,animalName))
        ii = ii+1;
        session_names{ii,1} = temp(i).name;
    end
end

session_paths = cellfun(@(x) fullfile(datafolder,x), session_names, 'uniform',0);
session_paths_files = cellfun(@(x) fullfile(datafolder,x, 'clean_lfp.mat'), session_names, 'uniform',0);
output_paths = cellfun(@(x) fullfile(datafolder, x,'ERP_LFP','per_loc'),session_names, 'uniform',0);

% default time vector (if needed)
L = linspace(-1,0.4,1401);

%% Configuration

zscore_run   = 0;      % 1 to zscore sessions and save zscored data
first_runperm= 1;      % 1 to compute trial perm matrices and send slurm jobs
permut_n     = 1000;   % number of permutations
plotting     = 1;
plotERPs     = 1;

%%  Create permutation matrices (per session) and compute real ERPs per location

if first_runperm
    
    for isess = 1:length(session_names)
        fprintf('Session %d / %d\n', isess, length(session_names))
        
        % load LFP data
        ESIload(session_paths_files{isess});
        lfpTrials = clean_data; clear clean_data
        
        % collect unique locations from column 16
        locations = unique(lfpTrials.trialinfo(:,16));
        locations(isnan(locations)) = [];
        
        % prepare output folder
        if ~isdir(output_paths{isess})
            mkdir(output_paths{isess})
        end
        
        % Make the structure
        trial_perm_ind = [];
        trial_perm_ind.sessname = session_names{isess};
        trial_perm_ind.trialinfo = lfpTrials.trialinfo;
        trial_perm_ind.locations = locations;
        trial_perm_ind.loc = struct([]);
        
        for iloc = 1:length(locations)
            loc = locations(iloc);
            trial_perm_ind.loc(iloc).locationID = loc;
            
            % find hits and misses at this location
            hit_idx  = find(lfpTrials.trialinfo(:,16)==loc & lfpTrials.trialinfo(:,20)==1);
            miss_idx = find(lfpTrials.trialinfo(:,16)==loc & lfpTrials.trialinfo(:,20)==5);
            
            trial_perm_ind.loc(iloc).hit = hit_idx;
            trial_perm_ind.loc(iloc).miss = miss_idx;
            
            % create combined trial vector for permutation only if both exist
            trial_ind_vec = [hit_idx; miss_idx];
            if isempty(trial_ind_vec)
                trial_perm_ind.loc(iloc).rand_matrix = [];
                continue
            end
            
            % shuffle many times
            rand_matrix = nan(permut_n, length(trial_ind_vec));
            for iter = 1:permut_n
                rand_matrix(iter,:) = permutate(trial_ind_vec');
            end
            trial_perm_ind.loc(iloc).rand_matrix = rand_matrix;
        end
        
        % save trial perm structure
        cd(output_paths{isess})
        save('trial_perm_ind','trial_perm_ind')
        
        % zscore (or load zscored data)
        if zscore_run
            disp(sprintf('session %d: running LFP z-scoring', isess))
            zlfptrials = fun_zscore_session(lfpTrials);
            cd(session_paths{isess})
            save('zlfptrials','zlfptrials')
            clear zlfptrials
        else
            % try to load precomputed z-scored data; if not present, use original lfpTrials
            cd(session_paths{isess})
            if exist('zlfptrials.mat','file')
                load('zlfptrials.mat')
            else
                zlfptrials = lfpTrials;
            end
        end
        
        % compute real ERPs per location
        for iloc = 1:length(trial_perm_ind.loc)
            locid = trial_perm_ind.loc(iloc).locationID;
            hit_trials = trial_perm_ind.loc(iloc).hit;
            miss_trials = trial_perm_ind.loc(iloc).miss;
            
            % skip if no trials
            if isempty(hit_trials) || isempty(miss_trials)
                fprintf('Session %s - location %g: skipping (no hits or no misses)\n', session_names{isess}, locid)
                continue
            end
            
            locERPdir = fullfile(output_paths{isess},sprintf('loc%d',locid),'ERP_real');
            if ~isdir(locERPdir), mkdir(locERPdir); end
            
            % average hit response
            cfg = [];
            cfg.trials = hit_trials';
            cfg.outputfile = fullfile(locERPdir,'hit_lfp_avg.mat');
            ft_timelockanalysis(cfg, zlfptrials)
            
            % average miss response
            cfg = [];
            cfg.trials = miss_trials';
            cfg.outputfile = fullfile(locERPdir,'miss_lfp_avg.mat');
            ft_timelockanalysis(cfg, zlfptrials)
            
            % baseline period average (use both hit+miss)
            cfg = [];
            cfg.latency = [-0.4 0];
            cfg.trials = [hit_trials' miss_trials'];
            cfg.outputfile = fullfile(locERPdir,'bsl_lfp_avg.mat');
            ft_timelockanalysis(cfg, zlfptrials)
        end
        
        clear lfpTrials zlfptrials trial_perm_ind
    end
    
    %%  Submit ft_timelockanalysis jobs for permuted averages (per session x location)

    for isess = 1:length(session_names)
        fprintf('Submitting perm jobs for session %d/%d\n', isess, length(session_names))
        cd(output_paths{isess})
        load('trial_perm_ind.mat')  
        
        % load zlfptrials path presence check
        zfile = fullfile(session_paths{isess},'zlfptrials.mat');
        if ~exist(zfile,'file')
            warning('zlfptrials.mat not found in %s; using clean_data (if available)', session_paths{isess})
        end
        
        for iloc = 1:length(trial_perm_ind.loc)
            locid = trial_perm_ind.loc(iloc).locationID;
            hit_trials = trial_perm_ind.loc(iloc).hit;
            miss_trials = trial_perm_ind.loc(iloc).miss;
            rand_matrix = trial_perm_ind.loc(iloc).rand_matrix;
            
            if isempty(rand_matrix)
                fprintf('Session %s loc %g: no rand_matrix -> skipping permutations\n', session_names{isess}, locid)
                continue
            end
            
            % for hits (first half of rand matrix columns)
            cfg = cell(1,permut_n);
            for iperm = 1:permut_n
                outdir = fullfile(output_paths{isess}, sprintf('loc%d',locid), num2str(iperm));
                if ~isdir(outdir), mkdir(outdir); end
                cfg{iperm}.trials = rand_matrix(iperm, 1:length(hit_trials));
                cfg{iperm}.inputfile = fullfile(session_paths{isess},'zlfptrials.mat');
                cfg{iperm}.outputfile = fullfile(outdir,'hit_lfp_avg.mat');
            end
            slurmfun(@ft_timelockanalysis, cfg, 'partition', '8GB', 'stopOnError', false, 'useUserPath', true);
            
            % for misses (second half)
            cfg = cell(1,permut_n);
            for iperm = 1:permut_n
                outdir = fullfile(output_paths{isess}, sprintf('loc%d',locid), num2str(iperm));
                if ~isdir(outdir), mkdir(outdir); end
                cfg{iperm}.trials = rand_matrix(iperm, length(hit_trials)+1:end);
                cfg{iperm}.inputfile = fullfile(session_paths{isess},'zlfptrials.mat');
                cfg{iperm}.outputfile = fullfile(outdir,'miss_lfp_avg.mat');
            end
            slurmfun(@ft_timelockanalysis, cfg, 'partition', '8GB', 'stopOnError', false, 'useUserPath', true);
        end
    end
end % first_runperm

%% Baseline normalization of real data (per session x location)

for isess = 1:length(session_names)
    fprintf('Baseline normalization session %d/%d\n', isess, length(session_names))
    cd(fullfile(output_paths{isess}))
    if ~exist('trial_perm_ind.mat','file')
        warning('trial_perm_ind.mat missing for session %s - skipping', session_names{isess})
        continue
    end
    load('trial_perm_ind.mat')  % trial_perm_ind
    if ~isfield(trial_perm_ind,'loc'), continue; end
    
    for iloc = 1:length(trial_perm_ind.loc)
        locid = trial_perm_ind.loc(iloc).locationID;
        locERPdir = fullfile(output_paths{isess}, sprintf('loc%d',locid), 'ERP_real');
        if ~isdir(locERPdir), continue; end
        
        % check required files
        hitfile = fullfile(locERPdir,'hit_lfp_avg.mat');
        missfile = fullfile(locERPdir,'miss_lfp_avg.mat');
        bslfile = fullfile(locERPdir,'bsl_lfp_avg.mat');
        if ~exist(hitfile,'file') || ~exist(missfile,'file') || ~exist(bslfile,'file')
            warning('Missing avg files for session %s loc %g - skipping norm', session_names{isess}, locid)
            continue
        end
        
        cd(locERPdir)
        hit_timelock = load('hit_lfp_avg.mat');      % saved as struct with .timelock
        miss_timelock = load('miss_lfp_avg.mat');
        bsl_timelock = load('bsl_lfp_avg.mat');
        
        bsl_avg = mean(bsl_timelock.timelock.avg, 2);
        save('bsl_avg','bsl_avg')
        
        norm_hit_timelock = hit_timelock.timelock;
        norm_hit_timelock.avg = hit_timelock.timelock.avg - bsl_avg;
        save('norm_hit_timelock','norm_hit_timelock')
        
        norm_miss_timelock = miss_timelock.timelock;
        norm_miss_timelock.avg = miss_timelock.timelock.avg - bsl_avg;
        save('norm_miss_timelock','norm_miss_timelock')
    end
end


%% Compute permuted condition differences per session x location (get perm_cond_max/min)

for isess = 1:length(session_names)
    fprintf('Permuted cond diffs session %d/%d\n', isess, length(session_names))
    cd(fullfile(output_paths{isess}))
    if ~exist('trial_perm_ind.mat','file'), continue; end
    load('trial_perm_ind.mat')
    
    for iloc = 1:length(trial_perm_ind.loc)
        locid = trial_perm_ind.loc(iloc).locationID;
        rand_matrix = trial_perm_ind.loc(iloc).rand_matrix;
        hit_trials = trial_perm_ind.loc(iloc).hit;
        miss_trials = trial_perm_ind.loc(iloc).miss;
        
        if isempty(rand_matrix)
            fprintf('Session %s loc %g: no rand_matrix -> skipping cond diffs\n', session_names{isess}, locid)
            continue
        end 
        
        % load baseline average
        locERPdir = fullfile(output_paths{isess}, sprintf('loc%d',locid), 'ERP_real');
        if ~isdir(locERPdir), continue; end
        cd(locERPdir)
        if ~exist('bsl_avg.mat','file'), continue; end
        load('bsl_avg')  % bsl_avg (channels x 1)
        
        % prepare storage
        nchan = length(bsl_avg);
        perm_cond_max = nan(nchan, permut_n);
        perm_cond_min = nan(nchan, permut_n);
        
        for iperm = 1:permut_n
            permdir = fullfile(output_paths{isess}, sprintf('loc%d',locid), num2str(iperm));
            if ~isdir(permdir)
                continue
            end
            cd(permdir)
            if ~exist('hit_lfp_avg.mat','file') || ~exist('miss_lfp_avg.mat','file')
                continue
            end
            hit_perm = load('hit_lfp_avg.mat');   % struct with .timelock
            miss_perm = load('miss_lfp_avg.mat');
            norm_hit = hit_perm.timelock.avg - bsl_avg;
            norm_miss = miss_perm.timelock.avg - bsl_avg;
            
            perm_diff = norm_hit - norm_miss;    % channels x time
            perm_cond_max(:,iperm) = max(perm_diff,[],2);
            perm_cond_min(:,iperm) = min(perm_diff,[],2);
        end
        
        % save location-specific perm distributions
        cd(fullfile(output_paths{isess}))
        locsavefile = fullfile(output_paths{isess}, sprintf('perm_cond_loc%d.mat', locid));
        save(locsavefile, 'perm_cond_max','perm_cond_min')
        clear perm_cond_max perm_cond_min
    end
end

%% Plot real data differences & percentiles per session x location

for isess = 1:length(session_names)
    fprintf('Plotting per-session per-location for session %d/%d\n', isess, length(session_names))
    cd(fullfile(output_paths{isess}))
    if ~exist('trial_perm_ind.mat','file'), continue; end
    load('trial_perm_ind.mat')
    
    for iloc = 1:length(trial_perm_ind.loc)
        locid = trial_perm_ind.loc(iloc).locationID;
        locERPdir = fullfile(output_paths{isess}, sprintf('loc%d',locid), 'ERP_real');
        if ~isdir(locERPdir), continue; end
        
        % load saved perm distributions
        permfile = fullfile(output_paths{isess}, sprintf('perm_cond_loc%d.mat',locid));
        if ~exist(permfile,'file'), continue; end
        load(permfile,'perm_cond_max','perm_cond_min')
        
        % load normalized real data
        cd(locERPdir)
        if ~exist('norm_hit_timelock.mat','file') || ~exist('norm_miss_timelock.mat','file'), continue; end
        hit = load('norm_hit_timelock.mat'); hit = hit.norm_hit_timelock;
        miss = load('norm_miss_timelock.mat'); miss = miss.norm_miss_timelock;
        
        real_cond_diff = hit.avg - miss.avg;   % channels x time
        % significance limits (per channel)
        limit_max = quantile([perm_cond_max perm_cond_min], 0.975, 2);
        limit_min = quantile([perm_cond_max perm_cond_min], 0.025, 2);
        
        if plotting && plotERPs
            fig = figure('Units','normalized','Position',[0 0 1 1]);
            nchan = length(hit.label);
            for ichan = 1:nchan
                sig_time = (real_cond_diff(ichan,:) >= limit_max(ichan,:) ) | (real_cond_diff(ichan,:) <= limit_min(ichan,:));
                sigonset  = find(conv(sig_time,[1 -1])==1);
                sigoffset = find(conv(sig_time,[1 -1])==-1)-1;
                if isempty(sigonset) || isempty(sigoffset)
                    durplot = [];
                else
                    durplot = [L(sigonset); L(sigoffset)];
                end
                
                subplot(8,8,ichan)
                plot(hit.time, hit.avg(ichan,:),'k'); hold on
                plot(miss.time, miss.avg(ichan,:),'r')
                if ~isempty(durplot)
                    for i = 1:size(durplot,2)
                        v = [durplot(1,i) -1.5; durplot(2,i) -1.5; durplot(2,i) 1; durplot(1,i) 1];
                        patch('Faces',[1 2 3 4],'Vertices',v,'FaceColor','blue','FaceAlpha',.2,'EdgeColor','blue','EdgeAlpha',.2)
                    end
                end
                xlim([-0.1 0.15])
                ylim([-1.5 1])
                title(sprintf('Ch %s', char(hit.label(ichan))))
            end
            % save figure
            locid = trial_perm_ind.loc(iloc).locationID;
            baseDir = '/mnt/hpc/projects/MWSampling/4Shivangi/Plots/Hitvsmiss/erpdiff_hermes/lfp/per_loc';
            saveDir = fullfile(baseDir, num2str(locid));
            
            if ~isfolder(saveDir)
                mkdir(saveDir);
            end
            
            savefig_filename = sprintf('%s_loc%d.pdf', session_names{isess}, locid);
            savefig_path = fullfile(saveDir, savefig_filename);
            
            set(fig, 'PaperPositionMode', 'auto')
            print(fig, savefig_path, '-dpdf', '-fillpage');
            close(fig);
            
        end
    end
end

%% Pooling across sessions per location (channel-wise thresholds and ERPs)

% Remove unwanted sessions if needed 
remove_sessions = [1, 3, 21];
keep_sessions = setdiff(1:length(session_names), remove_sessions);
session_namesp = session_names(keep_sessions);
output_pathsp = output_paths(keep_sessions);

% load channel list used previously
cd('/mnt/hpc/projects/MWSampling/4Shivangi/results_hermes/critical_time')
if exist('all_channels.mat','file')
    load('all_channels')  % loads all_channels
else
    error('Please provide all_channels.mat in the critical_time folder')
end
num_channels = length(all_channels);
num_permutations = permut_n;
num_sessions = length(session_namesp);
time_len = length(L);

% For each location, collect perm distributions across sessions
% First find union of locations across sessions
all_locs = [];
for isess = 1:length(session_namesp)
    cd(fullfile(output_pathsp{isess}))
    if exist('trial_perm_ind.mat','file')
        load('trial_perm_ind.mat')
        if isfield(trial_perm_ind,'locations')
            all_locs = union(all_locs, trial_perm_ind.locations(:));
        end
    end
end

% iterate locations
for iloc_all = 1:length(all_locs)
    locid = all_locs(iloc_all);
    fprintf('Pooling location %g across sessions\n', locid)
    
    % initialize storage
    t_dist_max = nan(num_channels, num_permutations, num_sessions);
    t_dist_min = nan(num_channels, num_permutations, num_sessions);
    
    ERPdiff_r = nan(num_sessions, num_channels, time_len);
    ERP_hit = nan(num_channels, time_len, num_sessions);
    ERP_miss = nan(num_channels, time_len, num_sessions);
    
    sess_counter = 0;
    for isess = 1:length(session_namesp)
        sessidx = isess;
        outp = output_pathsp{isess};
        permfile = fullfile(outp, sprintf('perm_cond_loc%d.mat',locid));
        locERPdir = fullfile(outp,sprintf('loc%d',locid), 'ERP_real');
        if ~exist(permfile,'file') || ~isdir(locERPdir)
            continue
        end
        
        sess_counter = sess_counter + 1;
        % load perm distributions if present
        load(permfile,'perm_cond_max','perm_cond_min')  % channel x perm
        % load baseline normalized real data to get channel labels and time
        cd(locERPdir)
        if ~exist('norm_hit_timelock.mat','file') || ~exist('norm_miss_timelock.mat','file')
            continue
        end
        nh = load('norm_hit_timelock.mat'); nh = nh.norm_hit_timelock;
        nm = load('norm_miss_timelock.mat'); nm = nm.norm_miss_timelock;
        timelock = nh;  % use to get labels/time
        
        % map permuted distributions to all_channels using timelock.label
        for ichan = 1:length(timelock.label)
            chan_name = timelock.label{ichan};
            idx = find(strcmp(all_channels, chan_name));
            if ~isempty(idx)
                t_dist_max(idx, :, sess_counter) = perm_cond_max(ichan, :);
                t_dist_min(idx, :, sess_counter) = perm_cond_min(ichan, :);
            end
        end
        
        % map normalized avg to global arrays
        hit_avg_mapped = nan(num_channels, time_len);
        miss_avg_mapped = nan(num_channels, time_len);
        for ichan = 1:length(nh.label)
            chan_name = nh.label{ichan};
            idx = find(strcmp(all_channels, chan_name));
            if ~isempty(idx)
                hit_avg_mapped(idx, :) = nh.avg(ichan, :);
                miss_avg_mapped(idx, :) = nm.avg(ichan, :);
            end
        end
        
        % compute channel-wise difference handling NaNs
        diff_mapped = nan(num_channels, time_len);
        for c = 1:num_channels
            h = hit_avg_mapped(c,:);
            m = miss_avg_mapped(c,:);
            valid = ~isnan(h) & ~isnan(m);
            d = nan(1,time_len);
            d(valid) = h(valid) - m(valid);
            diff_mapped(c,:) = d;
        end
        
        ERPdiff_r(sess_counter, :, :) = diff_mapped;
        ERP_hit(:, :, sess_counter) = hit_avg_mapped;
        ERP_miss(:, :, sess_counter) = miss_avg_mapped;
    end % sessions loop
    
    % if no valid sessions, skip
    if sess_counter == 0
        fprintf('No valid sessions for location %g - skipping pooling\n', locid)
        continue
    end
    
    % reshape distributions and compute channel-wise percentiles
    t_dist = cat(2, t_dist_max, t_dist_min);  % channels x (perms*2) x sess
    % reshape into channels x (perms*2*sessions)
    t_reshaped = reshape(t_dist, num_channels, []);
    limit_max = zeros(num_channels,1);
    limit_min = zeros(num_channels,1);
    for c = 1:num_channels
        valid_data = t_reshaped(c, ~isnan(t_reshaped(c,:)));
        if isempty(valid_data)
            limit_max(c) = NaN; limit_min(c) = NaN;
        else
            limit_max(c) = quantile(valid_data, 0.975);
            limit_min(c) = quantile(valid_data, 0.025);
        end
    end
    
    % average difference across sessions (handling NaNs)
    ERPdiff_avg = squeeze(nanmean(ERPdiff_r(1:sess_counter,:,:),1)); % channels x time
    
    % plotting per-channel pooled ERPs for this location
    if plotting && plotERPs
        fig = figure('Units','normalized','Position',[0 0 1 1]);
        for ichan = 1:num_channels
            sig_time = (ERPdiff_avg(ichan,:) >= limit_max(ichan)) | (ERPdiff_avg(ichan,:) <= limit_min(ichan));
            sigonset  = find(conv(double(sig_time), [1 -1])==1);
            sigoffset = find(conv(double(sig_time), [1 -1])==-1)-1;
            if isempty(sigonset) || isempty(sigoffset)
                durplot = [];
            else
                % clip indices to valid range
                sigonset(sigonset < 1) = [];
                sigoffset(sigoffset < 1) = [];
                durplot = [L(sigonset); L(sigoffset)];
            end
            
            subplot(8,8,ichan)
            plot(nh.time, nanmean(ERP_hit(ichan,:,:),3),'k'); hold on
            plot(nh.time, nanmean(ERP_miss(ichan,:,:),3),'r')
            if ~isempty(durplot)
                for i = 1:size(durplot,2)
                    v = [durplot(1,i) -1; durplot(2,i) -1; durplot(2,i) 1; durplot(1,i) 1];
                    patch('Faces',[1 2 3 4],'Vertices',v,'FaceColor','blue','FaceAlpha',.2,'EdgeColor','blue','EdgeAlpha',.2)
                end
            end
            xlim([-0.1 0.15]); ylim([-1 1])
            title(sprintf('Ch %d', ichan))
        end
        
        % Save figure for this location
        saveDir = '/mnt/hpc/projects/MWSampling/4Shivangi/Plots/Hitvsmiss/erpdiff_hermes/lfp/per_loc';
        if ~isdir(saveDir), mkdir(saveDir); end
        cd(saveDir)
        fname = sprintf('ERPpooled_loc%d', locid);
        savefig(fname)
        set(fig,'PaperPositionMode','auto')
        print(fname, '-dpdf', '-fillpage')
        close(fig)
    end
end

%% Grand average across all channels & sessions per location (pooled scalar thresholds)

for iloc_all = 1:length(all_locs)
    locid = all_locs(iloc_all);
    fprintf('Grand pooled plot for location %g\n', locid)
    
    % collect t_dist for this location across all sessions (flatten)
    all_tvals = [];
    ERP_hit_all = [];
    ERP_miss_all = [];
    timevec = [];
    for isess = 1:length(session_namesp)
        outp = output_pathsp{isess};
        permfile = fullfile(outp, sprintf('perm_cond_loc%d.mat',locid));
        locERPdir = fullfile(outp, sprintf('loc%d',locid), 'ERP_real');
        if ~exist(permfile,'file') || ~isdir(locERPdir)
            continue
        end
        load(permfile,'perm_cond_max','perm_cond_min')
        % collect non-nan values
        all_tvals = [all_tvals; perm_cond_max(:); perm_cond_min(:)]; %#ok<AGROW>
        
        cd(locERPdir)
        if ~exist('norm_hit_timelock.mat','file') || ~exist('norm_miss_timelock.mat','file')
            continue
        end
        nh = load('norm_hit_timelock.mat'); nh = nh.norm_hit_timelock;
        nm = load('norm_miss_timelock.mat'); nm = nm.norm_miss_timelock;
        ERP_hit_all = cat(3, ERP_hit_all, nh.avg); % channels x time x session
        ERP_miss_all = cat(3, ERP_miss_all, nm.avg);
        timevec = nh.time;
    end
    
    if isempty(all_tvals) || isempty(ERP_hit_all)
        fprintf('No data to make grand pooled plot for location %g\n', locid)
        continue
    end
    
    limit_max_all = quantile(all_tvals(~isnan(all_tvals)), 0.975);
    limit_min_all = quantile(all_tvals(~isnan(all_tvals)), 0.025);
    
    % compute grand average across channels & sessions
    % collapse channels then sessions
    hit_mean_over_channels = squeeze(nanmean(nanmean(ERP_hit_all,1),3)); % time x 1
    miss_mean_over_channels = squeeze(nanmean(nanmean(ERP_miss_all,1),3));
    ERPdiff_all_avg = hit_mean_over_channels - miss_mean_over_channels; % time x 1
    
    % find significant timepoints
    sig_time = (ERPdiff_all_avg >= limit_max_all) | (ERPdiff_all_avg <= limit_min_all);
    sig_time = double(sig_time(:)');
    sigonset  = find(conv(sig_time, [1 -1]) == 1);
    sigoffset = find(conv(sig_time, [1 -1]) == -1) - 1;
    if isempty(sigonset) || isempty(sigoffset)
        durplot = [];
    else
        durplot = [timevec(sigonset); timevec(sigoffset)];
    end
    
    % Plot grand pooled
    fig = figure('Name', sprintf('Pooled ERP - Loc %d', locid), 'Units', 'normalized', 'Position', [0.2, 0.2, 0.6, 0.5]);
    plot(timevec, hit_mean_over_channels, 'k', 'LineWidth', 1.7); hold on;
    plot(timevec, miss_mean_over_channels, 'r', 'LineWidth', 1.7);
    for i = 1:size(durplot, 2)
        v = [durplot(1, i) -1; durplot(2, i) -1; durplot(2, i) 1; durplot(1, i) 1];
        patch('Faces', [1 2 3 4], 'Vertices', v, 'FaceColor', 'blue', 'FaceAlpha', .2, 'EdgeColor', 'none')
    end
    xlabel('Time (s)');
    ylabel('Normalized ERP');
    title(sprintf('Pooled ERP Across All Channels & Sessions - Loc %d', locid));
    legend('Hit', 'Miss');
    xlim([-0.1 0.15]);
    ylim([-1 1]);
    grid on;
    
    saveDir = '/mnt/hpc/projects/MWSampling/4Shivangi/Plots/Hitvsmiss/erpdiff_hermes/lfp/per_loc';
    if ~isdir(saveDir), mkdir(saveDir); end
    cd(saveDir)
    fname = sprintf('ERPpooled_GrandAvg_loc%d', locid);
    savefig(fname)
    print(fname, '-dpdf', '-fillpage');
    close(fig)
end

disp('All done.')
