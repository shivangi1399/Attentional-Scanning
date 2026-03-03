clear all
close all
clc

%% Specify paths

addpath /opt/fieldtrip_github/
ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/coherence_analysis
addpath /mnt/hpc/projects/MWSampling/4Shivangi
clc

%% Paths and parameters

datafolder   = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks';
output_folder = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks/coherence';
animalName = 'klecks';
permut_n = 1000;
alpha = 0.05;

% Session names
temp = dir(datafolder);
session_names = {};
for i = 1:length(temp)
    if contains(temp(i).name, animalName)
        session_names{end+1,1} = temp(i).name;
    end
end

session_paths = cellfun(@(x) fullfile(datafolder,x), session_names, 'uniform',0);
output_paths = cellfun(@(x) fullfile(datafolder, x,'Coherence'),session_names, 'uniform',0);

%% Compute coherence per session

all_coh = struct();
all_freq = [];

for isess = 1%:length(session_names)
    
    cd(session_paths{isess});
    load('clean_lfp.mat'); load('clean_mua.mat');
    lfpTrials = clean_data;
    muaTrials = clean_mua;
    
    % Add prefixes
    lfpTrials.label = cellfun(@(x) ['lfp_' x], lfpTrials.label, 'UniformOutput', false);
    muaTrials.label = cellfun(@(x) ['mua_' x], muaTrials.label, 'UniformOutput', false);

    % Trim trials
    t1 = 1001; %starting point
    for t = 1:min(numel(lfpTrials.trial), numel(muaTrials.trial))
        nMin = 1200 %min(numel(lfpTrials.time{t}), numel(muaTrials.time{t}));
        lfpTrials.trial{t} = lfpTrials.trial{t}(:,t1:nMin);
        muaTrials.trial{t} = muaTrials.trial{t}(:,t1:nMin);
        lfpTrials.time{t} = lfpTrials.time{t}(t1:nMin);
        muaTrials.time{t} = muaTrials.time{t}(t1:nMin);
    end
    
    % Combine datasets
    cfg = []; cfg.keepsampleinfo = 'no';
    data = ft_appenddata(cfg, lfpTrials, muaTrials);
    
    % Frequency analysis
    cfg_freq = [];
    cfg_freq.output     = 'fourier';
    cfg_freq.method     = 'mtmfft';
    cfg_freq.foilim     = [2 100];
    cfg_freq.tapsmofrq  = 5;
    cfg_freq.keeptrials = 'yes';
    cfg_freq.channel    = {'lfp_*' 'mua_*'};
    freq = ft_freqanalysis(cfg_freq, data);
    
    % Coherence
    cfg_coh = [];
    cfg_coh.method     = 'coh';
    cfg_coh.channelcmb = {'lfp_*' 'mua_*'};
    fd = ft_connectivityanalysis(cfg_coh, freq);
    
    % Identify strictly matching pairs: mua_1-lfp_1, mua_2-lfp_2, etc.
    nPairs = min(length(lfpTrials.label), length(muaTrials.label));
    same_pairs_idx = zeros(nPairs,1);
    for p = 1:nPairs
        pair_label = [lfpTrials.label{p} ' - ' muaTrials.label{p}];
        idx = find(strcmp(fd.labelcmb(:,1), lfpTrials.label{p}) & strcmp(fd.labelcmb(:,2), muaTrials.label{p}));
        if ~isempty(idx)
            same_pairs_idx(p) = idx;
        end
    end
    same_pairs_idx = same_pairs_idx(same_pairs_idx>0); % remove zeros
    
    % Store
    all_coh(isess).cohspctrm = fd.cohspctrm;
    all_coh(isess).labelcmb  = fd.labelcmb;
    all_coh(isess).freq      = fd.freq;
    all_coh(isess).same_pairs_idx = same_pairs_idx;
    
    all_freq = fd.freq;
end

if ~isfolder(output_folder)
    mkdir(output_folder)
end
cd(fullfile(output_folder))
save all_coh all_coh

%% plotting single session to check

fd = all_coh(isess); 
freq = fd.freq;
coh = fd.cohspctrm;
labelcmb = fd.labelcmb;
same_pairs_idx = fd.same_pairs_idx;

figure('Color', 'w');
nPairs = length(same_pairs_idx);

for i = 1:nPairs
    subplot(ceil(sqrt(nPairs)), ceil(sqrt(nPairs)), i);
    plot(freq, coh(same_pairs_idx(i), :), 'k', 'LineWidth', 1.5);
    ylim([0 1]);
    xlabel('Frequency (Hz)');
    ylabel('Coherence');
    
    % Extract channel number/name without 'lfp_' or 'mua_'
    lfp_name = erase(labelcmb{same_pairs_idx(i),1}, 'lfp_');
    mua_name = erase(labelcmb{same_pairs_idx(i),2}, 'mua_');
    
    % Title with cleaned names
    title(sprintf('%s', lfp_name), 'Interpreter', 'none');
    
    grid on;
end

sgtitle(sprintf('Session %d: LFP - MUA Coherence', isess));

%% Permutation test (shuffle MUA trials)

perm_cfg_all = cell(length(session_names),1);

for isess = 1:length(session_names)
    
    cd(session_paths{isess});
    load('clean_lfp.mat'); load('clean_mua.mat');
    lfpTrials = clean_data;
    muaTrials = clean_mua;
    
    % Prefix
    lfpTrials.label = cellfun(@(x) ['lfp_' x], lfpTrials.label, 'UniformOutput', false);
    muaTrials.label = cellfun(@(x) ['mua_' x], muaTrials.label, 'UniformOutput', false);
    
    % Trim trials
    for t = 1:min(numel(lfpTrials.trial), numel(muaTrials.trial))
        nMin = min(numel(lfpTrials.time{t}), numel(muaTrials.time{t}));
        lfpTrials.trial{t} = lfpTrials.trial{t}(:,1:nMin);
        muaTrials.trial{t} = muaTrials.trial{t}(:,1:nMin);
        lfpTrials.time{t} = lfpTrials.time{t}(1:nMin);
        muaTrials.time{t} = muaTrials.time{t}(1:nMin);
    end
    
    nTrials = numel(muaTrials.trial);
    trial_perm_ind.rand_matrix = zeros(permut_n, nTrials);
    for iperm = 1:permut_n
        trial_perm_ind.rand_matrix(iperm,:) = randperm(nTrials);
    end
    
    cfg = cell(1,permut_n);
    
    for iperm = 1:permut_n
        perm_folder = fullfile(output_paths{isess}, num2str(iperm));
        if ~isfolder(perm_folder)
            mkdir(perm_folder)
        end
        
        % Shuffle MUA trials
        mua_shuff = muaTrials;
        mua_shuff.trial = mua_shuff.trial(trial_perm_ind.rand_matrix(iperm,:));
        mua_shuff.time  = mua_shuff.time(trial_perm_ind.rand_matrix(iperm,:));
        
        % Combine with LFP
        cfg_append = []; cfg_append.keepsampleinfo = 'no';
        data_shuff = ft_appenddata(cfg_append, lfpTrials, mua_shuff);
        
        cfg{iperm}.data        = data_shuff;
        cfg{iperm}.foilim      = [2 100];
        cfg{iperm}.tapsmofrq   = 5;
        cfg{iperm}.channel     = {'lfp_*' 'mua_*'};
        cfg{iperm}.same_pairs_idx = all_coh(isess).same_pairs_idx;
        cfg{iperm}.outputfile  = fullfile(perm_folder,'fd_shuff.mat');
    end
    
    perm_cfg_all{isess} = cfg;
    
end

for isess = 1:length(session_names)
    slurmfun(@perm_mua_lfp_coherence, perm_cfg_all{isess}, ...
        'partition', '8GB', 'stopOnError', false, 'useUserPath', true);
end

%% Combine coherence across sessions for each matching pair

% Collect all matching pairs across sessions
matching_pairs = {};
for s = 1:length(all_coh)
    matching_pairs = [matching_pairs; all_coh(s).labelcmb(all_coh(s).same_pairs_idx,:)];
end

pair_strings = strcat(matching_pairs(:,1), '_', matching_pairs(:,2));
[~, ia] = unique(pair_strings);
matching_pairs = matching_pairs(ia,:);  % unique rows

nPairs = size(matching_pairs,1);
nFreqs = length(all_freq);
combined_coh = nan(nPairs, nFreqs, length(all_coh));

for s = 1:length(all_coh)
    session_pairs = all_coh(s).labelcmb(all_coh(s).same_pairs_idx,:);
    for p = 1:nPairs
        idx = find(strcmp(session_pairs(:,1), matching_pairs{p,1}) & ...
                   strcmp(session_pairs(:,2), matching_pairs{p,2}));
        if ~isempty(idx)
            combined_coh(p,:,s) = all_coh(s).cohspctrm(all_coh(s).same_pairs_idx(idx),:);
        end
    end
end

mean_coh = nanmean(combined_coh,3);  % average across sessions

%% Compute permutation significance per pair across sessions

max_null = zeros(permut_n, nFreqs, nPairs);
min_null = zeros(permut_n, nFreqs, nPairs);

session_count = zeros(1,nPairs);

for s = 1:length(session_names)
    for iperm = 1:permut_n
        load(fullfile(output_paths{s}, num2str(iperm),'fd_shuff.mat'));
        session_pairs = fd_shuff.labelcmb;
        for p = 1:nPairs
            idx = find(strcmp(session_pairs(:,1), matching_pairs{p,1}) & ...
                       strcmp(session_pairs(:,2), matching_pairs{p,2}));
            if ~isempty(idx)
                max_null(iperm,:,p) = max_null(iperm,:,p) + fd_shuff.cohspctrm(idx,:);
                min_null(iperm,:,p) = min_null(iperm,:,p) + fd_shuff.cohspctrm(idx,:);
                session_count(p) = session_count(p) + 1;
            end
        end
    end
end

% Divide by number of sessions that actually contributed for each pair
for p = 1:nPairs
    max_null(:,:,p) = max_null(:,:,p) / session_count(p);
    min_null(:,:,p) = min_null(:,:,p) / session_count(p);
end

upper_thr = prctile(max_null, 100-(alpha/2*100),1);
lower_thr = prctile(min_null, (alpha/2*100),1);

%% Plot coherence per matching pair with significance

sig_mask = (mean_coh > upper_thr) | (mean_coh < lower_thr);

for p = 1:nPairs
    figure('Color','w'); hold on;
    plot(all_freq, mean_coh(p,:), 'k', 'LineWidth', 1.5);
    plot(all_freq, upper_thr(:,:,p), 'r--', 'LineWidth', 1);
    scatter(all_freq(sig_mask(p,:)), mean_coh(p,sig_mask(p,:)), 25, 'r', 'filled');

    xlabel('Frequency (Hz)');
    ylabel('Coherence');
    title(sprintf('Average LFPMUA Coherence (Pair %s - %s)', ...
        matching_pairs{p,1}, matching_pairs{p,2}));
    ylim([0 1]); grid on;
end



