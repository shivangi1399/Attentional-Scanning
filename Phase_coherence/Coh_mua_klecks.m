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

%% Create data paths

datafolder   = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks';

cd(datafolder),
animalName = 'klecks';
temp = dir;
session_names = [];
ii = 0;
for i = 1:length(temp)
    if strfind(temp(i).name,animalName)
        ii = ii+1;
        session_names{ii,1} = temp(i).name;
    end
end

session_paths_files = [];
session_paths_files = cellfun(@(x) fullfile(datafolder,x, 'clean_mua.mat'), session_names, 'uniform',0);

phase_paths = cellfun(@(x) fullfile(datafolder, x,'Phase_analysis/hit_miss'),session_names, 'uniform',0);
data_folder = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks/multi_lin_reg/cp10_till_100';
output_folder = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks/phase_coherence/cp10_till_100';
permut_n = 1000;

%% coherence all locations and difficulty levels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(data_folder)
load('ph_all_sess.mat')

% Real data
nCh = 64;
for ichan = 1:nCh
    ichan
    phase = ph_comb.phase_all(:,:,ichan);
    erp_amp = ph_comb.MUA_ERP_ampl_all(:,ichan);

    [coh, phase_spec] = phase_coherence(phase, erp_amp);

    chan_folder = fullfile(output_folder,'mua','all_loc_difflev', num2str(ichan));
    if ~exist(chan_folder,'dir'), mkdir(chan_folder); end
    save(fullfile(chan_folder,'coherence.mat'),'coh','phase_spec');
end

% Permutation test
nTrials  = length(ph_comb.trialinfo);
perm_indices = arrayfun(@(x) randperm(nTrials), 1:permut_n, 'UniformOutput', false);
trial_idx = 1:size(ph_comb.trialinfo, 1);

cfg = cell(1,nCh);
for ichan = 1:nCh
    cfg{ichan}.ichan        = ichan;
    cfg{ichan}.permut_n     = permut_n;
    cfg{ichan}.infile       = fullfile(data_folder);
    cfg{ichan}.outfile      = fullfile(output_folder, 'mua',...
        'all_loc_difflev');
    cfg{ichan}.perm_indices = perm_indices;
    cfg{ichan}.trial_idx    = trial_idx;
end

% Launch jobs
slurmfun(@phase_coherence_perm_mua, cfg, ...
    'partition',   '8GB', ...
    'stopOnError', false, ...
    'useUserPath', true);

%% coherence particular locations and difficulty levels - keep all Dlev %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(data_folder)
load('ph_all_sess.mat')

% Find unique locations and difficulty levels
targ_loc = unique(ph_comb.trialinfo(:,16));
diff_levels = unique(ph_comb.trialinfo(:,18));

min_difflev = [];
max_difflev = [];

for iloc = 1:length(targ_loc)
    hits = find(ph_comb.trialinfo(:,20) == 1 & ph_comb.trialinfo(:,16) == targ_loc(iloc));
    miss = find(ph_comb.trialinfo(:,20) == 5 & ph_comb.trialinfo(:,16) == targ_loc(iloc));
    de = ph_comb.trialinfo([hits; miss], 18);
    
    min_difflev = [min_difflev; min(de)];
    max_difflev = [max_difflev; max(de)];
end

%% coherence (real data)

for iloc = 1:length(targ_loc)
    iloc
    loc = targ_loc(iloc);
    min_th = min_difflev(iloc);
    max_th = max_difflev(iloc);
    
    % Select trials for this location and difficulty
    trial_idx = find(ph_comb.trialinfo(:,16) == loc & ...
                     ph_comb.trialinfo(:,18) >= min_th & ...
                     ph_comb.trialinfo(:,18) <= max_th);
    
    phase_data = ph_comb.phase_all(trial_idx,:,:);  % trials x freq x channels
    erp_amp    = ph_comb.MUA_ERP_ampl_all(trial_idx,:); % trials x channels
    
    for ichan = 1:64
        phase = squeeze(phase_data(:,:,ichan));  % trials x freq
        erp   = erp_amp(:,ichan);
        [coh, phase_spec] = phase_coherence(phase, erp);
        
        chan_folder = fullfile(output_folder, 'mua',...
            'loc_difflev_all', ...
            sprintf('loc%d', loc), ...
            sprintf('%d_%d', min_th, max_th), ...
            num2str(ichan));
        if ~exist(chan_folder, 'dir'), mkdir(chan_folder); end
        save(fullfile(chan_folder,'coherence.mat'),'coh','phase_spec');
    end
    
    %% Permutation
    
    permut_n = 1000;
    nTrials  = length(trial_idx);
    perm_indices = arrayfun(@(x) randperm(nTrials), 1:permut_n, 'UniformOutput', false);
    
    cfg = cell(1,64);
    for ichan = 1:64
        cfg{ichan}.ichan        = ichan;
        cfg{ichan}.permut_n     = permut_n;
        cfg{ichan}.infile       = fullfile(data_folder);
        cfg{ichan}.outfile      = fullfile(output_folder, 'mua',...
            'loc_difflev_all', ...
            sprintf('loc%d', loc), ...
            sprintf('%d_%d', min_th, max_th));
        cfg{ichan}.perm_indices = perm_indices;
        cfg{ichan}.trial_idx    = trial_idx;
    end
    
    % Launch jobs
    slurmfun(@phase_coherence_perm_mua, cfg, ...
        'partition',   '8GB', ...
        'stopOnError', false, ...
        'useUserPath', true);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting all locations and difficulty levels

save_root = '/mnt/hpc/projects/MWSampling/4Shivangi/Plots/phase_coherence/klecks/cp10_till_100/mua/all_loc_difflev';
if ~exist(save_root,'dir'), mkdir(save_root); end

cd(output_folder)
load('frequency.mat')
freq = frequency;
nCh = 64;

limit_maxc = nan(1, nCh);   % coherence thresholds
limit_maxp = nan(1, nCh);   % phase spec thresholds
Coh_chan   = false(nCh, numel(freq));
Phase_chan = false(nCh, numel(freq));

for ch = 1:nCh
    ch_folder = fullfile(output_folder,'mua','all_loc_difflev', num2str(ch));
    if ~exist(ch_folder, 'dir')
        warning(['Skipping channel ' num2str(ch) ' (folder missing)']);
        continue
    end
    
    cd(ch_folder);
    
    if ~exist('coherence.mat','file') || ~exist('coh_perm.mat','file') || ...
       ~exist('coherence.mat','file') || ~exist('phase_spec_perm.mat','file')
        warning(['Skipping channel ' num2str(ch) ' (missing required files)']);
        continue
    end
    
    load coherence
    load coh_perm
    load phase_spec_perm
    
    % Skip if data has NaN
    if any(isnan(coh)) || any(isnan(coh_perm(:))) || ...
       any(isnan(phase_spec)) || any(isnan(phase_spec_perm(:)))
        warning(['Skipping channel ' num2str(ch) ' (NaN values present)']);
        continue
    end
    
    % Thresholds
    tmaxc = max(coh_perm, [], 2);
    limit_maxc(ch) = quantile(tmaxc, 0.95);
    
    tmaxp = max(phase_spec_perm, [], 2);
    limit_maxp(ch) = quantile(tmaxp, 0.95);
    
    if isnan(limit_maxc(ch)) || isnan(limit_maxp(ch))
        warning(['Skipping channel ' num2str(ch) ' (NaN thresholds)']);
        continue
    end
    
    % Plot per channel
    f1 = figure(1);
    subplot(8, 8, ch);
    plot_sig(freq, coh, limit_maxc(ch), 'Frequency', 'Coherence');
    title(['Ch ' num2str(ch)])

    f2 = figure(2);
    subplot(8, 8, ch);
    plot_sig(freq, phase_spec, limit_maxp(ch),'Frequency', 'Phase spec');
    title(['Ch ' num2str(ch)])
    
    % Store significant frequencies
    Coh_chan(ch,:)   = coh   >= limit_maxc(ch);
    Phase_chan(ch,:) = phase_spec  >= limit_maxp(ch);
end

set(f1, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
set(f1, 'PaperPositionMode', 'auto');
set(f1, 'Renderer', 'opengl');
print(f1, fullfile(save_root, 'all_channels_coherence.png'),'-dpng', '-r0'); 
set(f2, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
set(f2, 'PaperPositionMode', 'auto');
set(f2, 'Renderer', 'opengl');
print(f2, fullfile(save_root, 'all_channels_phase.png'), '-dpng', '-r0');

% Summary heatmaps
f3 = figure;
imagesc(freq, 1:nCh, Coh_chan);
set(gca, 'YDir', 'normal');
xlabel('Frequency (Hz)'); ylabel('Channels');
title('Significant Coherence per Channel');
caxis([0 1]); colorbar;
saveas(f3, fullfile(save_root, 'summary_coherence.png'));

f4 = figure;
imagesc(freq, 1:nCh, Phase_chan);
set(gca, 'YDir', 'normal');
xlabel('Frequency (Hz)'); ylabel('Channels');
title('Significant Phase Spectrum per Channel');
caxis([0 1]); colorbar;
saveas(f4, fullfile(save_root, 'summary_phase.png'));

% Combine Across Channels
valid_idx = ~isnan(limit_maxc) & ~isnan(limit_maxp);
coh_all = [];
phase_all = [];
coh_perm_all = [];
phase_perm_all = [];

for ch = find(valid_idx)
    cd(fullfile(output_folder,'mua','all_loc_difflev', num2str(ch)));
    load coherence
    load coh_perm
    load phase_spec_perm
    
    coh_all = [coh_all; coh];
    phase_all = [phase_all; phase_spec];
    
    coh_perm_all   = cat(3, coh_perm_all, coh_perm);
    phase_perm_all = cat(3, phase_perm_all, phase_spec_perm);
end

if isempty(coh_all) || isempty(phase_all)
    warning('No valid channels for combined plot');
else
    % Coherence combined
    coh_avg = nanmean(coh_all, 1);
    coh_perm_avg = nanmean(coh_perm_all, 3);
    tmax_all = nanmax(coh_perm_avg, [], 2);
    limit_avg = quantile(tmax_all, 0.95);
    
    f5 = figure;
    plot_sig(freq, coh_avg, limit_avg,'Frequency', 'Coherence');
    title('All Channels Combined - Coherence');
    saveas(f5, fullfile(save_root, 'combined_coherence.png'));
    
    % Phase spec combined
    phase_avg = nanmean(phase_all, 1);
    phase_perm_avg = nanmean(phase_perm_all, 3);
    tmaxp_all = nanmax(phase_perm_avg, [], 2);
    limit_avgp = quantile(tmaxp_all, 0.95);
    
    f6 = figure;
    plot_sig(freq, phase_avg, limit_avgp,'Frequency', 'Phase spec');
    title('All Channels Combined - Phase Spec');
    saveas(f6, fullfile(save_root, 'combined_phase.png'));
end

%% Plotting combinations of difficulty level and locations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_folder = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks/phase_coherence/cp10_till_100';
save_root = '/mnt/hpc/projects/MWSampling/4Shivangi/Plots/coherence/klecks/cp10_till_100/mua/loc_difflev_all';
if ~exist(save_root,'dir'), mkdir(save_root); end

cd(output_folder)
load('frequency.mat')
freq = frequency;
nCh = 64;

% Load target locations and difficulty levels 
cd(data_folder)
load('ph_all_sess.mat')
targ_loc = unique(ph_comb.trialinfo(:,16));
diff_levels = unique(ph_comb.trialinfo(:,18));

min_difflev = [];
max_difflev = [];

for iloc = 1:length(targ_loc)
    hits = find(ph_comb.trialinfo(:,20) == 1 & ph_comb.trialinfo(:,16) == targ_loc(iloc));
    miss = find(ph_comb.trialinfo(:,20) == 5 & ph_comb.trialinfo(:,16) == targ_loc(iloc));
    de = ph_comb.trialinfo([hits; miss], 18);
    min_difflev = [min_difflev; min(de)];
    max_difflev = [max_difflev; max(de)];
end

% Plotting coherence + phase spec per location and difficulty level
for iloc = 1:length(targ_loc)
    loc = targ_loc(iloc);
    min_th = min_difflev(iloc);
    max_th = max_difflev(iloc);

    fprintf('Plotting coherence + phase spec for Loc%d Dlev%d_%d\n', loc, min_th, max_th);

    limit_maxc = nan(1,nCh);   % coherence thresholds
    limit_maxp = nan(1,nCh);   % phase spec thresholds
    Coh_chan   = false(nCh, numel(freq));
    Phase_chan = false(nCh, numel(freq));

    % PER CHANNEL
    f1 = figure(1); clf;
    f2 = figure(2); clf;
    for ch = 1:nCh
        ch_folder = fullfile(output_folder,'mua','loc_difflev_all',...
            sprintf('loc%d',loc), sprintf('%d_%d',min_th,max_th), num2str(ch));
        if ~exist(ch_folder,'dir'), continue; end
        cd(ch_folder);

        if ~exist('coherence.mat','file') || ~exist('coh_perm.mat','file') || ...
           ~exist('coherence.mat','file') || ~exist('phase_spec_perm.mat','file')
            continue; 
        end
        load coherence.mat coh phase_spec
        load coh_perm.mat coh_perm
        load phase_spec_perm.mat phase_spec_perm

        if any(isnan(coh)) || any(isnan(coh_perm(:))) || ...
           any(isnan(phase_spec)) || any(isnan(phase_spec_perm(:)))
            continue
        end

        % thresholds
        tmaxc = max(coh_perm,[],2);
        limit_maxc(ch) = quantile(tmaxc,0.95);

        tmaxp = max(phase_spec_perm,[],2);
        limit_maxp(ch) = quantile(tmaxp,0.95);

        % coherence subplot
        figure(f1);
        subplot(8,8,ch);
        plot_sig(freq, coh, limit_maxc(ch), 'Frequency','Coherence');
        title(['Ch ' num2str(ch)])
        Coh_chan(ch,:) = coh >= limit_maxc(ch);
        
        % phase spec subplot
        figure(f2);
        subplot(8,8,ch);
        plot_sig(freq, phase_spec, limit_maxp(ch), 'Frequency','Phase Spec');
        title(['Ch ' num2str(ch)])
        Phase_chan(ch,:) = phase_spec >= limit_maxp(ch);
    end
    
    % Save all channels grid
    set(f1, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    set(f1, 'PaperPositionMode', 'auto');
    set(f1, 'Renderer', 'opengl');
    print(f1, fullfile(save_root, sprintf('Loc%d_Dlev%d_%d_all_channels_coh.png',loc,min_th,max_th)),'-dpng', '-r0');
    set(f2, 'Units', 'normalized', 'OuterPosition', [0 0 1 1]);
    set(f2, 'PaperPositionMode', 'auto');
    set(f2, 'Renderer', 'opengl');
    print(f2, fullfile(save_root, sprintf('Loc%d_Dlev%d_%d_all_channels_phase.png',loc,min_th,max_th)),'-dpng', '-r0');
    
    % SUMMARY HEATMAPS
    f3 = figure;
    imagesc(freq,1:nCh,Coh_chan);
    set(gca,'YDir','normal');
    xlabel('Frequency (Hz)'); ylabel('Channels');
    title(sprintf('Significant Coherence Loc%d (%d_%d)',loc,min_th,max_th));
    caxis([0 1]); colorbar;
    saveas(f3, fullfile(save_root, sprintf('Loc%d_Dlev%d_%d_summary_coh.png',loc,min_th,max_th)));

    f4 = figure;
    imagesc(freq,1:nCh,Phase_chan);
    set(gca,'YDir','normal');
    xlabel('Frequency (Hz)'); ylabel('Channels');
    title(sprintf('Significant Phase Spec Loc%d (%d_%d)',loc,min_th,max_th));
    caxis([0 1]); colorbar;
    saveas(f4, fullfile(save_root, sprintf('Loc%d_Dlev%d_%d_summary_phase.png',loc,min_th,max_th)));
    close(f4);
    
    % COMBINED ACROSS CHANNELS
    valid_idx = ~isnan(limit_maxc) & ~isnan(limit_maxp);
    coh_all = []; phase_all = [];
    coh_perm_all = []; phase_perm_all = [];
    for ch = find(valid_idx)
        ch_folder = fullfile(output_folder,'mua','loc_difflev_all',...
            sprintf('loc%d',loc), sprintf('%d_%d',min_th,max_th), num2str(ch));
        cd(ch_folder);
        load coherence.mat coh phase_spec
        load coh_perm.mat coh_perm
        load phase_spec_perm.mat phase_spec_perm
        coh_all = [coh_all; coh];
        phase_all = [phase_all; phase_spec];
        coh_perm_all   = cat(3,coh_perm_all,coh_perm);
        phase_perm_all = cat(3,phase_perm_all,phase_spec_perm);
    end

    if ~isempty(coh_all)
        % coherence combined
        coh_avg = mean(coh_all,1);
        coh_perm_avg = mean(coh_perm_all,3);
        tmax_all = max(coh_perm_avg,[],2);
        limit_avg = quantile(tmax_all,0.95);
        f5 = figure;
        plot_sig(freq, coh_avg, limit_avg,'Frequency','Coherence');
        title(sprintf('Combined Loc%d Dlev%d_%d - Coherence',loc,min_th,max_th));
        saveas(f5, fullfile(save_root, sprintf('Loc%d_Dlev%d_%d_combined_coh.png',loc,min_th,max_th)));
        close(f5);
        
        % phase spec combined
        phase_avg = mean(phase_all,1);
        phase_perm_avg = mean(phase_perm_all,3);
        tmaxp_all = max(phase_perm_avg,[],2);
        limit_avgp = quantile(tmaxp_all,0.95);
        f6 = figure;
        plot_sig(freq, phase_avg, limit_avgp,'Frequency','Phase Spec');
        title(sprintf('Combined Loc%d Dlev%d_%d - Phase Spec',loc,min_th,max_th));
        saveas(f6, fullfile(save_root, sprintf('Loc%d_Dlev%d_%d_combined_phase.png',loc,min_th,max_th)));
        close(f6);
    else
        warning('No valid channels for Loc%d (%d-%d)',loc,min_th,max_th)
    end
    close(f1);
    close(f2);
end
