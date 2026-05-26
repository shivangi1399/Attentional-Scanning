% =====================================================================
% comparing_monkeys.m  — Merged-positions rerun for hermes Reg R^2 LFP
%
% Tests whether merging hermes' 16 fine-grained positions into 8 paired
% bins reveals the Reg R^2 LFP low-freq cluster that klecks shows but
% raw-16-position hermes does not.
%
% Scope: hermes only, Reg R^2 LFP only. Klecks unchanged (its 9 native
% positions are already coarse-spaced).
%
% Merge pairs (hermes positions -> merged group 1..8):
%   [46,52]  [58,64]  [70,76]   [82,89]
%   [95,101] [107,113][119,125] [131,144]
%
% Self-orchestrating in three phases (auto-detected by file presence):
%   1. Observed:  compute merged H2 observed R^2_phase per channel,
%                 cache (along with X/Y per (pos,freq) per channel needed
%                 for perm dispatch).
%   2. Dispatch:  if perm files not yet on disk, dispatch 1000 SLURM
%                 perms per channel via regress_perm_R_pos.
%   3. Collect:   once all perm files exist, collect per-channel and
%                 channel-average nulls, compute paired Jensen-corrected
%                 thresholds against existing H1 (complex) perms, and
%                 render the 3-panel comparison (hermes unmerged vs
%                 hermes merged vs klecks).
%
% Outputs:
%   <base>/results_hermes/multi_lin_reg/abs_per_pos_merged_LFP/
%       cp10_till_100/observed_cache.mat
%       cp10_till_100/perm_R_pos/LFP_ERP_ampl_all/<ch>/perm_NNNN.mat
%   <base>/Plots/sampling_compare/hypotheses/comparing_monkeys/
%       merged_vs_unmerged_hermes_LFP.pdf
%       merged_vs_unmerged_hermes_LFP_chan_avg.pdf
%
% Invoked from compare_hypotheses_per_chan.m when
% COMPARE_MERGED_HERMES_LFP_REG is true. Rerun the calling script after
% SLURM completes to enter the collect phase.
% =====================================================================

%% Settings
cm_base   = '/mnt/hpc/projects/MWSampling/4Shivangi';
cm_dv     = 'LFP_ERP_ampl_all';
cm_nPerm  = 1000;
cm_alpha  = 0.05;
cm_merge_pairs = {[46 52], [58 64], [70 76], [82 89], ...
                  [95 101], [107 113], [119 125], [131 144]};

cm_out_data = fullfile(cm_base, 'results_hermes', 'multi_lin_reg', ...
    'abs_per_pos_merged_LFP', 'cp10_till_100');
cm_out_plot = fullfile(cm_base, 'Plots','sampling_compare', ...
    'hypotheses','comparing_monkeys');
if ~exist(cm_out_data,'dir'), mkdir(cm_out_data); end
if ~exist(cm_out_plot,'dir'), mkdir(cm_out_plot); end

cm_perm_root = fullfile(cm_out_data, 'perm_R_pos', cm_dv);

%% Helpers on path
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/multiple_linear_reg/functions
addpath /opt/ESIsoftware/matlab/slurmfun/

%% Phase 1: observed + cache X/Y per channel for perm dispatch
cm_obs_cache = fullfile(cm_out_data, 'observed_cache.mat');
if isfile(cm_obs_cache)
    fprintf('[comparing_monkeys] Loading cached observed\n');
    load(cm_obs_cache, 'cm_R2_phase_merged', 'cm_nPos', 'cm_nFreq', 'cm_nCh');
    cm_X_per_ch = []; cm_Y_per_ch = [];  % only rebuilt if dispatch needed
else
    fprintf('[comparing_monkeys] Computing merged observed + caching X/Y\n');
    [cm_R2_phase_merged, cm_X_per_ch, cm_Y_per_ch, cm_nPos, cm_nFreq, cm_nCh] = ...
        cm_compute_observed(cm_base, cm_dv, cm_merge_pairs);
    save(cm_obs_cache, 'cm_R2_phase_merged', 'cm_nPos', 'cm_nFreq', 'cm_nCh', ...
        'cm_merge_pairs', '-v7.3');
    fprintf('[comparing_monkeys] Cached observed to %s\n', cm_obs_cache);
end

%% Phase 2: dispatch SLURM perms if not all present
cm_missing = cm_count_missing_perms(cm_perm_root, cm_nCh, cm_nPerm);
if cm_missing > 0
    fprintf('[comparing_monkeys] %d perm files missing — dispatching SLURM\n', cm_missing);
    if isempty(cm_X_per_ch)
        fprintf('[comparing_monkeys] Rebuilding X/Y per channel for dispatch\n');
        [~, cm_X_per_ch, cm_Y_per_ch, ~, ~, ~] = ...
            cm_compute_observed(cm_base, cm_dv, cm_merge_pairs);
    end
    cm_dispatch_slurm(cm_X_per_ch, cm_Y_per_ch, ...
        cm_nPos, cm_nFreq, cm_nCh, cm_nPerm, cm_perm_root);
    fprintf(['[comparing_monkeys] SLURM dispatched. If slurmfun returned\n' ...
             '  non-blocking, rerun compare_hypotheses_per_chan.m after\n' ...
             '  jobs finish to enter the collect phase.\n']);
    cm_remain = cm_count_missing_perms(cm_perm_root, cm_nCh, cm_nPerm);
    if cm_remain > 0
        fprintf('[comparing_monkeys] %d files still missing — exiting.\n', cm_remain);
        return
    end
end

%% Phase 3: collect perms, compute thresholds, plot
fprintf('[comparing_monkeys] All perm files present — collecting\n');

% Per-channel null_R_phase: nCh × nPerm × nFreq
cm_null_R_phase = nan(cm_nCh, cm_nPerm, cm_nFreq);
for ch = 1:cm_nCh
    for perm = 1:cm_nPerm
        f = fullfile(cm_perm_root, num2str(ch), sprintf('perm_%04d.mat', perm));
        if ~isfile(f), continue; end
        R = load(f, 'results');
        cm_null_R_phase(ch, perm, :) = R.results.null_R_phase(:)';
    end
end

% Channel-average null: mean across channels per (perm, freq)
cm_null_avg_R_phase = squeeze(mean(cm_null_R_phase, 1, 'omitnan'));  % nPerm × nFreq

%% Load H1 (complex) observed + perms for hermes & klecks, and H2 unmerged
fprintf('[comparing_monkeys] Loading H1 / unmerged H2 / klecks for pairing\n');

H2h = load(fullfile(cm_base,'results_hermes','multi_lin_reg', ...
    'abs_per_pos','cp10_till_100', ...
    'multi_regression_channelwise_R2_abs_per_pos.mat'), 'reg_results');
H1h = load(fullfile(cm_base,'results_hermes','multi_lin_reg', ...
    'complex','cp10_till_100', ...
    'multi_regression_channelwise_R2.mat'), 'reg_results');
H2k = load(fullfile(cm_base,'results_klecks','multi_lin_reg', ...
    'abs_per_pos','cp10_till_100', ...
    'multi_regression_channelwise_R2_abs_per_pos.mat'), 'reg_results');
H1k = load(fullfile(cm_base,'results_klecks','multi_lin_reg', ...
    'complex','cp10_till_100', ...
    'multi_regression_channelwise_R2.mat'), 'reg_results');

cm_hermes_H2_unm = H2h.reg_results.(cm_dv).R2_phase;
cm_hermes_H1     = H1h.reg_results.(cm_dv).R2_phase;
cm_klecks_H2     = H2k.reg_results.(cm_dv).R2_phase;
cm_klecks_H1     = H1k.reg_results.(cm_dv).R2_phase;

% H1 channel-avg perms (already paired with H2 unmerged via rng seed)
H1h_chan_avg = load(fullfile(cm_base,'results_hermes','multi_lin_reg', ...
    'complex','cp10_till_100','perm_R',cm_dv,'channel_avg_results.mat'), ...
    'null_avg_R_phase_freq');
cm_H1_perm_chan_avg_hermes = H1h_chan_avg.null_avg_R_phase_freq;

% H2 unmerged channel-avg perms (for the existing-comparison threshold)
H2h_chan_avg = load(fullfile(cm_base,'results_hermes','multi_lin_reg', ...
    'abs_per_pos','cp10_till_100','perm_R_pos',cm_dv,'channel_avg_results.mat'), ...
    'null_avg_R_phase_freq');
cm_H2unm_perm_chan_avg_hermes = H2h_chan_avg.null_avg_R_phase_freq;

% Klecks chan-avg perms (for reference threshold)
H2k_chan_avg = load(fullfile(cm_base,'results_klecks','multi_lin_reg', ...
    'abs_per_pos','cp10_till_100','perm_R_pos',cm_dv,'channel_avg_results.mat'), ...
    'null_avg_R_phase_freq');
H1k_chan_avg = load(fullfile(cm_base,'results_klecks','multi_lin_reg', ...
    'complex','cp10_till_100','perm_R',cm_dv,'channel_avg_results.mat'), ...
    'null_avg_R_phase_freq');
cm_H2_perm_chan_avg_klecks = H2k_chan_avg.null_avg_R_phase_freq;
cm_H1_perm_chan_avg_klecks = H1k_chan_avg.null_avg_R_phase_freq;

%% Paired Jensen-corrected thresholds (channel average)
cm_n_for_paired = min([cm_nPerm, size(cm_H1_perm_chan_avg_hermes,1)]);

% Hermes merged: paired diff_null = H2_merged_perm - H1_perm
cm_diff_null_merged = cm_null_avg_R_phase(1:cm_n_for_paired,:) ...
                    - cm_H1_perm_chan_avg_hermes(1:cm_n_for_paired,:);
cm_thr_merged_chan_avg = quantile(max(cm_diff_null_merged, [], 2), 1 - cm_alpha);
cm_obs_diff_merged_chan_avg   = mean(cm_R2_phase_merged, 1, 'omitnan') ...
                              - mean(cm_hermes_H1,       1, 'omitnan');

% Hermes unmerged paired chan-avg
cm_n_for_unm = min([size(cm_H2unm_perm_chan_avg_hermes,1), size(cm_H1_perm_chan_avg_hermes,1)]);
cm_diff_null_unm = cm_H2unm_perm_chan_avg_hermes(1:cm_n_for_unm,:) ...
                 - cm_H1_perm_chan_avg_hermes(1:cm_n_for_unm,:);
cm_thr_unm_chan_avg = quantile(max(cm_diff_null_unm, [], 2), 1 - cm_alpha);
cm_obs_diff_unm_chan_avg = mean(cm_hermes_H2_unm, 1, 'omitnan') ...
                         - mean(cm_hermes_H1,     1, 'omitnan');

% Klecks paired chan-avg
cm_n_for_klecks = min([size(cm_H2_perm_chan_avg_klecks,1), size(cm_H1_perm_chan_avg_klecks,1)]);
cm_diff_null_klecks = cm_H2_perm_chan_avg_klecks(1:cm_n_for_klecks,:) ...
                    - cm_H1_perm_chan_avg_klecks(1:cm_n_for_klecks,:);
cm_thr_klecks_chan_avg = quantile(max(cm_diff_null_klecks, [], 2), 1 - cm_alpha);
cm_obs_diff_klecks_chan_avg = mean(cm_klecks_H2, 1, 'omitnan') ...
                            - mean(cm_klecks_H1, 1, 'omitnan');

%% Per-channel paired thresholds for hermes merged
fprintf('[comparing_monkeys] Computing per-channel paired thresholds (merged)\n');
cm_pc_thr_merged = nan(cm_nCh, 1);
cm_pc_sig_merged = false(cm_nCh, cm_nFreq);
cm_pc_obs_merged = cm_R2_phase_merged - cm_hermes_H1;   % nCh × nFreq

for ch = 1:cm_nCh
    h1_pc_file = fullfile(cm_base, 'results_hermes', 'multi_lin_reg', ...
        'complex','cp10_till_100','perm_R', cm_dv, num2str(ch), 'per_channel_null.mat');
    if ~isfile(h1_pc_file), continue; end
    H1pc = load(h1_pc_file, 'null_R2_phase');
    if ~isfield(H1pc, 'null_R2_phase'), continue; end
    h1_perm_ch = H1pc.null_R2_phase;
    h2_perm_ch = squeeze(cm_null_R_phase(ch,:,:));
    n_pair = min([size(h1_perm_ch,1), size(h2_perm_ch,1)]);
    if n_pair < 1, continue; end
    diff_perm = h2_perm_ch(1:n_pair,:) - h1_perm_ch(1:n_pair,:);
    tmax = max(diff_perm, [], 2);
    cm_pc_thr_merged(ch) = quantile(tmax, 1 - cm_alpha);
    if isfinite(cm_pc_thr_merged(ch))
        cm_pc_sig_merged(ch,:) = cm_pc_obs_merged(ch,:) >= cm_pc_thr_merged(ch);
    end
end

%% Per-channel obs/sig for hermes unmerged and klecks (existing perms)
[cm_pc_obs_unm,    cm_pc_sig_unm]    = cm_load_existing_per_chan( ...
    cm_base, 'hermes', cm_dv, cm_hermes_H2_unm, cm_hermes_H1, cm_alpha);
[cm_pc_obs_klecks, cm_pc_sig_klecks] = cm_load_existing_per_chan( ...
    cm_base, 'klecks', cm_dv, cm_klecks_H2, cm_klecks_H1, cm_alpha);

%% Frequency axes
fS = load(fullfile(cm_base,'results_hermes','multi_lin_reg','cp10_till_100','frequency.mat'));
cm_freq_h = fS.frequency;
fS = load(fullfile(cm_base,'results_klecks','multi_lin_reg','cp10_till_100','frequency.mat'));
cm_freq_k = fS.frequency;

%% 3-panel heatmap with significance fade
cm_plot_heatmap_grid(cm_pc_obs_unm,    cm_pc_sig_unm,    cm_freq_h, ...
                     cm_pc_obs_merged, cm_pc_sig_merged, cm_freq_h, ...
                     cm_pc_obs_klecks, cm_pc_sig_klecks, cm_freq_k, ...
                     cm_nPos, cm_out_plot);

%% Channel-average comparison plot
cm_plot_chan_avg(cm_freq_h, cm_obs_diff_unm_chan_avg,    cm_thr_unm_chan_avg, ...
                 cm_freq_h, cm_obs_diff_merged_chan_avg, cm_thr_merged_chan_avg, ...
                 cm_freq_k, cm_obs_diff_klecks_chan_avg, cm_thr_klecks_chan_avg, ...
                 cm_nPos, cm_out_plot);

fprintf('[comparing_monkeys] Done. Plots in %s\n', cm_out_plot);

%% =====================================================================
%% Local helpers
%% =====================================================================

function [R2_phase, X_per_ch, Y_per_ch, nPos, nFreq, nCh] = ...
        cm_compute_observed(base, dv, merge_pairs)
% Compute hermes merged-position observed R^2_phase. Returns the cleaned
% regressors per (channel, position, freq) for downstream perm dispatch.
S = load(fullfile(base, 'results_hermes', 'multi_lin_reg', ...
    'cp10_till_100', 'ph_all_sess.mat'), 'ph_comb');
ph = S.ph_comb;

ph_all  = ph.phase_all;
amp_all = ph.amp_all;
pup_bsl = ph.pup_baseline;
mua_bsl = ph.MUA_baseline;
Y_orig  = ph.(dv);
ti      = ph.LFP_ERP_trialinfo;

pos_raw = ti(:,16);
pos     = nan(size(pos_raw));
for g = 1:numel(merge_pairs)
    pos(ismember(pos_raw, merge_pairs{g})) = g;
end
ti(:,16) = pos;
positions = (1:numel(merge_pairs))';
nPos = numel(positions);

[~, nFreq, nCh] = size(ph_all);
R2_phase  = nan(nCh, nFreq);
X_per_ch  = cell(nCh, 1);
Y_per_ch  = cell(nCh, 1);

for ch = 1:nCh
    if mod(ch,8)==1, fprintf('  observed ch %d/%d\n', ch, nCh); end
    Y_all = Y_orig(:, ch);
    X_freq_pos = cell(nPos, nFreq);
    Y_freq_pos = cell(nPos, nFreq);
    if all(isnan(Y_all))
        X_per_ch{ch} = X_freq_pos; Y_per_ch{ch} = Y_freq_pos;
        continue
    end
    for f = 1:nFreq
        X = [pup_bsl(:,ch), mua_bsl(:,ch), amp_all(:,f,ch), ...
             sin(ph_all(:,f,ch)), cos(ph_all(:,f,ch))];
        R2_pos = nan(nPos, 1);
        for p = 1:nPos
            mask   = ti(:,16) == positions(p);
            nanIdx = any(isnan(X),2) | isnan(Y_all) | ~mask;
            Xc = X(~nanIdx,:);  Yc = Y_all(~nanIdx);
            X_freq_pos{p,f} = Xc;
            Y_freq_pos{p,f} = Yc;
            if length(Yc) < size(Xc,2)+1, continue; end
            X_full   = [ones(size(Xc,1),1), Xc];
            b_full   = regress(Yc, X_full);
            RSS_full = sum((Yc - X_full*b_full).^2);
            RSS_null = sum((Yc - mean(Yc)).^2);
            X_red   = [ones(size(Xc,1),1), Xc(:,1:3)];
            b_red   = regress(Yc, X_red);
            RSS_red = sum((Yc - X_red*b_red).^2);
            R2_pos(p) = max(0, (RSS_red - RSS_full)/RSS_null);
        end
        R2_phase(ch, f) = mean(R2_pos, 'omitnan');
    end
    X_per_ch{ch} = X_freq_pos;
    Y_per_ch{ch} = Y_freq_pos;
end
end

function n_missing = cm_count_missing_perms(perm_root, nCh, nPerm)
n_missing = 0;
for ch = 1:nCh
    for p = 1:nPerm
        f = fullfile(perm_root, num2str(ch), sprintf('perm_%04d.mat', p));
        if ~isfile(f), n_missing = n_missing + 1; end
    end
end
end

function cm_dispatch_slurm(X_per_ch, Y_per_ch, nPos, nFreq, nCh, nPerm, perm_root)
% Build one slurmfun job per (channel × 10-perm bundle), matching the
% upstream regress_stats_R2_abs_per_pos.m dispatch pattern.
nJobs_per_ch = ceil(nPerm/10);
cfg_array = cell(nCh * nJobs_per_ch, 1);
idx = 0;
for ch = 1:nCh
    ch_dir = fullfile(perm_root, num2str(ch));
    if ~exist(ch_dir,'dir'), mkdir(ch_dir); end
    for j = 1:nJobs_per_ch
        ps = (j-1)*10 + 1;
        pe = min(j*10, nPerm);
        if all(arrayfun(@(pp) isfile(fullfile(ch_dir, sprintf('perm_%04d.mat', pp))), ps:pe)), continue; end
        idx = idx + 1;
        cfg_array{idx}.X_all_freq_pos = X_per_ch{ch};
        cfg_array{idx}.Y_all_freq_pos = Y_per_ch{ch};
        cfg_array{idx}.numFreq        = nFreq;
        cfg_array{idx}.nPos           = nPos;
        cfg_array{idx}.isLogistic     = false;
        cfg_array{idx}.perm_idx       = ps:pe;
        cfg_array{idx}.output_dir     = ch_dir;
    end
end
cfg_array = cfg_array(1:idx);
fprintf('  dispatching %d jobs (%d channels × %d bundles)\n', ...
    numel(cfg_array), nCh, nJobs_per_ch);
slurmfun(@regress_perm_R_pos, cfg_array, ...
    'partition','8GB','stopOnError',false,'useUserPath',true, ...
    'waitForToolboxes',{'statistics_toolbox'});
end

function [obs, sig] = cm_load_existing_per_chan(base, animal, dv, H2, H1, alpha)
% Per-channel paired sig using existing per_channel_null.mat files for
% both H1 and H2. Used for the unmerged hermes and klecks panels.
obs = H2 - H1;
[nCh, nFreq] = size(obs);
sig = false(nCh, nFreq);
h1_root = fullfile(base, ['results_' animal], 'multi_lin_reg', ...
    'complex','cp10_till_100','perm_R', dv);
h2_root = fullfile(base, ['results_' animal], 'multi_lin_reg', ...
    'abs_per_pos','cp10_till_100','perm_R_pos', dv);
for ch = 1:nCh
    h1f = fullfile(h1_root, num2str(ch), 'per_channel_null.mat');
    h2f = fullfile(h2_root, num2str(ch), 'per_channel_null.mat');
    if ~isfile(h1f) || ~isfile(h2f), continue; end
    H1pc = load(h1f, 'null_R2_phase');
    H2pc = load(h2f, 'null_R2_phase');
    if ~isfield(H1pc,'null_R2_phase') || ~isfield(H2pc,'null_R2_phase'), continue; end
    p1 = H1pc.null_R2_phase; p2 = H2pc.null_R2_phase;
    n_pair = min([size(p1,1), size(p2,1)]);
    if n_pair < 1, continue; end
    if size(p1,2) ~= nFreq || size(p2,2) ~= nFreq, continue; end
    diff_perm = p2(1:n_pair,:) - p1(1:n_pair,:);
    thr = quantile(max(diff_perm,[],2), 1 - alpha);
    if isfinite(thr), sig(ch,:) = obs(ch,:) >= thr; end
end
end

function cm_plot_heatmap_grid(obs_unm, sig_unm, freq_h_unm, ...
                              obs_mer, sig_mer, freq_h_mer, ...
                              obs_kl,  sig_kl,  freq_k, nPos, out_plot)
fig = figure('Visible','off','Units','centimeters','Position',[1 1 30 10]);
set(fig,'PaperUnits','centimeters','PaperSize',fig.Position(3:4), ...
    'PaperPosition',[0 0 fig.Position(3:4)]);

mats   = {obs_unm, obs_mer, obs_kl};
sigs   = {sig_unm, sig_mer, sig_kl};
freqs  = {freq_h_unm, freq_h_mer, freq_k};
titles = {'hermes unmerged (16 pos)', sprintf('hermes merged (%d pos)', nPos), 'klecks (9 pos)'};

all_vals = cellfun(@(M) M(~isnan(M)), mats, 'UniformOutput', false);
all_vals = vertcat(all_vals{:});
if isempty(all_vals), mx = 1;
else, mx = max(abs(prctile(all_vals,[2 98]))); if mx==0, mx = max(abs(all_vals))+eps; end
end

alpha_sig    = 1.0;
alpha_nonsig = 0.20;

for k = 1:3
    ax = subplot(1,3,k);
    img = imagesc(ax, freqs{k}, 1:size(mats{k},1), mats{k});
    set(ax,'YDir','normal','Color',[1 1 1]);
    colormap(ax, parula);
    caxis(ax, [-mx mx]);
    A = alpha_nonsig * ones(size(mats{k}));
    if ~isempty(sigs{k}), A(sigs{k}) = alpha_sig; end
    A(isnan(mats{k})) = 0;
    set(img, 'AlphaData', A);
    cb = colorbar(ax); cb.FontSize = 6;
    xlabel(ax,'Freq (Hz)'); ylabel(ax,'Channel');
    title(ax, titles{k}, 'FontSize',9);
    set(ax,'FontSize',7,'Box','on');
    xlim(ax,[min(freqs{k}) max(freqs{k})]);
end
sgtitle({'Reg R^2 LFP :: observed H2-H1 (channel x freq)', ...
         'shared color scale; faded = n.s. by per-channel paired test'}, ...
        'FontSize',10,'FontWeight','bold');

ps_file  = fullfile(out_plot, 'merged_vs_unmerged_hermes_LFP.ps');
pdf_file = fullfile(out_plot, 'merged_vs_unmerged_hermes_LFP.pdf');
if isfile(ps_file),  delete(ps_file);  end
if isfile(pdf_file), delete(pdf_file); end
print(fig,'-dpsc','-painters',ps_file);
close(fig);
[stat,msg] = system(sprintf('ps2pdf "%s" "%s"', ps_file, pdf_file));
if stat==0 && isfile(pdf_file), delete(ps_file);
else, warning('[comparing_monkeys] ps2pdf failed (%d): %s', stat, strtrim(msg));
end
end

function cm_plot_chan_avg(freq_h, d_unm, thr_unm, ...
                          freq_h2, d_mer, thr_mer, ...
                          freq_k, d_kl, thr_kl, nPos, out_plot)
fig = figure('Visible','off','Units','centimeters','Position',[1 1 30 8]);
set(fig,'PaperUnits','centimeters','PaperSize',fig.Position(3:4), ...
    'PaperPosition',[0 0 fig.Position(3:4)]);
titles = {'hermes unmerged (16 pos)', sprintf('hermes merged (%d pos)', nPos), 'klecks (9 pos)'};
datas  = {d_unm, d_mer, d_kl};
thrs   = {thr_unm, thr_mer, thr_kl};
freqs  = {freq_h, freq_h2, freq_k};
col    = [0.20 0.65 0.40];
for k = 1:3
    ax = subplot(1,3,k); hold(ax,'on');
    yline(ax, 0, 'k-','LineWidth',0.4);
    plot(ax, freqs{k}, datas{k}, 'Color', col, 'LineWidth',1.5);
    if isfinite(thrs{k})
        yline(ax, thrs{k}, ':','Color',col,'LineWidth',1);
    end
    xlabel(ax,'Freq (Hz)','FontSize',8);
    ylabel(ax,'\Delta R^2_{phase}','FontSize',8);
    title(ax, titles{k}, 'FontSize',9);
    set(ax,'FontSize',7,'Box','on');
end
sgtitle('Reg R^2 LFP :: channel-average H2-H1 (Jensen-corrected, paired)', ...
    'FontSize',10,'FontWeight','bold');

ps_file  = fullfile(out_plot, 'merged_vs_unmerged_hermes_LFP_chan_avg.ps');
pdf_file = fullfile(out_plot, 'merged_vs_unmerged_hermes_LFP_chan_avg.pdf');
if isfile(ps_file),  delete(ps_file);  end
if isfile(pdf_file), delete(pdf_file); end
print(fig,'-dpsc','-painters',ps_file);
close(fig);
[stat,msg] = system(sprintf('ps2pdf "%s" "%s"', ps_file, pdf_file));
if stat==0 && isfile(pdf_file), delete(ps_file);
else, warning('[comparing_monkeys] ps2pdf failed (%d): %s', stat, strtrim(msg));
end
end
