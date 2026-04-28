% Overlay H1 (complex / pooled), H2 (abs_per_pos), and H3 (abs_per_pos_diff)
% magnitudes on the same axes for coherence, correlation, and regression R^2,
% across the four DVs (MUA, LFP, RT, Hit/Miss). Both per-animal and
% monkey-average panels are produced.
%
% Loads from the channel-average .mat files saved by each hypothesis script:
%   H1: results_<animal>/phase_coherence/complex/cp10_till_100/...
%   H2: results_<animal>/phase_coherence/abs_per_pos/cp10_till_100/...
%   H3: results_<animal>/phase_coherence/abs_per_pos_diff/cp10_till_100/...
%   (and the parallel phase_correlation / multi_lin_reg trees)

clearvars; close all; clc

%% Settings
animal = 'klecks';   % 'hermes' or 'klecks'

addpath /opt/fieldtrip_github/; ft_defaults

base_results = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ['results_' animal]);
results_combined = '/mnt/hpc/projects/MWSampling/4Shivangi/results_combined';

save_root = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/sampling_compare/hypotheses', animal);
if ~exist(save_root, 'dir'), mkdir(save_root); end
monkey_save_root = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/sampling_compare/hypotheses/monkey_avg');
if ~exist(monkey_save_root, 'dir'), mkdir(monkey_save_root); end

% Hypothesis subfolder names + line styles
hyp_keys   = {'complex',     'abs_per_pos',  'abs_per_pos_diff'};
hyp_labels = {'H1 (pooled)', 'H2 (per pos)', 'H3 (per pos x diff)'};
hyp_colors = [0.20 0.40 0.80; 0.20 0.65 0.40; 0.80 0.30 0.30];

row_labels = {'MUA','LFP','RT','Hit/Miss'};

% Coherence channel-avg files: one per (DV, hypothesis)
coh_files = @(root) {
    fullfile(root, 'phase_coherence',   'mua', 'all_loc_difflev', 'channel_avg_results.mat');
    fullfile(root, 'phase_coherence',   'lfp', 'all_loc_difflev', 'channel_avg_results.mat');
    fullfile(root, 'phase_coherence',   'RT',  'all_loc_difflev', 'channel_avg_results.mat');
    fullfile(root, 'phase_correlation', 'hit_miss', 'all_loc_difflev', 'channel_avg_results_itc.mat');
    };

corr_files = @(root) {
    fullfile(root, 'phase_correlation', 'mua',      'all_loc_difflev', 'channel_avg_results.mat');
    fullfile(root, 'phase_correlation', 'lfp',      'all_loc_difflev', 'channel_avg_results.mat');
    fullfile(root, 'phase_correlation', 'RT',       'all_loc_difflev', 'channel_avg_results.mat');
    fullfile(root, 'phase_correlation', 'hit_miss', 'all_loc_difflev', 'channel_avg_results_pos.mat');
    };

%% Helper: extract value/threshold from a channel-avg file (handles all 3 H's)

    function [val, thr, freq] = pick_coh(file)
    val = []; thr = NaN; freq = [];
    if ~isfile(file), return; end
    s = load(file);
    % H1 ITC files store itc_chan_avg / thresh_chan_avg_itc, others store coh_*
    if isfield(s,'coh_chan_avg'),       val = s.coh_chan_avg; end
    if isfield(s,'itc_chan_avg'),       val = s.itc_chan_avg; end
    if isfield(s,'thresh_chan_avg'),    thr = s.thresh_chan_avg; end
    if isfield(s,'thresh_chan_avg_itc'),thr = s.thresh_chan_avg_itc; end
    if isfield(s,'freq'),               freq = s.freq; end
    end

    function [val, thr, freq] = pick_corr(file)
    val = []; thr = NaN; freq = [];
    if ~isfile(file), return; end
    s = load(file);
    if isfield(s,'corr_chan_avg'),         val = s.corr_chan_avg; end
    if isfield(s,'pos_chan_avg'),          val = s.pos_chan_avg; end
    if isfield(s,'thresh_chan_avg'),       thr = s.thresh_chan_avg; end
    if isfield(s,'thresh_chan_avg_pos'),   thr = s.thresh_chan_avg_pos; end
    if isfield(s,'freq'),                  freq = s.freq; end
    end

    function [val, thr, freq] = pick_monkey_coh(file)
    val = []; thr = NaN; freq = [];
    if ~isfile(file), return; end
    s = load(file);
    if isfield(s,'coh_monkey_avg'),     val = s.coh_monkey_avg; end
    if isfield(s,'itc_monkey_avg'),     val = s.itc_monkey_avg; end
    if isfield(s,'thresh_monkey_avg'),  thr = s.thresh_monkey_avg; end
    if isfield(s,'thresh_monkey_avg_itc'), thr = s.thresh_monkey_avg_itc; end
    if isfield(s,'freq'),               freq = s.freq; end
    end

    function [val, thr, freq] = pick_monkey_corr(file)
    val = []; thr = NaN; freq = [];
    if ~isfile(file), return; end
    s = load(file);
    if isfield(s,'corr_monkey_avg'),    val = s.corr_monkey_avg; end
    if isfield(s,'pos_monkey_avg'),     val = s.pos_monkey_avg; end
    if isfield(s,'thresh_monkey_avg'),  thr = s.thresh_monkey_avg; end
    if isfield(s,'thresh_monkey_avg_pos'), thr = s.thresh_monkey_avg_pos; end
    if isfield(s,'freq'),               freq = s.freq; end
    end

%% Load per-animal data

coh_curves = cell(3,4); coh_thrs = NaN(3,4);
corr_curves = cell(3,4); corr_thrs = NaN(3,4);
freq_axis = [];

for h = 1:numel(hyp_keys)
    root_h_coh  = fullfile(base_results, 'phase_coherence');
    root_h_corr = fullfile(base_results, 'phase_correlation');
    % insert hypothesis folder + cp10_till_100
    coh_root  = fullfile(root_h_coh,  hyp_keys{h}, 'cp10_till_100');
    corr_root = fullfile(root_h_corr, hyp_keys{h}, 'cp10_till_100');
    cf = coh_files(fullfile(base_results, 'placeholder')); %#ok<NASGU>
    % rebuild per hypothesis (placeholder above unused)
    cfs = {
        fullfile(coh_root,  'mua', 'all_loc_difflev', 'channel_avg_results.mat');
        fullfile(coh_root,  'lfp', 'all_loc_difflev', 'channel_avg_results.mat');
        fullfile(coh_root,  'RT',  'all_loc_difflev', 'channel_avg_results.mat');
        fullfile(corr_root, 'hit_miss', 'all_loc_difflev', 'channel_avg_results_itc.mat');
        };
    crs = {
        fullfile(corr_root, 'mua',      'all_loc_difflev', 'channel_avg_results.mat');
        fullfile(corr_root, 'lfp',      'all_loc_difflev', 'channel_avg_results.mat');
        fullfile(corr_root, 'RT',       'all_loc_difflev', 'channel_avg_results.mat');
        fullfile(corr_root, 'hit_miss', 'all_loc_difflev', 'channel_avg_results_pos.mat');
        };
    for r = 1:4
        [v, t, f] = pick_coh(cfs{r});
        coh_curves{h,r} = v; coh_thrs(h,r) = t;
        if ~isempty(f) && isempty(freq_axis), freq_axis = f; end

        [v2, t2, ~] = pick_corr(crs{r});
        corr_curves{h,r} = v2; corr_thrs(h,r) = t2;
    end
end

%% Load regression per animal (channel_avg_R.phase + channel_avg_thresh.phase)

reg_curves = cell(3,4); reg_thrs = NaN(3,4); freqs_reg = [];
reg_dvs = {'MUA_ERP_ampl_all','LFP_ERP_ampl_all','RT','hit_miss'};

for h = 1:numel(hyp_keys)
    if strcmp(hyp_keys{h},'complex')
        rfile = fullfile(base_results,'multi_lin_reg','complex','cp10_till_100', ...
            'multi_regression_channelwise_R2.mat');
    elseif strcmp(hyp_keys{h},'abs_per_pos')
        rfile = fullfile(base_results,'multi_lin_reg','abs_per_pos','cp10_till_100', ...
            'multi_regression_channelwise_R2_abs_per_pos.mat');
    else
        rfile = fullfile(base_results,'multi_lin_reg','abs_per_pos_diff','cp10_till_100', ...
            'multi_regression_channelwise_R2_abs_per_pos_diff.mat');
    end
    if ~isfile(rfile), continue; end
    s = load(rfile,'reg_results'); rr = s.reg_results;
    for r = 1:4
        dv = reg_dvs{r};
        if isfield(rr,dv) && isfield(rr.(dv),'channel_avg_R') ...
                && isfield(rr.(dv).channel_avg_R,'phase')
            reg_curves{h,r} = rr.(dv).channel_avg_R.phase;
            reg_thrs(h,r)   = rr.(dv).channel_avg_thresh.phase;
        end
    end
    if isempty(freqs_reg)
        sd = dir(fullfile(base_results,[animal '_*']));
        if ~isempty(sd)
            ff = fullfile(base_results, sd(1).name, ...
                'Phase_analysis','hit_miss','100iter_cut@cp_m10','1','freqpow.mat');
            if exist(ff,'file'), tmp = load(ff); freqs_reg = tmp.freqpow.freq; end
        end
        if isempty(freqs_reg), freqs_reg = freq_axis; end
    end
end

%% FIGURE 1 — Per-animal H1/H2/H3 overlay

f1 = figure('Name', sprintf('H1/H2/H3 overlay - %s', animal), ...
    'Units','centimeters','Position',[1 1 48 38]);
set(f1,'PaperUnits','centimeters','PaperSize',[48 38],'PaperPosition',[0 0 48 38]);

for r = 1:4
    % Coherence
    subplot(4,3,(r-1)*3+1); hold on;
    for h = 1:3
        v = coh_curves{h,r};
        if ~isempty(v)
            plot(freq_axis, v, 'Color', hyp_colors(h,:), 'LineWidth', 2);
            if ~isnan(coh_thrs(h,r))
                yline(coh_thrs(h,r), '--', 'Color', hyp_colors(h,:), 'LineWidth', 1);
            end
        end
    end
    title(sprintf('Coh: %s', row_labels{r}), 'FontSize', 9);
    xlabel('Freq (Hz)'); ylabel('Coh / ITC');
    if r == 1, legend(hyp_labels, 'Location', 'best', 'FontSize', 7); end
    set(gca,'FontSize',8,'Box','on');

    % Correlation
    subplot(4,3,(r-1)*3+2); hold on;
    for h = 1:3
        v = corr_curves{h,r};
        if ~isempty(v)
            plot(freq_axis, v, 'Color', hyp_colors(h,:), 'LineWidth', 2);
            if ~isnan(corr_thrs(h,r))
                yline(corr_thrs(h,r), '--', 'Color', hyp_colors(h,:), 'LineWidth', 1);
            end
        end
    end
    title(sprintf('Corr: %s', row_labels{r}), 'FontSize', 9);
    xlabel('Freq (Hz)'); ylabel('Corr / POS');
    set(gca,'FontSize',8,'Box','on');

    % Regression
    subplot(4,3,(r-1)*3+3); hold on;
    for h = 1:3
        v = reg_curves{h,r};
        if ~isempty(v)
            plot(freqs_reg, v, 'Color', hyp_colors(h,:), 'LineWidth', 2);
            if ~isnan(reg_thrs(h,r))
                yline(reg_thrs(h,r), '--', 'Color', hyp_colors(h,:), 'LineWidth', 1);
            end
        end
    end
    title(sprintf('Reg R^2: %s', row_labels{r}), 'FontSize', 9);
    xlabel('Freq (Hz)'); ylabel('R^2');
    set(gca,'FontSize',8,'Box','on');
end

sgtitle(sprintf('H1 vs H2 vs H3 — channel-avg, %s', animal), 'FontSize', 13, 'FontWeight', 'bold');
print(f1, fullfile(save_root, 'compare_hypotheses.pdf'), '-dpdf');
fprintf('Per-animal H1/H2/H3 overlay saved: %s\n', fullfile(save_root,'compare_hypotheses.pdf'));

%% Load monkey-average data

mk_coh_curves = cell(3,4); mk_coh_thrs = NaN(3,4);
mk_corr_curves = cell(3,4); mk_corr_thrs = NaN(3,4);
freq_mk = [];

for h = 1:numel(hyp_keys)
    coh_root  = fullfile(results_combined, 'phase_coherence',   hyp_keys{h}, 'cp10_till_100');
    corr_root = fullfile(results_combined, 'phase_correlation', hyp_keys{h}, 'cp10_till_100');
    cfs = {
        fullfile(coh_root,  'mua', 'all_loc_difflev', 'monkey_avg_results.mat');
        fullfile(coh_root,  'lfp', 'all_loc_difflev', 'monkey_avg_results.mat');
        fullfile(coh_root,  'RT',  'all_loc_difflev', 'monkey_avg_results.mat');
        fullfile(corr_root, 'hit_miss_itc', 'all_loc_difflev', 'monkey_avg_results_itc.mat');
        };
    crs = {
        fullfile(corr_root, 'mua',      'all_loc_difflev', 'monkey_avg_results.mat');
        fullfile(corr_root, 'lfp',      'all_loc_difflev', 'monkey_avg_results.mat');
        fullfile(corr_root, 'RT',       'all_loc_difflev', 'monkey_avg_results.mat');
        fullfile(corr_root, 'hit_miss', 'all_loc_difflev', 'monkey_avg_results_pos.mat');
        };
    for r = 1:4
        [v, t, f] = pick_monkey_coh(cfs{r});
        mk_coh_curves{h,r} = v; mk_coh_thrs(h,r) = t;
        if ~isempty(f) && isempty(freq_mk), freq_mk = f; end
        [v2, t2, ~] = pick_monkey_corr(crs{r});
        mk_corr_curves{h,r} = v2; mk_corr_thrs(h,r) = t2;
    end
end

% Regression monkey avg per hypothesis
mk_reg_curves = cell(3,4); mk_reg_thrs = NaN(3,4);
for h = 1:numel(hyp_keys)
    reg_root = fullfile(results_combined,'multi_lin_reg', hyp_keys{h}, 'cp10_till_100');
    for r = 1:4
        rf = fullfile(reg_root, reg_dvs{r}, 'monkey_avg_results.mat');
        if ~isfile(rf), continue; end
        s = load(rf);
        if isfield(s,'monkey_avg_obs')
            mk_reg_curves{h,r} = s.monkey_avg_obs.phase;
            mk_reg_thrs(h,r)   = s.thresh_monkey.phase;
        end
    end
end

%% FIGURE 2 — Monkey-average H1/H2/H3 overlay

f2 = figure('Name','H1/H2/H3 overlay - Monkey average', ...
    'Units','centimeters','Position',[1 1 48 38]);
set(f2,'PaperUnits','centimeters','PaperSize',[48 38],'PaperPosition',[0 0 48 38]);

if isempty(freq_mk), freq_mk = freq_axis; end

for r = 1:4
    subplot(4,3,(r-1)*3+1); hold on;
    for h = 1:3
        v = mk_coh_curves{h,r};
        if ~isempty(v)
            plot(freq_mk, v, 'Color', hyp_colors(h,:), 'LineWidth', 2);
            if ~isnan(mk_coh_thrs(h,r))
                yline(mk_coh_thrs(h,r), '--', 'Color', hyp_colors(h,:), 'LineWidth', 1);
            end
        end
    end
    title(sprintf('Coh: %s', row_labels{r}),'FontSize',9);
    xlabel('Freq (Hz)'); ylabel('Coh / ITC');
    if r == 1, legend(hyp_labels,'Location','best','FontSize',7); end
    set(gca,'FontSize',8,'Box','on');

    subplot(4,3,(r-1)*3+2); hold on;
    for h = 1:3
        v = mk_corr_curves{h,r};
        if ~isempty(v)
            plot(freq_mk, v, 'Color', hyp_colors(h,:), 'LineWidth', 2);
            if ~isnan(mk_corr_thrs(h,r))
                yline(mk_corr_thrs(h,r), '--', 'Color', hyp_colors(h,:), 'LineWidth', 1);
            end
        end
    end
    title(sprintf('Corr: %s', row_labels{r}),'FontSize',9);
    xlabel('Freq (Hz)'); ylabel('Corr / POS');
    set(gca,'FontSize',8,'Box','on');

    subplot(4,3,(r-1)*3+3); hold on;
    for h = 1:3
        v = mk_reg_curves{h,r};
        if ~isempty(v)
            plot(freqs_reg, v, 'Color', hyp_colors(h,:), 'LineWidth', 2);
            if ~isnan(mk_reg_thrs(h,r))
                yline(mk_reg_thrs(h,r), '--', 'Color', hyp_colors(h,:), 'LineWidth', 1);
            end
        end
    end
    title(sprintf('Reg R^2: %s', row_labels{r}),'FontSize',9);
    xlabel('Freq (Hz)'); ylabel('R^2');
    set(gca,'FontSize',8,'Box','on');
end

sgtitle('H1 vs H2 vs H3 — monkey-average channel curves','FontSize',13,'FontWeight','bold');
print(f2, fullfile(monkey_save_root, 'compare_hypotheses_monkey_avg.pdf'), '-dpdf');
fprintf('Monkey-avg H1/H2/H3 overlay saved: %s\n', fullfile(monkey_save_root,'compare_hypotheses_monkey_avg.pdf'));
