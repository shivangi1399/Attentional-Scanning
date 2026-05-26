% =====================================================================
% Per-animal + per-channel paired H2 − H1 test (Jensen-corrected)
%
% Question: "Do positions disagree on preferred phase?"
%
% This script applies the paired Jensen-corrected max-stat test from
% compare_hypotheses.m at two finer levels than the monkey-average one:
%
%   1. Per-animal channel-average — uses each animal's existing
%      channel_avg_results.mat. No new permutations needed because the
%      upstream pipelines seed with rng(2025) so H1 and H2 perms are
%      already matched.
%
%   2. Per-channel within each animal — uses each channel's own H1 and
%      H2 per-perm null files (matched by the same seed). Paired
%      diff_obs vs Jensen-corrected max-stat threshold from
%      null_diff = H2_perm_ch − H1_perm_ch.
%
% Pipelines × DVs covered:
%   - Coherence:          mua, lfp, RT             (chan-avg + per-channel)
%   - Correlation:        mua, lfp, RT, hit_miss   (chan-avg + per-channel)
%   - Regression R²:      mua, lfp, RT, hit_miss   (chan-avg + per-channel*)
%   - Regression R_phase: mua, lfp, RT, hit_miss   (chan-avg + per-channel*)
%
%   *Per-channel regression requires the per-channel null shards to be
%    pre-aggregated into one mat per channel. Run
%    aggregate_regression_per_channel_nulls.m once before this script
%    to produce <perm_R[_pos]>/<dv>/<ch>/per_channel_null.mat. If those
%    files don't exist, the regression per-channel test logs a warning
%    and leaves the per-channel panels empty (channel-avg still works).
%
% Outputs to Plots/sampling_compare/hypotheses/<animal>/:
%   channel_avg_H2-H1.pdf           4 DV × 4 metric grid
%   n_sig_channels_vs_freq.pdf      4 DV × 4 metric grid (all pipelines)
%   per_channel_<pipe>_<dv>.pdf     8×8 channel grid per pipeline × dv
% =====================================================================

clearvars; close all; clc;

%% Settings
animals    = {'hermes','klecks'};
pipelines  = {'coherence','correlation','reg_R2','reg_Rphase'};
pipe_label = {'Coherence','Correlation','Reg R^2','Reg R_{phase}'};
dvs        = {'mua','lfp','RT','hit_miss'};
dv_label   = {'MUA','LFP','RT','Hit/Miss'};
alpha      = 0.05;
nCh        = 64;
base       = '/mnt/hpc/projects/MWSampling/4Shivangi';
out_root   = fullfile(base, 'Plots','sampling_compare','hypotheses');

% Toggle: after the main loop, also run comparing_monkeys.m which
% recomputes hermes Reg R^2 LFP H2 with paired positions (16 -> 8 bins)
% and plots merged-hermes vs unmerged-hermes vs klecks side by side.
COMPARE_MERGED_HERMES_LFP_REG = false;

% Color used for H2−H1 (positions). Matches compare_hypotheses.m green.
COL_POS = [0.20 0.65 0.40];

nA = numel(animals); nP = numel(pipelines); nD = numel(dvs);

%% Per-animal main loop
ca_results       = cell(nA, nP, nD);
nsig_ch_per_freq = cell(nA, nP, nD);
pc_data          = cell(nA, nP, nD);
freq_per_animal  = cell(nA, 1);

for a = 1:nA
    animalName = animals{a};
    fprintf('\n=== %s ===\n', animalName);

    out_dir = fullfile(out_root, animalName);
    if ~exist(out_dir,'dir'), mkdir(out_dir); end

    % Common frequency axis
    freq_file = fullfile(base, ['results_' animalName], ...
        'multi_lin_reg','cp10_till_100','frequency.mat');
    freq = [];
    if isfile(freq_file)
        fS   = load(freq_file);
        freq = fS.frequency;
    end
    freq_per_animal{a} = freq;

    for p = 1:nP
        pipeline = pipelines{p};
        for d = 1:nD
            dv = dvs{d};

            % ── Channel-avg paired test ──
            [obsH1, obsH2, permH1, permH2, fr_l] = ...
                load_chan_avg(animalName, pipeline, dv, base);
            if isempty(obsH1) || isempty(obsH2) ...
                    || isempty(permH1) || isempty(permH2)
                continue
            end
            if isempty(freq) && ~isempty(fr_l), freq = fr_l; end
            if isempty(freq), freq = 1:numel(obsH1); end

            obs_d  = obsH2(:)' - obsH1(:)';                  % 1 × nFreq
            null_d = permH2 - permH1;                         % nPerm × nFreq
            tmax   = max(null_d, [], 2);
            thr    = quantile(tmax, 1-alpha);
            sig    = obs_d >= thr;

            ca_results{a,p,d} = struct( ...
                'animal',animalName,'pipeline',pipeline,'dv',dv, ...
                'obs',obs_d,'null_mean',mean(null_d,1,'omitnan'), ...
                'thr',thr,'sig',sig,'freq',freq);

            % ── Per-channel paired test (coh/corr only) ──
            if ~has_per_chan(pipeline), continue; end
            [pc_obs, pc_null_mean, pc_thr, pc_sig, freq_pc] = ...
                paired_per_channel(animalName, pipeline, dv, nCh, base);
            if isempty(pc_obs), continue; end
            if isempty(freq_pc), freq_pc = freq; end

            nsig_ch_per_freq{a,p,d} = sum(pc_sig, 1, 'omitnan');
            pc_data{a,p,d} = struct( ...
                'obs', pc_obs, 'sig', pc_sig, 'freq', freq_pc);

            % 8×8 per-channel grid
            fig = fig_per_channel_grid( ...
                pc_obs, pc_null_mean, pc_thr, pc_sig, freq_pc, COL_POS, ...
                sprintf('%s — %s — %s :: per-channel H2−H1', ...
                    animalName, pipe_label{p}, dv_label{d}));
            ps_file  = fullfile(out_dir, sprintf('per_channel_%s_%s.ps',  pipeline, dv));
            pdf_file = fullfile(out_dir, sprintf('per_channel_%s_%s.pdf', pipeline, dv));
            save_pdf(fig, ps_file, pdf_file);
            fprintf('  Saved %s\n', pdf_file);
        end
    end

    % Channel-avg grid: 4 DV × 4 pipeline
    fig = fig_chan_avg_grid( ...
        ca_results(a,:,:), pipe_label, dv_label, COL_POS, ...
        sprintf('%s :: channel-average H2−H1 (positions disagree?)', animalName));
    save_pdf(fig, ...
        fullfile(out_dir,'channel_avg_H2-H1.ps'), ...
        fullfile(out_dir,'channel_avg_H2-H1.pdf'));
    fprintf('  Saved %s\n', fullfile(out_dir,'channel_avg_H2-H1.pdf'));

    % Per-channel heatmap grid: 4 DV × 4 pipeline, each panel = channel × freq
    % heatmap of obs H2-H1, with non-significant cells faded.
    fig = fig_per_channel_heatmap_grid( ...
        pc_data(a,:,:), pipe_label, dv_label, ...
        sprintf('%s :: per-channel H2−H1 (channel × freq; faded = n.s.)', animalName));
    save_pdf(fig, ...
        fullfile(out_dir,'per_channel_heatmap.ps'), ...
        fullfile(out_dir,'per_channel_heatmap.pdf'));
    fprintf('  Saved %s\n', fullfile(out_dir,'per_channel_heatmap.pdf'));

    % N sig channels vs freq — all 4 pipelines (regression panels will
    % be empty unless aggregate_regression_per_channel_nulls.m has been
    % run first; fig_nsig_grid renders a "(no data)" placeholder for
    % empty slots).
    fig = fig_nsig_grid( ...
        nsig_ch_per_freq(a, :, :), freq, pipe_label, dv_label, COL_POS, ...
        sprintf('%s :: # channels with sig H2−H1 per freq', animalName));
    save_pdf(fig, ...
        fullfile(out_dir,'n_sig_channels_vs_freq.ps'), ...
        fullfile(out_dir,'n_sig_channels_vs_freq.pdf'));
    fprintf('  Saved %s\n', fullfile(out_dir,'n_sig_channels_vs_freq.pdf'));

    % Diagnostics: usable channels per (pipeline × DV) and trial counts
    % per (position × difficulty) per DV. Helps adjudicate whether across-
    % monkey differences in sig clusters reflect biology vs. methodological
    % asymmetries in channel attrition and trial counts.
    diag_file = fullfile(out_dir, 'diagnostics.txt');
    write_diagnostics(diag_file, animalName, dvs, pipe_label, dv_label, ...
        squeeze(pc_data(a,:,:)), base);
    fprintf('  Saved %s\n', diag_file);
end

fprintf('\nDone.\n');

%% Optional: merged-positions comparison for hermes Reg R^2 LFP
if COMPARE_MERGED_HERMES_LFP_REG
    this_dir = fileparts(mfilename('fullpath'));
    run(fullfile(this_dir, 'comparing_monkeys.m'));
end

%% ───────────────────────────────────────────────────────────────────────
%% Local helpers
%% ───────────────────────────────────────────────────────────────────────

function tf = has_per_chan(pipeline)
% Coherence / correlation always have per-channel data on disk in a
% single nPerm × nFreq mat per channel. Regression (R² and R_phase)
% have per-channel data only after running
% aggregate_regression_per_channel_nulls.m, which consolidates the
% 1000 perm shards into a single per_channel_null.mat per channel.
tf = any(strcmp(pipeline, {'coherence','correlation','reg_R2','reg_Rphase'}));
end

function save_pdf(fig, ps_file, pdf_file)
% Render fig via PS + ps2pdf (R2020-safe; matches plot_scanning_figures).
if isfile(ps_file),  delete(ps_file);  end
if isfile(pdf_file), delete(pdf_file); end
print(fig, '-dpsc', '-painters', ps_file);
close(fig);
[stat, msg] = system(sprintf('ps2pdf "%s" "%s"', ps_file, pdf_file));
if stat == 0 && isfile(pdf_file)
    delete(ps_file);
else
    warning('ps2pdf failed for %s (status %d): %s\nLeaving PS at %s', ...
        pdf_file, stat, strtrim(msg), ps_file);
end
end

%% ── Channel-average loaders ───────────────────────────────────────────
function [obsH1, obsH2, permH1, permH2, freq] = load_chan_avg(animal, pipeline, dv, base)
obsH1 = []; obsH2 = []; permH1 = []; permH2 = []; freq = [];
switch pipeline
    case 'coherence'
        if strcmp(dv,'hit_miss'), return; end   % no per-animal coh hit_miss
        H1 = fullfile(base, ['results_' animal], ...
            'phase_coherence','complex','cp10_till_100', dv, ...
            'all_loc_difflev', 'channel_avg_results.mat');
        H2 = fullfile(base, ['results_' animal], ...
            'phase_coherence','abs_per_pos','cp10_till_100', dv, ...
            'all_loc_difflev', 'channel_avg_results.mat');
        if ~isfile(H1) || ~isfile(H2), return; end
        S1 = load(H1, 'coh_chan_avg','coh_perm_chan_avg','freq');
        S2 = load(H2, 'coh_chan_avg','coh_perm_chan_avg');
        obsH1  = S1.coh_chan_avg;       obsH2  = S2.coh_chan_avg;
        permH1 = S1.coh_perm_chan_avg;  permH2 = S2.coh_perm_chan_avg;
        if isfield(S1,'freq'), freq = S1.freq; end

    case 'correlation'
        if strcmp(dv,'hit_miss')
            H1 = fullfile(base, ['results_' animal], ...
                'phase_correlation','complex','cp10_till_100', ...
                'hit_miss','all_loc_difflev','channel_avg_results_pos.mat');
            H2 = fullfile(base, ['results_' animal], ...
                'phase_correlation','abs_per_pos','cp10_till_100', ...
                'hit_miss','all_loc_difflev','channel_avg_results_pos.mat');
            if ~isfile(H1) || ~isfile(H2), return; end
            S1 = load(H1, 'pos_chan_avg','pos_perm_chan_avg','freq');
            S2 = load(H2, 'pos_chan_avg','pos_perm_chan_avg');
            obsH1  = S1.pos_chan_avg;       obsH2  = S2.pos_chan_avg;
            permH1 = S1.pos_perm_chan_avg;  permH2 = S2.pos_perm_chan_avg;
        else
            H1 = fullfile(base, ['results_' animal], ...
                'phase_correlation','complex','cp10_till_100', dv, ...
                'all_loc_difflev', 'channel_avg_results.mat');
            H2 = fullfile(base, ['results_' animal], ...
                'phase_correlation','abs_per_pos','cp10_till_100', dv, ...
                'all_loc_difflev', 'channel_avg_results.mat');
            if ~isfile(H1) || ~isfile(H2), return; end
            S1 = load(H1, 'corr_chan_avg','corr_perm_chan_avg','freq');
            S2 = load(H2, 'corr_chan_avg','corr_perm_chan_avg');
            obsH1  = S1.corr_chan_avg;       obsH2  = S2.corr_chan_avg;
            permH1 = S1.corr_perm_chan_avg;  permH2 = S2.corr_perm_chan_avg;
        end
        if isfield(S1,'freq'), freq = S1.freq; end

    case {'reg_R2','reg_Rphase'}
        dv_key = reg_dv_key(dv);
        H1 = fullfile(base, ['results_' animal], ...
            'multi_lin_reg','complex','cp10_till_100', ...
            'perm_R', dv_key, 'channel_avg_results.mat');
        H2 = fullfile(base, ['results_' animal], ...
            'multi_lin_reg','abs_per_pos','cp10_till_100', ...
            'perm_R_pos', dv_key, 'channel_avg_results.mat');
        if ~isfile(H1) || ~isfile(H2), return; end
        S1 = load(H1); S2 = load(H2);
        switch pipeline
            case 'reg_R2'
                if isfield(S1,'obs_avg') && isfield(S1.obs_avg,'phase')
                    obsH1 = S1.obs_avg.phase;
                end
                if isfield(S2,'obs_avg') && isfield(S2.obs_avg,'phase')
                    obsH2 = S2.obs_avg.phase;
                end
                if isfield(S1,'null_avg_R_phase_freq'), permH1 = S1.null_avg_R_phase_freq; end
                if isfield(S2,'null_avg_R_phase_freq'), permH2 = S2.null_avg_R_phase_freq; end
            case 'reg_Rphase'
                if isfield(S1,'obs_avg') && isfield(S1.obs_avg,'R_phase')
                    obsH1 = S1.obs_avg.R_phase;
                end
                if isfield(S2,'obs_avg') && isfield(S2.obs_avg,'R_phase')
                    obsH2 = S2.obs_avg.R_phase;
                end
                if isfield(S1,'null_avg_R_phase_mag_freq'), permH1 = S1.null_avg_R_phase_mag_freq; end
                if isfield(S2,'null_avg_R_phase_mag_freq'), permH2 = S2.null_avg_R_phase_mag_freq; end
        end
end
end

function k = reg_dv_key(dv)
switch dv
    case 'mua',      k = 'MUA_ERP_ampl_all';
    case 'lfp',      k = 'LFP_ERP_ampl_all';
    case 'RT',       k = 'RT';
    case 'hit_miss', k = 'hit_miss';
end
end

%% ── Per-channel paired test (dispatcher) ─────────────────────────────
function [pc_obs, pc_null_mean, pc_thr, pc_sig, freq] = ...
        paired_per_channel(animal, pipeline, dv, nCh, base)
switch pipeline
    case {'coherence','correlation'}
        [pc_obs, pc_null_mean, pc_thr, pc_sig, freq] = ...
            paired_per_channel_cohcorr(animal, pipeline, dv, nCh, base);
    case {'reg_R2','reg_Rphase'}
        [pc_obs, pc_null_mean, pc_thr, pc_sig, freq] = ...
            paired_per_channel_reg(animal, pipeline, dv, nCh, base);
    otherwise
        pc_obs = []; pc_null_mean = []; pc_thr = []; pc_sig = []; freq = [];
end
end

%% ── Coherence / correlation per-channel test ─────────────────────────
function [pc_obs, pc_null_mean, pc_thr, pc_sig, freq] = ...
        paired_per_channel_cohcorr(animal, pipeline, dv, nCh, base)
pc_obs = []; pc_null_mean = []; pc_thr = []; pc_sig = []; freq = [];

[h1_root, h2_root, obs_field, obs_filename, ...
    h1_perm_file, h1_perm_field, ...
    h2_perm_file, h2_perm_field] = perchan_paths(pipeline, dv, animal, base);
if isempty(h1_root) || isempty(h2_root), return; end

% Frequency axis
freq_file = fullfile(base, ['results_' animal], ...
    'multi_lin_reg','cp10_till_100','frequency.mat');
if isfile(freq_file), fS = load(freq_file); freq = fS.frequency; end
if isempty(freq), return; end
nFreq = numel(freq);

pc_obs       = nan(nCh, nFreq);
pc_null_mean = nan(nCh, nFreq);
pc_thr       = nan(nCh, 1);
pc_sig       = false(nCh, nFreq);

for ch = 1:nCh
    H1_obs_file  = fullfile(h1_root, num2str(ch), obs_filename);
    H2_obs_file  = fullfile(h2_root, num2str(ch), obs_filename);
    H1_perm_full = fullfile(h1_root, num2str(ch), h1_perm_file);
    H2_perm_full = fullfile(h2_root, num2str(ch), h2_perm_file);
    if ~isfile(H1_obs_file) || ~isfile(H2_obs_file), continue; end
    if ~isfile(H1_perm_full) || ~isfile(H2_perm_full), continue; end

    S1 = load(H1_obs_file,  obs_field);
    S2 = load(H2_obs_file,  obs_field);
    P1 = load(H1_perm_full, h1_perm_field);
    P2 = load(H2_perm_full, h2_perm_field);
    obs1 = S1.(obs_field)(:)';
    obs2 = S2.(obs_field)(:)';
    p1   = P1.(h1_perm_field);
    p2   = P2.(h2_perm_field);

    if numel(obs1) ~= nFreq || numel(obs2) ~= nFreq, continue; end
    if ~isequal(size(p1), size(p2)), continue; end
    % Normalize to nPerm × nFreq if stored the other way around.
    if size(p1, 2) ~= nFreq
        if size(p1, 1) == nFreq
            p1 = p1.'; p2 = p2.';
        else
            continue
        end
    end
    if any(isnan(obs1)) || any(isnan(obs2)) ...
            || any(isnan(p1(:))) || any(isnan(p2(:)))
        continue
    end

    diff_obs  = obs2 - obs1;          % 1 × nFreq
    diff_null = p2 - p1;               % nPerm × nFreq
    tmax      = max(diff_null, [], 2);
    thr       = quantile(tmax, 0.95);

    pc_obs(ch, :)       = diff_obs;
    pc_null_mean(ch, :) = mean(diff_null, 1, 'omitnan');
    pc_thr(ch)          = thr;
    pc_sig(ch, :)       = diff_obs >= thr;
end
end

function [h1_root, h2_root, obs_field, obs_filename, ...
        h1_perm_file, h1_perm_field, h2_perm_file, h2_perm_field] = ...
        perchan_paths(pipeline, dv, animal, base)
h1_root = ''; h2_root = ''; obs_field = ''; obs_filename = '';
h1_perm_file = ''; h1_perm_field = '';
h2_perm_file = ''; h2_perm_field = '';
switch pipeline
    case 'coherence'
        if strcmp(dv,'hit_miss'), return; end
        h1_root = fullfile(base, ['results_' animal], ...
            'phase_coherence','complex','cp10_till_100', dv, 'all_loc_difflev');
        h2_root = fullfile(base, ['results_' animal], ...
            'phase_coherence','abs_per_pos','cp10_till_100', dv, 'all_loc_difflev');
        obs_field    = 'coh';
        obs_filename = 'coherence.mat';
        h1_perm_file = 'coh_perm.mat';     h1_perm_field = 'coh_perm';
        h2_perm_file = 'coh_perm_pos.mat'; h2_perm_field = 'coh_perm_pos';
    case 'correlation'
        h1_root = fullfile(base, ['results_' animal], ...
            'phase_correlation','complex','cp10_till_100', dv, 'all_loc_difflev');
        h2_root = fullfile(base, ['results_' animal], ...
            'phase_correlation','abs_per_pos','cp10_till_100', dv, 'all_loc_difflev');
        if strcmp(dv,'hit_miss')
            obs_field    = 'pos';
            obs_filename = 'pos.mat';
            h1_perm_file = 'pos_perm.mat';     h1_perm_field = 'pos_perm';
            h2_perm_file = 'pos_perm_pos.mat'; h2_perm_field = 'pos_perm_pos';
        else
            obs_field    = 'correlation';
            obs_filename = 'correlation.mat';
            h1_perm_file = 'corr_perm.mat';     h1_perm_field = 'corr_perm';
            h2_perm_file = 'corr_perm_pos.mat'; h2_perm_field = 'corr_perm_pos';
        end
end
end

%% ── Regression per-channel test ──────────────────────────────────────
function [pc_obs, pc_null_mean, pc_thr, pc_sig, freq] = ...
        paired_per_channel_reg(animal, pipeline, dv, nCh, base)
% Per-channel paired H2-H1 test for regression R² or R_phase.
% Requires aggregate_regression_per_channel_nulls.m to have been run
% (which produces per_channel_null.mat per channel under perm_R[_pos]).
pc_obs = []; pc_null_mean = []; pc_thr = []; pc_sig = []; freq = [];

dv_key = reg_dv_key(dv);

% Observed: one file per (animal, hypothesis) holding all 64 channels.
H1_obs_file = fullfile(base, ['results_' animal], ...
    'multi_lin_reg','complex','cp10_till_100', ...
    'multi_regression_channelwise_R2.mat');
H2_obs_file = fullfile(base, ['results_' animal], ...
    'multi_lin_reg','abs_per_pos','cp10_till_100', ...
    'multi_regression_channelwise_R2_abs_per_pos.mat');
if ~isfile(H1_obs_file) || ~isfile(H2_obs_file), return; end

S1 = load(H1_obs_file, 'reg_results');
S2 = load(H2_obs_file, 'reg_results');
if ~isfield(S1, 'reg_results') || ~isfield(S2, 'reg_results'), return; end
if ~isfield(S1.reg_results, dv_key) || ~isfield(S2.reg_results, dv_key), return; end

switch pipeline
    case 'reg_R2'
        obs_field_name = 'R2_phase';
        null_field     = 'null_R2_phase';
    case 'reg_Rphase'
        obs_field_name = 'R_phase';
        null_field     = 'null_R_phase_mag';
end

obs1 = S1.reg_results.(dv_key).(obs_field_name);   % nCh × nFreq
obs2 = S2.reg_results.(dv_key).(obs_field_name);
if isempty(obs1) || isempty(obs2), return; end
[nCh_obs, nFreq] = size(obs1);
nCh = min(nCh, nCh_obs);

% Frequency axis
freq_file = fullfile(base, ['results_' animal], ...
    'multi_lin_reg','cp10_till_100','frequency.mat');
if isfile(freq_file), fS = load(freq_file); freq = fS.frequency; end

h1_root = fullfile(base, ['results_' animal], ...
    'multi_lin_reg','complex','cp10_till_100','perm_R', dv_key);
h2_root = fullfile(base, ['results_' animal], ...
    'multi_lin_reg','abs_per_pos','cp10_till_100','perm_R_pos', dv_key);

pc_obs       = nan(nCh, nFreq);
pc_null_mean = nan(nCh, nFreq);
pc_thr       = nan(nCh, 1);
pc_sig       = false(nCh, nFreq);

missing_aggregate = 0;
for ch = 1:nCh
    H1_null_file = fullfile(h1_root, num2str(ch), 'per_channel_null.mat');
    H2_null_file = fullfile(h2_root, num2str(ch), 'per_channel_null.mat');
    if ~isfile(H1_null_file) || ~isfile(H2_null_file)
        missing_aggregate = missing_aggregate + 1;
        continue
    end
    P1 = load(H1_null_file, null_field);
    P2 = load(H2_null_file, null_field);
    if ~isfield(P1, null_field) || ~isfield(P2, null_field), continue; end
    p1 = P1.(null_field);    % nPerm × nFreq
    p2 = P2.(null_field);
    if ~isequal(size(p1), size(p2)) || size(p1, 2) ~= nFreq, continue; end

    o1 = obs1(ch, :); o2 = obs2(ch, :);
    if all(isnan(o1)) || all(isnan(o2)), continue; end
    if any(isnan(p1(:))) && all(isnan(p1(:))), continue; end

    diff_obs  = o2 - o1;
    diff_null = p2 - p1;
    tmax      = max(diff_null, [], 2, 'omitnan');
    thr       = quantile(tmax, 0.95);

    pc_obs(ch, :)       = diff_obs;
    pc_null_mean(ch, :) = mean(diff_null, 1, 'omitnan');
    pc_thr(ch)          = thr;
    pc_sig(ch, :)       = diff_obs >= thr;
end

if missing_aggregate == nCh
    % No per-channel aggregates exist yet → return empty to skip plot
    pc_obs = []; pc_null_mean = []; pc_thr = []; pc_sig = [];
    warning(['Per-channel regression aggregates not found for ' ...
        '%s / %s / %s. Run aggregate_regression_per_channel_nulls.m first.'], ...
        animal, pipeline, dv);
elseif missing_aggregate > 0
    fprintf('    note: %s/%s/%s — %d/%d channels missing per_channel_null.mat\n', ...
        animal, pipeline, dv, missing_aggregate, nCh);
end
end

%% ── Figure builders ──────────────────────────────────────────────────
function fig = fig_chan_avg_grid(R, pipe_label, dv_label, col, sgtitle_str)
nP = numel(pipe_label); nD = numel(dv_label);
fig = figure('Visible','off','Units','centimeters','Position',[1 1 36 24]);
set(fig,'PaperUnits','centimeters','PaperSize',fig.Position(3:4), ...
    'PaperPosition',[0 0 fig.Position(3:4)]);
for d = 1:nD
    for p = 1:nP
        ax = subplot(nD, nP, (d-1)*nP + p); hold(ax,'on');
        S = R{1, p, d};
        if isempty(S)
            title(ax, sprintf('%s — %s (no data)', pipe_label{p}, dv_label{d}), 'FontSize',7);
            set(ax,'XTick',[],'YTick',[]); continue
        end
        yline(ax, 0, 'k-','LineWidth',0.4,'HandleVisibility','off');
        plot(ax, S.freq, S.obs,       'Color',col,'LineWidth',1.5, ...
            'DisplayName','obs H2−H1');
        plot(ax, S.freq, S.null_mean, '--','Color',col,'LineWidth',0.8, ...
            'DisplayName','null mean (Jensen)');
        if isfinite(S.thr)
            yline(ax, S.thr, ':','Color',col,'LineWidth',1, ...
                'DisplayName','95% threshold');
        end
        yl = ylim(ax); span = yl(2)-yl(1);
        if any(S.sig)
            scatter(ax, S.freq(S.sig), repmat(yl(2)-0.05*span,1,sum(S.sig)), ...
                8, col, 'filled','HandleVisibility','off');
        end
        ylim(ax, yl);
        if d == nD, xlabel(ax,'Freq (Hz)','FontSize',7); end
        if p == 1,  ylabel(ax,sprintf('%s\n\\Delta', dv_label{d}),'FontSize',7); end
        if d == 1, title(ax, pipe_label{p}, 'FontSize',8); end
        if d == 1 && p == 1, legend(ax,'Location','best','FontSize',5); end
        set(ax,'FontSize',6,'Box','on');
    end
end
sgtitle(sgtitle_str, 'FontSize',10,'FontWeight','bold');
end

function fig = fig_nsig_grid(NS, freq, pipe_label, dv_label, col, sgtitle_str)
nP = numel(pipe_label); nD = numel(dv_label);
fig = figure('Visible','off','Units','centimeters','Position',[1 1 36 22]);
set(fig,'PaperUnits','centimeters','PaperSize',fig.Position(3:4), ...
    'PaperPosition',[0 0 fig.Position(3:4)]);
for d = 1:nD
    for p = 1:nP
        ax = subplot(nD, nP, (d-1)*nP + p); hold(ax,'on');
        v = NS{1, p, d};
        if isempty(v) || isempty(freq)
            title(ax, sprintf('%s — %s (no data)', pipe_label{p}, dv_label{d}), 'FontSize',7);
            set(ax,'XTick',[],'YTick',[]); continue
        end
        bar(ax, freq, v, 'FaceColor', col, 'EdgeColor','none');
        ylim(ax, [0 max(max(v(:)),1)*1.1]);
        if d == nD, xlabel(ax,'Freq (Hz)','FontSize',7); end
        if p == 1,  ylabel(ax,sprintf('%s\n# sig ch', dv_label{d}),'FontSize',7); end
        if d == 1, title(ax, pipe_label{p}, 'FontSize',8); end
        set(ax,'FontSize',6,'Box','on');
    end
end
sgtitle(sgtitle_str, 'FontSize',10,'FontWeight','bold');
end

function fig = fig_per_channel_heatmap_grid(D, pipe_label, dv_label, sgtitle_str)
% 4×4 grid of channel × freq heatmaps. Cells where the paired test was
% not significant are rendered with reduced alpha so significant
% (channel, freq) bins stand out without losing the underlying H2-H1
% magnitudes.
nP = numel(pipe_label); nD = numel(dv_label);
fig = figure('Visible','off','Units','centimeters','Position',[1 1 36 24]);
set(fig,'PaperUnits','centimeters','PaperSize',fig.Position(3:4), ...
    'PaperPosition',[0 0 fig.Position(3:4)]);

alpha_sig    = 1.0;
alpha_nonsig = 0.20;

for d = 1:nD
    for p = 1:nP
        ax = subplot(nD, nP, (d-1)*nP + p);
        S = D{1, p, d};
        if isempty(S) || ~isfield(S,'obs') || isempty(S.obs) ...
                || isempty(S.freq)
            title(ax, sprintf('%s — %s (no data)', pipe_label{p}, dv_label{d}), 'FontSize',7);
            set(ax,'XTick',[],'YTick',[]); continue
        end

        nCh = size(S.obs, 1);
        img = imagesc(ax, S.freq, 1:nCh, S.obs);
        set(ax,'YDir','normal','Color',[1 1 1]);
        colormap(ax, parula);

        v = S.obs(~isnan(S.obs));
        if ~isempty(v)
            mx = max(abs(prctile(v,[2 98])));
            if mx > 0, caxis(ax, [-mx mx]); end
        end

        A = alpha_nonsig * ones(size(S.obs));
        if ~isempty(S.sig), A(S.sig) = alpha_sig; end
        A(isnan(S.obs)) = 0;
        set(img, 'AlphaData', A);

        cb = colorbar(ax); cb.FontSize = 5;
        if d == nD, xlabel(ax,'Freq (Hz)','FontSize',7); end
        if p == 1,  ylabel(ax,sprintf('%s\nChannel', dv_label{d}),'FontSize',7); end
        if d == 1,  title(ax, pipe_label{p}, 'FontSize',8); end
        set(ax,'FontSize',6,'Box','on');
        xlim(ax, [min(S.freq) max(S.freq)]);
        ylim(ax, [0.5 nCh+0.5]);
    end
end
sgtitle(sgtitle_str, 'FontSize',10,'FontWeight','bold');
end

function fig = fig_per_channel_grid(pc_obs, pc_null_mean, pc_thr, pc_sig, freq, col, sgtitle_str)
nCh  = size(pc_obs, 1);
cols = ceil(sqrt(nCh)); rows = ceil(nCh/cols);
fig  = figure('Visible','off','Units','centimeters', ...
    'Position',[1 1 max(36, 2.4*cols) max(24, 2.4*rows)]);
set(fig,'PaperUnits','centimeters','PaperSize',fig.Position(3:4), ...
    'PaperPosition',[0 0 fig.Position(3:4)]);

all_vals = [pc_obs(:); pc_thr(:); 0];
all_vals = all_vals(~isnan(all_vals));
if isempty(all_vals), all_vals = [0 1]; end
ymin = min(all_vals); ymax = max(all_vals);
if ymin == ymax, ymax = ymin + 1; end
pad  = 0.05*(ymax-ymin);
yl   = [ymin-pad ymax+pad];

for ch = 1:nCh
    ax = subplot(rows, cols, ch); hold(ax,'on');
    obs = pc_obs(ch,:); nm = pc_null_mean(ch,:);
    t = pc_thr(ch);     sig = pc_sig(ch,:);
    if all(isnan(obs))
        title(ax, sprintf('ch%d (no data)', ch), 'FontSize',6);
        set(ax,'XTick',[],'YTick',[]); continue
    end
    if any(sig), shade_sig(ax, freq, sig, yl, col); end
    yline(ax, 0, 'k-','LineWidth',0.3,'HandleVisibility','off');
    plot(ax, freq, obs, 'Color', col, 'LineWidth',1);
    plot(ax, freq, nm,  '--', 'Color', col, 'LineWidth',0.5);
    if ~isnan(t)
        yline(ax, t, ':','Color',col,'LineWidth',0.5);
    end
    ylim(ax, yl); xlim(ax, [min(freq) max(freq)]);
    title(ax, sprintf('ch%d  n_{sig}=%d', ch, sum(sig)), 'FontSize',6);
    ax.FontSize = 5;
    if mod(ch-1, cols) ~= 0, ax.YTickLabel = {}; end
    if ch <= (rows-1)*cols,  ax.XTickLabel = {}; end
end
sgtitle(sgtitle_str, 'FontSize',9, 'FontWeight','bold');
end

%% ── Diagnostics: usable channels + (pos × diff) trial counts ─────────
function write_diagnostics(fname, animalName, dvs, pipe_label, dv_label, pc_data_animal, base)
% pc_data_animal: nP × nD cell of structs with .obs (nCh × nFreq).
% Writes two tables:
%   (1) usable per-channel rows per (pipeline × DV) — a channel is
%       counted as usable if any freq bin in its row is non-NaN.
%   (2) trial counts per (position × difficulty bin) per DV, using the
%       same within-position quantile binning the regression pipeline
%       uses (nDiffBins = 4, col 18 = difficulty, col 16 = position).
fid = fopen(fname, 'w');
if fid < 0, return; end
cleanup = onCleanup(@() fclose(fid));

nP = numel(pipe_label); nD = numel(dv_label);

fprintf(fid, '=== %s: usable per-channel rows in paired H2-H1 test ===\n', animalName);
fprintf(fid, '(channel counted as usable if any freq bin has non-NaN obs)\n\n');
fprintf(fid, '%-14s', 'pipeline\dv');
for d = 1:nD, fprintf(fid, '%-12s', dv_label{d}); end
fprintf(fid, '\n');
for p = 1:nP
    fprintf(fid, '%-14s', pipe_label{p});
    for d = 1:nD
        S = pc_data_animal{p, d};
        if isempty(S) || ~isfield(S,'obs') || isempty(S.obs)
            fprintf(fid, '%-12s', '—');
        else
            usable = sum(~all(isnan(S.obs), 2));
            fprintf(fid, '%-12s', sprintf('%d/%d', usable, size(S.obs,1)));
        end
    end
    fprintf(fid, '\n');
end

ph_file = fullfile(base, ['results_' animalName], ...
    'multi_lin_reg','cp10_till_100','ph_all_sess.mat');
if ~isfile(ph_file)
    fprintf(fid, '\n(ph_all_sess.mat not found — skipping trial counts)\n');
    return
end
S = load(ph_file, 'ph_comb');
if ~isfield(S, 'ph_comb')
    fprintf(fid, '\n(ph_comb not in ph_all_sess.mat — skipping trial counts)\n');
    return
end
ph_comb   = S.ph_comb;
nDiffBins = 4;

dv_to_ti = containers.Map( ...
    {'mua','lfp','RT','hit_miss'}, ...
    {'MUA_ERP_trialinfo','LFP_ERP_trialinfo','RT_trialinfo','trialinfo'});

fprintf(fid, '\n=== %s: trial counts per (position × difficulty bin) ===\n', animalName);
fprintf(fid, '(difficulty col 18 binned into %d within-position quantiles)\n', nDiffBins);

for d = 1:nD
    field = dv_to_ti(dvs{d});
    if ~isfield(ph_comb, field) || isempty(ph_comb.(field))
        fprintf(fid, '\n-- %s: %s not in ph_comb --\n', dv_label{d}, field);
        continue
    end
    ti = ph_comb.(field);
    if size(ti,2) < 18
        fprintf(fid, '\n-- %s: trialinfo has only %d cols, need >=18 --\n', ...
            dv_label{d}, size(ti,2));
        continue
    end
    pos       = ti(:,16);
    diff_vals = ti(:,18);
    positions = unique(pos(~isnan(pos)));
    nPos      = numel(positions);
    bin_idx   = bin_difficulty_inline(diff_vals, pos, positions, nDiffBins);

    M = zeros(nPos, nDiffBins);
    for ip = 1:nPos
        for ib = 1:nDiffBins
            M(ip, ib) = sum(pos == positions(ip) & bin_idx == ib);
        end
    end

    fprintf(fid, '\n-- %s (n_trials = %d, n_pos = %d) --\n', ...
        dv_label{d}, size(ti,1), nPos);
    fprintf(fid, '%-8s', 'pos\bin');
    for ib = 1:nDiffBins, fprintf(fid, '%-8s', sprintf('d%d', ib)); end
    fprintf(fid, '%-8s\n', 'total');
    for ip = 1:nPos
        fprintf(fid, '%-8s', sprintf('p%g', positions(ip)));
        for ib = 1:nDiffBins, fprintf(fid, '%-8d', M(ip, ib)); end
        fprintf(fid, '%-8d\n', sum(M(ip,:)));
    end
    fprintf(fid, '%-8s', 'min');
    for ib = 1:nDiffBins, fprintf(fid, '%-8d', min(M(:,ib))); end
    fprintf(fid, '%-8d\n', min(sum(M,2)));
end
end

function bin_idx = bin_difficulty_inline(diff_vals, pos_labels, positions, nDiffBins)
% Inline copy of multiple_linear_reg/functions/bin_difficulty_per_pos.m so
% this script does not depend on that path being added.
bin_idx = nan(size(diff_vals));
for p = 1:numel(positions)
    pos_mask = pos_labels == positions(p);
    de = diff_vals(pos_mask); valid = ~isnan(de);
    if sum(valid) < nDiffBins, continue; end
    edges      = quantile(de(valid), linspace(0, 1, nDiffBins + 1));
    edges(1)   = -Inf;
    edges(end) =  Inf;
    [~,~,bin_local] = histcounts(de, edges);
    bin_local(bin_local == 0) = NaN;
    pos_idx = find(pos_mask);
    bin_idx(pos_idx) = bin_local;
end
end

function shade_sig(ax, freq, sig, yl, col)
sig    = logical(sig(:)');
starts = find(diff([0 sig]) ==  1);
ends   = find(diff([sig 0]) == -1);
for k = 1:numel(starts)
    xp = [freq(starts(k)) freq(ends(k)) freq(ends(k)) freq(starts(k))];
    yp = [yl(1) yl(1) yl(2) yl(2)];
    fill(ax, xp, yp, col, 'FaceAlpha',0.15,'EdgeColor','none', ...
        'HandleVisibility','off');
end
end
