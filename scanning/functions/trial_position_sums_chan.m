function out = trial_position_sums_chan(cfg_fun)
% Per-channel worker for the TRIAL-LEVEL de-rotation estimator.
%
% The estimator asks the de-rotation question on the raw trials instead of on
% the per-position preferred phases:
%     c(f,v) = ( 1/W ) * SUM over ALL TRIALS  y_t * exp( i*( phi_tf - k*d_p(t) ) )
% i.e. every trial is rotated by the wave model's prediction for the location
% that trial had, and one phase coherence is computed over the whole trial set.
% At k=0 this is exactly the phase_coherence/ measure with all locations pooled.
%
% The de-rotation depends on a trial only through its location, so the whole
% estimator collapses onto per-location complex sums:
%     S(p,f) = SUM over trials t at location p of  y_t * exp(i*phi_tf)
%     c(f,v) = ( 1/W ) * SUM over p  S(p,f) * exp(-i*k*d_p)
% so this worker needs no speed grid, geometry or mode: it produces S(p,f) for
% the observed labels and every permutation, and the caller sweeps (frequency x
% speed) and the three modes cheaply from S. That is what makes the trial-level
% version affordable.
%
% Against the phase-progression (level-2) estimator: that one uses the
% per-location MEAN c_p = S(p,:)/n_p and weights each location by |c_p|, which
% goes as 1/sqrt(n_p) under noise, so a location with FEWER trials gets a
% LARGER weight. The raw sums here make every trial count exactly once.
%
% cfg_fun fields:
%   ichan          - channel index
%   nPerm          - number of position-label shuffles
%   dv             - 'mua' | 'lfp' | 'RT' | 'hit_miss'
%   infile         - folder containing ph_all_sess.mat
%   outfile        - folder; a <ichan>/ subdir is created
%   perm_seed_base - base seed. Permutation k is seeded rng(perm_seed_base+k)
%                    in EVERY channel, so the same trial->location relabelling
%                    is applied across channels. This keeps the cross-channel
%                    dependence of array LFP in the null, which the channel
%                    average and the coherent mode both rely on. Default 2025.
%
% Saves <outfile>/<ichan>/trial_position_sums.mat with:
%   S_obs   [nPos x nFreq]          complex, observed labels
%   S_perm  [nPos x nFreq x nPerm]  complex single, permuted labels
%   W       scalar                  sum |y| over the valid trials (normaliser)
%   n_p     [1 x nPos]              trials per location
%   positions, nValid, seed

ichan   = cfg_fun.ichan;
nPerm   = cfg_fun.nPerm;
dv      = cfg_fun.dv;
infile  = cfg_fun.infile;
outfile = cfg_fun.outfile;
if isfield(cfg_fun,'perm_seed_base'), perm_seed_base = cfg_fun.perm_seed_base;
else,                                 perm_seed_base = 2025; end

S = load(fullfile(infile, 'ph_all_sess.mat'));
ph_comb = S.ph_comb;

[Y_full, trlInfo, keepIdx, isPerCh] = local_get_dv_data(ph_comb, dv);

positions = unique(trlInfo(keepIdx, 16));
nPos      = numel(positions);
nTrials   = size(ph_comb.phase_all, 1);
nFreq     = size(ph_comb.phase_all, 2);

pos_idx_all = zeros(nTrials, 1);
for p = 1:nPos
    pos_idx_all(trlInfo(:,16) == positions(p) & keepIdx) = p;
end

if isPerCh, y_ch = Y_full(:, ichan);
else,       y_ch = Y_full;
end
phs_ch = ph_comb.phase_all(:, :, ichan);
valid  = ~isnan(y_ch) & ~any(isnan(phs_ch), 2) & (pos_idx_all > 0);

S_obs  = complex(nan(nPos, nFreq));
S_perm = complex(nan(nPos, nFreq, nPerm, 'single'));
W      = 0; n_p = zeros(1, nPos); nValid = sum(valid);

if nValid >= 2*nPos
    y_v   = y_ch(valid);
    phs_v = phs_ch(valid, :);
    pos_v = pos_idx_all(valid);
    nV    = numel(y_v);

    % weight(t,f) = y_t * exp(i*phi_tf) — one row per surviving trial
    weight = bsxfun(@times, exp(1i .* phs_v), y_v);
    W      = sum(abs(y_v));                 % normaliser: |c| <= 1 always
    n_p    = accumarray(pos_v, 1, [nPos 1]).';

    % Per-location sums via a sparse indicator matrix: S = P' * weight.
    P     = sparse(1:nV, pos_v, 1, nV, nPos);
    S_obs = full(P.' * weight);

    % Synchronised permutation null: seed by PERMUTATION INDEX, not channel, so
    % perm k is the same relabelling in every channel (see header).
    base_trials = find(pos_idx_all > 0);
    base_labels = pos_idx_all(base_trials);
    for k = 1:nPerm
        rng(perm_seed_base + k);
        pos_perm              = zeros(nTrials, 1);
        pos_perm(base_trials) = base_labels(randperm(numel(base_labels)));
        pos_k                 = pos_perm(valid);
        Pk                    = sparse(1:nV, pos_k, 1, nV, nPos);
        S_perm(:,:,k)         = single(full(Pk.' * weight));
    end
end

chan_dir = fullfile(outfile, num2str(ichan));
if ~exist(chan_dir, 'dir'), mkdir(chan_dir); end
save(fullfile(chan_dir, 'trial_position_sums.mat'), ...
    'S_obs','S_perm','W','n_p','positions','nValid','-v7.3');

out = 1;
end

% =====================================================================
% Local helper — duplicated from phase_progression_chan.m so the SLURM
% worker stays self-contained.
% =====================================================================
function [Y, trlInfo, keepIdx, isPerCh] = local_get_dv_data(ph_comb, dv)
switch dv
    case 'mua'
        Y       = ph_comb.MUA_ERP_ampl_all;
        trlInfo = ph_comb.trialinfo;
        keepIdx = true(size(Y,1), 1);
        isPerCh = true;
    case 'lfp'
        Y       = ph_comb.LFP_ERP_ampl_all;
        trlInfo = ph_comb.trialinfo;
        keepIdx = true(size(Y,1), 1);
        isPerCh = true;
    case 'RT'
        Y       = ph_comb.RT;
        trlInfo = ph_comb.RT_trialinfo;
        keepIdx = trlInfo(:,20) == 1;
        isPerCh = true;
    case 'hit_miss'
        Y           = ph_comb.trialinfo(:,20);
        Y(Y == 5)   = 0;
        trlInfo     = ph_comb.trialinfo;
        keepIdx     = true(size(Y,1), 1);
        isPerCh     = false;
    otherwise
        error('Unknown DV: %s', dv);
end
end
