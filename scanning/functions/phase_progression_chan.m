function out = phase_progression_chan(cfg_fun)
% Per-channel phase-progression worker for SLURM parallelisation.
% Mirrors the per-channel block of phase_progression_per_chan.m so that
% the systematicity stats can be computed for one channel × DV at a time
% on a single SLURM job.
%
% cfg_fun fields:
%   ichan   — channel index
%   nPerm   — number of label-shuffle permutations
%   dv      — DV type: 'mua' | 'lfp' | 'RT' | 'hit_miss'
%   infile  — folder containing ph_all_sess.mat
%   outfile — folder where a <ichan>/ subdir will be created
%   seed    — RNG seed for this channel (per-channel reproducibility)
%   perm_seed_base — (optional) base seed for the SYNCHRONISED permutation
%                    null; perm k is seeded rng(perm_seed_base + k) in EVERY
%                    channel so the channel-average null keeps cross-channel
%                    dependence. Default 2025 (matches regress_perm_R_pos.m).
%
% Saves <outfile>/<ichan>/phase_progression_chan.mat with:
%   R_obs (1×nFreq), R_null (nFreq×nPerm), pref_phase (nFreq×nPos),
%   coh_mag (nFreq×nPos), step_phase (nFreq×nPos-1), mean_step (1×nFreq),
%   n_p (1×nPos), p_val (1×nFreq), positions (nPos×1), seed.

ichan   = cfg_fun.ichan;
nPerm   = cfg_fun.nPerm;
dv      = cfg_fun.dv;
infile  = cfg_fun.infile;
outfile = cfg_fun.outfile;
seed    = cfg_fun.seed;
if isfield(cfg_fun,'perm_seed_base'), perm_seed_base = cfg_fun.perm_seed_base;
else,                                 perm_seed_base = 2025; end

rng(seed);

S = load(fullfile(infile, 'ph_all_sess.mat'));
ph_comb = S.ph_comb;

[Y_full, trlInfo, keepIdx, isPerCh] = local_get_dv_data(ph_comb, dv);

positions = unique(trlInfo(keepIdx, 16));
nPos      = numel(positions);

nTrials = size(ph_comb.phase_all, 1);
nFreq   = size(ph_comb.phase_all, 2);

pos_idx_all = zeros(nTrials, 1);
for p = 1:nPos
    pos_idx_all(trlInfo(:,16) == positions(p) & keepIdx) = p;
end

if isPerCh, y_ch = Y_full(:, ichan);
else,       y_ch = Y_full;
end
phs_ch = ph_comb.phase_all(:, :, ichan);

valid = ~isnan(y_ch) & ~any(isnan(phs_ch), 2) & (pos_idx_all > 0);

R_obs      = nan(1, nFreq);
R_null     = nan(nFreq, nPerm);
pref_phase = nan(nFreq, nPos);
coh_mag    = nan(nFreq, nPos);
step_phase = nan(nFreq, max(nPos-1, 0));
mean_step  = nan(1, nFreq);
n_p        = zeros(1, nPos);
p_val      = nan(1, nFreq);

if nPos >= 2 && sum(valid) >= 2*nPos
    y_v   = y_ch(valid);
    phs_v = phs_ch(valid, :);
    pos_v = pos_idx_all(valid);

    % weight(t,f) = y(t) · exp(i·phase(t,f))
    weight = bsxfun(@times, exp(1i .* phs_v), y_v);

    [R_o, phi_pos, mag_pos, n_p, dphi_o, mstep_o] = ...
        local_systematicity(weight, pos_v, nPos);
    R_obs      = R_o;
    pref_phase = phi_pos;
    coh_mag    = mag_pos;
    step_phase = dphi_o;
    mean_step  = mstep_o;

    % --- Synchronised permutation null (matches the sampling pipeline) ---
    % Seed each shuffle by the PERMUTATION INDEX (not the channel) and apply
    % the SAME position-label relabelling to every channel. base_trials and
    % base_labels are identical across channels for a given DV, and the RNG
    % is seeded by perm index, so perm k is the same trial->position
    % reassignment in every channel. Averaging this null across channels
    % therefore preserves the spatial (cross-channel) dependence of the
    % array LFP. Without this — i.e. shuffling each channel independently —
    % the channel-average null is far too tight (its spread shrinks like
    % 1/sqrt(nCh) assuming channel independence, which is false for array
    % LFP) and the channel-average significance is anti-conservative.
    % (cf. multiple_linear_reg/functions/regress_perm_R_pos.m: rng(2025+perm_idx).)
    base_trials = find(pos_idx_all > 0);     % identical across channels for this DV
    base_labels = pos_idx_all(base_trials);
    for k = 1:nPerm
        rng(perm_seed_base + k);             % SHARED across channels
        pos_perm              = zeros(nTrials, 1);
        pos_perm(base_trials) = base_labels(randperm(numel(base_labels)));
        pos_k                 = pos_perm(valid);   % valid is a subset of base_trials -> labels in 1..nPos
        R_null(:, k) = local_systematicity(weight, pos_k, nPos).';
    end

    p_val = mean(bsxfun(@ge, R_null, R_obs(:)), 2).';
end

chan_dir = fullfile(outfile, num2str(ichan));
if ~exist(chan_dir, 'dir'), mkdir(chan_dir); end
save(fullfile(chan_dir, 'phase_progression_chan.mat'), ...
    'R_obs','R_null','pref_phase','coh_mag','step_phase','mean_step', ...
    'n_p','p_val','positions','seed','-v7.3');

out = 1;
end

% =====================================================================
% Local helpers — duplicated from phase_progression_per_chan.m so the
% SLURM worker is self-contained.
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

function [R, phi_pos, mag_pos, n_p, dphi_wrapped, mean_step] = ...
        local_systematicity(weight, pos_k, nPos)
nFreq   = size(weight, 2);
phi_pos = nan(nFreq, nPos);
mag_pos = nan(nFreq, nPos);
n_p     = zeros(1, nPos);

for p = 1:nPos
    m = pos_k == p;
    n_p(p) = sum(m);
    if n_p(p) < 2, continue; end
    c = mean(weight(m,:), 1);
    phi_pos(:,p) = angle(c).';
    mag_pos(:,p) = abs(c).';
end

x   = (1:nPos)';
x_c = x - mean(x);
x_s = sqrt(mean(x_c.^2));

S = sin(phi_pos);
C = cos(phi_pos);
S_c = bsxfun(@minus, S, mean(S, 2));
C_c = bsxfun(@minus, C, mean(C, 2));
S_s = sqrt(mean(S_c.^2, 2));
C_s = sqrt(mean(C_c.^2, 2));

rxs = (S_c * x_c) ./ (nPos .* S_s .* x_s);
rxc = (C_c * x_c) ./ (nPos .* C_s .* x_s);
rcs = mean(S_c .* C_c, 2) ./ (S_s .* C_s);

denom = 1 - rcs.^2;
inner = (rxc.^2 + rxs.^2 - 2 .* rxc .* rxs .* rcs) ./ denom;
inner(denom <= eps | isnan(denom)) = NaN;
R = sqrt(max(0, inner)).';

dphi         = diff(phi_pos, 1, 2);
dphi_wrapped = angle(exp(1i * dphi));
zbar         = mean(exp(1i * dphi), 2);
mean_step    = angle(zbar).';
end
