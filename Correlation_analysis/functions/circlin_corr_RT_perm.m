function corr_perm = circlin_corr_RT_perm(cfg_fun)
% circlin_corr_RT_perm  Permutation test for circular-linear correlation
%                       between phase and reaction time (RT).
%
% Shuffles RT across trials, computes circ_corrcl(phase, RT_perm) per freq.

ichan = cfg_fun.ichan;
permut_n = cfg_fun.permut_n;
perm_indices = cfg_fun.perm_indices;  % shared permutations
trial_idx = cfg_fun.trial_idx;        % hit trials with valid RT

cd(cfg_fun.infile)
load('ph_all_sess.mat')

phase  = ph_comb.phase_all(trial_idx, :, ichan);  % [nHits x nFreq]
rt     = ph_comb.RT(trial_idx, ichan);             % [nHits x 1]

% Remove NaN RT trials (safety check)
valid = ~isnan(rt);
phase = phase(valid, :);
rt    = rt(valid);

if isempty(rt)
    fprintf('Channel %d skipped (no valid RT)\n', ichan);
    return
end

nFreq = size(phase, 2);
corr_perm   = nan(permut_n, nFreq);
pvalue_perm = nan(permut_n, nFreq);

for perm = 1:permut_n
    rt_perm = rt(perm_indices{perm});  % shuffle RT
    for foi = 1:nFreq
        [corr_perm(perm, foi), pvalue_perm(perm, foi)] = circ_corrcl(phase(:, foi), rt_perm);
    end
end

cd(cfg_fun.outfile)
if ~exist(num2str(ichan), 'dir')
    mkdir(num2str(ichan))
end
cd(num2str(ichan))

ESIsave corr_perm corr_perm
ESIsave pvalue_perm pvalue_perm
ESIsave perm_indices perm_indices
end
