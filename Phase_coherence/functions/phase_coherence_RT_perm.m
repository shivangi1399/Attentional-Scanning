function phase_coherence_RT_perm(cfg_fun)
% phase_coherence_RT_perm  Permutation test for phase coherence with RT.
%
% Shuffles RT across trials, computes phase coherence per freq.

ichan = cfg_fun.ichan;
permut_n = cfg_fun.permut_n;
perm_indices = cfg_fun.perm_indices;
trial_idx = cfg_fun.trial_idx;  % hit trials with valid RT

cd(cfg_fun.infile)
load('ph_all_sess.mat')

phase = ph_comb.phase_all(trial_idx, :, ichan);  % [nHits x nFreq]
rt    = ph_comb.RT(trial_idx, ichan);             % [nHits x 1]

% Remove NaN RT trials (safety check)
valid = ~isnan(rt);
phase = phase(valid, :);
rt    = rt(valid);

% Skip channel if no valid RT values
if isempty(rt)
    fprintf('Channel %d skipped (no valid RT)\n', ichan);
    continue
end

nFreq = size(phase, 2);
coh_perm        = nan(permut_n, nFreq);
phase_spec_perm = nan(permut_n, nFreq);

for perm = 1:permut_n
    rt_perm = rt(perm_indices{perm});  % shuffle RT

    for foi = 1:nFreq
        vec = exp(1i * phase(:, foi)) .* rt_perm(:);
        cavg = mean(vec);
        coh_perm(perm, foi)        = abs(cavg);
        phase_spec_perm(perm, foi) = angle(cavg);
    end
end

cd(cfg_fun.outfile)
if ~exist(num2str(ichan), 'dir')
    mkdir(num2str(ichan))
end
cd(num2str(ichan))

ESIsave coh_perm coh_perm
ESIsave phase_spec_perm phase_spec_perm
ESIsave perm_indices perm_indices
end
