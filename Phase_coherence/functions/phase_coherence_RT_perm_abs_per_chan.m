function phase_coherence_RT_perm_abs_per_chan(cfg_fun)
% H1+H4 permutation for phase coherence with RT.
% Pool all hit trials in complex space within a channel, then take abs()
% per channel. Cross-channel and cross-animal averaging is on magnitudes.
%
% Saved variable: coh_perm  [permut_n x nFreq], real-valued magnitudes.

ichan        = cfg_fun.ichan;
permut_n     = cfg_fun.permut_n;
perm_indices = cfg_fun.perm_indices;   % shared permutations
trial_idx    = cfg_fun.trial_idx;      % hit trials with valid RT

cd(cfg_fun.infile)
load('ph_all_sess.mat')

phase = ph_comb.phase_all(trial_idx, :, ichan);   % [nHits x nFreq]
rt    = ph_comb.RT(trial_idx, ichan);              % [nHits x 1]

if isempty(rt)
    fprintf('Channel %d skipped (no valid RT)\n', ichan);
    return
end

nFreq    = size(phase, 2);
coh_perm = nan(permut_n, nFreq);

for perm = 1:permut_n
    rt_perm    = rt(perm_indices{perm});
    valid_perm = ~isnan(rt_perm);   % shuffled NaNs may land in valid positions
    rt_clean   = rt_perm(valid_perm);
    ph_clean   = phase(valid_perm, :);

    for foi = 1:nFreq
        vec = exp(1i * ph_clean(:,foi)) .* rt_clean(:);
        coh_perm(perm,foi) = abs(mean(vec));   % abs at channel level (H4)
    end
end

cd(cfg_fun.outfile)
if ~exist(num2str(ichan), 'dir'), mkdir(num2str(ichan)); end
cd(num2str(ichan))

ESIsave coh_perm coh_perm
ESIsave perm_indices perm_indices
end
