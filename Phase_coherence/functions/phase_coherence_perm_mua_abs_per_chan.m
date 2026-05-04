function phase_coherence_perm_mua_abs_per_chan(cfg_fun)
% H1+H4 permutation: pool all trials in complex space within a channel,
% then take abs() per channel. The cross-channel and cross-animal averages
% are taken on these magnitudes downstream (arithmetic mean of |coh|).
%
% Saved variable: coh_perm  [permut_n x nFreq], real-valued magnitudes.

ichan        = cfg_fun.ichan;
permut_n     = cfg_fun.permut_n;
perm_indices = cfg_fun.perm_indices;   % shared permutations
trial_idx    = cfg_fun.trial_idx;

cd(cfg_fun.infile)
load('ph_all_sess.mat')

phase   = ph_comb.phase_all(trial_idx, :, ichan);       % [trials x freq]
erp_amp = ph_comb.MUA_ERP_ampl_all(trial_idx, ichan);   % [trials x 1]

nFreq    = size(phase, 2);
coh_perm = nan(permut_n, nFreq);

for perm = 1:permut_n
    erp_p = erp_amp(perm_indices{perm});
    for foi = 1:nFreq
        vec = exp(1i * phase(:,foi)) .* erp_p(:);
        coh_perm(perm,foi) = abs(mean(vec));   % abs at channel level (H4)
    end
end

cd(cfg_fun.outfile)
if ~exist(num2str(ichan), 'dir'), mkdir(num2str(ichan)); end
cd(num2str(ichan))

ESIsave coh_perm coh_perm
ESIsave perm_indices perm_indices
end
