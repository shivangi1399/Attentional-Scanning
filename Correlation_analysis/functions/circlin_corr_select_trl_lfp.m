function [corr_perm] = circlin_corr_select_trl_lfp(cfg_fun)

ichan = cfg_fun.ichan;
permut_n = cfg_fun.permut_n;
perm_indices = cfg_fun.perm_indices;  % shared permutations
trial_idx = cfg_fun.trial_idx; %selected trials

cd(cfg_fun.infile)
load('ph_all_sess.mat')

phase = ph_comb.phase_all(trial_idx,:,ichan);    % [trials × freq]
erp_amp = ph_comb.LFP_ERP_ampl_all(trial_idx,ichan); % [trials × 1]

nFreq = size(phase, 2);

corr_perm = nan(permut_n, nFreq);
pvalue_perm = nan(permut_n, nFreq);

for perm = 1:permut_n
    erp_perm = erp_amp(perm_indices{perm}); % same shuffle for all channels
    for foi = 1:nFreq
        [corr_perm(perm, foi), pvalue_perm(perm, foi)] = circ_corrcl(phase(:, foi), erp_perm);
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
