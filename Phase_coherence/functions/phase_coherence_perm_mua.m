function coh_perm = phase_coherence_perm_mua(cfg_fun)
% phase_coherence_permutation computes coherence under shuffled ERP amplitudes

ichan = cfg_fun.ichan;
permut_n = cfg_fun.permut_n;
perm_indices = cfg_fun.perm_indices;  % shared permutations
trial_idx = cfg_fun.trial_idx; %selected trials

cd(cfg_fun.infile)
load('ph_all_sess.mat')

phase = ph_comb.phase_all(trial_idx,:,ichan);    % [trials × freq]
erp_amp = ph_comb.MUA_ERP_ampl_all(trial_idx,ichan); % [trials × 1]

nFreq = size(phase, 2);
coh_perm = nan(permut_n, nFreq);
phase_spec_perm = nan(permut_n, nFreq);

for perm = 1:permut_n
    erp_perm = erp_amp(perm_indices{perm}); 
    
    for foi = 1:nFreq
        vec = exp(1i * phase(:,foi)) .* erp_perm(:);
        cavg = mean(vec);
        coh_perm(perm,foi) = abs(cavg);
        phase_spec_perm(perm,foi) = angle(cavg);
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
