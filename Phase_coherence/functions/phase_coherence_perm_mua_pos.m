function phase_coherence_perm_mua_pos(cfg_fun)
% H2 permutation: shuffle MUA amplitude globally, then compute coherence
% within each stimulus position, abs() per position, average across positions.

ichan        = cfg_fun.ichan;
permut_n     = cfg_fun.permut_n;
perm_indices = cfg_fun.perm_indices;
trial_idx    = cfg_fun.trial_idx;
positions    = cfg_fun.positions;

cd(cfg_fun.infile)
load('ph_all_sess.mat')

phase      = ph_comb.phase_all(trial_idx, :, ichan);
erp_amp    = ph_comb.MUA_ERP_ampl_all(trial_idx, ichan);
pos_labels = ph_comb.trialinfo(trial_idx, 16);

nFreq = size(phase, 2);
nPos  = numel(positions);
coh_perm_pos = nan(permut_n, nFreq);

for perm = 1:permut_n
    erp_perm = erp_amp(perm_indices{perm});   % shuffle amplitude globally
    coh_pos  = nan(nPos, nFreq);
    for p = 1:nPos
        mask = pos_labels == positions(p);
        if sum(mask) < 2, continue; end
        ph_p  = phase(mask, :);
        amp_p = erp_perm(mask);
        for foi = 1:nFreq
            vec = exp(1i * ph_p(:,foi)) .* amp_p;
            coh_pos(p,foi) = abs(mean(vec));   % abs per position
        end
    end
    coh_perm_pos(perm,:) = mean(coh_pos, 1, 'omitnan');
end

cd(cfg_fun.outfile)
if ~exist(num2str(ichan), 'dir'), mkdir(num2str(ichan)); end
cd(num2str(ichan))
ESIsave coh_perm_pos coh_perm_pos
end
