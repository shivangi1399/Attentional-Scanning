function circlin_correlation_mua_pos(cfg_fun)
% H2 permutation: shuffle MUA amplitude globally, then compute circ_corrcl
% within each stimulus position, average correlations across positions.

ichan        = cfg_fun.ichan;
permut_n     = cfg_fun.permut_n;
perm_indices = cfg_fun.perm_indices;
positions    = cfg_fun.positions;

cd(cfg_fun.infile)
load('ph_all_sess.mat')

phase      = ph_comb.phase_all(:, :, ichan);
erp_amp    = ph_comb.MUA_ERP_ampl_all(:, ichan);
pos_labels = ph_comb.trialinfo(:, 16);

nFreq = size(phase, 2);
nPos  = numel(positions);
corr_perm_pos = nan(permut_n, nFreq);

for perm = 1:permut_n
    erp_perm  = erp_amp(perm_indices{perm});
    corr_pos  = nan(nPos, nFreq);
    for p = 1:nPos
        mask = pos_labels == positions(p);
        if sum(mask) < 2, continue; end
        ph_p  = phase(mask, :);
        amp_p = erp_perm(mask);
        for foi = 1:nFreq
            corr_pos(p,foi) = circ_corrcl(ph_p(:,foi), amp_p);
        end
    end
    corr_perm_pos(perm,:) = mean(corr_pos, 1, 'omitnan');
end

cd(cfg_fun.outfile)
if ~exist(num2str(ichan), 'dir'), mkdir(num2str(ichan)); end
cd(num2str(ichan))
ESIsave corr_perm_pos corr_perm_pos
end
