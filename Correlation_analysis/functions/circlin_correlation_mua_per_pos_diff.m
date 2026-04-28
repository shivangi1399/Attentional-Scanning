function circlin_correlation_mua_per_pos_diff(cfg_fun)
% H3 permutation: shuffle MUA amplitude globally, then compute circ_corrcl
% within each (position x difficulty bin) cell, average correlations across
% cells.

ichan        = cfg_fun.ichan;
permut_n     = cfg_fun.permut_n;
perm_indices = cfg_fun.perm_indices;
cell_pos     = cfg_fun.cell_pos;
cell_dif     = cfg_fun.cell_dif;
diff_bin     = cfg_fun.diff_bin;

cd(cfg_fun.infile)
load('ph_all_sess.mat')

phase      = ph_comb.phase_all(:, :, ichan);
erp_amp    = ph_comb.MUA_ERP_ampl_all(:, ichan);
pos_labels = ph_comb.trialinfo(:, 16);

nFreq = size(phase, 2);
nCell = numel(cell_pos);
corr_perm_pos_diff = nan(permut_n, nFreq);

for perm = 1:permut_n
    erp_perm  = erp_amp(perm_indices{perm});
    corr_cell = nan(nCell, nFreq);
    for c = 1:nCell
        mask = (pos_labels == cell_pos(c)) & (diff_bin == cell_dif(c));
        if sum(mask) < 2, continue; end
        ph_c  = phase(mask, :);
        amp_c = erp_perm(mask);
        for foi = 1:nFreq
            corr_cell(c,foi) = circ_corrcl(ph_c(:,foi), amp_c);
        end
    end
    corr_perm_pos_diff(perm,:) = mean(corr_cell, 1, 'omitnan');
end

cd(cfg_fun.outfile)
if ~exist(num2str(ichan), 'dir'), mkdir(num2str(ichan)); end
cd(num2str(ichan))
ESIsave corr_perm_pos_diff corr_perm_pos_diff
end
