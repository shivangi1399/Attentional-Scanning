function circlin_corr_hitmiss_perm_per_pos_diff(cfg_fun)
% H3 permutation for phase vs hit/miss.
% Shuffles hit/miss labels globally, then within each (position x difficulty
% bin) cell:
%   POS: ITC_hits + ITC_misses  (abs already inherent in ITC)
%   ITC: abs(mean(exp(i*phase_inverted)))  — abs per cell before averaging
% Averages both measures across cells.

ichan        = cfg_fun.ichan;
permut_n     = cfg_fun.permut_n;
perm_indices = cfg_fun.perm_indices;
trial_idx    = cfg_fun.trial_idx;     % all trials (hits + misses)
hit_labels   = cfg_fun.hit_labels;     % logical: true=hit, false=miss
cell_pos     = cfg_fun.cell_pos;
cell_dif     = cfg_fun.cell_dif;
diff_bin     = cfg_fun.diff_bin;       % bin per trial in trial_idx

cd(cfg_fun.infile)
load('ph_all_sess.mat')

phase      = ph_comb.phase_all(trial_idx, :, ichan);
pos_labels = ph_comb.trialinfo(trial_idx, 16);

nFreq = size(phase, 2);
nCell = numel(cell_pos);
pos_perm_pos_diff = nan(permut_n, nFreq);
itc_perm_pos_diff = nan(permut_n, nFreq);

for perm = 1:permut_n
    labels_perm = hit_labels(perm_indices{perm});

    pos_c_mat = nan(nCell, nFreq);
    itc_c_mat = nan(nCell, nFreq);

    for c = 1:nCell
        mask     = (pos_labels == cell_pos(c)) & (diff_bin == cell_dif(c));
        if sum(mask) < 2, continue; end
        ph_c     = phase(mask, :);
        hits_c   = labels_perm(mask);
        misses_c = ~hits_c;
        if sum(hits_c) < 1 || sum(misses_c) < 1, continue; end

        for foi = 1:nFreq
            itc_h = abs(mean(exp(1i * ph_c(hits_c,   foi))));
            itc_m = abs(mean(exp(1i * ph_c(misses_c, foi))));
            pos_c_mat(c,foi) = itc_h + itc_m;
        end

        phase_inv_c = ph_c;
        phase_inv_c(misses_c, :) = mod(ph_c(misses_c, :) + pi, 2*pi) - pi;
        for foi = 1:nFreq
            itc_c_mat(c,foi) = abs(mean(exp(1i * phase_inv_c(:,foi))));
        end
    end

    pos_perm_pos_diff(perm,:) = mean(pos_c_mat, 1, 'omitnan');
    itc_perm_pos_diff(perm,:) = mean(itc_c_mat, 1, 'omitnan');
end

cd(cfg_fun.outfile)
if ~exist(num2str(ichan), 'dir'), mkdir(num2str(ichan)); end
cd(num2str(ichan))
ESIsave pos_perm_pos_diff pos_perm_pos_diff
ESIsave itc_perm_pos_diff itc_perm_pos_diff
end
