function phase_coherence_RT_perm_per_pos_diff(cfg_fun)
% H3 permutation: shuffle RT globally across hit trials, then compute coherence
% within each (position x difficulty bin) cell, abs() per cell, average
% magnitudes across cells.

ichan        = cfg_fun.ichan;
permut_n     = cfg_fun.permut_n;
perm_indices = cfg_fun.perm_indices;
trial_idx    = cfg_fun.trial_idx;     % hit trial indices into RT_trialinfo
cell_pos     = cfg_fun.cell_pos;
cell_dif     = cfg_fun.cell_dif;
diff_bin     = cfg_fun.diff_bin;       % bin per trial in trial_idx

cd(cfg_fun.infile)
load('ph_all_sess.mat')

phase      = ph_comb.phase_all(trial_idx, :, ichan);
rt         = ph_comb.RT(trial_idx, ichan);
pos_labels = ph_comb.RT_trialinfo(trial_idx, 16);

if isempty(rt)
    fprintf('Channel %d skipped (no valid RT)\n', ichan);
    return
end

nFreq = size(phase, 2);
nCell = numel(cell_pos);
coh_perm_pos_diff = nan(permut_n, nFreq);

for perm = 1:permut_n
    rt_perm  = rt(perm_indices{perm});
    coh_cell = nan(nCell, nFreq);
    for c = 1:nCell
        mask = (pos_labels == cell_pos(c)) & (diff_bin == cell_dif(c));
        if sum(mask) < 2, continue; end
        ph_c  = phase(mask, :);
        rt_c  = rt_perm(mask);
        valid = ~isnan(rt_c);
        if sum(valid) < 2, continue; end
        ph_c  = ph_c(valid, :);
        rt_c  = rt_c(valid);
        for foi = 1:nFreq
            vec = exp(1i * ph_c(:,foi)) .* rt_c;
            coh_cell(c,foi) = abs(mean(vec));
        end
    end
    coh_perm_pos_diff(perm,:) = mean(coh_cell, 1, 'omitnan');
end

cd(cfg_fun.outfile)
if ~exist(num2str(ichan), 'dir'), mkdir(num2str(ichan)); end
cd(num2str(ichan))
ESIsave coh_perm_pos_diff coh_perm_pos_diff
end
