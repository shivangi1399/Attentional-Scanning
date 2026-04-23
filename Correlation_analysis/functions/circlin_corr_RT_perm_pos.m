function circlin_corr_RT_perm_pos(cfg_fun)
% H2 permutation: shuffle RT globally across hit trials, then compute circ_corrcl
% within each stimulus position, average correlations across positions.

ichan        = cfg_fun.ichan;
permut_n     = cfg_fun.permut_n;
perm_indices = cfg_fun.perm_indices;
trial_idx    = cfg_fun.trial_idx;   % hit trial indices into RT_trialinfo
positions    = cfg_fun.positions;

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
nPos  = numel(positions);
corr_perm_pos = nan(permut_n, nFreq);

for perm = 1:permut_n
    rt_perm  = rt(perm_indices{perm});
    corr_pos = nan(nPos, nFreq);
    for p = 1:nPos
        mask = pos_labels == positions(p);
        if sum(mask) < 2, continue; end
        ph_p  = phase(mask, :);
        rt_p  = rt_perm(mask);
        valid = ~isnan(rt_p);
        if sum(valid) < 2, continue; end
        ph_p  = ph_p(valid, :);
        rt_p  = rt_p(valid);
        for foi = 1:nFreq
            corr_pos(p,foi) = circ_corrcl(ph_p(:,foi), rt_p);
        end
    end
    corr_perm_pos(perm,:) = mean(corr_pos, 1, 'omitnan');
end

cd(cfg_fun.outfile)
if ~exist(num2str(ichan), 'dir'), mkdir(num2str(ichan)); end
cd(num2str(ichan))
ESIsave corr_perm_pos corr_perm_pos
end
