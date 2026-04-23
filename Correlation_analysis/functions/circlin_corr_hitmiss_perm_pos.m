function circlin_corr_hitmiss_perm_pos(cfg_fun)
% H2 permutation for phase vs hit/miss.
% Shuffles hit/miss labels globally, then within each stimulus position:
%   POS: ITC_hits + ITC_misses  (abs already inherent in ITC)
%   ITC: abs(mean(exp(i*phase_inverted)))  — abs per position before averaging
% Averages both measures across positions.

ichan        = cfg_fun.ichan;
permut_n     = cfg_fun.permut_n;
perm_indices = cfg_fun.perm_indices;
trial_idx    = cfg_fun.trial_idx;    % all trials (hits + misses)
hit_labels   = cfg_fun.hit_labels;   % logical: true=hit, false=miss
positions    = cfg_fun.positions;

cd(cfg_fun.infile)
load('ph_all_sess.mat')

phase      = ph_comb.phase_all(trial_idx, :, ichan);
pos_labels = ph_comb.trialinfo(trial_idx, 16);

nFreq = size(phase, 2);
nPos  = numel(positions);
pos_perm_pos = nan(permut_n, nFreq);
itc_perm_pos = nan(permut_n, nFreq);

for perm = 1:permut_n
    labels_perm = hit_labels(perm_indices{perm});   % shuffle hit/miss labels globally

    pos_p_mat = nan(nPos, nFreq);
    itc_p_mat = nan(nPos, nFreq);

    for p = 1:nPos
        mask     = pos_labels == positions(p);
        if sum(mask) < 2, continue; end
        ph_p     = phase(mask, :);
        hits_p   = labels_perm(mask);
        misses_p = ~hits_p;
        if sum(hits_p) < 1 || sum(misses_p) < 1, continue; end

        % POS per position
        for foi = 1:nFreq
            itc_h = abs(mean(exp(1i * ph_p(hits_p,   foi))));
            itc_m = abs(mean(exp(1i * ph_p(misses_p, foi))));
            pos_p_mat(p,foi) = itc_h + itc_m;
        end

        % ITC with inverted miss phases — abs per position
        phase_inv_p = ph_p;
        phase_inv_p(misses_p, :) = mod(ph_p(misses_p, :) + pi, 2*pi) - pi;
        for foi = 1:nFreq
            itc_p_mat(p,foi) = abs(mean(exp(1i * phase_inv_p(:,foi))));
        end
    end

    pos_perm_pos(perm,:) = mean(pos_p_mat, 1, 'omitnan');
    itc_perm_pos(perm,:) = mean(itc_p_mat, 1, 'omitnan');
end

cd(cfg_fun.outfile)
if ~exist(num2str(ichan), 'dir'), mkdir(num2str(ichan)); end
cd(num2str(ichan))
ESIsave pos_perm_pos pos_perm_pos
ESIsave itc_perm_pos itc_perm_pos
end
