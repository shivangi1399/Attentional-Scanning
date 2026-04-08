function circlin_corr_hitmiss_perm(cfg_fun)
% circlin_corr_hitmiss_perm  Permutation test for phase vs hit/miss.
%
% Two measures computed per permutation:
%   1) POS (Phase Opposition Sum) = ITC_hits + ITC_misses  -> pos_perm
%   2) ITC with inverted miss phases (+pi)                 -> itc_perm
%
% Shuffles hit/miss labels across trials.

ichan = cfg_fun.ichan;
permut_n = cfg_fun.permut_n;
perm_indices = cfg_fun.perm_indices;  % shared permutations
trial_idx = cfg_fun.trial_idx;        % all trials (hits + misses)
hit_labels = cfg_fun.hit_labels;      % logical: true=hit, false=miss

cd(cfg_fun.infile)
load('ph_all_sess.mat')

phase = ph_comb.phase_all(trial_idx, :, ichan);  % [nTrials x nFreq]

nFreq = size(phase, 2);
pos_perm         = nan(permut_n, nFreq);
itc_perm_complex = complex(nan(permut_n, nFreq));

for perm = 1:permut_n
    % Shuffle hit/miss labels
    labels_perm = hit_labels(perm_indices{perm});
    hit_idx  = labels_perm;
    miss_idx = ~labels_perm;

    % 1) POS: ITC_hits + ITC_misses
    for foi = 1:nFreq
        itc_h = abs(mean(exp(1i * phase(hit_idx, foi))));
        itc_m = abs(mean(exp(1i * phase(miss_idx, foi))));
        pos_perm(perm, foi) = itc_h + itc_m;
    end

    % 2) ITC with inverted miss phases
    phase_combined = phase;
    phase_combined(miss_idx, :) = mod(phase(miss_idx, :) + pi, 2*pi) - pi;

    for foi = 1:nFreq
        itc_perm_complex(perm, foi) = mean(exp(1i * phase_combined(:, foi)));
    end
end
itc_perm = abs(itc_perm_complex);

cd(cfg_fun.outfile)
if ~exist(num2str(ichan), 'dir')
    mkdir(num2str(ichan))
end
cd(num2str(ichan))

ESIsave pos_perm pos_perm
ESIsave itc_perm itc_perm itc_perm_complex
ESIsave perm_indices perm_indices
end
