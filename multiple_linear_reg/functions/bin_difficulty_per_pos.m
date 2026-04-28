function bin_idx = bin_difficulty_per_pos(diff_vals, pos_labels, positions, nDiffBins)
% bin_difficulty_per_pos  Within-position quantile bins of difficulty (dE00).
%
%   bin_idx = bin_difficulty_per_pos(diff_vals, pos_labels, positions, nDiffBins)
%
% INPUTS
%   diff_vals   nTrials x 1, continuous difficulty (e.g. trialinfo(:,18))
%   pos_labels  nTrials x 1, stimulus position per trial (trialinfo(:,16))
%   positions   nPos x 1, list of unique positions to bin against
%   nDiffBins   scalar, number of within-position quantile bins (e.g. 4)
%
% OUTPUT
%   bin_idx     nTrials x 1, integer bin index in 1..nDiffBins, NaN if
%               position is missing or difficulty value is NaN, or if the
%               position has fewer trials than nDiffBins.
%
% Used by H3 (per position x difficulty) analyses to group trials into
% (position, difficulty) cells with roughly equal counts within each
% position. Bin edges are the within-position quantiles of the valid
% difficulty values, with the outer edges extended to +/- Inf so that
% extreme values are still captured.

bin_idx = nan(size(diff_vals));

for p = 1:numel(positions)
    pos_mask = pos_labels == positions(p);
    de       = diff_vals(pos_mask);
    valid    = ~isnan(de);
    if sum(valid) < nDiffBins, continue; end

    edges      = quantile(de(valid), linspace(0, 1, nDiffBins + 1));
    edges(1)   = -Inf;
    edges(end) =  Inf;

    [~,~,bin_local]      = histcounts(de, edges);
    bin_local(bin_local == 0) = NaN;

    pos_idx           = find(pos_mask);
    bin_idx(pos_idx)  = bin_local;
end
end
