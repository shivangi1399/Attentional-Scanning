function [h, crit_p, adj_p, sorted_p] = fdr_bh(pvals, q, method, report)
% Benjamini & Hochberg (1995) FDR correction for multiple comparisons
%
% INPUTS:
% pvals  - vector of p-values
% q      - desired false discovery rate (default 0.05)
% method - 'pdep' (independent/positive dependence) or 'dep' (general dependence), default 'pdep'
% report - 1 to print summary, 0 (default) no print
%
% OUTPUTS:
% h      - hypothesis test results (1=reject null, 0=accept)
% crit_p - critical p-value threshold
% adj_p  - adjusted p-values
% sorted_p - sorted input p-values
%
% Usage: [h, crit_p, adj_p, sorted_p] = fdr_bh(pvals, 0.05, 'pdep', 1);

if nargin<2 || isempty(q)
    q = 0.05;
end
if nargin<3 || isempty(method)
    method = 'pdep';
end
if nargin<4
    report = 0;
end

pvals = pvals(:);
m = length(pvals);
[sorted_p, sort_ids] = sort(pvals);
[~, unsort_ids] = sort(sort_ids);

if strcmpi(method,'pdep')
    % Benjamini & Hochberg (1995)
    thresh = (1:m)'/m * q;
elseif strcmpi(method,'dep')
    % Benjamini & Yekutieli (2001)
    denom = sum(1./(1:m));
    thresh = (1:m)'/m * q / denom;
else
    error('Method must be ''pdep'' or ''dep''');
end

w = find(sorted_p <= thresh, 1, 'last'); % largest p that satisfies p <= thresh
if isempty(w)
    crit_p = 0;
    h = zeros(m,1);
else
    crit_p = sorted_p(w);
    h = pvals <= crit_p;
end

% Adjusted p-values
adj_p = zeros(m,1);
for i = m:-1:1
    if i == m
        adj_p(i) = sorted_p(i);
    else
        adj_p(i) = min(sorted_p(i)*m/i, adj_p(i+1));
    end
end
adj_p = adj_p(unsort_ids);

if report
    fprintf('FDR correction results:\n');
    fprintf('Number of tests: %d\n', m);
    fprintf('Number of rejected hypotheses: %d\n', sum(h));
    fprintf('Critical p-value: %.4g\n', crit_p);
end

end
