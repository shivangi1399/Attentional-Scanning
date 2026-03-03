function regress_perm_R(cfg)
% REGRESS_PERM_R - Single permutation iteration for nested model R²/pseudo-R²
%
% Input:
%   cfg - structure with fields:
%       .X_all_freq     - cell array of predictor matrices per frequency
%       .Y_all_freq     - cell array of response vectors per frequency
%       .numFreq        - number of frequencies
%       .isLogistic     - logical, true for logistic regression
%       .perm_idx       - vector of permutation indices (for saving)
%       .output_dir     - directory to save results

% Extract configuration
X_all_freq = cfg.X_all_freq;
Y_all_freq = cfg.Y_all_freq;
numFreq    = cfg.numFreq;
isLogistic = cfg.isLogistic;
output_dir = cfg.output_dir;

if ~isvector(cfg.perm_idx)
    error('cfg.perm_idx must be a vector of permutation indices');
end
perm_indices = cfg.perm_idx(:)';

for perm_idx = perm_indices

    % Initialize null distributions for this permutation
    null_R_phase    = zeros(numFreq, 1);
    null_R_MUA      = zeros(numFreq, 1);
    null_R_Amp      = zeros(numFreq, 1);
    null_R_AmpPhase = zeros(numFreq, 1);
    null_R_any      = zeros(numFreq, 1);

    for f = 1:numFreq
        X_clean = X_all_freq{f};
        Y_clean = Y_all_freq{f};

        if isempty(Y_clean), continue; end

        % Permute response variable
        Y_perm = Y_clean(randperm(length(Y_clean)));

        if ~isLogistic
            % LINEAR REGRESSION
            X_full = [ones(size(X_clean,1),1), X_clean];
            b_full = regress(Y_perm, X_full);
            RSS_full = sum((Y_perm - X_full*b_full).^2);

            % Null model (intercept only) for normalisation
            RSS_null = sum((Y_perm - mean(Y_perm)).^2);

            % Any predictor
            null_R_any(f) = max(0, (RSS_null - RSS_full) / RSS_null);

            % Phase (drops sin/cos, cols 4-5)
            X_red = [ones(size(X_clean,1),1), X_clean(:,1:3)];
            b_red = regress(Y_perm, X_red);
            RSS_red = sum((Y_perm - X_red*b_red).^2);
            null_R_phase(f) = max(0, (RSS_red - RSS_full) / RSS_null);

            % MUA (drops col 2)
            X_red = [ones(size(X_clean,1),1), X_clean(:,[1 3 4 5])];
            b_red = regress(Y_perm, X_red);
            RSS_red = sum((Y_perm - X_red*b_red).^2);
            null_R_MUA(f) = max(0, (RSS_red - RSS_full) / RSS_null);

            % Amp (drops col 3)
            X_red = [ones(size(X_clean,1),1), X_clean(:,[1 2 4 5])];
            b_red = regress(Y_perm, X_red);
            RSS_red = sum((Y_perm - X_red*b_red).^2);
            null_R_Amp(f) = max(0, (RSS_red - RSS_full) / RSS_null);

            % Amp+Phase (drops cols 3-5)
            X_red = [ones(size(X_clean,1),1), X_clean(:,1:2)];
            b_red = regress(Y_perm, X_red);
            RSS_red = sum((Y_perm - X_red*b_red).^2);
            null_R_AmpPhase(f) = max(0, (RSS_red - RSS_full) / RSS_null);

        else
            % LOGISTIC REGRESSION
            b_full = glmfit(X_clean, Y_perm, 'binomial', 'link', 'logit');
            p_full = glmval(b_full, X_clean, 'logit');
            LL_full = sum(Y_perm.*log(p_full+eps) + (1-Y_perm).*log(1-p_full+eps));

            % Null model
            b_null = glmfit(ones(size(Y_perm,1),1), Y_perm, 'binomial', 'link', 'logit');
            p_null = glmval(b_null, ones(size(Y_perm,1),1), 'logit');
            LL_null = sum(Y_perm.*log(p_null+eps) + (1-Y_perm).*log(1-p_null+eps));

            % Any predictor
            null_R_any(f) = 1 - (LL_full / LL_null);

            % Phase
            b_red = glmfit(X_clean(:,1:3), Y_perm, 'binomial', 'link', 'logit');
            p_red = glmval(b_red, X_clean(:,1:3), 'logit');
            LL_red = sum(Y_perm.*log(p_red+eps) + (1-Y_perm).*log(1-p_red+eps));
            null_R_phase(f) = 1 - (LL_full / LL_red);

            % MUA
            b_red = glmfit(X_clean(:,[1 3 4 5]), Y_perm, 'binomial', 'link', 'logit');
            p_red = glmval(b_red, X_clean(:,[1 3 4 5]), 'logit');
            LL_red = sum(Y_perm.*log(p_red+eps) + (1-Y_perm).*log(1-p_red+eps));
            null_R_MUA(f) = 1 - (LL_full / LL_red);

            % Amp
            b_red = glmfit(X_clean(:,[1 2 4 5]), Y_perm, 'binomial', 'link', 'logit');
            p_red = glmval(b_red, X_clean(:,[1 2 4 5]), 'logit');
            LL_red = sum(Y_perm.*log(p_red+eps) + (1-Y_perm).*log(1-p_red+eps));
            null_R_Amp(f) = 1 - (LL_full / LL_red);

            % Amp+Phase
            b_red = glmfit(X_clean(:,1:2), Y_perm, 'binomial', 'link', 'logit');
            p_red = glmval(b_red, X_clean(:,1:2), 'logit');
            LL_red = sum(Y_perm.*log(p_red+eps) + (1-Y_perm).*log(1-p_red+eps));
            null_R_AmpPhase(f) = 1 - (LL_full / LL_red);
        end
    end

    % Max-stat across frequencies for multiple comparison correction
    results.null_max_phase    = max(null_R_phase);
    results.null_max_MUA      = max(null_R_MUA);
    results.null_max_Amp      = max(null_R_Amp);
    results.null_max_AmpPhase = max(null_R_AmpPhase);
    results.null_max_any      = max(null_R_any);

    % Store frequency-wise results
    results.null_R_phase    = null_R_phase;
    results.null_R_MUA      = null_R_MUA;
    results.null_R_Amp      = null_R_Amp;
    results.null_R_AmpPhase = null_R_AmpPhase;
    results.null_R_any      = null_R_any;

    % Save results
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    save(fullfile(output_dir, sprintf('perm_%04d.mat', perm_idx)), 'results');

end