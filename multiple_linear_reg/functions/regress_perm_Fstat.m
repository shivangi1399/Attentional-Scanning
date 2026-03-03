function regress_perm_Fstat(cfg)
% REGRESS_PERM_FSTAT - Single permutation iteration for nested model F-tests
%
% Input:
%   cfg - structure with fields:
%       .X_all_freq     - cell array of predictor matrices per frequency
%       .Y_all_freq     - cell array of response vectors per frequency
%       .numFreq        - number of frequencies
%       .isLogistic     - logical, true for logistic regression
%       .perm_idx       - permutation index (for saving)
%       .output_dir     - directory to save results

% Extract configuration
X_all_freq = cfg.X_all_freq;
Y_all_freq = cfg.Y_all_freq;
numFreq = cfg.numFreq;
isLogistic = cfg.isLogistic;
output_dir = cfg.output_dir;

if ~isvector(cfg.perm_idx)
    error('cfg.perm_idx must be a vector of permutation indices');
end
perm_indices = cfg.perm_idx(:)';

for perm_idx = perm_indices
    
    % Initialize null distributions for this permutation
    null_F_phase    = zeros(numFreq, 1);
    null_F_MUA      = zeros(numFreq, 1);
    null_F_Amp      = zeros(numFreq, 1);
    null_F_AmpPhase = zeros(numFreq, 1);
    null_F_any      = zeros(numFreq, 1);
    
    for f = 1:numFreq
        X_clean = X_all_freq{f};
        Y_clean = Y_all_freq{f};
        
        if isempty(Y_clean), continue; end
        
        % Permute response variable
        Y_perm = Y_clean(randperm(length(Y_clean)));
        
        if ~isLogistic
            % LINEAR REGRESSION
            X_full = [ones(size(X_clean, 1), 1), X_clean];
            b_full_perm = regress(Y_perm, X_full);
            RSS_full = sum((Y_perm - X_full * b_full_perm).^2);
            df2 = size(X_clean, 1) - size(X_clean, 2) - 1;
            
            % Phase
            X_red = [ones(size(X_clean, 1), 1), X_clean(:, 1:3)];
            b_red_perm = regress(Y_perm, X_red);
            RSS_red = sum((Y_perm - X_red * b_red_perm).^2);
            df1 = 2;
            null_F_phase(f) = ((RSS_red - RSS_full) / df1) / (RSS_full / df2);
            
            % MUA
            X_red = [ones(size(X_clean, 1), 1), X_clean(:, [1 3 4 5])];
            b_red = regress(Y_perm, X_red);
            RSS_red = sum((Y_perm - X_red * b_red).^2);
            df1 = 1;
            null_F_MUA(f) = ((RSS_red - RSS_full) / df1) / (RSS_full / df2);
            
            % Amp
            X_red = [ones(size(X_clean, 1), 1), X_clean(:, 3)];
            b_red = regress(Y_perm, X_red);
            RSS_red = sum((Y_perm - X_red * b_red).^2);
            df1 = 1;
            null_F_Amp(f) = ((RSS_red - RSS_full) / df1) / (RSS_full / df2);
            
            % Amp+Phase
            X_red = [ones(size(X_clean, 1), 1), X_clean(:, 1:2)];
            b_red = regress(Y_perm, X_red);
            RSS_red = sum((Y_perm - X_red * b_red).^2);
            df1 = 3;
            null_F_AmpPhase(f) = ((RSS_red - RSS_full) / df1) / (RSS_full / df2);
            
            % Any predictor
            X_null = ones(size(X_clean, 1), 1);
            b_null = regress(Y_perm, X_null);
            RSS_null = sum((Y_perm - X_null * b_null).^2);
            df1 = size(X_clean, 2);
            null_F_any(f) = ((RSS_null - RSS_full) / df1) / (RSS_full / df2);
            
        else
            % LOGISTIC REGRESSION
            b_full_perm = glmfit(X_clean, Y_perm, 'binomial', 'link', 'logit');
            y_hat_full = glmval(b_full_perm, X_clean, 'logit');
            dev_full = -2 * sum(Y_perm .* log(y_hat_full + eps) + ...
                (1 - Y_perm) .* log(1 - y_hat_full + eps));
            
            % Phase
            b_red = glmfit(X_clean(:, 1:3), Y_perm, 'binomial', 'link', 'logit');
            y_hat_red = glmval(b_red, X_clean(:, 1:3), 'logit');
            dev_red = -2 * sum(Y_perm .* log(y_hat_red + eps) + ...
                (1 - Y_perm) .* log(1 - y_hat_red + eps));
            null_F_phase(f) = dev_red - dev_full;
            
            % MUA
            b_red = glmfit(X_clean(:, [1 3:5]), Y_perm, 'binomial', 'link', 'logit');
            y_hat_red = glmval(b_red, X_clean(:, [1 3:5]), 'logit');
            dev_red = -2 * sum(Y_perm .* log(y_hat_red + eps) + ...
                (1 - Y_perm) .* log(1 - y_hat_red + eps));
            null_F_MUA(f) = dev_red - dev_full;
            
            % Amp
            b_red = glmfit(X_clean(:, 3), Y_perm, 'binomial', 'link', 'logit');
            y_hat_red = glmval(b_red, X_clean(:, 3), 'logit');
            dev_red = -2 * sum(Y_perm .* log(y_hat_red + eps) + ...
                (1 - Y_perm) .* log(1 - y_hat_red + eps));
            null_F_Amp(f) = dev_red - dev_full;
            
            % Amp+Phase
            b_red = glmfit(X_clean(:, 1:2), Y_perm, 'binomial', 'link', 'logit');
            y_hat_red = glmval(b_red, X_clean(:, 1:2), 'logit');
            dev_red = -2 * sum(Y_perm .* log(y_hat_red + eps) + ...
                (1 - Y_perm) .* log(1 - y_hat_red + eps));
            null_F_AmpPhase(f) = dev_red - dev_full;
            
            % Any predictor
            b_null = glmfit(ones(size(Y_perm, 1), 1), Y_perm, 'binomial', 'link', 'logit');
            y_hat_null = glmval(b_null, ones(size(Y_perm)), 'logit');
            dev_null = -2 * sum(Y_perm .* log(y_hat_null + eps) + ...
                (1 - Y_perm) .* log(1 - y_hat_null + eps));
            null_F_any(f) = dev_null - dev_full;
        end
    end
    
    % Max-stat across frequencies for multiple comparison correction
    results.null_max_phase    = max(null_F_phase);
    results.null_max_MUA      = max(null_F_MUA);
    results.null_max_Amp      = max(null_F_Amp);
    results.null_max_AmpPhase = max(null_F_AmpPhase);
    results.null_max_any      = max(null_F_any);
    
    % Store frequency-wise results
    results.null_F_phase    = null_F_phase;
    results.null_F_MUA      = null_F_MUA;
    results.null_F_Amp      = null_F_Amp;
    results.null_F_AmpPhase = null_F_AmpPhase;
    results.null_F_any      = null_F_any;
    
    % Save results
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end
    save(fullfile(cfg.output_dir, sprintf('perm_%04d.mat', perm_idx)), 'results');
    
end

