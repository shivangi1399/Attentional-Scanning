function regress_perm_R_pos(cfg)
% H2 permutation for regression R².
% Within each permutation: shuffle Y within each stimulus position independently,
% run regression per position, average R² across positions.
% This mirrors the real-data H2 computation: within-position regression → average.

X_all_freq_pos = cfg.X_all_freq_pos;   % cell(nPos, numFreq)
Y_all_freq_pos = cfg.Y_all_freq_pos;   % cell(nPos, numFreq)
numFreq    = cfg.numFreq;
nPos       = cfg.nPos;
isLogistic = cfg.isLogistic;
output_dir = cfg.output_dir;
perm_indices = cfg.perm_idx(:)';

for perm_idx = perm_indices

    rng(2025 + perm_idx)

    null_R_phase    = zeros(numFreq, 1);
    null_R_MUA      = zeros(numFreq, 1);
    null_R_Amp      = zeros(numFreq, 1);
    null_R_AmpPhase = zeros(numFreq, 1);
    null_R_any      = zeros(numFreq, 1);

    % Per-position full-model sin/cos betas — for paired R_phase test.
    null_b_sin_pos = nan(nPos, numFreq);
    null_b_cos_pos = nan(nPos, numFreq);

    for f = 1:numFreq
        R_phase_pos    = nan(nPos,1);
        R_MUA_pos      = nan(nPos,1);
        R_Amp_pos      = nan(nPos,1);
        R_AmpPhase_pos = nan(nPos,1);
        R_any_pos      = nan(nPos,1);

        for p = 1:nPos
            X_clean = X_all_freq_pos{p,f};
            Y_clean = Y_all_freq_pos{p,f};
            if isempty(Y_clean) || length(Y_clean) < size(X_clean,2)+1, continue; end

            % Shuffle Y within this position
            Y_perm = Y_clean(randperm(length(Y_clean)));

            if ~isLogistic
                X_full   = [ones(size(X_clean,1),1), X_clean];
                b_full   = regress(Y_perm, X_full);
                RSS_full = sum((Y_perm - X_full*b_full).^2);
                RSS_null = sum((Y_perm - mean(Y_perm)).^2);

                R_any_pos(p) = max(0,(RSS_null - RSS_full)/RSS_null);

                X_red = [ones(size(X_clean,1),1), X_clean(:,1:3)];
                b_red = regress(Y_perm,X_red); RSS_red = sum((Y_perm-X_red*b_red).^2);
                R_phase_pos(p) = max(0,(RSS_red - RSS_full)/RSS_null);

                X_red = [ones(size(X_clean,1),1), X_clean(:,[1 3 4 5])];
                b_red = regress(Y_perm,X_red); RSS_red = sum((Y_perm-X_red*b_red).^2);
                R_MUA_pos(p) = max(0,(RSS_red - RSS_full)/RSS_null);

                X_red = [ones(size(X_clean,1),1), X_clean(:,[1 2 4 5])];
                b_red = regress(Y_perm,X_red); RSS_red = sum((Y_perm-X_red*b_red).^2);
                R_Amp_pos(p) = max(0,(RSS_red - RSS_full)/RSS_null);

                X_red = [ones(size(X_clean,1),1), X_clean(:,1:2)];
                b_red = regress(Y_perm,X_red); RSS_red = sum((Y_perm-X_red*b_red).^2);
                R_AmpPhase_pos(p) = max(0,(RSS_red - RSS_full)/RSS_null);

                null_b_sin_pos(p,f) = b_full(end-1);
                null_b_cos_pos(p,f) = b_full(end);

            else
                b_full = glmfit(X_clean,Y_perm,'binomial','link','logit');
                p_full = glmval(b_full,X_clean,'logit');
                LL_full = sum(Y_perm.*log(p_full+eps)+(1-Y_perm).*log(1-p_full+eps));

                b_null = glmfit(ones(size(Y_perm,1),1),Y_perm,'binomial','link','logit');
                p_null = glmval(b_null,ones(size(Y_perm,1),1),'logit');
                LL_null = sum(Y_perm.*log(p_null+eps)+(1-Y_perm).*log(1-p_null+eps));
                R_any_pos(p) = 1-(LL_full/LL_null);

                b_red = glmfit(X_clean(:,1:3),Y_perm,'binomial','link','logit');
                p_red = glmval(b_red,X_clean(:,1:3),'logit');
                LL_red = sum(Y_perm.*log(p_red+eps)+(1-Y_perm).*log(1-p_red+eps));
                R_phase_pos(p) = 1-(LL_full/LL_red);

                b_red = glmfit(X_clean(:,[1 3 4 5]),Y_perm,'binomial','link','logit');
                p_red = glmval(b_red,X_clean(:,[1 3 4 5]),'logit');
                LL_red = sum(Y_perm.*log(p_red+eps)+(1-Y_perm).*log(1-p_red+eps));
                R_MUA_pos(p) = 1-(LL_full/LL_red);

                b_red = glmfit(X_clean(:,[1 2 4 5]),Y_perm,'binomial','link','logit');
                p_red = glmval(b_red,X_clean(:,[1 2 4 5]),'logit');
                LL_red = sum(Y_perm.*log(p_red+eps)+(1-Y_perm).*log(1-p_red+eps));
                R_Amp_pos(p) = 1-(LL_full/LL_red);

                b_red = glmfit(X_clean(:,1:2),Y_perm,'binomial','link','logit');
                p_red = glmval(b_red,X_clean(:,1:2),'logit');
                LL_red = sum(Y_perm.*log(p_red+eps)+(1-Y_perm).*log(1-p_red+eps));
                R_AmpPhase_pos(p) = 1-(LL_full/LL_red);

                null_b_sin_pos(p,f) = b_full(end-1);
                null_b_cos_pos(p,f) = b_full(end);
            end
        end

        % Average R² across positions for this frequency
        null_R_phase(f)    = mean(R_phase_pos,    'omitnan');
        null_R_MUA(f)      = mean(R_MUA_pos,      'omitnan');
        null_R_Amp(f)      = mean(R_Amp_pos,      'omitnan');
        null_R_AmpPhase(f) = mean(R_AmpPhase_pos, 'omitnan');
        null_R_any(f)      = mean(R_any_pos,       'omitnan');
    end

    results.null_max_phase    = max(null_R_phase);
    results.null_max_MUA      = max(null_R_MUA);
    results.null_max_Amp      = max(null_R_Amp);
    results.null_max_AmpPhase = max(null_R_AmpPhase);
    results.null_max_any      = max(null_R_any);
    results.null_R_phase      = null_R_phase;
    results.null_R_MUA        = null_R_MUA;
    results.null_R_Amp        = null_R_Amp;
    results.null_R_AmpPhase   = null_R_AmpPhase;
    results.null_R_any        = null_R_any;

    % Per-position full-model sin/cos betas (nPos × numFreq) — used to
    % build paired R_phase = |complex β| nulls under H2 aggregation
    % (Way 2: |β_pos| per pos, then mean across positions).
    results.null_b_sin_pos = null_b_sin_pos;
    results.null_b_cos_pos = null_b_cos_pos;

    if ~exist(output_dir, 'dir'), mkdir(output_dir); end
    save(fullfile(output_dir, sprintf('perm_%04d.mat', perm_idx)), 'results');
end
