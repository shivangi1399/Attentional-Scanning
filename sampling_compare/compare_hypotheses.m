% =====================================================================
% compare_hypotheses.m  —  H1 / H2 / H3 / H4 comparison + paired tests
%
% Overlays hypothesis magnitude curves (monkey-average) and tests
% statistically where relaxing each level (position, difficulty, channel)
% finds more effect than chance allows.
%
% STATISTICAL FRAMEWORK
% ---------------------
% By Jensen's inequality (mean(|x_i|) >= |mean(x_i)|), a "relaxed"
% hypothesis (H2/H3/H4 — abs taken earlier) is structurally larger than
% H1 even when phase and DV are independent. So comparing raw H_n to
% raw H_{n-1} is meaningless on its own; some Jensen advantage is
% baked in by the recipe.
%
% Solution: PAIRED permutation test on the difference H_n - H_{n-1}.
% All analysis scripts now seed with rng(2025) before generating perm
% indices, so for each permutation index i:
%     H1_perm(i), H2_perm(i), H3_perm(i), H4_perm(i)
% are computed from the SAME shuffled DV. The null distribution of
% (H_n,perm - H_{n-1},perm) therefore captures exactly the Jensen
% advantage that relaxing the level provides under the null. The
% observed difference is significant only if it exceeds that advantage.
%
% Test logic per pair (one-sided, max-stat across frequencies):
%     null_diff(perm, f) = H_n,perm(perm, f) - H_{n-1},perm(perm, f)
%     obs_diff(f)        = H_n,obs(f)        - H_{n-1},obs(f)
%     tmax(perm)         = max_f null_diff(perm, f)
%     threshold          = quantile(tmax, 0.95)
%     significant freqs  = freqs where obs_diff(f) >= threshold
%
% Interpretation of a significant pair:
%     H2 - H1   ->  POSITIONS disagree on preferred phase
%     H3 - H2   ->  DIFFICULTY levels disagree (within position)
%     H4 - H1   ->  CHANNELS disagree on preferred phase
%
% Figures produced (monkey-average; hypotheses are claims about the
% pooled dataset, not individual animals):
%   1. compare_hypotheses_monkey_avg.pdf
%      H1-H4 overlay with per-hypothesis significance shading + ticks
%   2. compare_hypotheses_differences.pdf
%      Pairwise differences with PAIRED null mean (Jensen reference)
%      and 95% max-stat threshold; significant freqs marked
%   3. compare_hypotheses_ratios.pdf
%      Pairwise ratios H_n / H_{n-1} (visualization only, not tested)
%   4. compare_hypotheses_sig_pattern.pdf
%      Significance heatmap (frequency x hypothesis), per-hypothesis test
%   5. compare_hypotheses_nsig.pdf
%      Bar: n significant freqs per hypothesis (per-hypothesis test)
%   6. compare_hypotheses_paired_nsig.pdf
%      Bar: n significant freqs per pair (paired difference test)
%
% Pipelines compared (4 columns per figure):
%   1. Coherence              — |mean(complex coherence)| / mean(|.|)
%   2. Correlation            — circ_corrcl-style magnitude
%   3. Regression R²          — partial R² for sin+cos predictors
%                                = (RSS_reduced - RSS_full) / RSS_null
%                                "How much variance does phase explain
%                                after controlling for MUA and Amp?"
%                                Scale: [0, 1]. Threshold from perm null
%                                of the same partial-R² statistic.
%   4. Regression R_phase     — magnitude of full-model complex β
%                                = sqrt(b_sin^2 + b_cos^2)
%                                "How strong is the phase tuning, and
%                                where is the preferred phase?" This is
%                                the direct analogue of the coherence /
%                                correlation magnitude — same H1-H4
%                                Jensen logic applies (|mean(complex β)|
%                                vs mean(|complex β|) across strata).
%                                Scale: regression-coefficient units.
%                                NO threshold currently saved; observed
%                                curves shown without significance test.
%
% Notes:
% - Coherence, correlation, regression R², regression R_phase: paired
%   test fully wired. Each pipeline persists [nPerm × nFreq] null curves
%   under monkey_avg_results.mat -> perm_monkey_avg, which the loader
%   below populates into perm{p,h,d}.
% - Regression: requires the updated regress_perm_R*.m (saves full-model
%   sin/cos betas) AND the updated stats scripts (aggregate betas into
%   per-perm R_phase null curves) to have been re-run. Old result files
%   missing perm_monkey_avg fall back to "no paired test" gracefully.
% - Correlation H4-H1 is identically zero by construction: circ_corrcl
%   returns a non-negative magnitude per channel, so Way-1 and Way-2
%   across channels collapse to the same arithmetic mean. Both observed
%   and null differences are flat zero — the H4-H1 panel is structural,
%   not a result.
% =====================================================================
clearvars; close all; clc

%% Settings
addpath /opt/fieldtrip_github/; ft_defaults

base_comb = '/mnt/hpc/projects/MWSampling/4Shivangi/results_combined';
save_dir  = '/mnt/hpc/projects/MWSampling/4Shivangi/Plots/sampling_compare/hypotheses/monkey_avg';
if ~exist(save_dir,'dir'), mkdir(save_dir); end

hyp_keys   = {'complex','abs_per_pos','abs_per_pos_diff','abs_per_chan'};
hyp_labels = {'H1 (pooled)','H2 (per pos)','H3 (per pos×diff)','H4 (per chan)'};
hyp_colors = [0.20 0.40 0.80;
              0.20 0.65 0.40;
              0.80 0.30 0.30;
              0.70 0.40 0.80];
nH = numel(hyp_keys);

dv_labels   = {'MUA','LFP','RT','Hit/Miss'};
reg_dv_keys = {'MUA_ERP_ampl_all','LFP_ERP_ampl_all','RT','hit_miss'};
% Regression appears twice: column 3 = partial R² (variance share),
% column 4 = R_phase = |complex β| (tuning strength, directly analogous
% to coherence/correlation magnitudes).
pipe_labels = {'Coherence','Correlation','Regression R²','Regression R_phase'};
nDV = 4;  nP = 4;

% Pairwise comparisons — one per axis being tested
cmp_hi  = [2 3 4];  cmp_lo  = [1 2 1];
cmp_lab = {'H2−H1 (position)','H3−H2 (difficulty)','H4−H1 (channel)'};
cmp_col = [0.20 0.65 0.40; 0.80 0.30 0.30; 0.70 0.40 0.80];
nC = 3;

%% Data loading (monkey-average only)
% obs{p,h,d}   — observed curve [1 x nFreq]
% thr(p,h,d)   — per-hypothesis permutation threshold (scalar)
% perm{p,h,d}  — per-perm null [nPerm x nFreq] (for paired test; empty if missing)
obs  = cell(nP,nH,nDV);
thr  = nan(nP,nH,nDV);
perm = cell(nP,nH,nDV);
freq_axis = [];

% ── Coherence (pipe 1) ─
coh_sub  = {'mua','lfp','RT','hit_miss_itc'};
coh_file = {'monkey_avg_results.mat','monkey_avg_results.mat', ...
            'monkey_avg_results.mat','monkey_avg_results_itc.mat'};
coh_vf   = {'coh_monkey_avg','coh_monkey_avg','coh_monkey_avg','itc_monkey_avg'};
coh_tf   = {'thresh_monkey_avg','thresh_monkey_avg','thresh_monkey_avg','thresh_monkey_avg_itc'};
coh_pf   = {'perm_monkey_avg','perm_monkey_avg','perm_monkey_avg','perm_monkey_avg_itc'};

for h = 1:nH
    root = fullfile(base_comb, 'phase_coherence', hyp_keys{h}, 'cp10_till_100');
    for d = 1:nDV
        f = fullfile(root, coh_sub{d}, 'all_loc_difflev', coh_file{d});
        [obs{1,h,d}, thr(1,h,d), perm{1,h,d}, freq_axis] = ...
            load_field(f, coh_vf{d}, coh_tf{d}, coh_pf{d}, freq_axis);
    end
end

% ── Correlation (pipe 2) ─
corr_sub  = {'mua','lfp','RT','hit_miss'};
corr_file = {'monkey_avg_results.mat','monkey_avg_results.mat', ...
             'monkey_avg_results.mat','monkey_avg_results_pos.mat'};
corr_vf   = {'corr_monkey_avg','corr_monkey_avg','corr_monkey_avg','pos_monkey_avg'};
corr_tf   = {'thresh_monkey_avg','thresh_monkey_avg','thresh_monkey_avg','thresh_monkey_avg_pos'};
corr_pf   = {'perm_monkey_avg','perm_monkey_avg','perm_monkey_avg','perm_monkey_avg_pos'};

for h = 1:nH
    root = fullfile(base_comb, 'phase_correlation', hyp_keys{h}, 'cp10_till_100');
    for d = 1:nDV
        f = fullfile(root, corr_sub{d}, 'all_loc_difflev', corr_file{d});
        [obs{2,h,d}, thr(2,h,d), perm{2,h,d}, freq_axis] = ...
            load_field(f, corr_vf{d}, corr_tf{d}, corr_pf{d}, freq_axis);
    end
end

% ── Regression R² (pipe 3) ─ partial R² for phase: variance share that
%   sin+cos uniquely explain after MUA/Amp are accounted for.
%   Loaded from monkey_avg_obs.phase; threshold from same partial-R² null;
%   per-perm null curve from perm_monkey_avg.phase (paired test).
for h = 1:nH
    root = fullfile(base_comb, 'multi_lin_reg', hyp_keys{h}, 'cp10_till_100');
    for d = 1:nDV
        rf = fullfile(root, reg_dv_keys{d}, 'monkey_avg_results.mat');
        if ~isfile(rf), continue; end
        s = load(rf);
        if isfield(s,'monkey_avg_obs') && isfield(s.monkey_avg_obs,'phase')
            obs{3,h,d} = s.monkey_avg_obs.phase;
        end
        if isfield(s,'thresh_monkey') && isfield(s.thresh_monkey,'phase')
            thr(3,h,d) = s.thresh_monkey.phase;
        end
        if isfield(s,'perm_monkey_avg') && isfield(s.perm_monkey_avg,'phase')
            perm{3,h,d} = s.perm_monkey_avg.phase;
        end
    end
end

% ── Regression R_phase (pipe 4) ─ |complex β| tuning strength: magnitude
%   of full-model sin/cos coefficients, aggregated via complex mean
%   (Way 1 / H1) or mean of magnitudes (Way 2 / H4). Paired-null curve
%   from perm_monkey_avg.R_phase (built by stats scripts from per-perm
%   full-model betas saved by regress_perm_R*.m).
for h = 1:nH
    root = fullfile(base_comb, 'multi_lin_reg', hyp_keys{h}, 'cp10_till_100');
    for d = 1:nDV
        rf = fullfile(root, reg_dv_keys{d}, 'monkey_avg_results.mat');
        if ~isfile(rf), continue; end
        s = load(rf);
        if isfield(s,'monkey_avg_obs') && isfield(s.monkey_avg_obs,'R_phase')
            obs{4,h,d} = s.monkey_avg_obs.R_phase;
        end
        if isfield(s,'perm_monkey_avg') && isfield(s.perm_monkey_avg,'R_phase')
            perm{4,h,d} = s.perm_monkey_avg.R_phase;
            % No per-hypothesis threshold for R_phase (different scale than
            % partial R²); paired test handles significance for differences.
        end
    end
end

nFreq = numel(freq_axis);
if nFreq == 0
    error('No frequency axis found — check that at least one result file exists.');
end

%% Paired permutation test on differences
% diff_obs{p,c,d}        — observed difference curve [1 x nFreq]
% diff_thr(p,c,d)        — paired-null max-stat 95th percentile (scalar)
% diff_sig{p,c,d}        — logical mask [1 x nFreq] of significant freqs
% diff_null_mean{p,c,d}  — per-freq null mean (visualised Jensen advantage)
diff_obs       = cell(nP,nC,nDV);
diff_thr       = nan(nP,nC,nDV);
diff_sig       = cell(nP,nC,nDV);
diff_null_mean = cell(nP,nC,nDV);

for p = 1:nP
    for c = 1:nC
        for d = 1:nDV
            v_hi = obs{p,cmp_hi(c),d};
            v_lo = obs{p,cmp_lo(c),d};
            P_hi = perm{p,cmp_hi(c),d};
            P_lo = perm{p,cmp_lo(c),d};
            if isempty(v_hi)||isempty(v_lo), continue; end

            % Observed difference always available
            diff_obs{p,c,d} = v_hi - v_lo;

            % Paired test only if both per-perm null arrays are present
            if isempty(P_hi)||isempty(P_lo) || ~isequal(size(P_hi), size(P_lo))
                continue
            end

            null_d = P_hi - P_lo;                  % [nPerm x nFreq]
            tmax   = max(null_d, [], 2);            % one-sided max-stat
            t_d    = quantile(tmax, 0.95);

            diff_thr(p,c,d)       = t_d;
            diff_sig{p,c,d}       = diff_obs{p,c,d}(:)' >= t_d;
            diff_null_mean{p,c,d} = mean(null_d, 1, 'omitnan');
        end
    end
end

%% Figure 1: Monkey-average overlay with per-hypothesis significance shading
f1 = figure('Name','H1-H4 overlay — monkey avg', ...
    'Units','centimeters','Position',[1 1 70 42]);
set(f1,'PaperUnits','centimeters','PaperSize',[70 42],'PaperPosition',[0 0 70 42]);

for d = 1:nDV
    for p = 1:nP
        subplot(nDV,nP,(d-1)*nP+p); hold on;
        for h = 1:nH
            v = obs{p,h,d};
            if isempty(v), continue; end
            plot(freq_axis, v, 'Color',hyp_colors(h,:), 'LineWidth',2, ...
                'DisplayName',hyp_labels{h});
            t = thr(p,h,d);
            if ~isnan(t)
                yline(t,'--','Color',hyp_colors(h,:),'LineWidth',0.8, ...
                    'Alpha',0.7,'HandleVisibility','off');
            end
        end
        for h = 1:nH
            do_shade(freq_axis, obs{p,h,d}, thr(p,h,d), hyp_colors(h,:));
        end
        % Significance tick marks near top (one row per hypothesis)
        yl = ylim; span = yl(2)-yl(1);
        for h = 1:nH
            v = obs{p,h,d}; t = thr(p,h,d);
            if isempty(v)||isnan(t), continue; end
            sig_f = freq_axis(v >= t);
            if ~isempty(sig_f)
                scatter(sig_f, repmat(yl(2)-(h-1)*0.04*span, 1,numel(sig_f)), ...
                    8, hyp_colors(h,:), 'filled','s','HandleVisibility','off');
            end
        end
        if d==1, title(sprintf('%s — %s',pipe_labels{p},dv_labels{d}),'FontSize',8);
        else,    title(dv_labels{d},'FontSize',8); end
        if p==1 && d==1, legend('Location','best','FontSize',6); end
        xlabel('Freq (Hz)','FontSize',7);
        set(gca,'FontSize',7,'Box','on');
    end
end
sgtitle('H1–H4 overlay — monkey average  (shading + ticks = per-hypothesis significant)', ...
    'FontSize',12,'FontWeight','bold');
print(f1, fullfile(save_dir,'compare_hypotheses_monkey_avg.pdf'),'-dpdf');
fprintf('Saved: %s\n', fullfile(save_dir,'compare_hypotheses_monkey_avg.pdf'));

%% Figure 2: Pairwise differences with PAIRED null and significance
% Solid coloured = observed difference H_n - H_{n-1}
% Dashed coloured = mean of paired null (Jensen advantage under H0)
% Dotted coloured = 95% max-stat threshold (one-sided)
% Dots near top = frequencies where observed exceeds threshold

f2 = figure('Name','Pairwise differences (paired test)', ...
    'Units','centimeters','Position',[1 1 56 38]);
set(f2,'PaperUnits','centimeters','PaperSize',[56 38],'PaperPosition',[0 0 56 38]);

for d = 1:nDV
    for p = 1:nP
        ax2 = subplot(nDV,nP,(d-1)*nP+p); hold on;
        yline(0,'k-','LineWidth',0.5,'HandleVisibility','off');
        any_line = false;
        for c = 1:nC
            od  = diff_obs{p,c,d};
            t_d = diff_thr(p,c,d);
            nm  = diff_null_mean{p,c,d};
            if isempty(od), continue; end

            % Observed difference
            plot(freq_axis, od, 'Color',cmp_col(c,:),'LineWidth',2, ...
                'DisplayName',cmp_lab{c});
            any_line = true;

            % Null mean (Jensen advantage reference)
            if ~isempty(nm)
                plot(freq_axis, nm, '--', 'Color',cmp_col(c,:),'LineWidth',0.8, ...
                    'HandleVisibility','off');
            end

            % Threshold line
            if ~isnan(t_d)
                yline(t_d,':','Color',cmp_col(c,:),'LineWidth',1, ...
                    'HandleVisibility','off');
            end
        end

        % Significance dots near top (after axis stabilises)
        yl = ylim; span = yl(2)-yl(1);
        for c = 1:nC
            sig = diff_sig{p,c,d};
            if isempty(sig)||~any(sig), continue; end
            y_marker = yl(2) - (c-1)*0.05*span;
            scatter(freq_axis(sig), repmat(y_marker, 1,sum(sig)), ...
                10, cmp_col(c,:),'filled','HandleVisibility','off');
        end
        ylim(yl);

        if d==1, title(sprintf('%s — %s',pipe_labels{p},dv_labels{d}),'FontSize',8);
        else,    title(dv_labels{d},'FontSize',8); end
        xlabel('Freq (Hz)','FontSize',7); ylabel('\Delta','FontSize',7);
        if p==1 && d==1 && any_line, legend('Location','best','FontSize',6); end
        set(ax2,'FontSize',7,'Box','on');
    end
end
sgtitle({'Paired differences (monkey avg)', ...
    'solid = observed | dashed = null mean (Jensen advantage) | dotted = 95% max-stat threshold | dots = significant'}, ...
    'FontSize',10,'FontWeight','bold');
print(f2, fullfile(save_dir,'compare_hypotheses_differences.pdf'),'-dpdf');
fprintf('Saved: %s\n', fullfile(save_dir,'compare_hypotheses_differences.pdf'));

%% Figure 3: Pairwise ratio curves (visualization only, not tested)
% H_n / H_{n-1}: >1 = relaxing finds more effect. Y clipped at [0 3].
f3 = figure('Name','Pairwise ratios — monkey avg', ...
    'Units','centimeters','Position',[1 1 54 36]);
set(f3,'PaperUnits','centimeters','PaperSize',[54 36],'PaperPosition',[0 0 54 36]);

for d = 1:nDV
    for p = 1:nP
        ax3 = subplot(nDV,nP,(d-1)*nP+p); hold on;
        yline(1,'k--','LineWidth',1,'HandleVisibility','off');
        any_line = false;
        for c = 1:nC
            vhi = obs{p,cmp_hi(c),d};
            vlo = obs{p,cmp_lo(c),d};
            if isempty(vhi)||isempty(vlo), continue; end
            plot(freq_axis, vhi./max(vlo,1e-10), 'Color',cmp_col(c,:),'LineWidth',2, ...
                'DisplayName',cmp_lab{c});
            any_line = true;
        end
        ylim([0 3]);
        if d==1, title(sprintf('%s — %s',pipe_labels{p},dv_labels{d}),'FontSize',8);
        else,    title(dv_labels{d},'FontSize',8); end
        xlabel('Freq (Hz)','FontSize',7); ylabel('Ratio','FontSize',7);
        if p==1 && d==1 && any_line, legend('Location','best','FontSize',6); end
        set(ax3,'FontSize',7,'Box','on');
    end
end
sgtitle({'Pairwise ratios (monkey avg) — visualization only, no formal test', ...
    '>1 = relaxing finds more effect   (y clipped at 3)'}, ...
    'FontSize',11,'FontWeight','bold');
print(f3, fullfile(save_dir,'compare_hypotheses_ratios.pdf'),'-dpdf');
fprintf('Saved: %s\n', fullfile(save_dir,'compare_hypotheses_ratios.pdf'));

%% Figure 4: Significance pattern heatmap (per-hypothesis test)
f4 = figure('Name','Significance pattern — monkey avg', ...
    'Units','centimeters','Position',[1 1 70 42]);
set(f4,'PaperUnits','centimeters','PaperSize',[70 42],'PaperPosition',[0 0 70 42]);

for d = 1:nDV
    for p = 1:nP
        ax4 = subplot(nDV,nP,(d-1)*nP+p);
        img = ones(nFreq, nH, 3);   % white background
        for h = 1:nH
            v = obs{p,h,d}; t = thr(p,h,d);
            if isempty(v)||isnan(t)||numel(v)~=nFreq, continue; end
            mask = v(:) >= t;
            for ch = 1:3
                col_ch = img(:,h,ch);
                col_ch(mask) = hyp_colors(h,ch);
                img(:,h,ch) = col_ch;
            end
        end
        imagesc(1:nH, freq_axis, img);
        set(ax4,'XTick',1:nH,'XTickLabel',{'H1','H2','H3','H4'}, ...
            'XTickLabelRotation',30,'YDir','normal','FontSize',7);
        ylabel('Freq (Hz)','FontSize',7);
        if d==1, title(sprintf('%s — %s',pipe_labels{p},dv_labels{d}),'FontSize',8);
        else,    title(dv_labels{d},'FontSize',8); end
    end
end
sgtitle({'Significance pattern (monkey avg) — per-hypothesis test', ...
    'H1 blue | H2 green | H3 red | H4 purple   (white = not significant)'}, ...
    'FontSize',11,'FontWeight','bold');
print(f4, fullfile(save_dir,'compare_hypotheses_sig_pattern.pdf'),'-dpdf');
fprintf('Saved: %s\n', fullfile(save_dir,'compare_hypotheses_sig_pattern.pdf'));

%% Figure 5: N significant freqs per hypothesis (per-hypothesis test)
n_sig = zeros(nP,nH,nDV);
for p = 1:nP; for h = 1:nH; for d = 1:nDV
    v = obs{p,h,d}; t = thr(p,h,d);
    if ~isempty(v)&&~isnan(t), n_sig(p,h,d) = sum(v>=t); end
end; end; end

f5 = figure('Name','N sig frequencies per hypothesis', ...
    'Units','centimeters','Position',[1 1 54 28]);
set(f5,'PaperUnits','centimeters','PaperSize',[54 28],'PaperPosition',[0 0 54 28]);

for d = 1:nDV
    for p = 1:nP
        ax5 = subplot(nDV,nP,(d-1)*nP+p);
        vals = squeeze(n_sig(p,:,d));
        b = bar(1:nH, vals, 'FaceColor','flat');
        for h = 1:nH, b.CData(h,:) = hyp_colors(h,:); end
        set(ax5,'XTick',1:nH,'XTickLabel',{'H1','H2','H3','H4'}, ...
            'XTickLabelRotation',30,'FontSize',7);
        ylabel('n sig freqs','FontSize',7);
        if d==1, title(sprintf('%s — %s',pipe_labels{p},dv_labels{d}),'FontSize',8);
        else,    title(dv_labels{d},'FontSize',8); end
    end
end
sgtitle('N significant frequencies per hypothesis (per-hypothesis test, monkey avg)', ...
    'FontSize',11,'FontWeight','bold');
print(f5, fullfile(save_dir,'compare_hypotheses_nsig.pdf'),'-dpdf');
fprintf('Saved: %s\n', fullfile(save_dir,'compare_hypotheses_nsig.pdf'));

%% Figure 6: N significant freqs per pair (PAIRED test)
n_sig_pair = zeros(nP,nC,nDV);
for p = 1:nP; for c = 1:nC; for d = 1:nDV
    s = diff_sig{p,c,d};
    if ~isempty(s), n_sig_pair(p,c,d) = sum(s); end
end; end; end

f6 = figure('Name','N sig freqs per paired comparison', ...
    'Units','centimeters','Position',[1 1 54 28]);
set(f6,'PaperUnits','centimeters','PaperSize',[54 28],'PaperPosition',[0 0 54 28]);

for d = 1:nDV
    for p = 1:nP
        ax6 = subplot(nDV,nP,(d-1)*nP+p);
        vals = squeeze(n_sig_pair(p,:,d));
        b = bar(1:nC, vals, 'FaceColor','flat');
        for c = 1:nC, b.CData(c,:) = cmp_col(c,:); end
        set(ax6,'XTick',1:nC,'XTickLabel',{'H2−H1','H3−H2','H4−H1'}, ...
            'XTickLabelRotation',30,'FontSize',7);
        ylabel('n sig freqs','FontSize',7);
        if d==1, title(sprintf('%s — %s',pipe_labels{p},dv_labels{d}),'FontSize',8);
        else,    title(dv_labels{d},'FontSize',8); end
    end
end
sgtitle({'N significant frequencies per paired comparison (monkey avg)', ...
    'tests positions / difficulty / channels — significant only if exceeds Jensen advantage'}, ...
    'FontSize',11,'FontWeight','bold');
print(f6, fullfile(save_dir,'compare_hypotheses_paired_nsig.pdf'),'-dpdf');
fprintf('Saved: %s\n', fullfile(save_dir,'compare_hypotheses_paired_nsig.pdf'));

fprintf('\nAll figures saved to: %s\n', save_dir);

%% ── Local functions ──────────────────────────────────────────────────

function [val, thr, perm, freq_out] = load_field(filepath, val_field, thr_field, perm_field, freq_in)
% Load observed curve, threshold, and per-perm null array (if present).
val = []; thr = NaN; perm = []; freq_out = freq_in;
if ~isfile(filepath), return; end
s = load(filepath);
if isfield(s, val_field), val = s.(val_field); end
if isfield(s, thr_field), thr = s.(thr_field); end
if ~isempty(perm_field) && isfield(s, perm_field), perm = s.(perm_field); end
if isfield(s,'freq') && isempty(freq_in), freq_out = s.freq; end
end

function do_shade(freq_axis, val, thr, col)
% Shade frequency regions where val >= thr with semi-transparent fill.
if isempty(val)||isnan(thr)||numel(val)~=numel(freq_axis), return; end
sig = val(:)' >= thr;
if ~any(sig), return; end
yl = ylim;
starts = find(diff([0 sig]) ==  1);
ends   = find(diff([sig 0]) == -1);
for k = 1:numel(starts)
    xp = [freq_axis(starts(k)) freq_axis(ends(k)) freq_axis(ends(k)) freq_axis(starts(k))];
    yp = [yl(1) yl(1) yl(2) yl(2)];
    fill(xp, yp, col, 'FaceAlpha',0.12,'EdgeColor','none','HandleVisibility','off');
end
ylim(yl);
end
