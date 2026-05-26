% =====================================================================
% Phase progression across stimulus positions, per channel
%
% Question: for each recording channel, does the preferred phase
% (angle of complex coherence) depend *systematically* on stimulus
% position?
%
% Per (channel × frequency):
%   - Compute complex coherence c_pos = mean(y · exp(i·phase)) across
%     trials within each stimulus position. y is the DV weight (MUA /
%     LFP amplitude, RT, or binary hit_miss). Angle of c_pos = preferred
%     phase at that position; abs(c_pos) = coherence magnitude.
%   - Systematicity statistic: **circular-linear correlation**
%     r_cl between position index (1..nPos, peripheral→foveal) and the
%     preferred-phase vector (φ_1,…,φ_nPos), per Mardia/Jupp:
%         r_cl = sqrt[(rxc² + rxs² − 2·rxc·rxs·rcs) / (1 − rcs²)]
%     where rxc = corr(pos, cos φ), rxs = corr(pos, sin φ),
%           rcs = corr(cos φ, sin φ).
%     r_cl ∈ [0, 1], unsigned magnitude of the position→phase
%     association.
%   - Direction of progression (advance vs retreat with position) is
%     captured separately by `mean_step` = circular mean of consecutive
%     Δφ across positions. `step_phase` (per-step Δφ values, wrapped)
%     is also saved for the polar-histogram plots.
%
% Null: shuffle the position label across trials (within channel),
% recompute r_cl. Repeat nPerm times → per-channel p-value at each
% freq. Output variables retain the names R_obs / R_null so downstream
% plotting code is unchanged; they now hold the r_cl values.
%
% NOTE on position ordering: positions are taken from unique(trialinfo
% col 16) in sorted numerical order. r_cl is unsigned (doesn't care
% which direction phase advances) — but `mean_step` does. If the
% spatial layout isn't monotonic in the col 16 coding, reorder
% `positions` by hand to keep `mean_step`'s sign meaningful.
% =====================================================================
clearvars; close all; clc

%% Dependencies
addpath /opt/fieldtrip_github/; ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/Phase_coherence/functions
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/scanning/functions
addpath /mnt/hpc/projects/MWSampling/4Shivangi
clc

%% Settings
% Dependent variables used as the weight in the complex coherence:
%   mua      — MUA ERP amplitude (per trial × channel)
%   lfp      — LFP ERP amplitude (per trial × channel)
%   RT       — reaction time, hit trials only (per trial × channel,
%              hits selected via RT_trialinfo col 20 == 1)
%   hit_miss — binary hit (1) / miss (0) outcome (per trial); for binary
%              y, mean(y · exp(i·φ)) = (n_hits/n_total) · ITC_hits,
%              so the angle is the preferred phase that predicts hits.
dv_types = {'mua','lfp','RT','hit_miss'};
animals  = {'hermes','klecks'};
nPerm    = 1000;
alpha    = 0.05;
rng(2025);

%% Loop over animals
for a = 1:numel(animals)
    animalName = animals{a};
    fprintf('\n=== %s ===\n', animalName);

    data_load_folder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ...
        ['results_' animalName], 'multi_lin_reg', 'cp10_till_100');
    output_folder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ...
        ['results_' animalName], 'scanning', 'phase_progression', 'cp10_till_100');
    if ~exist(output_folder,'dir'), mkdir(output_folder); end

    cd(data_load_folder); load('ph_all_sess.mat'); load('frequency.mat');
    freq = frequency;

    [nTrials, nFreq, nCh] = size(ph_comb.phase_all);

    for d = 1:numel(dv_types)
        dv = dv_types{d};
        [Y_full, trlInfo, keepIdx, isPerCh] = get_dv_data(ph_comb, dv);

        fprintf('--- %s (%s) ---\n', animalName, upper(dv));

        % Positions taken from the trials this DV actually uses.
        positions = unique(trlInfo(keepIdx, 16));
        nPos      = numel(positions);
        if nPos < 2
            warning('%s/%s: <2 positions, skipping.', animalName, dv); continue
        end

        % Position index per trial (1..nPos in sorted-position order);
        % trials excluded by keepIdx get 0 and are filtered out below.
        pos_idx_all = zeros(nTrials, 1);
        for p = 1:nPos
            pos_idx_all(trlInfo(:,16) == positions(p) & keepIdx) = p;
        end

        dv_dir = fullfile(output_folder, dv);
        if ~exist(dv_dir,'dir'), mkdir(dv_dir); end

        % ── Submit one SLURM job per channel ──────────────────────────
        cfg = cell(1, nCh);
        for ch = 1:nCh
            cfg{ch}.ichan   = ch;
            cfg{ch}.nPerm   = nPerm;
            cfg{ch}.dv      = dv;
            cfg{ch}.infile  = data_load_folder;
            cfg{ch}.outfile = dv_dir;
            cfg{ch}.seed    = 2025 + ch;
        end

        slurmfun(@phase_progression_chan, cfg, ...
            'partition',     '8GB', ...
            'stopOnError',   false, ...
            'useUserPath',   true);

        % ── Aggregate per-channel results from disk ───────────────────
        pref_phase = nan(nCh, nFreq, nPos);     % φ(ch, f, pos)
        coh_mag    = nan(nCh, nFreq, nPos);     % |c|(ch, f, pos)
        step_phase = nan(nCh, nFreq, nPos-1);   % Δφ between consecutive positions
                                                %   (peripheral → foveal direction);
                                                %   wrapped to (-π, π], in radians.
        mean_step  = nan(nCh, nFreq);           % typical Δφ per channel × freq
                                                %   = angle(mean exp(i·Δφ_k))
        R_obs      = nan(nCh, nFreq);           % observed systematicity
        R_null     = nan(nCh, nFreq, nPerm);    % permutation null
        p_val      = nan(nCh, nFreq);           % per-freq p-value
        n_pos      = zeros(nCh, nPos);          % trials per (ch, pos)

        for ch = 1:nCh
            chan_file = fullfile(dv_dir, num2str(ch), 'phase_progression_chan.mat');
            if ~isfile(chan_file)
                warning('Missing per-channel file for %s/%s ch %d.', animalName, dv, ch);
                continue
            end
            tmp = load(chan_file);
            R_obs(ch, :)        = tmp.R_obs;
            R_null(ch, :, :)    = tmp.R_null;
            pref_phase(ch, :, :)= reshape(tmp.pref_phase, [1 nFreq nPos]);
            coh_mag(ch, :, :)   = reshape(tmp.coh_mag,    [1 nFreq nPos]);
            step_phase(ch, :, :)= reshape(tmp.step_phase, [1 nFreq nPos-1]);
            mean_step(ch, :)    = tmp.mean_step;
            p_val(ch, :)        = tmp.p_val;
            n_pos(ch, :)        = tmp.n_p;
        end

        save(fullfile(dv_dir, 'phase_progression.mat'), ...
            'pref_phase','coh_mag','step_phase','mean_step', ...
            'R_obs','R_null','p_val','positions','freq','n_pos','-v7.3');

        % ── Channel-average + per-channel significance summary ─────────
        % Per-channel max-stat correction across frequencies:
        % for each channel, take max across freq per perm → one threshold.
        tmax_per_chan   = squeeze(max(R_null, [], 2));      % nCh × nPerm
        thresh_per_chan = quantile(tmax_per_chan, 1-alpha, 2); % nCh × 1

        % Channel-average uses the same max-stat logic as the existing
        % coherence pipeline: average null curves across channels first,
        % then take the freq-wise max for a single threshold.
        R_null_chan_avg = mean(R_null, 1, 'omitnan');     % 1 × nFreq × nPerm
        R_null_chan_avg = squeeze(R_null_chan_avg);       % nFreq × nPerm
        R_obs_chan_avg  = mean(R_obs, 1, 'omitnan');      % 1 × nFreq

        tmax_chan_avg     = max(R_null_chan_avg, [], 1);  % 1 × nPerm
        thresh_chan_avg   = quantile(tmax_chan_avg, 1-alpha);

        save(fullfile(dv_dir, 'channel_avg_results.mat'), ...
            'R_obs','R_obs_chan_avg','R_null_chan_avg','tmax_chan_avg', ...
            'thresh_chan_avg','thresh_per_chan','freq','positions');

        fprintf('Saved %s (%s); channel-avg threshold = %.4f\n', ...
            upper(dv), animalName, thresh_chan_avg);
    end
end

%% Monkey-average

nAnimals = numel(animals);
for d = 1:numel(dv_types)
    dv = dv_types{d};
    R_animals = []; null_animals = []; thresh_animals = nan(1, nAnimals);
    for a = 1:nAnimals
        avg_file = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi', ...
            ['results_' animals{a}], 'scanning','phase_progression','cp10_till_100', ...
            dv, 'channel_avg_results.mat');
        if ~isfile(avg_file)
            warning('Channel-avg not found for %s (%s).', dv, animals{a}); continue
        end
        tmp = load(avg_file);
        R_animals    = cat(1, R_animals,    tmp.R_obs_chan_avg);
        null_animals = cat(3, null_animals, tmp.R_null_chan_avg);  % nFreq × nPerm × nAnimals
        thresh_animals(a) = tmp.thresh_chan_avg;
    end
    if size(R_animals,1) < nAnimals
        warning('Skipping monkey-avg for %s (not all animals present).', dv); continue
    end

    R_monkey_avg     = mean(R_animals, 1);
    null_monkey      = mean(null_animals, 3);            % nFreq × nPerm
    tmax_monkey      = max(null_monkey, [], 1);
    thresh_monkey    = quantile(tmax_monkey, 1-alpha);
    sig_freqs        = tmp.freq(R_monkey_avg >= thresh_monkey);

    save_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/results_combined', ...
        'scanning','phase_progression','cp10_till_100', dv);
    if ~exist(save_dir,'dir'), mkdir(save_dir); end
    save(fullfile(save_dir, 'monkey_avg_results.mat'), ...
        'R_monkey_avg','null_monkey','tmax_monkey','thresh_monkey', ...
        'R_animals','animals','sig_freqs','-v7.3');

    fprintf('%s monkey-avg threshold = %.4f, n sig freqs = %d\n', ...
        upper(dv), thresh_monkey, numel(sig_freqs));

    %% Quick plots: R_obs vs freq with threshold, per-channel heatmap
    plot_dir = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots', ...
        'scanning','phase_progression','cp10_till_100', dv);
    if ~exist(plot_dir,'dir'), mkdir(plot_dir); end

    f1 = figure('Visible','off','Position',[100 100 700 380]); hold on;
    plot(tmp.freq, R_monkey_avg, 'k', 'LineWidth', 2, 'DisplayName','Monkey avg');
    cols = lines(nAnimals);
    for a = 1:nAnimals
        plot(tmp.freq, R_animals(a,:), 'Color', cols(a,:), 'LineWidth', 1, ...
            'DisplayName', animals{a});
    end
    % Monkey-average max-stat threshold (solid red)
    yline(thresh_monkey,'-','Color',[0.85 0.15 0.15],'LineWidth',1.4, ...
        'DisplayName','Monkey 95% threshold');
    % Per-animal max-stat thresholds (dashed, in each animal's colour)
    for a = 1:nAnimals
        if isnan(thresh_animals(a)), continue; end
        yline(thresh_animals(a),'--','Color',cols(a,:),'LineWidth',1.0, ...
            'DisplayName',sprintf('%s 95%% threshold',animals{a}));
    end
    xlabel('Frequency (Hz)'); ylabel('r_{cl} (position–phase circular-linear correlation)');
    title(sprintf('%s — phase progression systematicity (channel avg)', upper(dv)));
    legend('Location','best','FontSize',7); grid on;
    saveas(f1, fullfile(plot_dir,'systematicity_vs_freq.pdf'));
    close(f1);
end

fprintf('\nDone. Per-channel results in results_<animal>/scanning/phase_progression/.\n');

%% =====================================================================
%% Local functions
%% =====================================================================
function [Y, trlInfo, keepIdx, isPerCh] = get_dv_data(ph_comb, dv)
% Return the DV values, the trialinfo aligned with phase_all, a
% keep-mask, and whether Y is per-channel.
%   Y       — either nTrials × nCh (mua, lfp, RT) or nTrials × 1 (hit_miss)
%   trlInfo — nTrials × nCols, gives stimulus position in col 16
%   keepIdx — logical(nTrials × 1), trials this DV is defined on
%             (e.g. RT keeps only hit trials)
%   isPerCh — true if Y has a per-channel column
switch dv
    case 'mua'
        Y       = ph_comb.MUA_ERP_ampl_all;
        trlInfo = ph_comb.trialinfo;
        keepIdx = true(size(Y,1), 1);
        isPerCh = true;
    case 'lfp'
        Y       = ph_comb.LFP_ERP_ampl_all;
        trlInfo = ph_comb.trialinfo;
        keepIdx = true(size(Y,1), 1);
        isPerCh = true;
    case 'RT'
        Y       = ph_comb.RT;
        trlInfo = ph_comb.RT_trialinfo;
        keepIdx = trlInfo(:,20) == 1;       % hits only
        isPerCh = true;
    case 'hit_miss'
        Y           = ph_comb.trialinfo(:,20);
        Y(Y == 5)   = 0;                    % recode miss code 5 → 0
        trlInfo     = ph_comb.trialinfo;
        keepIdx     = true(size(Y,1), 1);
        isPerCh     = false;
    otherwise
        error('Unknown DV: %s', dv);
end
end

