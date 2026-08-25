clear all
close all
clc

% Code description:
% -----------------
% Pools klecks and hermes into a SINGLE result for the TFR and the pre-stimulus
% power analyses, and keeps the per-animal results alongside it. Reads what
% TFRdiff_hitmiss.m and powdiff_hitmiss.m already wrote, so nothing has to be
% recomputed from the raw data.
%
% Why channel-AVERAGED and not per channel:
% -----------------------------------------
% The two arrays share no channel identity - klecks is V4- 65 .. V4-128 and
% hermes is V4-  1 .. V4- 63, different electrodes in different brains. A
% pooled per-channel map would therefore be 125 channels each of which exists
% in exactly one animal, which is a patchwork, not a pooled map. Every pooled
% quantity below is the average over that session's own channels first, so a
% session contributes one map (TFR) or one spectrum (POW), and the pooling is
% across the 48 sessions of both animals.
%
% Unit of observation (differs by analysis, and is forced, not chosen):
% ---------------------------------------------------------------------
%   TFR  - TRIAL. fun_tfr_perm_session relabels trials within a session, so the
%          null is a trial-level exchange stratified by session. Sessions are
%          then averaged and the max-pixel correction runs on that average.
%   POW  - SESSION. fooof has to be fit on an averaged spectrum, so trials are
%          already gone by the time the statistic exists. The group test is a
%          sign-flip across sessions, i.e. a random-effects test.
%
% No stratification of hit/miss trial counts:
% -------------------------------------------
% Both analyses use the geometric mean, mean_i(log10 P_i), whose expectation
% does not depend on the number of trials. Matching hit and miss counts would
% discard roughly half the hit trials and buy nothing. (This is NOT true of the
% magnitude of an average, e.g. an ERP RMS, which is inflated by ~sigma^2/N.)
%
% Reporting:
% ----------
% For every pooled number the script also prints how many of the 48 sessions
% fall on the same side as the pooled mean. That count is the honest companion
% to the p value: a pooled effect carried by one animal shows up as a session
% count near chance even when p is small.
%
% Outputs:
%   results_combined/group_TFR/<run_name>/tfr_hitmiss_pooled_both.mat
%   results_combined/group_POW/<pow_run>/pow_hitmiss_pooled_both.mat
%   Plots/Hitvsmiss/pooled_both/*.{fig,pdf}

%% Specify paths

addpath /opt/fieldtrip_github/
ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/LFP_MUA_data/functions
clc

%% Define important variables

do_pow       = 1;              % pool the pre-stimulus power analysis
do_tfr       = 1;              % pool the TFR analysis (reads ~318 MB per session)
do_erp_lfp   = 1;              % pool the LFP evoked difference
do_erp_mua   = 1;              % pool the MUA evoked difference
plotting     = 1;

% ERP options (apply to both do_erp_lfp and do_erp_mua)
erp_signed   = 1;              % pool the signed channel-averaged difference
erp_rect     = 1;              % also pool the rectified difference, mean|d| over
                               % channels. Non-negative, so it cannot cancel
                               % across animals, at the cost of losing the sign
                               % and of a one-sided test. See the note below.
erp_test_win = [];             % latency (s) the max-statistic is corrected over.
                               % [] = the whole epoch, matching erpdiff_*.m.
                               % Narrowing it (e.g. [-0.1 0.4]) lowers the
                               % threshold and gains sensitivity in the response,
                               % but then the pre-stimulus effect is untested.
erp_recache  = 0;              % 1 = rebuild the cached per-session permutation
                               % waveforms even if they already exist

nperm        = 1000;
alpha_level  = 0.05;
min_sess_frac= 0.8;            % a bin is tested only if this fraction of the
                               % pooled sessions contributes a value to it

tfr_run      = 'unified_4to80Hz_2cyc';
pow_run      = 'prestim_2to80Hz';

base         = '/mnt/hpc/projects/MWSampling/4Shivangi';
combdir      = fullfile(base,'results_combined');
plotfolder   = fullfile(base,'Plots','Hitvsmiss','pooled_both');
if ~isdir(plotfolder), mkdir(plotfolder), end

animals(1).name            = 'klecks';
animals(1).resultsfolder   = fullfile(base,'results_klecks');
animals(1).remove_sessions = [2 15 16 18];

animals(2).name            = 'hermes';
animals(2).resultsfolder   = fullfile(base,'results_hermes');
animals(2).remove_sessions = [1 3 21];

%% ====================== POW: pool the two animals ====================== %%
% powdiff_hitmiss.m already saves the per-session, per-channel fooof output in
% its pooled file (FH/FM periodic, PH/PM raw, OH/OM offset, EH/EM exponent,
% all nSess x nChan [x nFreq]). Averaging over channels and concatenating the
% session dimension is all the pooling that is needed.

if do_pow

    F = []; P = []; O = []; E = [];      % hit-miss, sessions x (freq)
    FHa = []; FMa = [];                  % periodic spectra per condition
    PHa = []; PMa = [];                  % raw log power per condition
    animal_of = [];                      % which animal each session came from
    freqs = [];
    per_animal = struct();

    for ianimal = 1:length(animals)

        f = fullfile(animals(ianimal).resultsfolder,'group_POW',pow_run, ...
            'pow_hitmiss_pooled.mat');
        if ~exist(f,'file')
            error('%s: %s not found - run powdiff_hitmiss.m first', ...
                animals(ianimal).name, f)
        end
        S = load(f);

        if isempty(freqs)
            freqs = S.freqs(:)';
        elseif ~isequal(S.freqs(:)', freqs)
            error('the two animals have different frequency grids')
        end

        % average over channels, per session. The per-condition spectra are
        % kept as well as their difference, because the figure has to show
        % what the fooof-corrected spectra actually look like, not only the
        % contrast between them.
        fh = squeeze(nanmean(S.FH, 2));  fm = squeeze(nanmean(S.FM, 2));  % periodic
        ph = squeeze(nanmean(S.PH, 2));  pm = squeeze(nanmean(S.PM, 2));  % raw log power
        dF = fh - fm;                               % nSess x nFreq
        dP = ph - pm;                               % nSess x nFreq
        dO = nanmean(S.OH - S.OM, 2);               % nSess x 1
        dE = nanmean(S.EH - S.EM, 2);               % nSess x 1

        F = [F; dF]; P = [P; dP]; O = [O; dO]; E = [E; dE]; %#ok<AGROW>
        FHa = [FHa; fh]; FMa = [FMa; fm];                    %#ok<AGROW>
        PHa = [PHa; ph]; PMa = [PMa; pm];                    %#ok<AGROW>
        animal_of = [animal_of; ianimal*ones(size(dF,1),1)]; %#ok<AGROW>

        per_animal(ianimal).name = animals(ianimal).name;
        per_animal(ianimal).dF   = dF;
        per_animal(ianimal).FH   = fh;
        per_animal(ianimal).FM   = fm;
        per_animal(ianimal).dO   = dO;
        per_animal(ianimal).dE   = dE;
        per_animal(ianimal).n    = size(dF,1);

        fprintf('POW %-7s: %d sessions\n', animals(ianimal).name, size(dF,1))
    end

    nS = size(F,1);
    fprintf('POW pooled: %d sessions\n', nS)

    [stat_flat, thr_flat] = signflip_test(F, nperm, alpha_level);
    [stat_raw , thr_raw ] = signflip_test(P, nperm, alpha_level);
    [stat_off , thr_off ] = signflip_test(O, nperm, alpha_level);
    [stat_exp , thr_exp ] = signflip_test(E, nperm, alpha_level);

    fprintf('\n--- POW, both animals pooled (%d sessions) ---\n', nS)
    report_pooled('aperiodic offset  ', stat_off, O, animal_of, animals, [])
    report_pooled('aperiodic exponent', stat_exp, E, animal_of, animals, [])
    report_pooled('periodic power    ', stat_flat, F, animal_of, animals, freqs)

    outdir = fullfile(combdir,'group_POW',pow_run);
    if ~isdir(outdir), mkdir(outdir), end
    save(fullfile(outdir,'pow_hitmiss_pooled_both.mat'), ...
        'stat_flat','stat_raw','stat_off','stat_exp', ...
        'thr_flat','thr_raw','thr_off','thr_exp', ...
        'F','P','O','E','FHa','FMa','PHa','PMa', ...
        'animal_of','per_animal','freqs','nS','nperm','-v7.3')
    fprintf('  saved to %s\n', outdir)

    if plotting
        fig = figure('Name','POW hit-miss, both animals', ...
            'Units','normalized','Position',[0.03 0.08 0.94 0.82]);
        xl = [min(freqs) max(freqs)];

        % (a) raw log power spectra per condition, as in powdiff_hitmiss.m
        subplot(2,3,1)
        hh = plot(freqs, nanmean(PHa,1), 'k','LineWidth',1.6); hold on
        hm = plot(freqs, nanmean(PMa,1), 'r','LineWidth',1.6);
        legend([hh hm], {'hit','miss'}, 'Location','best'), legend boxoff
        xlabel('Frequency (Hz)'), ylabel('log_{10} power')
        title(sprintf('(a) raw spectrum, %d sessions', nS))
        xlim(xl), box on

        % (b) THE fooof-corrected spectra per condition. This is the panel that
        % shows what the periodic power actually looks like once the 1/f
        % background is removed; the difference in (c) is the contrast between
        % these two curves.
        subplot(2,3,2)
        hh = plot(freqs, nanmean(FHa,1), 'k','LineWidth',1.6); hold on
        hm = plot(freqs, nanmean(FMa,1), 'r','LineWidth',1.6);
        plot(freqs, zeros(size(freqs)), 'Color',[.6 .6 .6])
        shade_sig(freqs, stat_flat.mask)
        legend([hh hm], {'hit','miss'}, 'Location','best'), legend boxoff
        xlabel('Frequency (Hz)'), ylabel('log_{10} periodic power')
        title('(b) 1/f removed (fooof), pooled')
        xlim(xl), box on

        % (c) their difference, with SEM over sessions and the corrected mask
        subplot(2,3,3)
        d  = nanmean(F,1);
        se = nanstd(F,0,1)./sqrt(sum(~isnan(F),1));
        plot(freqs, d, 'k','LineWidth',1.8), hold on
        plot(freqs, d+se, 'k:', freqs, d-se, 'k:')
        plot(freqs, zeros(size(freqs)), 'Color',[.6 .6 .6])
        shade_sig(freqs, stat_flat.mask)
        xlabel('Frequency (Hz)'), ylabel('hit - miss')
        title('(c) periodic difference \pm SEM')
        xlim(xl), box on

        % (d) fooof-corrected spectra per animal, hit and miss
        subplot(2,3,4)
        cols = {'b','r'};
        h = gobjects(1,2*length(animals)); lbl = {};
        for ianimal = 1:length(animals)
            h(2*ianimal-1) = plot(freqs, nanmean(per_animal(ianimal).FH,1), ...
                'Color', cols{ianimal}, 'LineWidth',1.5); hold on
            h(2*ianimal)   = plot(freqs, nanmean(per_animal(ianimal).FM,1), ...
                'Color', cols{ianimal}, 'LineWidth',1.5, 'LineStyle','--');
            lbl = [lbl, {[per_animal(ianimal).name ' hit'], ...
                         [per_animal(ianimal).name ' miss']}]; %#ok<AGROW>
        end
        plot(freqs, zeros(size(freqs)), 'Color',[.6 .6 .6])
        legend(h, lbl, 'Location','best'), legend boxoff
        xlabel('Frequency (Hz)'), ylabel('log_{10} periodic power')
        title('(d) 1/f removed, per animal')
        xlim(xl), box on

        % (e) per-animal difference - where a disagreement between them shows
        subplot(2,3,5)
        h = gobjects(1,length(animals));
        for ianimal = 1:length(animals)
            h(ianimal) = plot(freqs, nanmean(per_animal(ianimal).dF,1), ...
                cols{ianimal}, 'LineWidth',1.6); hold on
        end
        plot(freqs, zeros(size(freqs)), 'Color',[.6 .6 .6])
        legend(h, {per_animal.name}, 'Location','best'), legend boxoff
        xlabel('Frequency (Hz)'), ylabel('hit - miss')
        title('(e) difference per animal')
        xlim(xl), box on

        % (f) session consistency per frequency. A pooled effect the animals
        % share sits well above 50%; one carried by a single animal does not.
        subplot(2,3,6)
        agree = nan(size(freqs));
        for k = 1:length(freqs)
            s = sign(F(:,k)); s = s(~isnan(s) & s~=0);
            if isempty(s), continue, end
            agree(k) = max(sum(s>0), sum(s<0)) / numel(s);
        end
        plot(freqs, 100*agree, 'k','LineWidth',1.6), hold on
        plot(freqs, 50*ones(size(freqs)), 'Color',[.6 .6 .6])
        xlabel('Frequency (Hz)'), ylabel('% sessions on the majority side')
        title('(f) session consistency (50% = chance)')
        ylim([40 100]), xlim(xl), box on

        sgtitle(sprintf('Pre-stimulus power, hit - miss, %d sessions pooled over both animals', nS))

        cd(plotfolder)
        savefig(fig,'POWdiff_pooled_both')
        set(fig,'PaperPositionMode','auto')
        print(fig, 'POWdiff_pooled_both','-dpdf','-fillpage')
    end
end

%% ====================== TFR: pool the two animals ====================== %%
% Reads the per-session tfr_perm.mat files. Each one is ~318 MB, so the loop
% averages over channels immediately and keeps only the freq x time result,
% which is a few MB per session.

if do_tfr

    real_sum = []; perm_sum = []; count = [];
    sess_real = []; animal_of_tfr = []; nSess = 0;
    freqs = []; times = [];
    per_animal_tfr = struct();

    for ianimal = 1:length(animals)

        animalName = animals(ianimal).name;
        datafolder = animals(ianimal).resultsfolder;
        fprintf('\n=============== TFR %s ===============\n', animalName)

        temp = dir(datafolder);
        session_names = {};
        ii = 0;
        for i = 1:length(temp)
            if temp(i).isdir && ~isempty(strfind(temp(i).name, animalName)) %#ok<STREMP>
                ii = ii+1;
                session_names{ii,1} = temp(i).name; %#ok<SAGROW>
            end
        end
        keep_sessions = setdiff(1:length(session_names), animals(ianimal).remove_sessions);
        output_paths  = cellfun(@(x) fullfile(datafolder,x,'TFR',tfr_run), ...
            session_names, 'uniform',0);

        n_this = 0;
        for k = 1:length(keep_sessions)
            isess = keep_sessions(k);
            f = fullfile(output_paths{isess},'tfr_perm.mat');
            if ~exist(f,'file')
                fprintf('  %s: no tfr_perm.mat, skipped\n', session_names{isess})
                continue
            end
            S = load(f);
            if isfield(S,'skipped') && S.skipped
                fprintf('  %s: skipped by the job (%s)\n', session_names{isess}, S.reason)
                continue
            end

            if isempty(real_sum)
                freqs = S.freq; times = S.time;
                real_sum = zeros(length(freqs), length(times), 'single');
                perm_sum = zeros(nperm, length(freqs), length(times), 'single');
                count    = zeros(length(freqs), length(times), 'single');
            end
            if ~isequal(S.freq, freqs) || ~isequal(S.time, times)
                fprintf('  %s: different time-frequency grid, skipped\n', session_names{isess})
                continue
            end

            % ---- average over this session's own channels, then pool.
            % squeeze over the channel dimension: diff_real is chan x freq x
            % time, diff_perm is nperm x chan x freq x time.
            dr = squeeze(nanmean(S.diff_real, 1));            % freq x time
            dp = squeeze(nanmean(S.diff_perm, 2));            % nperm x freq x time

            valid = ~isnan(dr);
            dr(~valid) = 0;
            dp(isnan(dp)) = 0;

            real_sum = real_sum + single(dr);
            perm_sum = perm_sum + single(dp);
            count    = count    + single(valid);

            nSess  = nSess  + 1;
            n_this = n_this + 1;
            sess_real(:,:,nSess) = single(dr); %#ok<SAGROW>
            animal_of_tfr(nSess) = ianimal;    %#ok<SAGROW>

            fprintf('  %s: %d hits, %d misses\n', session_names{isess}, S.n_hit, S.n_miss)
            clear S dr dp valid
        end

        per_animal_tfr(ianimal).name = animalName;
        per_animal_tfr(ianimal).n    = n_this;
    end

    if nSess < 2, error('fewer than 2 sessions pooled'), end

    % ---- session average, per bin
    cnt = count;
    cnt(cnt < max(1, min_sess_frac*nSess)) = NaN;
    fprintf('\n  %d / %d bins have at least %d of %d sessions\n', ...
        sum(~isnan(cnt(:))), numel(cnt), ceil(max(1,min_sess_frac*nSess)), nSess)

    real_avg = real_sum ./ cnt;
    perm_avg = perm_sum ./ reshape(cnt,[1 size(cnt)]);
    clear real_sum perm_sum

    % ---- normalise each bin against its own permutation distribution, then
    % correct with the max pixel over the whole freq x time map. Same
    % construction as the per-animal test in TFRdiff_hitmiss.m, applied to the
    % channel-averaged map because that is what is pooled and plotted.
    mu = nanmean(perm_avg,1); sd = nanstd(perm_avg,0,1);
    sd(sd == 0) = NaN;
    z_real   = (real_avg - squeeze(mu)) ./ squeeze(sd);
    perm_avg = (perm_avg - mu) ./ sd;

    maxdist = max(abs(reshape(perm_avg, nperm, [])), [], 2);
    thr     = quantile(maxdist, 1 - alpha_level);
    mask    = abs(z_real) > thr;
    p_corr  = reshape(mean(maxdist >= abs(z_real(:))', 1), size(z_real));
    p_corr(isnan(z_real)) = NaN;
    clear perm_avg

    fprintf('\n--- TFR, both animals pooled ---\n')
    fprintf('  %d sessions (%s), %d permutations\n', nSess, ...
        strjoin(arrayfun(@(a) sprintf('%s %d', a.name, a.n), per_animal_tfr, ...
        'uniform',0), ', '), nperm)
    fprintf('  |z| threshold at p<%g corrected: %.2f\n', alpha_level, thr)
    fprintf('  %d / %d usable bins significant\n', sum(mask(:)), sum(~isnan(z_real(:))))
    if any(mask(:))
        [fi, ti] = find(mask);
        fprintf('  extent: %g-%g Hz, %g-%g s\n', ...
            min(freqs(fi)), max(freqs(fi)), min(times(ti)), max(times(ti)))
        fprintf('  smallest corrected p: %.4f\n', min(p_corr(~isnan(p_corr))))

        % session consistency at the strongest bin
        [~, ib] = max(abs(z_real(:)));
        [fb, tb] = ind2sub(size(z_real), ib);
        v = squeeze(sess_real(fb,tb,:));
        s = sign(v); s = s(~isnan(s) & s~=0);
        fprintf('  strongest bin %g Hz %g s: %d/%d sessions on the same side', ...
            freqs(fb), times(tb), max(sum(s>0),sum(s<0)), numel(s))
        for ianimal = 1:length(animals)
            vi = v(animal_of_tfr == ianimal);
            fprintf(' | %s %+.4f', animals(ianimal).name, nanmean(vi));
        end
        fprintf('\n')
    end

    outdir = fullfile(combdir,'group_TFR',tfr_run);
    if ~isdir(outdir), mkdir(outdir), end
    save(fullfile(outdir,'tfr_hitmiss_pooled_both.mat'), ...
        'real_avg','z_real','mask','p_corr','thr','maxdist', ...
        'sess_real','animal_of_tfr','count','freqs','times', ...
        'nSess','nperm','per_animal_tfr','-v7.3')
    fprintf('  saved to %s\n', outdir)

    if plotting
        fig = figure('Name','TFR hit-miss, both animals', ...
            'Units','normalized','Position',[0.05 0.15 0.9 0.6]);
        clim = quantile(abs(real_avg(~isnan(real_avg))), 0.99);

        subplot(1,3,1)
        plot_map(times, freqs, real_avg, mask, clim)
        title(sprintf('pooled, %d sessions', nSess))
        ylabel('Frequency (Hz)')

        for ianimal = 1:length(animals)
            subplot(1,3,1+ianimal)
            m = nanmean(sess_real(:,:,animal_of_tfr == ianimal), 3);
            plot_map(times, freqs, m, false(size(m)), clim)
            title(sprintf('%s, %d sessions', animals(ianimal).name, ...
                per_animal_tfr(ianimal).n))
        end

        cd(plotfolder)
        savefig(fig,'TFRdiff_pooled_both')
        set(fig,'PaperPositionMode','auto')
        print(fig, 'TFRdiff_pooled_both','-dpdf','-fillpage')
    end
end

%% ================ ERP (LFP / MUA): pool the two animals ================ %%
% Same session-level construction as POW and TFR: average over that session's
% own channels first, then pool the sessions of both animals.
%
% Two variants, because they answer different questions:
%
%   SIGNED     mean over channels of (hit - miss). Keeps the direction of the
%              effect. Measured on the observed data, MUA pools essentially
%              perfectly this way (cross-animal r = +0.89, retention 0.996,
%              both animals peaking near 92-95 ms with the same sign), while
%              the LFP does not (r = +0.30, retention 0.754) because klecks
%              peaks -0.73 at 93 ms and hermes +0.77 at 132 ms.
%
%   RECTIFIED  mean over channels of |hit - miss|. Non-negative, so animals
%              cannot cancel: it lifts LFP retention from 0.754 to 0.859. The
%              cost is real - the sign is gone, so the direction of the effect
%              is no longer readable, and the two animals' peaks still sit at
%              different latencies, so a pooled rectified curve shows both as
%              separate humps rather than one effect.
%
% The rectified null MUST be rectified at the same stage as the observed data,
% which is why the permutation waveforms below are stored in both forms. Its
% test is ONE-SIDED (the statistic cannot be negative); the signed test stays
% two-sided, matching erpdiff_*.m.
%
% NOTE ON WHAT IS BEING TESTED: erpdiff_lfp.m and erpdiff_mua.m take the max
% over the whole epoch, so the test covers the pre-stimulus period as well as
% the response. For MUA the hit-miss difference is already positive a full
% second before the stimulus (hermes 28/28 sessions), which is a real sustained
% effect and not a baselining error - the common baseline cancels in the
% difference and cannot remove it. erp_test_win exists so that can be scoped
% deliberately rather than by accident.

erp_jobs = {};
if do_erp_lfp
    erp_jobs{end+1} = struct('sig','lfp','sub','ERP_LFP', ...
        'hit','hit_lfp_avg.mat','miss','miss_lfp_avg.mat'); %#ok<SAGROW>
end
if do_erp_mua
    erp_jobs{end+1} = struct('sig','mua','sub','ERP_MUA', ...
        'hit','hit_mua_avg.mat','miss','miss_mua_avg.mat'); %#ok<SAGROW>
end

for ijob = 1:numel(erp_jobs)

    J = erp_jobs{ijob};
    fprintf('\n=============== ERP %s, both animals ===============\n', upper(J.sig))

    obs_sg = []; obs_rc = [];          % nSess x nTime, observed per session
    sum_sg = []; sum_rc = [];          % nperm x nTime, running sums of the null
    animal_of_erp = []; tvec = []; nS = 0;
    per_animal_erp = struct();

    for ianimal = 1:length(animals)

        animalName = animals(ianimal).name;
        datafolder = animals(ianimal).resultsfolder;

        temp = dir(datafolder);
        session_names = {};
        ii = 0;
        for i = 1:length(temp)
            if temp(i).isdir && ~isempty(strfind(temp(i).name, animalName)) %#ok<STREMP>
                ii = ii+1;
                session_names{ii,1} = temp(i).name; %#ok<SAGROW>
            end
        end
        keep_sessions = setdiff(1:length(session_names), animals(ianimal).remove_sessions);

        A = load(fullfile(datafolder,'critical_time','all_channels.mat'));
        all_channels = A.all_channels(:);

        n_this = 0;
        for k = 1:length(keep_sessions)
            isess = keep_sessions(k);
            op = fullfile(datafolder, session_names{isess}, J.sub);

            fh = fullfile(op,'ERP_real','norm_hit_timelock.mat');
            fm = fullfile(op,'ERP_real','norm_miss_timelock.mat');
            if ~exist(fh,'file') || ~exist(fm,'file')
                fprintf('  %s: no ERP_real, skipped\n', session_names{isess})
                continue
            end
            H = load(fh); H = H.norm_hit_timelock;
            M = load(fm); M = M.norm_miss_timelock;

            keepch = ismember(H.label, all_channels);
            if ~any(keepch)
                fprintf('  %s: no channel in all_channels, skipped\n', session_names{isess})
                continue
            end

            if isempty(tvec)
                tvec  = H.time(:)';
                ntime = length(tvec);
                sum_sg = zeros(nperm, ntime, 'single');
                sum_rc = zeros(nperm, ntime, 'single');
            end
            if length(H.time) ~= ntime
                fprintf('  %s: different time axis, skipped\n', session_names{isess})
                continue
            end

            d = H.avg(keepch,:) - M.avg(keepch,:);

            nS     = nS     + 1;
            n_this = n_this + 1;
            obs_sg(nS,:) = nanmean(d, 1);
            obs_rc(nS,:) = nanmean(abs(d), 1);
            animal_of_erp(nS) = ianimal; %#ok<SAGROW>

            % permuted waveforms for this session, cached on disk
            Pp = erp_perm_pooled(op, keepch, nperm, ntime, J, erp_recache);
            sum_sg = sum_sg + Pp.sg;
            sum_rc = sum_rc + Pp.rc;

            fprintf('  %s: %d channels\n', session_names{isess}, sum(keepch))
            clear H M d Pp
        end

        per_animal_erp(ianimal).name = animalName;
        per_animal_erp(ianimal).n    = n_this;
    end

    if nS < 2, error('%s: fewer than 2 sessions pooled', J.sig), end

    % ---- session averages, observed and null
    obs_sg_m = nanmean(obs_sg, 1);
    obs_rc_m = nanmean(obs_rc, 1);
    null_sg  = double(sum_sg) / nS;
    null_rc  = double(sum_rc) / nS;

    % ---- window the max-statistic is corrected over
    if isempty(erp_test_win)
        wtest = true(1, ntime);
    else
        wtest = tvec >= erp_test_win(1) & tvec <= erp_test_win(2);
    end
    fprintf('  %d sessions, correcting over %g to %g s (%d of %d samples)\n', ...
        nS, tvec(find(wtest,1)), tvec(find(wtest,1,'last')), sum(wtest), ntime)

    erp = struct('sig', J.sig, 'time', tvec, 'nS', nS, ...
        'animal_of', animal_of_erp, 'test_win', erp_test_win, ...
        'per_animal', per_animal_erp);

    if erp_signed
        % two-sided: pool the per-permutation max AND min over the window
        dist = [max(null_sg(:,wtest),[],2); min(null_sg(:,wtest),[],2)];
        lim_hi = quantile(dist, 1 - alpha_level/2);
        lim_lo = quantile(dist,     alpha_level/2);
        mask   = obs_sg_m >= lim_hi | obs_sg_m <= lim_lo;
        erp.signed = struct('obs',obs_sg,'mean',obs_sg_m, ...
            'lim_hi',lim_hi,'lim_lo',lim_lo,'mask',mask);
        report_erp('signed   ', obs_sg, obs_sg_m, mask, tvec, animal_of_erp, animals)
    end

    if erp_rect
        % one-sided: the statistic is non-negative, so only the upper tail
        dist   = max(null_rc(:,wtest),[],2);
        lim_hi = quantile(dist, 1 - alpha_level);
        mask   = obs_rc_m >= lim_hi;
        erp.rect = struct('obs',obs_rc,'mean',obs_rc_m, ...
            'lim_hi',lim_hi,'lim_lo',NaN,'mask',mask);
        report_erp('rectified', obs_rc, obs_rc_m, mask, tvec, animal_of_erp, animals)
    end

    outdir = fullfile(combdir,'group_ERP',J.sig);
    if ~isdir(outdir), mkdir(outdir), end
    save(fullfile(outdir,'erp_hitmiss_pooled_both.mat'),'erp','-v7.3')
    fprintf('  saved to %s\n', outdir)

    if plotting
        vars = {}; if erp_signed, vars{end+1}='signed'; end %#ok<SAGROW>
        if erp_rect, vars{end+1}='rect'; end %#ok<SAGROW>

        fig = figure('Name',sprintf('ERP %s hit-miss, both animals',upper(J.sig)), ...
            'Units','normalized','Position',[0.05 0.15 0.9 0.75]);
        cols = {'b','r'};

        for iv = 1:numel(vars)
            V = erp.(vars{iv});

            % pooled
            subplot(numel(vars), 2, 2*iv-1)
            plot(tvec, V.mean, 'k','LineWidth',1.8), hold on
            plot(tvec([1 end]), [V.lim_hi V.lim_hi], 'Color',[.5 .5 .5])
            if ~isnan(V.lim_lo)
                plot(tvec([1 end]), [V.lim_lo V.lim_lo], 'Color',[.5 .5 .5])
            end
            plot(tvec, zeros(size(tvec)), 'Color',[.8 .8 .8])
            shade_sig(tvec, V.mask)
            xlabel('Time (s)'), ylabel('hit - miss')
            title(sprintf('%s %s, pooled (%d sessions)', upper(J.sig), vars{iv}, nS))
            xlim([tvec(1) tvec(end)]), box on

            % per animal
            subplot(numel(vars), 2, 2*iv)
            hh = gobjects(1,length(animals));
            for ianimal = 1:length(animals)
                hh(ianimal) = plot(tvec, nanmean(V.obs(animal_of_erp==ianimal,:),1), ...
                    cols{ianimal}, 'LineWidth',1.5); hold on
            end
            plot(tvec, zeros(size(tvec)), 'Color',[.8 .8 .8])
            legend(hh, {animals.name}, 'Location','best'), legend boxoff
            xlabel('Time (s)'), ylabel('hit - miss')
            title(sprintf('%s %s, per animal', upper(J.sig), vars{iv}))
            xlim([tvec(1) tvec(end)]), box on
        end

        cd(plotfolder)
        savefig(fig, sprintf('ERPdiff_%s_pooled_both', J.sig))
        set(fig,'PaperPositionMode','auto')
        print(fig, sprintf('ERPdiff_%s_pooled_both', J.sig), '-dpdf','-fillpage')
    end
end

%% ============================== helpers ================================ %%

function [stat, thr] = signflip_test(D, nrand, alpha_level)
% Paired test across sessions with a sign-flip permutation. Same construction
% as the local function of the same name in powdiff_hitmiss.m; copied so this
% script stands alone.

sz = size(D); n = sz(1);
Dm = reshape(D, n, []);
cnt = sum(~isnan(Dm),1);
tobs = nanmean(Dm,1) ./ (nanstd(Dm,0,1) ./ sqrt(cnt));

signs   = sign(randn(nrand, n));
maxdist = zeros(nrand,1);
for p = 1:nrand
    X  = Dm .* signs(p,:)';
    tp = nanmean(X,1) ./ (nanstd(X,0,1) ./ sqrt(cnt));
    maxdist(p) = max(abs(tp));
end

thr = quantile(maxdist, 1 - alpha_level);
pv  = arrayfun(@(v) mean(maxdist >= abs(v)), tobs);
pv(~isfinite(tobs)) = NaN;

stat         = [];
stat.t       = reshape(tobs, [sz(2:end) 1]);
stat.d       = reshape(nanmean(Dm,1), [sz(2:end) 1]);
stat.mask    = abs(stat.t) > thr;
stat.p       = reshape(pv, [sz(2:end) 1]);
stat.n       = n;
stat.maxdist = maxdist;
end


function report_pooled(name, stat, D, animal_of, animals, freqs)
% Pooled effect, its p value, how many sessions agree with it, and the
% per-animal split. The session count is what separates an effect the animals
% share from one that a single animal carries.

[~, k] = max(abs(stat.t(:)));           % report at the strongest element
col = D(:,k);
s   = sign(col); s = s(~isnan(s) & s~=0);
agree = max(sum(s>0), sum(s<0));

fprintf('  %s: d=%+.4f  p=%.4f  %d/%d sessions agree', ...
    name, stat.d(k), stat.p(k), agree, numel(s));
for ianimal = 1:length(animals)
    vi = col(animal_of == ianimal);
    fprintf(' | %s %+.4f (%d/%d)', animals(ianimal).name, nanmean(vi), ...
        max(sum(sign(vi)>0), sum(sign(vi)<0)), sum(~isnan(vi)));
end
if ~isempty(freqs) && numel(stat.d) == numel(freqs)
    fprintf(' @ %g Hz', freqs(k));
end
fprintf('\n');
end


function plot_map(times, freqs, M, mask, clim)
% TFR map with the non-significant part faded, matching the style of
% TFRdiff_hitmiss.m.

imagesc(times, freqs, M), axis xy, hold on
caxis([-clim clim]), colormap(jet), colorbar
if any(mask(:))
    contour(times, freqs, double(mask), [.5 .5], 'k', 'LineWidth', 1.2)
end
xlabel('Time (s)')
end


function shade_sig(x, mask)
% Blue patches over the contiguous runs where mask is true. Copied from
% powdiff_hitmiss.m so this script stands alone.

mask = double(mask(:))';
if ~any(mask), return, end

onset  = find(conv(mask,[1 -1]) == 1);
offset = find(conv(mask,[1 -1]) == -1) - 1;
offset = offset(1:length(onset));

yl = ylim;
for i = 1:length(onset)
    x1 = x(onset(i)); x2 = x(offset(i));
    if x2 == x1
        w = (x(min(end,2)) - x(1))/2;
        x1 = x1 - w; x2 = x2 + w;
    end
    patch('Faces',[1 2 3 4], ...
        'Vertices',[x1 yl(1); x2 yl(1); x2 yl(2); x1 yl(2)], ...
        'FaceColor','blue','FaceAlpha',0.2, ...
        'EdgeColor','blue','EdgeAlpha',0.2)
end
ylim(yl)
uistack(findobj(gca,'Type','patch'),'bottom')
end


function P = erp_perm_pooled(op, keepch, nperm, ntime, J, recache)
% Channel-averaged permuted hit-miss waveforms for one session, in both the
% signed and the rectified form, cached next to the permutation folders.
%
% The rectified version has to be built HERE, at the same stage as the observed
% data: rectifying the observed difference but not the null would make every
% sample significant. Both forms are stored so the two tests stay matched.
%
% The common baseline cancels in the difference, (hit-bsl)-(miss-bsl) =
% hit-miss, so bsl_avg is not needed - same as pooled_perm_limits in
% erpdiff_lfp.m.
%
% Reading 1000 permutation folders per session is the slow part of this whole
% script, so the result is cached. Delete perm_pooled_<sig>.mat, or set
% erp_recache = 1, to rebuild.

cf = fullfile(op, ['perm_pooled_' J.sig '.mat']);
if exist(cf,'file') && ~recache
    C = load(cf);
    if isequal(C.P.keepch, keepch) && C.P.nperm == nperm && size(C.P.sg,2) == ntime
        P = C.P;
        return
    end
end

P = struct('keepch', keepch, 'nperm', nperm, ...
    'sg', zeros(nperm, ntime, 'single'), ...
    'rc', zeros(nperm, ntime, 'single'));

for iperm = 1:nperm
    h = load(fullfile(op, num2str(iperm), J.hit));
    m = load(fullfile(op, num2str(iperm), J.miss));
    d = h.timelock.avg(keepch,:) - m.timelock.avg(keepch,:);
    P.sg(iperm,:) = nanmean(d, 1);
    P.rc(iperm,:) = nanmean(abs(d), 1);
end

save(cf, 'P', '-v7.3')
end


function report_erp(name, obs, obsm, mask, tvec, animal_of, animals)
% Peak of the pooled waveform, how many sessions fall on the same side there,
% and the per-animal split - the same disclosure the POW and TFR sections print.

[~, ip] = max(abs(obsm));
v = obs(:, ip);
s = sign(v); s = s(~isnan(s) & s ~= 0);

fprintf('  %s: peak %+.4f at %.0f ms, %d/%d sessions agree', ...
    name, obsm(ip), tvec(ip)*1000, max(sum(s>0), sum(s<0)), numel(s));
for ianimal = 1:length(animals)
    fprintf(' | %s %+.4f', animals(ianimal).name, nanmean(v(animal_of == ianimal)));
end
fprintf('\n');

if any(mask)
    on  = find(conv(double(mask(:))',[1 -1]) ==  1);
    off = find(conv(double(mask(:))',[1 -1]) == -1) - 1;
    off = off(1:length(on));
    fprintf('      significant %d/%d samples, %d run(s), first %g to %g s\n', ...
        sum(mask), numel(mask), length(on), tvec(on(1)), tvec(off(1)));
else
    fprintf('      nothing significant\n');
end
end
