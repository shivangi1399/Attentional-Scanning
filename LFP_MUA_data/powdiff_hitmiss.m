
clear all
close all
clc

% Code description:
% -----------------
% Pre-stimulus power spectrum, hits vs misses, run separately for klecks and
% hermes. No time axis - this is the spectral counterpart of
% TFRdiff_hitmiss.m, built on pow_freq.m.
%
% The 1/f background is removed with fooof, and the two parts are tested
% separately, because they answer different questions:
%
%   aperiodic OFFSET   is the whole spectrum shifted up on hits? This is the
%                      broadband effect that dominates the TFR analysis.
%   aperiodic EXPONENT is the 1/f slope different? A slope change points at
%                      an arousal / excitation-inhibition difference rather
%                      than a simple gain change.
%   PERIODIC power     once the background is removed, is there a genuine
%                      oscillatory difference, and at which frequencies?
%
% Trials are combined as a GEOMETRIC mean (log10 on single trials, then
% average), matching fun_tfr_perm_session.m. See fun_powspec_session.m.
%
% Statistics:
% -----------------
% fooof has to be fit on averaged spectra, not single trials, so the SESSION
% is the unit of observation. Per session the hit and miss spectra are a
% matched pair, and the group test is a sign-flip permutation across sessions
% (exchanging the hit/miss label within a session), which makes this a
% random-effects test: it asks whether the difference replicates.
%
% Correction: max-statistic over the whole channel x frequency map for the
% periodic spectra, and over channels for offset and exponent.
%
% Outputs (per run_name, so runs never overwrite each other):
%   results_<animal>/<session>/POW/<run_name>/pow_spec.mat
%   results_<animal>/group_POW/<run_name>/pow_hitmiss_pooled.mat
%   Plots/Hitvsmiss/POWdiff/<run_name>/POWdiff_<animal>_<run_name>_*.{fig,pdf}

%% Specify paths

addpath /opt/fieldtrip_github/
ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/LFP_MUA_data/functions
clc

%% Define important variables

run_spec    = 1;               % 1 = submit the slurm jobs (one per session)
run_stats   = 1;               % 1 = group test across sessions
plotting    = 1;

min_trials  = 10;
min_sess_frac= 0.8;            % test an element only if this fraction of the
                               % pooled sessions contributes a value to it
alpha_level = 0.05;            % two-sided corrected alpha
nperm       = 5000;            % sign flips; 'all' is used when few sessions

run_name = 'prestim_2to80Hz';

switch run_name
    case 'prestim_2to80Hz'
        % the full pre-stimulus second, as in pow_freq.m.
        % mtmfft over 1 s gives 1 Hz resolution, so 2 Hz is usable here -
        % unlike in the TFR, where a 2 Hz sliding window does not fit
        toilim = [-1 0];
        foi    = 2:2:80;
    otherwise
        error('unknown run_name: %s', run_name)
end

fprintf('run "%s": %g to %g s, foi %g-%g Hz\n', run_name, toilim(1), toilim(2), foi(1), foi(end))

runinfo = struct('name',run_name,'toilim',toilim,'foi',foi, ...
    'nperm',nperm,'min_trials',min_trials,'alpha_level',alpha_level, ...
    'created',datestr(now));

plotfolder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/Hitvsmiss/POWdiff', run_name);
if ~isdir(plotfolder), mkdir(plotfolder), end

%% Animal configuration

animals(1).name            = 'klecks';
animals(1).resultsfolder   = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks';
animals(1).remove_sessions = [2 15 16 18];

animals(2).name            = 'hermes';
animals(2).resultsfolder   = '/mnt/hpc/projects/MWSampling/4Shivangi/results_hermes';
animals(2).remove_sessions = [1 3 21];

%% ========================= one animal at a time ======================== %%

for ianimal = 1:length(animals)

    animalName = animals(ianimal).name;
    datafolder = animals(ianimal).resultsfolder;
    fprintf('\n=============== %s ===============\n', animalName)

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
    fprintf('%d sessions found, %d kept\n', length(session_names), length(keep_sessions))

    session_paths = cellfun(@(x) fullfile(datafolder,x), session_names, 'uniform',0);
    output_paths  = cellfun(@(x) fullfile(datafolder,x,'POW',run_name), session_names, 'uniform',0);

    %% 1) Spectra on slurm, one job per session

    if run_spec
        cfg = cell(1,length(keep_sessions));
        for k = 1:length(keep_sessions)
            isess = keep_sessions(k);
            if ~isdir(output_paths{isess}), mkdir(output_paths{isess}), end
            cfg{k}.inputfile  = fullfile(session_paths{isess},'zlfptrials.mat');
            cfg{k}.outputfile = fullfile(output_paths{isess},'pow_spec.mat');
            cfg{k}.toilim     = toilim;
            cfg{k}.foi        = foi;
            cfg{k}.min_trials = min_trials;
            cfg{k}.run_name   = run_name;
        end
        slurmfun(@fun_powspec_session, cfg, ...
            'partition',   '8GB', ...
            'stopOnError', false, ...
            'useUserPath', true    );
    end

    %% 2) Group test across sessions

    if run_stats

        f = fullfile(datafolder,'critical_time','all_channels.mat');
        if exist(f,'file')
            A = load(f); chans = A.all_channels(:);
        else
            chans = {};
        end

        FH=[]; FM=[]; PH=[]; PM=[]; OH=[]; OM=[]; EH=[]; EM=[];
        sess_use = {}; freqs = [];

        for k = 1:length(keep_sessions)
            isess = keep_sessions(k);
            fpath = fullfile(output_paths{isess},'pow_spec.mat');
            if ~exist(fpath,'file')
                fprintf('  %s: no pow_spec.mat\n', session_names{isess}), continue
            end
            S = load(fpath);
            if isfield(S,'skipped') && S.skipped
                fprintf('  %s: skipped (%s)\n', session_names{isess}, S.reason), continue
            end
            if isempty(chans), chans = S.label(:); end
            if isempty(freqs)
                freqs = S.freq;
                nC = length(chans); nF = length(freqs);
                FH = nan(0,nC,nF); FM = FH; PH = FH; PM = FH;
                OH = nan(0,nC);    OM = OH; EH = OH; EM = OH;
            end
            if ~isequal(S.freq, freqs)
                fprintf('  %s: different frequency axis\n', session_names{isess}), continue
            end

            [tf, idx] = ismember(S.label, chans);
            ci = find(tf); co = idx(tf);

            fh = nan(1,nC,nF); fm = fh; ph = fh; pm = fh;
            oh = nan(1,nC); om = oh; eh = oh; em = oh;
            fh(1,co,:) = S.flat_hit(ci,:);   fm(1,co,:) = S.flat_miss(ci,:);
            ph(1,co,:) = S.pow_hit(ci,:);    pm(1,co,:) = S.pow_miss(ci,:);
            oh(1,co)   = S.off_hit(ci);      om(1,co)   = S.off_miss(ci);
            eh(1,co)   = S.exp_hit(ci);      em(1,co)   = S.exp_miss(ci);

            FH(end+1,:,:)=fh; FM(end+1,:,:)=fm; %#ok<SAGROW>
            PH(end+1,:,:)=ph; PM(end+1,:,:)=pm; %#ok<SAGROW>
            OH(end+1,:)=oh;   OM(end+1,:)=om;   %#ok<SAGROW>
            EH(end+1,:)=eh;   EM(end+1,:)=em;   %#ok<SAGROW>
            sess_use{end+1} = session_names{isess}; %#ok<SAGROW>
            fprintf('  %s: %d hits, %d misses\n', session_names{isess}, S.n_hit, S.n_miss)
        end

        nSess = size(FH,1);
        fprintf('\n[%s] group test on %d sessions\n', animalName, nSess)
        if nSess < 3, warning('too few sessions'), continue, end

        % ---- drop elements built from too few sessions. An SD estimated from
        % 2 sessions can collapse to ~0 and blow up t, and since the correction
        % takes the max |t|, one such element sets the threshold for the whole
        % family. Also removes all-NaN frequencies (the fooof edge at 80 Hz).
        min_sess = max(3, ceil(min_sess_frac*nSess));
        [dFlat, nFlat] = drop_sparse(FH - FM, min_sess);
        [dRaw , nRaw ] = drop_sparse(PH - PM, min_sess);
        [dOff , nOff ] = drop_sparse(OH - OM, min_sess);
        [dExp , nExp ] = drop_sparse(EH - EM, min_sess);
        fprintf('  elements needing >=%d of %d sessions; dropped: periodic %d, raw %d, offset %d, exponent %d\n', ...
            min_sess, nSess, nFlat, nRaw, nOff, nExp)

        % ---- the four tests, per channel
        [stat_flat, thr_flat] = signflip_test(dFlat, nperm, alpha_level);
        [stat_raw,  thr_raw ] = signflip_test(dRaw,  nperm, alpha_level);
        [stat_off,  thr_off ] = signflip_test(dOff,  nperm, alpha_level);
        [stat_exp,  thr_exp ] = signflip_test(dExp,  nperm, alpha_level);

        % ---- the same tests on the channel average, which is what the figures
        % plot: corrected over frequencies only
        [stat_flat_avg, thr_flat_avg] = signflip_test(squeeze(nanmean(dFlat,2)), nperm, alpha_level);
        [stat_raw_avg,  thr_raw_avg ] = signflip_test(squeeze(nanmean(dRaw ,2)), nperm, alpha_level);

        % offset and exponent averaged over channels: one number per session,
        % so nothing to correct over
        [stat_off_avg, thr_off_avg] = signflip_test(nanmean(dOff,2), nperm, alpha_level);
        [stat_exp_avg, thr_exp_avg] = signflip_test(nanmean(dExp,2), nperm, alpha_level);

        outdir = fullfile(datafolder,'group_POW',run_name);
        if ~isdir(outdir), mkdir(outdir), end
        save(fullfile(outdir,'pow_hitmiss_pooled.mat'), ...
            'stat_flat','stat_raw','stat_off','stat_exp', ...
            'thr_flat','thr_raw','thr_off','thr_exp', ...
            'stat_flat_avg','stat_raw_avg','thr_flat_avg','thr_raw_avg', ...
            'stat_off_avg','stat_exp_avg','thr_off_avg','thr_exp_avg', ...
            'FH','FM','PH','PM','OH','OM','EH','EM', ...
            'chans','freqs','sess_use','nSess','runinfo','-v7.3')

        %% 3) Report

        fprintf('\n--- %s, %d sessions, pre-stimulus %g to %g s ---\n', ...
            animalName, nSess, toilim(1), toilim(2))
        fprintf(' per channel:\n')
        report_one('aperiodic offset  ', stat_off, thr_off, alpha_level, [])
        report_one('aperiodic exponent', stat_exp, thr_exp, alpha_level, [])
        report_one('periodic power    ', stat_flat, thr_flat, alpha_level, freqs)
        report_one('raw log power     ', stat_raw,  thr_raw,  alpha_level, freqs)

        fprintf(' averaged over channels:\n')
        report_one('aperiodic offset  ', stat_off_avg, thr_off_avg, alpha_level, [])
        report_one('aperiodic exponent', stat_exp_avg, thr_exp_avg, alpha_level, [])
        report_one('periodic power    ', stat_flat_avg, thr_flat_avg, alpha_level, freqs)
        report_one('raw log power     ', stat_raw_avg,  thr_raw_avg,  alpha_level, freqs)

        %% 4) Plot

        if plotting
            fig = figure('Name',sprintf('POW hit-miss %s',animalName), ...
                'Units','centimeters','Position',[0 0 30 20]);

            % raw and periodic spectra, mean over channels and sessions.
            % Blue shading marks frequencies where the hit-miss difference of
            % the plotted (channel-averaged) quantity is significant, in the
            % same style as the ERP figures in erpdiff_lfp_*.m
            % (a) raw spectra. Sanity panel only - deliberately NOT shaded:
            % the raw test is dominated by the aperiodic rotation shown in (d),
            % so marking it here would report the same effect twice. The raw
            % test is still computed and printed, as a check that the
            % decomposition behaved.
            subplot(2,2,1), hold on
            hh = plot(freqs, squeeze(nanmean(nanmean(PH,2),1)), 'k', 'LineWidth',1.5);
            hm = plot(freqs, squeeze(nanmean(nanmean(PM,2),1)), 'r', 'LineWidth',1.5);
            xlabel('Frequency (Hz)'), ylabel('log_{10} power')
            legend([hh hm], {'hit','miss'}, 'Location','best')   % attach to the LINES,
            title('(a) raw spectra')                             % not the shading patches
            box on

            % (b) periodic spectra
            subplot(2,2,2), hold on
            hh = plot(freqs, squeeze(nanmean(nanmean(FH,2),1)), 'k', 'LineWidth',1.5);
            hm = plot(freqs, squeeze(nanmean(nanmean(FM,2),1)), 'r', 'LineWidth',1.5);
            shade_sig(freqs, stat_flat_avg.mask)
            xlabel('Frequency (Hz)'), ylabel('log_{10} periodic power')
            legend([hh hm], {'hit','miss'}, 'Location','best')
            title('(b) 1/f removed (fooof)'), box on

            % (c) periodic difference, with the peak annotated
            subplot(2,2,3), hold on
            d  = squeeze(nanmean(nanmean(FH-FM,2),1));
            se = squeeze(nanstd(nanmean(FH-FM,2),0,1))./sqrt(nSess);
            fill([freqs fliplr(freqs)],[d'+se' fliplr(d'-se')],[.8 .8 .8], ...
                'EdgeColor','none','FaceAlpha',.5)
            plot(freqs, d, 'k', 'LineWidth',1.5)
            yline(0,'k--')
            shade_sig(freqs, stat_flat_avg.mask)

            % label the largest significant deflection
            tv = stat_flat_avg.t(:); msk = logical(stat_flat_avg.mask(:));
            if any(msk)
                cand = find(msk);
                [~,k] = max(abs(d(cand))); ip = cand(k);
                plot(freqs(ip), d(ip), 'ro', 'MarkerSize',7, 'LineWidth',1.2)
                text(freqs(ip)+3, d(ip), sprintf(' %.1f Hz: %+.3f (%+.1f%%), t(%d) = %.2f', ...
                    freqs(ip), d(ip), 100*(10^d(ip)-1), nSess-1, tv(ip)), 'FontSize',8)
            end
            xlabel('Frequency (Hz)'), ylabel('hit - miss (log_{10})')
            title(sprintf('(c) periodic difference (blue: p<%g, corrected over %d freqs)', ...
                alpha_level, length(freqs)))
            box on

            % (d) aperiodic parameters, one point per channel
            subplot(2,2,4), hold on
            x = stat_off.d(:); y = stat_exp.d(:);
            ok = isfinite(x) & isfinite(y);
            plot(x(ok), y(ok), 'k.', 'MarkerSize',12)

            % regression line through the cloud; its slope gives the frequency
            % the 1/f component pivots about, since d(offset) = d(exp)*log10(f0)
            if sum(ok) > 2
                b  = polyfit(x(ok), y(ok), 1);
                rr = corrcoef(x(ok), y(ok)); r = rr(1,2);
                xf = linspace(min(x(ok)), max(x(ok)), 10);
                plot(xf, polyval(b,xf), 'b-', 'LineWidth',1)
                % Frequency the two aperiodic lines cross at. They differ by
                %   d_offset - d_exponent*log10(f),
                % so the crossing is at 10^(d_offset/d_exponent), using the
                % GROUP MEAN deltas. Not the slope of the scatter - that line
                % has an intercept, so its slope is not the ratio of the means.
                f0 = 10^( stat_off_avg.d / stat_exp_avg.d );
                if isfinite(f0) && f0 > 1 && f0 < 500
                    pivtxt = sprintf(', cross ~%.0f Hz', f0);
                else
                    pivtxt = '';
                end
            else
                r = NaN; pivtxt = '';
            end

            % keep zero inside the axes so the reader can see which side the
            % cloud sits on - otherwise the zero lines end up on the border
            xl = [min(0,min(x(ok))) max(0,max(x(ok)))]; px = 0.12*diff(xl) + eps;
            yl = [min(0,min(y(ok))) max(0,max(y(ok)))]; py = 0.12*diff(yl) + eps;
            xlim(xl + [-px px]); ylim(yl + [-py py])
            xline(0,'k--'), yline(0,'k--')

            text(0.03, 0.97, sprintf(['offset:   t(%d) = %.2f, p = %.3g (%d/%d ch)\n' ...
                'exponent: t(%d) = %.2f, p = %.3g (%d/%d ch)\nr = %.2f%s'], ...
                nSess-1, stat_off_avg.t, stat_off_avg.p, sum(stat_off.mask(:)), numel(stat_off.mask), ...
                nSess-1, stat_exp_avg.t, stat_exp_avg.p, sum(stat_exp.mask(:)), numel(stat_exp.mask), ...
                r, pivtxt), ...
                'Units','normalized', 'VerticalAlignment','top', 'FontSize',8)

            xlabel('\Delta aperiodic offset (hit - miss)')
            ylabel('\Delta aperiodic exponent (hit - miss)')
            title('(d) aperiodic parameters, one dot per channel'), box on

            sgtitle(sprintf('%s [%s]: pre-stimulus %g to %g s, %d sessions', ...
                animalName, run_name, toilim(1), toilim(2), nSess), 'Interpreter','none')
            fname = sprintf('POWdiff_%s_%s', animalName, run_name);
            savefig(fig, fullfile(plotfolder, [fname '.fig']))
            print(fig, fullfile(plotfolder, [fname '.pdf']), '-dpdf','-fillpage')
        end
    end
end

%% ============================== helpers ================================ %%

function [stat, thr] = signflip_test(D, nrand, alpha_level)
% Paired test across sessions with a sign-flip permutation.
% D is nSess x nChan (x nFreq) of hit-miss differences.
% The null exchanges the hit/miss label within a session, which for a paired
% difference is a sign flip. The max |t| over the whole map gives a
% family-wise threshold covering channels (and frequencies).

sz = size(D); n = sz(1);
Dm = reshape(D, n, []);
cnt = sum(~isnan(Dm),1);
tobs = nanmean(Dm,1) ./ (nanstd(Dm,0,1) ./ sqrt(cnt));

% random sign flips. With n sessions there are 2^n possible assignments;
% for n > ~15 that is far more than we would ever enumerate, so sample.
signs   = sign(randn(nrand, n));
maxdist = zeros(nrand,1);
for p = 1:nrand
    X  = Dm .* signs(p,:)';
    tp = nanmean(X,1) ./ (nanstd(X,0,1) ./ sqrt(cnt));
    maxdist(p) = max(abs(tp));
end

thr = quantile(maxdist, 1 - alpha_level);

pv = arrayfun(@(v) mean(maxdist >= abs(v)), tobs);
pv(~isfinite(tobs)) = NaN;      % dropped elements have no p value

stat        = [];
stat.t      = reshape(tobs, [sz(2:end) 1]);
stat.d      = reshape(nanmean(Dm,1), [sz(2:end) 1]);
stat.mask   = abs(stat.t) > thr;
stat.p      = reshape(pv, [sz(2:end) 1]);
stat.n      = n;
stat.maxdist= maxdist;
end


function [D, ndropped] = drop_sparse(D, min_sess)
% NaN out any element that fewer than min_sess sessions contribute to.
% D is nSess x ... ; the first dimension is sessions.

sz  = size(D);
Dm  = reshape(D, sz(1), []);
bad = sum(~isnan(Dm), 1) < min_sess;
ndropped = sum(bad);
Dm(:, bad) = NaN;
D = reshape(Dm, sz);
end


function shade_sig(x, mask)
% Blue patches over the contiguous runs where mask is true, spanning the full
% height of the axis. Same style as the significance shading in
% erpdiff_lfp_*.m, and it rescales with the axis rather than being pinned to
% hard-coded y limits.

mask = double(mask(:))';
if ~any(mask), return, end

onset  = find(conv(mask,[1 -1]) == 1);
offset = find(conv(mask,[1 -1]) == -1) - 1;
offset = offset(1:length(onset));

yl = ylim;
for i = 1:length(onset)
    x1 = x(onset(i)); x2 = x(offset(i));
    if x2 == x1                       % single sample: give it visible width
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


function report_one(name, stat, thr, alpha_level, freqs)
m = stat.mask;
fprintf('  %s: |t| threshold %.2f, %d/%d significant, mean diff %+.4f, min p %.4f\n', ...
    name, thr, sum(m(:)), numel(m), nanmean(stat.d(:)), min(stat.p(:)));
if any(m(:)) && ~isempty(freqs)
    if isvector(m)
        % channel-averaged: one value per frequency
        fi   = find(m);
        nchs = 1;
    else
        % per channel: rows are channels, columns are frequencies
        [ri, fi] = find(m);
        nchs = numel(unique(ri));
    end
    fprintf('      frequencies %g-%g Hz, %d channel(s)\n', ...
        min(freqs(fi)), max(freqs(fi)), nchs);
end
end
