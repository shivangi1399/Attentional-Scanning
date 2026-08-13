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

        % ---- the four tests, per channel
        [stat_flat, thr_flat] = signflip_test(FH - FM, nperm, alpha_level);
        [stat_raw,  thr_raw ] = signflip_test(PH - PM, nperm, alpha_level);
        [stat_off,  thr_off ] = signflip_test(OH - OM, nperm, alpha_level);
        [stat_exp,  thr_exp ] = signflip_test(EH - EM, nperm, alpha_level);

        % ---- the same tests on the channel average.
        % The plotted curves are averaged over channels, so the shading in the
        % figures has to come from a test on that same quantity - correcting
        % over frequencies only, since there is no longer a channel dimension.
        [stat_flat_avg, thr_flat_avg] = signflip_test(squeeze(nanmean(FH-FM,2)), nperm, alpha_level);
        [stat_raw_avg,  thr_raw_avg ] = signflip_test(squeeze(nanmean(PH-PM,2)), nperm, alpha_level);

        % offset and exponent averaged over channels: one number per session,
        % so this is a plain paired test with nothing to correct over
        [stat_off_avg, thr_off_avg] = signflip_test(nanmean(OH-OM,2), nperm, alpha_level);
        [stat_exp_avg, thr_exp_avg] = signflip_test(nanmean(EH-EM,2), nperm, alpha_level);

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
            subplot(2,2,1), hold on
            plot(freqs, squeeze(nanmean(nanmean(PH,2),1)), 'k', 'LineWidth',1.5)
            plot(freqs, squeeze(nanmean(nanmean(PM,2),1)), 'r', 'LineWidth',1.5)
            shade_sig(freqs, stat_raw_avg.mask)
            xlabel('Frequency (Hz)'), ylabel('log_{10} power')
            legend({'hit','miss'},'Location','best'), title('raw spectra'), box on

            subplot(2,2,2), hold on
            plot(freqs, squeeze(nanmean(nanmean(FH,2),1)), 'k', 'LineWidth',1.5)
            plot(freqs, squeeze(nanmean(nanmean(FM,2),1)), 'r', 'LineWidth',1.5)
            shade_sig(freqs, stat_flat_avg.mask)
            xlabel('Frequency (Hz)'), ylabel('log_{10} periodic power')
            title('1/f removed (fooof)'), box on

            % difference with the corrected threshold
            subplot(2,2,3), hold on
            d  = squeeze(nanmean(nanmean(FH-FM,2),1));
            se = squeeze(nanstd(nanmean(FH-FM,2),0,1))./sqrt(nSess);
            fill([freqs fliplr(freqs)],[d'+se' fliplr(d'-se')],[.8 .8 .8], ...
                'EdgeColor','none','FaceAlpha',.5)
            plot(freqs, d, 'k', 'LineWidth',1.5)
            yline(0,'k--')
            shade_sig(freqs, stat_flat_avg.mask)
            xlabel('Frequency (Hz)'), ylabel('hit - miss (log_{10})')
            title('periodic difference (blue = p<0.05 corrected)'), box on

            % offset and exponent per channel
            subplot(2,2,4), hold on
            plot(stat_off.d, stat_exp.d, 'k.', 'MarkerSize',12)
            xline(0,'k--'), yline(0,'k--')
            xlabel('\Delta aperiodic offset (hit - miss)')
            ylabel('\Delta aperiodic exponent (hit - miss)')
            title('aperiodic parameters, one dot per channel'), box on

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

stat        = [];
stat.t      = reshape(tobs, [sz(2:end) 1]);
stat.d      = reshape(nanmean(Dm,1), [sz(2:end) 1]);
stat.mask   = abs(stat.t) > thr;
stat.p      = reshape(arrayfun(@(v) mean(maxdist >= abs(v)), tobs), [sz(2:end) 1]);
stat.n      = n;
stat.maxdist= maxdist;
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
