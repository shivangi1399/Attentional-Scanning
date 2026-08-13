clear all
close all
clc

% Code description:
% -----------------
% Hit vs miss time-frequency contrast of the LFP, run separately for klecks
% and hermes.
%
% Data organisation follows pow_freq.m / erpdiff_lfp_*.m:
%   results_<animal>/<session>/zlfptrials.mat        -> z-scored LFP
%   results_<animal>/critical_time/all_channels.mat  -> channel frame
%   trialinfo(:,20): 1 = hit, 5 = miss
%
% Everything is written under a folder named after run_name (see below), so
% several parameter choices can coexist:
%   results_<animal>/<session>/TFR/<run_name>/tfr_perm.mat
%   results_<animal>/group_TFR/<run_name>/tfr_hitmiss_pooled.mat
%   Plots/Hitvsmiss/TFRdiff/<run_name>/TFRdiff_<animal>_<run_name>_*.{fig,pdf}
% The per-session files carry the settings that made them in 'settings', and
% the pooled file in 'runinfo'.
%
% Statistics:
% -----------------
% Trial is the unit of observation, sessions are pooled.
% Per session (one slurm job each) the single-trial TFR is computed once, and
% the hit-miss difference map is computed for the real labels and for nperm
% relabelings. Relabelling happens within a session only, keeping that
% session's hit and miss counts.
% The permuted maps are then averaged over sessions exactly like the real
% ones, so the null distribution describes the pooled statistic.
% Correction: max-pixel (max-statistic) over the whole channel x frequency x
% time map, on maps normalised bin by bin against their own permutation
% distribution.

%% Specify paths

addpath /opt/fieldtrip_github/
ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/code/LFP_MUA_data/functions
clc

%% Define important variables

run_perm     = 1;              % 1 = submit the slurm jobs (one per session)
run_stats    = 1;              % 1 = pool the sessions and run the max-pixel test
plotting     = 1;

nperm        = 1000;           % within-session relabelings
min_trials   = 10;             % skip a session with fewer hits or misses
min_sess_frac= 0.8;            % a time-frequency bin is tested only if at least
                               % this fraction of the pooled sessions has a
                               % valid (non-NaN) value there
alpha_level  = 0.05;           % two-sided corrected alpha

%% Analysis run
% foi and toi are constrained by the epoch. The epochs run from -1 to 0.4 s,
% and the sliding window at the LOWEST frequency is ncycles/foi(1) seconds
% long, so a time point is only usable if
%     toi +- ncycles/foi(1)/2   stays inside [-1 0.4].
% Low frequencies therefore cost post-target time, and there is no single
% setting that covers both ends. Each named run below is one choice, and
% everything it produces - per session files, pooled results, figures - is
% written into a folder of that name, so runs never overwrite each other.

run_name = 'unified_4to80Hz_2cyc';

switch run_name

    case 'pretarget_6to80Hz'
        % low frequencies, long pre-target baseline.
        % 3 cycles at 6 Hz = 0.5 s window -> toi can reach 0.15 s
        foi     = 6:2:80;
        ncycles = 3;
        toi     = -0.5:0.02:0.15;

    case 'poststim_10to80Hz'
        % post-target response.
        % 3 cycles at 10 Hz = 0.3 s window -> toi can reach 0.25 s
        foi     = 10:2:80;
        ncycles = 3;
        toi     = -0.4:0.02:0.24;

    case 'unified_4to80Hz_2cyc'
        % 2 cycles instead of 3: shorter windows, so low frequencies and
        % post-target time fit in ONE map.
        % Cost is frequency smearing - resolution is roughly foi/ncycles, so a
        % 20 Hz estimate integrates ~10 Hz rather than ~7 Hz.
        % 6 Hz and above are valid across the whole toi. The 4 Hz row (0.5 s
        % window) goes NaN after 0.15 s and drops out of the test there, which
        % is fine - NaN bins are handled per bin.
        % 2 Hz is deliberately NOT included: a 2-cycle window at 2 Hz is 1 s
        % wide, so it would be NaN after -0.1 s and every surviving time point
        % would be built from nearly the same samples.
        foi     = 4:2:80;
        ncycles = 2;
        toi     = -0.5:0.02:0.22;

    otherwise
        error('unknown run_name: %s', run_name)
end

% check the run against the epoch before anything is submitted
maxwin = max(ncycles ./ foi);
fprintf('run "%s": foi %g-%g Hz, %g cycles, toi %g to %g s, longest window %g s\n', ...
    run_name, foi(1), foi(end), ncycles, toi(1), toi(end), maxwin)
fprintf('  needs data from %g to %g s\n', toi(1)-maxwin/2, toi(end)+maxwin/2)

runinfo = struct('name',run_name,'foi',foi,'ncycles',ncycles,'toi',toi, ...
    'nperm',nperm,'min_trials',min_trials,'min_sess_frac',min_sess_frac, ...
    'alpha_level',alpha_level,'created',datestr(now));

plotfolder = fullfile('/mnt/hpc/projects/MWSampling/4Shivangi/Plots/Hitvsmiss/TFRdiff', run_name);
if ~isdir(plotfolder), mkdir(plotfolder), end

%% Animal configuration
% remove_sessions are the session indices excluded in erpdiff_lfp_*.m - check
% that these are still the ones you want to drop

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

    %% Find sessions

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
    output_paths  = cellfun(@(x) fullfile(datafolder,x,'TFR',run_name), session_names, 'uniform',0);

    %% 1) Permutations on slurm, one job per session

    if run_perm

        cfg = cell(1,length(keep_sessions));
        for k = 1:length(keep_sessions)
            isess = keep_sessions(k);
            if ~isdir(output_paths{isess}), mkdir(output_paths{isess}), end

            cfg{k}.inputfile  = fullfile(session_paths{isess},'zlfptrials.mat');
            cfg{k}.outputfile = fullfile(output_paths{isess},'tfr_perm.mat');
            cfg{k}.foi        = foi;
            cfg{k}.ncycles    = ncycles;
            cfg{k}.toi        = toi;
            cfg{k}.nperm      = nperm;
            cfg{k}.min_trials = min_trials;
            cfg{k}.seed       = isess;    % different relabelings per session
            cfg{k}.run_name   = run_name; % stored in the output for traceability
        end

        slurmfun(@fun_tfr_perm_session, cfg, ...
            'partition',   '8GB', ...
            'stopOnError', false, ...
            'useUserPath', true    );
    end

    %% 2) Pool the sessions and run the max-pixel test

    if run_stats

        % channel frame, so that sessions with different surviving channels
        % can still be pooled
        f = fullfile(datafolder,'critical_time','all_channels.mat');
        if exist(f,'file')
            A = load(f);
            chans = A.all_channels(:);
        else
            chans = {};
        end

        real_sum = []; perm_sum = []; count = []; nSess = 0;

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

            if isempty(chans), chans = S.label(:); end
            if isempty(real_sum)
                freqs    = S.freq;
                times    = S.time;
                nC       = length(chans);
                real_sum = zeros(nC, length(freqs), length(times), 'single');
                perm_sum = zeros(nperm, nC, length(freqs), length(times), 'single');
                count    = zeros(nC, length(freqs), length(times), 'single');
            end
            if ~isequal(S.freq, freqs) || ~isequal(S.time, times)
                fprintf('  %s: different time-frequency grid, skipped\n', session_names{isess})
                continue
            end

            [tf, idx] = ismember(S.label, chans);
            chan_in   = find(tf);
            chan_out  = idx(tf);

            % Accumulate bin by bin, treating NaN as "this session has nothing
            % to say about this bin" rather than poisoning it for every
            % session. The NaN pattern is identical in the real and the
            % permuted maps, so both stay consistent with the same count.
            d     = S.diff_real(chan_in,:,:);
            valid = ~isnan(d);
            d(~valid) = 0;
            dp = S.diff_perm(:,chan_in,:,:);
            dp(isnan(dp)) = 0;

            real_sum(chan_out,:,:)   = real_sum(chan_out,:,:)   + d;
            perm_sum(:,chan_out,:,:) = perm_sum(:,chan_out,:,:) + dp;
            count(chan_out,:,:)      = count(chan_out,:,:)      + single(valid);
            clear d dp valid
            nSess = nSess + 1;
            fprintf('  %s: %d hits, %d misses, %d channels, %.1f%% NaN bins\n', ...
                session_names{isess}, S.n_hit, S.n_miss, length(chan_in), 100*S.nan_frac)
            clear S
        end

        if nSess < 2
            warning('%s: fewer than 2 sessions pooled - skipping', animalName)
            continue
        end

        % ---- session averages, per bin.
        % A bin is only tested if enough sessions contributed to it, so that
        % the map is not a patchwork of averages over wildly different numbers
        % of sessions. Everything below the cutoff becomes NaN and drops out.
        cnt = count;
        cnt(cnt < max(1, min_sess_frac*nSess)) = NaN;
        fprintf('  %d / %d bins have at least %d of %d sessions\n', ...
            sum(~isnan(cnt(:))), numel(cnt), ceil(max(1,min_sess_frac*nSess)), nSess)

        real_avg = real_sum ./ cnt;
        perm_avg = perm_sum ./ reshape(cnt,[1 size(cnt)]);
        clear real_sum perm_sum

        % ---- channel-averaged map, tested in its own right.
        % The channel-average figure plots the mean over channels, so its
        % significance has to come from a test on that same quantity: average
        % the permuted maps over channels FIRST, then normalise and correct
        % over frequency x time only. Marking it with the per-channel result
        % would flag bins that do not correspond to what is drawn. Averaging
        % before testing is also more sensitive to an effect shared across the
        % array, for the same reason a pooled test beats correcting over every
        % channel separately.
        real_chan = squeeze(nanmean(real_avg, 1));            % freq x time
        perm_chan = squeeze(nanmean(perm_avg, 2));            % nperm x freq x time
        mu_c = nanmean(perm_chan,1); sd_c = nanstd(perm_chan,0,1);
        sd_c(sd_c == 0) = NaN;
        z_chan       = (real_chan - squeeze(mu_c)) ./ squeeze(sd_c);
        perm_chan    = (perm_chan - mu_c) ./ sd_c;
        maxdist_chan = max(abs(reshape(perm_chan, nperm, [])), [], 2);
        thr_chan     = quantile(maxdist_chan, 1 - alpha_level);
        mask_chan    = abs(z_chan) > thr_chan;
        p_chan       = reshape(mean(maxdist_chan >= abs(z_chan(:))', 1), size(z_chan));
        p_chan(isnan(z_chan)) = NaN;
        clear perm_chan

        % ---- normalise each bin against its own permutation distribution,
        % otherwise the max is set by whichever bins are noisiest
        mu = nanmean(perm_avg, 1);
        sd = nanstd(perm_avg, 0, 1);
        sd(sd == 0) = NaN;
        z_real   = (real_avg - squeeze(mu)) ./ squeeze(sd);
        perm_avg = (perm_avg - mu) ./ sd;

        % ---- max-pixel correction over the whole channel x freq x time map
        maxdist = max(abs(reshape(perm_avg, nperm, [])), [], 2);
        thr     = quantile(maxdist, 1 - alpha_level);
        mask    = abs(z_real) > thr;
        clear perm_avg

        % ---- corrected p value per bin
        p_corr = nan(size(z_real));
        for ichan = 1:size(z_real,1)
            zc = abs(squeeze(z_real(ichan,:,:)));
            pc = reshape(mean(maxdist >= zc(:)', 1), size(zc));
            pc(isnan(zc)) = NaN;          % unusable bins have no p value
            p_corr(ichan,:,:) = pc;
        end

        %% 3) Report and save

        fprintf('\n--- %s: pooled max-pixel test ---\n', animalName)
        fprintf('  %d sessions, %d channels, %d permutations\n', nSess, length(chans), nperm)
        fprintf('  |z| threshold at p<%g corrected: %.2f\n', alpha_level, thr)
        fprintf('  %d / %d usable bins significant\n', sum(mask(:)), sum(~isnan(z_real(:))))
        fprintf('  smallest corrected p: %.4f\n', min(p_corr(~isnan(p_corr))))

        fprintf('  -- average over channels --\n')
        fprintf('  |z| threshold at p<%g corrected: %.2f\n', alpha_level, thr_chan)
        fprintf('  %d / %d usable bins significant, smallest corrected p: %.4f\n', ...
            sum(mask_chan(:)), sum(~isnan(z_chan(:))), min(p_chan(~isnan(p_chan))))
        if any(mask_chan(:))
            [fi_c, ti_c] = find(mask_chan);
            fprintf('  extent: %g-%g Hz, %g-%g s\n', ...
                min(freqs(fi_c)), max(freqs(fi_c)), min(times(ti_c)), max(times(ti_c)))
        end
        if any(mask(:))
            [ci, fi, ti] = ind2sub(size(mask), find(mask));
            fprintf('  extent: %g-%g Hz, %g-%g s, %d channels\n', ...
                min(freqs(fi)), max(freqs(fi)), min(times(ti)), max(times(ti)), ...
                length(unique(ci)))
        end

        outdir = fullfile(datafolder,'group_TFR',run_name);
        if ~isdir(outdir), mkdir(outdir), end
        save(fullfile(outdir,'tfr_hitmiss_pooled.mat'), ...
            'real_avg','z_real','mask','p_corr','thr','maxdist', ...
            'real_chan','z_chan','mask_chan','p_chan','thr_chan','maxdist_chan', ...
            'chans','freqs','times','count','nSess','nperm','runinfo','-v7.3')
        fprintf('  saved to %s\n', outdir)

        %% 4) Plot

        if plotting

            nC   = length(chans);
            nrow = ceil(sqrt(nC));
            clim = quantile(abs(real_avg(~isnan(real_avg))), 0.99);

            % per channel
            fig = figure('Name',sprintf('TFR hit-miss %s',animalName), ...
                'Units','normalized','Position',[0 0 1 1]);
            for ichan = 1:nC
                subplot(nrow,nrow,ichan)
                imagesc(times, freqs, squeeze(real_avg(ichan,:,:)))
                set(gca,'YDir','normal'), caxis([-clim clim]), hold on
                m = squeeze(mask(ichan,:,:));
                if any(m(:))
                    contour(times, freqs, double(m), [0.5 0.5], 'k', 'LineWidth', 1)
                end
                xline(0,'k--')
                title(chans{ichan}, 'Interpreter','none', 'FontSize',6)
                set(gca,'FontSize',6)
            end
            colormap(bluewhitered_local)
            c = colorbar('south'); c.Label.String = 'log_{10}(hit/miss)';
            sgtitle(sprintf('%s [%s]: hit - miss log power, %d sessions pooled (outline = p<%g, max-pixel)', ...
                animalName, run_name, nSess, alpha_level), 'Interpreter','none')
            fname = sprintf('TFRdiff_%s_%s_perchannel', animalName, run_name);
            savefig(fig, fullfile(plotfolder, [fname '.fig']))
            print(fig, fullfile(plotfolder, [fname '.pdf']), '-dpdf','-fillpage')

            % channel average, with its own corrected significance
            fig = figure('Name',sprintf('TFR hit-miss %s (chan avg)',animalName), ...
                'Units','centimeters','Position',[0 0 18 12]);
            imagesc(times, freqs, real_chan)
            set(gca,'YDir','normal')
            caxis([-1 1]*quantile(abs(real_chan(~isnan(real_chan))), 0.99))
            hold on
            if any(mask_chan(:))
                contour(times, freqs, double(mask_chan), [0.5 0.5], 'k', 'LineWidth', 1.5)
            end
            xline(0,'k--','LineWidth',1)
            colormap(bluewhitered_local)
            c = colorbar; c.Label.String = 'log_{10}(hit/miss)';
            xlabel('Time from target (s)'), ylabel('Frequency (Hz)')
            title(sprintf('%s [%s]: hit - miss, average over %d channels, %d sessions (outline = p<%g)', ...
                animalName, run_name, nC, nSess, alpha_level), 'Interpreter','none')
            fname = sprintf('TFRdiff_%s_%s_chanavg', animalName, run_name);
            savefig(fig, fullfile(plotfolder, [fname '.fig']))
            print(fig, fullfile(plotfolder, [fname '.pdf']), '-dpdf','-fillpage')
        end

        clear real_avg z_real mask p_corr maxdist
        clear real_chan z_chan mask_chan p_chan maxdist_chan
    end
end

%% ============================== helper ================================= %%

function cmap = bluewhitered_local
% symmetric blue-white-red colormap, so that zero is white

n = 128;
b = [linspace(0,1,n)' linspace(0,1,n)' ones(n,1)];
r = [ones(n,1) linspace(1,0,n)' linspace(1,0,n)'];
cmap = [b; r];
end
