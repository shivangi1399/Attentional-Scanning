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
foi          = 4:2:80;         % Hz
ncycles      = 3;              % sliding window length, in cycles
toi          = -0.6:0.02:0.4;  % s, relative to target onset
min_trials   = 10;             % skip a session with fewer hits or misses
alpha_level  = 0.05;           % two-sided corrected alpha

plotfolder   = '/mnt/hpc/projects/MWSampling/4Shivangi/Plots/Hitvsmiss/TFRdiff';
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
    output_paths  = cellfun(@(x) fullfile(datafolder,x,'TFR'), session_names, 'uniform',0);

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

            if isempty(chans), chans = S.label(:); end
            if isempty(real_sum)
                freqs    = S.freq;
                times    = S.time;
                nC       = length(chans);
                real_sum = zeros(nC, length(freqs), length(times), 'single');
                perm_sum = zeros(nperm, nC, length(freqs), length(times), 'single');
                count    = zeros(nC, 1);
            end
            if ~isequal(S.freq, freqs) || ~isequal(S.time, times)
                fprintf('  %s: different time-frequency grid, skipped\n', session_names{isess})
                continue
            end

            [tf, idx] = ismember(S.label, chans);
            chan_in   = find(tf);
            chan_out  = idx(tf);

            real_sum(chan_out,:,:)   = real_sum(chan_out,:,:)   + S.diff_real(chan_in,:,:);
            perm_sum(:,chan_out,:,:) = perm_sum(:,chan_out,:,:) + S.diff_perm(:,chan_in,:,:);
            count(chan_out)          = count(chan_out) + 1;
            nSess = nSess + 1;
            fprintf('  %s: %d hits, %d misses, %d channels\n', ...
                session_names{isess}, S.n_hit, S.n_miss, length(chan_in))
            clear S
        end

        if nSess < 2
            warning('%s: fewer than 2 sessions pooled - skipping', animalName)
            continue
        end

        % ---- session averages
        cnt = single(count);
        cnt(cnt == 0) = NaN;
        real_avg = real_sum ./ cnt;
        perm_avg = perm_sum ./ reshape(cnt,[1 numel(cnt)]);
        clear real_sum perm_sum

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
            p_corr(ichan,:,:) = reshape(mean(maxdist >= zc(:)', 1), size(zc));
        end

        %% 3) Report and save

        fprintf('\n--- %s: pooled max-pixel test ---\n', animalName)
        fprintf('  %d sessions, %d channels, %d permutations\n', nSess, length(chans), nperm)
        fprintf('  |z| threshold at p<%g corrected: %.2f\n', alpha_level, thr)
        fprintf('  %d / %d bins significant\n', sum(mask(:)), sum(~isnan(z_real(:))))
        if any(mask(:))
            [ci, fi, ti] = ind2sub(size(mask), find(mask));
            fprintf('  extent: %g-%g Hz, %g-%g s, %d channels\n', ...
                min(freqs(fi)), max(freqs(fi)), min(times(ti)), max(times(ti)), ...
                length(unique(ci)))
            fprintf('  smallest corrected p: %.4f\n', min(p_corr(:)))
        end

        outdir = fullfile(datafolder,'group_TFR');
        if ~isdir(outdir), mkdir(outdir), end
        save(fullfile(outdir,'tfr_hitmiss_pooled.mat'), ...
            'real_avg','z_real','mask','p_corr','thr','maxdist', ...
            'chans','freqs','times','count','nSess','nperm','-v7.3')

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
            sgtitle(sprintf('%s: hit - miss log power, %d sessions pooled (outline = p<%g, max-pixel)', ...
                animalName, nSess, alpha_level))
            savefig(fig, fullfile(plotfolder, sprintf('TFRdiff_%s_perchannel.fig',animalName)))
            print(fig, fullfile(plotfolder, sprintf('TFRdiff_%s_perchannel.pdf',animalName)), ...
                '-dpdf','-fillpage')

            % channel average, descriptive
            fig = figure('Name',sprintf('TFR hit-miss %s (chan avg)',animalName), ...
                'Units','centimeters','Position',[0 0 18 12]);
            D = squeeze(nanmean(real_avg,1));
            imagesc(times, freqs, D)
            set(gca,'YDir','normal')
            caxis([-1 1]*quantile(abs(D(~isnan(D))), 0.99))
            hold on, xline(0,'k--','LineWidth',1)
            colormap(bluewhitered_local)
            c = colorbar; c.Label.String = 'log_{10}(hit/miss)';
            xlabel('Time from target (s)'), ylabel('Frequency (Hz)')
            title(sprintf('%s: hit - miss, average over %d channels, %d sessions', ...
                animalName, nC, nSess))
            savefig(fig, fullfile(plotfolder, sprintf('TFRdiff_%s_chanavg.fig',animalName)))
            print(fig, fullfile(plotfolder, sprintf('TFRdiff_%s_chanavg.pdf',animalName)), ...
                '-dpdf','-fillpage')
        end

        clear real_avg z_real mask p_corr maxdist
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
