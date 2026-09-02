% =====================================================================
% CHECK (not a main analysis): could a travelling wave be hidden by trial
% averaging?
%
% WHY THIS EXISTS
%   cortical_planar_wave_PGD.m and cortical_planar_wave_derotation.m both run
%   on pref_phase, which is a circular mean OVER TRIALS. If a wave existed on
%   single trials but pointed in a different direction on each one, those ramps
%   would cancel in the average and both analyses would report nothing.
%   THIS SCRIPT TESTS IT.
%
% THE KEY IDEA
%   You do not have to track where each wave went. A wave at speed v has
%   k = 2*pi*f/v on EVERY trial, whatever its direction, and k is a MAGNITUDE.
%   So the frequency-scaling test survives direction variability completely:
%       real wave at one speed  ->  log-log slope of k vs f = +1
%       fixed spatial offset    ->  log-log slope of k vs f =  0
%   Direction consistency is reported separately, as a description of how much
%   the trial average was actually losing.
%
% THE TRAP THIS SCRIPT AVOIDS
%   A plane fitted to ~60 noisy single-trial phases returns a non-zero k by
%   chance, and that noise floor is set by phase scatter and array geometry,
%   NOT by frequency -- so it is roughly flat in f. A floor-dominated fit
%   therefore produces slope ~0 automatically and would "confirm" the fixed
%   offset while measuring nothing. Every k below is reported as an EXCESS
%   over a per-frequency spatial-shuffle null, and the raw and excess slopes
%   are printed side by side so floor domination is visible.
%
% WHAT THIS IS NOT
%   Fits the raw single-trial phase across electrodes, NOT pref_phase. It is a
%   control on the averaging question, not a re-run of the main analysis.
%   Numbers are not expected to equal those from cortical_planar_wave_PGD.m.
%
% INPUT   results_<animal>/<session>/Phase_analysis/hit_miss/
%             100iter_cut@cp_m10/phase/<chan_folder>/phase.mat
%         (the 'hit_miss' folder holds the LFP phase estimate itself; it is
%          not DV-specific)
%
% OUTPUT  Plots/scanning/checks/check_single_trial_planes.pdf
%         results_combined/scanning/checks/check_single_trial_planes.mat
% =====================================================================

clearvars; close all; clc

%% --- Settings -------------------------------------------------------
base        = '/mnt/hpc/projects/MWSampling/4Shivangi';
animals     = {'klecks','hermes'};
SPACING_MM  = 0.4;
GRID        = 8;
K_MAX       = 2.5;          % rad/mm, top of the swept wavenumber grid
N_K         = 26;           % wavenumber steps
N_THETA     = 24;           % direction steps (15 deg each)
nPerm       = 1000;          % channel-position shuffles per frequency
MAX_SESS    = Inf;          % set to e.g. 3 for a quick run
MIN_TRIALS  = 40;
rng(0);

out_dir = fullfile(base,'Plots','scanning','checks');
res_dir = fullfile(base,'results_combined','scanning','checks');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
if ~exist(res_dir,'dir'), mkdir(res_dir); end

ch_col = ceil((1:GRID*GRID)'/GRID);
ch_row = GRID - mod((1:GRID*GRID)'-1, GRID);

R = struct();

%% --- Per animal -----------------------------------------------------
for ia = 1:numel(animals)
    animalName = animals{ia};
    sessions = list_sessions(base, animalName);
    if isempty(sessions)
        warning('No sessions found for %s -- skipping.', animalName); continue
    end
    nS = min(numel(sessions), MAX_SESS);
    fprintf('\n##### %s : %d sessions #####\n', animalName, nS);

    per_sess = {}; freq = []; nTot = 0; sess_used = {};

    for is = 1:nS
        [Z, freq, present] = load_session(base, animalName, sessions{is}, GRID);
        if isempty(Z), fprintf('  ! skip %s (unreadable)\n', sessions{is}); continue, end

        A = Z(:,:,present);
        good = all(all(isfinite(A),2),3);
        A = A(good,:,:);
        if size(A,1) < MIN_TRIALS
            fprintf('  . skip %s (only %d clean trials)\n', sessions{is}, size(A,1)); continue
        end
        nTot = nTot + size(A,1);

        x = (ch_col(present)-1)*SPACING_MM;
        y = (ch_row(present)-1)*SPACING_MM;
        [E, kmag, kdir] = build_grid(x, y, K_MAX, N_K, N_THETA);   % E: grid x chan

        rec = nan(numel(freq),3);          % [median k, median k_null, dir concentration]
        for f = 1:numel(freq)
            Mf = squeeze(A(:,f,:));                       % trials x chan
            sw = max(sum(abs(Mf),2), eps);
            Rg = abs(Mf * E.') ./ sw;                      % trials x grid
            [~, b] = max(Rg, [], 2);
            kn = nan(nPerm,1);
            for p = 1:nPerm
                pm = randperm(size(Mf,2));
                Rp = abs(Mf(:,pm) * E.') ./ sw;
                [~, bp] = max(Rp, [], 2);
                kn(p) = median(kmag(bp));
            end
            rec(f,:) = [median(kmag(b)), mean(kn), abs(mean(exp(1i*kdir(b))))];
        end
        per_sess{end+1} = rec; %#ok<AGROW>
        sess_used{end+1} = sessions{is}; %#ok<AGROW>
        exc = max(rec(:,1)-rec(:,2), 1e-6);
        fprintf('  . %-45s %5d trials   slope(excess) %+.2f   dirR %.2f\n', ...
            sessions{is}, size(A,1), loglog_slope(freq, exc), median(rec(:,3)));
    end

    if isempty(per_sess), continue, end

    M   = mean(cat(3, per_sess{:}), 3);          % freq x 3, mean over sessions
    exc = max(M(:,1)-M(:,2), 1e-6);
    sl_raw = loglog_slope(freq, M(:,1));
    sl_exc = loglog_slope(freq, exc);

    fprintf('\n  pooled over %d sessions (%d trials)\n', numel(per_sess), nTot);
    fprintf('    median k                     %.3f rad/mm\n', median(M(:,1)));
    fprintf('    median shuffled-null k       %.3f rad/mm   (signal/floor %.2f)\n', ...
        median(M(:,2)), median(M(:,1))/max(median(M(:,2)),eps));
    fprintf('    slope k~f, raw               %+.2f\n', sl_raw);
    fprintf('    slope k~f, EXCESS over null  %+.2f    <-- the number to read\n', sl_exc);
    fprintf('    freqs clamped at floor       %d of %d   (excess<=0; these bias the slope)\n', sum(M(:,1)-M(:,2)<=0), numel(freq));
    fprintf('    direction concentration      %.2f     (0 = random per trial, 1 = identical)\n', ...
        median(M(:,3)));
    if sl_exc > 0.6
        verdict = 'consistent with a travelling wave (slope near +1)';
    else
        verdict = 'NO wave even on single trials: k does not scale with f (slope far below +1)';
    end
    fprintf('  --> %s\n', verdict);

    R.(animalName) = struct('freq',freq,'M',M,'per_sess',{per_sess},'sessions',{sess_used}, ...
        'slope_raw',sl_raw,'slope_excess',sl_exc,'nTrials',nTot,'verdict',verdict);
end

%% --- Figure ---------------------------------------------------------
% Row 1: k vs frequency, observed / null / excess, with the two model lines.
%        A wave must climb parallel to the "+1" reference line.
% Row 2: how consistent the wave direction is across trials -- i.e. how much
%        the trial average was throwing away.
an = fieldnames(R);
if ~isempty(an)
    fg = new_fig(430*numel(an), 640);
    for k = 1:numel(an)
        A = R.(an{k}); fr = A.freq(:); M = A.M;
        exc = max(M(:,1)-M(:,2), 1e-6);

        ax = subplot(2, numel(an), k); hold(ax,'on');
        plot(ax, fr, M(:,1),  '-',  'Color',[0.15 0.35 0.75], 'LineWidth',1.6, 'DisplayName','observed k');
        plot(ax, fr, M(:,2),  ':',  'Color',[0.45 0.45 0.45], 'LineWidth',1.3, 'DisplayName','shuffled null');
        plot(ax, fr, exc,     '-',  'Color',[0.85 0.35 0.10], 'LineWidth',2.0, 'DisplayName','EXCESS (read this)');
        % model reference lines, anchored at the low-frequency excess
        ref = exc(1);
        plot(ax, fr, ref*(fr/fr(1)), '--', 'Color',[0.10 0.55 0.35], 'LineWidth',1.2, ...
             'DisplayName','wave: slope +1');
        plot(ax, fr, ref*ones(size(fr)), '--', 'Color',[0.85 0.55 0.10], 'LineWidth',1.2, ...
             'DisplayName','fixed offset: slope 0');
        set(ax,'XScale','log','YScale','log'); grid(ax,'on');
        xlabel(ax,'frequency (Hz)'); ylabel(ax,'wavenumber k (rad/mm)');
        title(ax, sprintf('%s -- excess slope %+.2f (%d trials)', an{k}, A.slope_excess, A.nTrials), ...
              'FontSize', 9);
        legend(ax,'Location','best','Box','off','FontSize',7);

        ax2 = subplot(2, numel(an), numel(an)+k); hold(ax2,'on');
        plot(ax2, fr, M(:,3), '-', 'Color',[0.35 0.20 0.55], 'LineWidth',1.8);
        yline(1, 'k--', 'identical direction every trial');
        set(ax2,'XScale','log'); ylim(ax2,[0 1]); grid(ax2,'on');
        xlabel(ax2,'frequency (Hz)'); ylabel(ax2,'direction concentration across trials');
        title(ax2, {'low = direction varies trial to trial', ...
                    'i.e. how much the trial average discards'}, 'FontSize', 8);
    end
    sgtitle({'CHECK: could a travelling wave be hidden by trial averaging?', ...
             'k is a magnitude, so the frequency-scaling test survives direction variability.', ...
             'Read the EXCESS curve against the two dashed model lines: +1 = wave, 0 = fixed offset.'}, ...
             'FontSize', 10);
    save_pdf(fg, fullfile(out_dir,'check_single_trial_planes.pdf'));
end

save(fullfile(res_dir,'check_single_trial_planes.mat'), 'R', 'K_MAX','N_K','N_THETA','nPerm', '-v7.3');
fprintf('\nSaved figure to %s\n', fullfile(out_dir,'check_single_trial_planes.pdf'));

%% ===================== local functions ==============================
function s = list_sessions(base, animalName)
d = dir(fullfile(base, ['results_' animalName], [animalName '_*']));
d = d([d.isdir]);
s = {};
for i = 1:numel(d)
    if exist(fullfile(d(i).folder, d(i).name, ...
             'Phase_analysis','hit_miss','100iter_cut@cp_m10','phase'), 'dir')
        s{end+1} = d(i).name; %#ok<AGROW>
    end
end
s = sort(s);
end

function [Z, freq, present] = load_session(base, animalName, sess, GRID)
% -> Z (trials x freq x 64 complex, NaN where a channel is absent)
% The i-th SORTED channel folder corresponds to the i-th entry of
% phase.label, and canonical channel = label number - 64. (Verified against
% channels_sessions.mat.)
Z = []; freq = []; present = [];
pdir = fullfile(base, ['results_' animalName], sess, ...
                'Phase_analysis','hit_miss','100iter_cut@cp_m10','phase');
d = dir(pdir); d = d([d.isdir]);
nm = {d.name}; nm = nm(~ismember(nm,{'.','..'}));
num = str2double(nm); keep = ~isnan(num);
[~,ord] = sort(num(keep)); nm = nm(keep); nm = nm(ord);
if isempty(nm), return, end

S = load(fullfile(pdir, nm{1}, 'phase.mat'));
lab = S.phase.label;
enum = nan(numel(lab),1);
for i = 1:numel(lab)
    t = regexp(lab{i}, '-\s*(\d+)', 'tokens', 'once');
    enum(i) = str2double(t{1});
end
% Electrode numbering differs BETWEEN ANIMALS: klecks labels run V4-65..128,
% hermes V4-1..64. Detect the offset instead of assuming one; both map onto
% canonical slots 1..64 (verified against channels_sessions.mat).
if max(enum) > 64, offset = 64; else, offset = 0; end
canon = enum - offset;
assert(all(canon >= 1 & canon <= GRID*GRID), ...
    'Channel labels in %s map outside 1..%d', pdir, GRID*GRID);
freq = S.phase.freq(:);
[nT, nF] = size(S.phase.ar_phase);
Z = nan(nT, nF, GRID*GRID);
for i = 1:numel(nm)
    P = load(fullfile(pdir, nm{i}, 'phase.mat'));
    Z(:,:,canon(i)) = P.phase.amp_vec .* exp(1i*P.phase.ar_phase);
end
present = sort(canon);
end

function [E, kmag, kdir] = build_grid(x, y, K_MAX, N_K, N_THETA)
% de-rotation design matrix over (wavenumber x direction), plus the k=0 cell
ks = linspace(K_MAX/N_K, K_MAX, N_K);
th = linspace(0, 2*pi, N_THETA+1); th(end) = [];
kx = 0; ky = 0; kmag = 0; kdir = 0;
for i = 1:numel(ks)
    kx   = [kx,   ks(i)*cos(th)]; %#ok<AGROW>
    ky   = [ky,   ks(i)*sin(th)]; %#ok<AGROW>
    kmag = [kmag, repmat(ks(i),1,N_THETA)]; %#ok<AGROW>
    kdir = [kdir, th]; %#ok<AGROW>
end
E = exp(-1i*(kx(:)*x(:).' + ky(:)*y(:).'));      % grid x chan
kmag = kmag(:); kdir = kdir(:);
end

function s = loglog_slope(freq, v)
p = polyfit(log(freq(:)), log(max(v(:), 1e-6)), 1);
s = p(1);
end

function f = new_fig(w, h)
ss = get(0,'ScreenSize');
w = min(w, ss(3)-80); h = min(h, ss(4)-120);
f = figure('Units','pixels','Position',[40 40 w h],'Color','w');
end

function save_pdf(fig, fname)
drawnow;
set(fig,'Color','w','InvertHardcopy','off');
set(fig,'Units','inches'); p = get(fig,'Position');
set(fig,'PaperUnits','inches','PaperSize',[p(3) p(4)], ...
        'PaperPosition',[0 0 p(3) p(4)],'PaperPositionMode','manual');
set(fig,'Units','pixels');
try
    exportgraphics(fig, fname, 'ContentType','vector','BackgroundColor','white');
catch
    print(fig, fname, '-dpdf', '-painters', '-r300');
end
end
