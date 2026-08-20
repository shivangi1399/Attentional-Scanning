% =====================================================================
% CHECK (not a main analysis): is the planar phase offset produced by the
% shared recording reference?
%
% WHY THIS EXISTS
%   cortical_planar_wave_PGD.m finds a planar phase gradient across the 8x8
%   array whose wavenumber does not scale with frequency -- a fixed spatial
%   phase offset rather than a travelling wave. A shared recording reference
%   is one candidate explanation for such an offset. THIS SCRIPT TESTS IT.
%
% THE IDEA
%   Every channel is recorded against one shared reference, so
%       Z_c = S_c - Z_R           for every electrode c.
%   The SAME vector Z_R is subtracted from source vectors of DIFFERENT
%   amplitude, so it rotates weak channels more than strong ones. Amplitude
%   varies smoothly across the array, so the induced rotation does too, and a
%   smooth phase ramp appears with nothing propagating. It is frequency-
%   independent, which is exactly the observed signature.
%
%   Common-average re-referencing removes it EXACTLY:
%       Z_c - mean_c(Z) = (S_c - Z_R) - (mean_c S - Z_R) = S_c - mean_c S
%   The reference cancels. So if the gradient is caused by the reference it
%   MUST disappear under CAR; if it survives CAR, the reference did not make it.
%
% WHY THERE IS NO LAPLACIAN HERE
%   A Laplacian re-reference is the other obvious control, but it cannot
%   answer this question: the Laplacian of a linear ramp is zero by
%   construction, so it removes a smooth spatial gradient whether that
%   gradient is neural or an artefact. A null Laplacian result would say
%   nothing about the reference, so it is not computed.
%
% WHAT IS COMPARED
%   Each referencing scheme is judged against ITS OWN spatially shuffled null.
%   This matters: re-referencing removes the large shared component, so the
%   residual is smaller and noisier and the fitted wavenumber inflates under
%   every scheme. Only the EXCESS over the matched null is interpretable.
%
% WHAT THIS IS NOT
%   Works on the raw between-channel phase relationship, NOT on pref_phase
%   (which is additionally weighted by the dependent variable). It is a
%   control on the referencing question, not a re-run of the main analysis.
%   Numbers are not expected to equal those from cortical_planar_wave_PGD.m.
%
% INPUT   results_<animal>/<session>/Phase_analysis/hit_miss/
%             100iter_cut@cp_m10/phase/<chan_folder>/phase.mat
%         Fields used: ar_phase (trials x freq), amp_vec (trials x freq),
%         label (session channel list).
%         NOTE 1: the 'hit_miss' folder holds the LFP phase estimate itself;
%                 it is not DV-specific.
%         NOTE 2: amp_vec is used because it is the RAW amplitude. The
%                 amplitude in ph_all_sess.mat is z-scored per channel, which
%                 destroys the across-channel scale that re-referencing needs.
%
% OUTPUT  Plots/scanning/checks/check_rereferencing.pdf
%         results_combined/scanning/checks/check_rereferencing.mat
% =====================================================================

clearvars; close all; clc

%% --- Settings -------------------------------------------------------
base        = '/mnt/hpc/projects/MWSampling/4Shivangi';
animals     = {'klecks','hermes'};
SPACING_MM  = 0.4;
GRID        = 8;
nPerm       = 200;          % spatial shuffles for the per-scheme null
MAX_SESS    = Inf;          % set to e.g. 3 for a quick run
MIN_TRIALS  = 40;           % skip sessions with fewer clean trials
schemes     = {'ORIGINAL','CAR'};
rng(1);

out_dir = fullfile(base,'Plots','scanning','checks');
res_dir = fullfile(base,'results_combined','scanning','checks');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
if ~exist(res_dir,'dir'), mkdir(res_dir); end

% canonical channel 1..64 -> grid position (same convention as the PGD script)
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

    acc = struct(); for s = 1:numel(schemes), acc.(schemes{s}) = []; end
    freq = []; nTot = 0;

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
        fprintf('  . %-45s %5d trials\n', sessions{is}, size(A,1));

        x = (ch_col(present)-1)*SPACING_MM;
        y = (ch_row(present)-1)*SPACING_MM;

        % --- the two referencing schemes ---
        A_car = A - mean(A,3,'omitnan');

        % One global phase origin per trial, shared by all schemes. It is a
        % single rotation applied to every channel alike, so it cannot create
        % or destroy SPATIAL structure; it only lets trials be averaged.
        theta = angle(sum(A,3,'omitnan'));                       % trials x freq

        M = {A, A_car};
        for s = 1:numel(schemes)
            W = squeeze(mean(M{s} .* exp(-1i*theta), 1, 'omitnan'));   % freq x chan
            rec = nan(numel(freq),4);
            for f = 1:numel(freq)
                a = angle(W(f,:)).';
                a = angle(exp(1i*(a - angle(mean(exp(1i*a))))));       % centre, no wrap issue
                [k,~]  = plane_fit(a, x, y);
                kn = nan(nPerm,1); pn = nan(nPerm,1);
                for b = 1:nPerm
                    p = randperm(numel(a));
                    kn(b) = plane_fit(a(p), x, y);
                    pn(b) = planarity(a(p), x, y);
                end
                rec(f,:) = [rad2deg(k*(GRID-1)*SPACING_MM), ...
                            rad2deg(mean(kn)*(GRID-1)*SPACING_MM), ...
                            planarity(a,x,y), mean(pn)];
            end
            acc.(schemes{s}) = [acc.(schemes{s}); rec];
        end
    end

    if isempty(acc.ORIGINAL), continue, end

    % EXCESS is the difference of the medians shown in the same row, so the
    % printed table is internally consistent (a reader can subtract the columns
    % and get the excess column). %>null is the fraction of frequency bins in
    % which the observed value beat its own null: 50% is chance, and it is the
    % most robust of the three because it does not depend on the size of the
    % difference, only on its sign.
    fprintf('\n  %-10s %10s %10s %10s %8s   %10s %10s %10s\n', ...
        'scheme','tilt(deg)','null(deg)','EXCESS','%>null','planarity','null','EXCESS');
    for s = 1:numel(schemes)
        Mv = acc.(schemes{s});
        fprintf('  %-10s %10.1f %10.1f %+10.1f %7.0f%%   %10.2f %10.2f %+10.2f\n', schemes{s}, ...
            median(Mv(:,1)), median(Mv(:,2)), median(Mv(:,1))-median(Mv(:,2)), ...
            100*mean(Mv(:,1) > Mv(:,2)), ...
            median(Mv(:,3)), median(Mv(:,4)), median(Mv(:,3)-Mv(:,4)));
    end

    % Three outcomes, in order. The first guard matters: if ORIGINAL shows no
    % gradient above its own null there is nothing to explain, and the check is
    % simply uninformative for this animal -- NOT evidence of an artefact.
    % The criterion is deliberately the SIGN test (fraction of frequency bins
    % above the matched null, chance = 50%) rather than the size of the tilt
    % difference, which is skewed across frequencies and can disagree in sign
    % with the difference of the medians.
    hasSig = @(Mv) mean(Mv(:,1) > Mv(:,2)) > 0.75 && median(Mv(:,3)-Mv(:,4)) > 0.02;
    if ~hasSig(acc.ORIGINAL)
        verdict = 'no gradient above the null even BEFORE re-referencing -> underpowered, check not informative';
    elseif hasSig(acc.CAR)
        verdict = 'SURVIVES common-average referencing -> the gradient is NOT a reference artefact';
    else
        verdict = 'does NOT survive CAR -> consistent with a reference artefact';
    end
    fprintf('  --> %s   [%d trials]\n', verdict, nTot);

    R.(animalName).acc = acc; R.(animalName).freq = freq;
    R.(animalName).nTrials = nTot; R.(animalName).verdict = verdict;
    R.(animalName).nSess = nS;
end

%% --- Figure ---------------------------------------------------------
% Row 1: tilt vs frequency, observed against its own shuffled null, per scheme.
% Row 2: the summary bar -- EXCESS over null. This is the panel to read.
an = fieldnames(R);
if ~isempty(an)
    fg = new_fig(420*numel(an), 640);
    cols = [0.15 0.35 0.75; 0.85 0.35 0.10];
    for k = 1:numel(an)
        A = R.(an{k}); fr = A.freq(:);
        % -- row 1: tilt vs frequency
        ax = subplot(2, numel(an), k); hold(ax,'on');
        for s = 1:numel(schemes)
            Mv = A.acc.(schemes{s});
            nF = numel(fr); nRep = size(Mv,1)/nF;
            obs = median(reshape(Mv(:,1), nF, nRep), 2);
            nul = median(reshape(Mv(:,2), nF, nRep), 2);
            plot(ax, fr, obs, '-',  'Color', cols(s,:), 'LineWidth', 1.6, 'DisplayName', schemes{s});
            plot(ax, fr, nul, ':',  'Color', cols(s,:), 'LineWidth', 1.1, 'HandleVisibility','off');
        end
        set(ax,'XScale','log'); grid(ax,'on');
        xlabel(ax,'frequency (Hz)'); ylabel(ax,'tilt across array (deg)');
        title(ax, sprintf('%s -- solid = observed, dotted = shuffled null', an{k}), 'FontSize', 9);
        legend(ax,'Location','best','Box','off','FontSize',7);

        % -- row 2: excess over null (the panel that answers the question)
        ax2 = subplot(2, numel(an), numel(an)+k); hold(ax2,'on');
        ex = nan(1,numel(schemes));
        for s = 1:numel(schemes)
            Mv = A.acc.(schemes{s});
            ex(s) = median(Mv(:,1)) - median(Mv(:,2));   % consistent with the printed table
        end
        for s = 1:numel(schemes)      % one bar at a time: CData indexing needs R2017b+
            bar(ax2, s, ex(s), 'FaceColor', cols(s,:), 'EdgeColor','k');
        end
        yline(0, 'k-');
        set(ax2,'XTick',1:numel(schemes),'XTickLabel',schemes,'XTickLabelRotation',20);
        ylabel(ax2,'EXCESS tilt over null (deg)'); grid(ax2,'on');
        title(ax2, {'positive = real gradient above chance', ...
                    sprintf('CAR is the decisive scheme (%d trials)', A.nTrials)}, 'FontSize', 8);
    end
    sgtitle({'CHECK: is the planar phase offset produced by the shared reference?', ...
             'CAR cancels the reference exactly -- if the gradient survives CAR, the reference did not make it.'}, ...
             'FontSize', 10);
    save_pdf(fg, fullfile(out_dir,'check_rereferencing.pdf'));
end

save(fullfile(res_dir,'check_rereferencing.mat'), 'R', 'schemes', 'nPerm', '-v7.3');
fprintf('\nSaved figure to %s\n', fullfile(out_dir,'check_rereferencing.pdf'));

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

function [k, dir] = plane_fit(a, x, y)
% least-squares plane through centred phases -> wavenumber (rad/mm), direction
b = [ones(numel(x),1) x y] \ a;
k = hypot(b(2), b(3)); dir = atan2(b(3), b(2));
end

function p = planarity(a, x, y)
% fraction of phase variance explained by the best plane, 0..1
[k, d] = plane_fit(a, x, y);
pred = k*(x*cos(d) + y*sin(d)); pred = pred - mean(pred);
p = max(0, 1 - var(a - pred)/max(var(a), eps));
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
