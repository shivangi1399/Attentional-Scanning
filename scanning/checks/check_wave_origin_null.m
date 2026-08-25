% =====================================================================
% CHECK (completes an existing analysis): a permutation null for the
% WAVE-ORIGIN alignment reported by cortical_planar_wave_PGD.m.
%
% WHY THIS EXISTS
%   cortical_planar_wave_PGD.m computes, per frequency and stimulus position,
%   the angle between the fitted propagation direction and the axis running
%   from the RF-DRIVEN patch to the centroid of the reliable channels, folded
%   to [0, 90] deg. Small means the wave radiates away from the patch the
%   stimulus drove -- which is what a scan starting at the stimulus would look
%   like. The number is reported, but with no null attached, so the Limitations
%   section of the paper has to describe it as untested rather than negative.
%   This script attaches the null. It changes nothing upstream.
%
% WHY THE OBVIOUS NULL (45 deg) IS WRONG
%   An axial angle folded to [0, 90] has a chance mean of 45 deg only if the
%   two axes are independent AND each is uniformly distributed. Neither holds
%   here. The fitted directions cluster tightly -- that is the planar gradient
%   the paper reports -- and the driven patches are retinotopically arranged,
%   so their source-to-centre axes are structured too. Two structured
%   distributions can produce a mean well below 45 deg with no position-
%   specific relationship whatever. Testing against 45 would manufacture a
%   result.
%
% THE NULL THAT MATCHES THE CLAIM  (primary, NULL_MODE = 'position')
%   The claim is that the wave direction at position p is aligned with the
%   patch driven by THAT stimulus. So hold both marginals fixed -- the set of
%   fitted directions and the set of source-to-centre axes, exactly as
%   observed, frequency by frequency -- and permute only which direction is
%   paired with which axis. Clustering of directions survives, retinotopic
%   arrangement of patches survives, and only the position-specific pairing is
%   destroyed. That isolates the claim and nothing else.
%
% THE SUPPLEMENTARY NULL  (NULL_MODE = 'phase')
%   Shuffle phases across electrode locations, refit the direction from the
%   shuffled map, and recompute the alignment. This is the same spatial
%   shuffle cortical_planar_wave_PGD.m already uses for PGD, and it asks the
%   weaker question of whether the alignment beats a random phase map. Run
%   both: 'position' is the one to report, 'phase' shows the alignment is not
%   an artefact of the map being unstructured.
%
% SELF-CHECK BEFORE ANY P-VALUE
%   The script reconstructs origin_align from the saved inputs and compares it
%   with the origin_align stored by cortical_planar_wave_PGD.m. If the two do
%   not agree to within RECON_TOL the reconstruction is wrong, the null would
%   be testing something other than the reported statistic, and the script
%   stops. Do not relax the tolerance to get past this -- fix the mismatch.
%
% INPUT
%   results_combined/scanning/planar_wave_existence/cp10_till_100/<dv>/
%       planar_wave_existence.mat   -> DIRp, sigp, ch_covers, pos_keep,
%                                      coh_sig, origin_align, freq
%   results_<animal>/scanning/phase_progression/cp10_till_100/<dv>/
%       phase_progression.mat       -> pref_phase, coh_mag
%   (the raw phase arrays are stripped from the saved A struct, exactly as the
%    "already in phase_progression.mat" comment there says, so they are
%    reloaded here the same way the figure code does)
%
% OUTPUT
%   Plots/scanning/checks/check_wave_origin_null.pdf
%   results_combined/scanning/checks/check_wave_origin_null.mat
%
% WHAT TO WRITE
%   A mean alignment BELOW the null with p < alpha in both animals means the
%   wave radiates from the stimulus-driven patch -- which would be a positive
%   for the scanning account and would need reconciling with the de-rotation
%   results. A mean alignment indistinguishable from the null means the
%   origin check is NEGATIVE rather than untested, and the Limitations
%   sentence can be deleted.
% =====================================================================

clearvars; close all; clc

%% --- Settings -------------------------------------------------------
animals    = {'hermes','klecks'};
base       = '/mnt/hpc/projects/MWSampling/4Shivangi';
dv         = 'lfp';
grid_rows  = 8; grid_cols = 8;
MIN_CH     = 8;                    % as in cortical_planar_wave_PGD.m
nPerm      = 10000;                % cheap here: no phase maps are refitted
alpha      = 0.05;
RECON_TOL  = 1e-6;                 % deg, for the self-check
NULL_MODES = {'position','phase'}; % primary first
nPermPhase = 1000;                 % for the supplementary null (refits maps)
rng(2025);

out_dir = fullfile(base,'Plots','scanning','checks');
res_dir = fullfile(base,'results_combined','scanning','checks');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
if ~exist(res_dir,'dir'), mkdir(res_dir); end

% Array coordinates -- IDENTICAL construction to cortical_planar_wave_PGD.m
nChTot = grid_rows*grid_cols;
ch_col = ceil((1:nChTot)' / grid_rows);
ch_row = grid_rows - mod((1:nChTot)' - 1, grid_rows);

fp_pgd = fullfile(base,'results_combined','scanning','planar_wave_existence', ...
    'cp10_till_100', dv, 'planar_wave_existence.mat');
assert(isfile(fp_pgd), 'planar_wave_existence.mat not found: %s', fp_pgd);
S = load(fp_pgd,'results'); AA = S.results.A;

R = struct();

%% --- Per animal -----------------------------------------------------
for ia = 1:numel(animals)
    animalName = animals{ia};
    ja = find(strcmp({AA.animal}, animalName), 1);
    if isempty(ja)
        warning('%s not in planar_wave_existence.mat -- skipped.', animalName); continue
    end
    a = AA(ja);
    if ~isfield(a,'origin_align') || all(isnan(a.origin_align(:)))
        warning('%s has no origin_align (was DO_PER_POSITION off?) -- skipped.', animalName); continue
    end

    % raw phase maps, reloaded and subset exactly as the PGD script did
    fp_pp = fullfile(base, ['results_' animalName], 'scanning', ...
        'phase_progression','cp10_till_100', dv, 'phase_progression.mat');
    assert(isfile(fp_pp), 'phase_progression.mat not found: %s', fp_pp);
    P = load(fp_pp,'pref_phase','coh_mag');
    phi_all = P.pref_phase(:, :, a.pos_keep);
    coh_all = P.coh_mag(:,    :, a.pos_keep);

    freq  = a.freq(:);
    nFreq = numel(freq);
    nPos  = numel(a.pos_keep);
    fprintf('\n##### %s : %d freqs x %d positions #####\n', animalName, nFreq, nPos);

    % ---- rebuild the source->centre axis and the alignment ------------
    axis_ang = nan(nFreq,nPos);        % angle of the driven->centre vector
    recon    = nan(nFreq,nPos);        % reconstructed origin_align, for the self-check
    for p = 1:nPos
        drv = find(a.ch_covers(:,p));
        if isempty(drv), continue, end
        src = [mean(ch_row(drv)), mean(ch_col(drv))];
        for f = 1:nFreq
            gp = build_grid_local(phi_all(:,f,p), coh_all(:,f,p), ch_row, ch_col, ...
                grid_rows, grid_cols, a.coh_sig(:,f));
            rc = ~isnan(gp);
            if sum(rc(:)) < MIN_CH, continue, end
            [rr,cc] = find(rc); ctr = [mean(rr) mean(cc)];
            outv = (ctr(2)-src(2)) + 1i*(ctr(1)-src(1));   % col + i*row, as upstream
            axis_ang(f,p) = angle(outv);
            if isfinite(a.DIRp(f,p))
                da = angle(exp(1i*(a.DIRp(f,p) - axis_ang(f,p))));
                recon(f,p) = rad2deg(min(abs(da), pi-abs(da)));
            end
        end
    end

    % ---- SELF-CHECK: does the reconstruction reproduce what was saved? --
    both = isfinite(recon) & isfinite(a.origin_align);
    dmax = max(abs(recon(both) - a.origin_align(both)));
    fprintf('  self-check: %d cells compared, max |recon - saved| = %.3g deg\n', sum(both(:)), dmax);
    assert(~isempty(dmax) && dmax < RECON_TOL, ...
        ['Reconstruction does not match the saved origin_align (max diff %.4g deg). ' ...
         'The null would not be testing the reported statistic. Fix this before reading any p-value.'], dmax);

    % ---- the cells the statistic is evaluated on ----------------------
    sel = logical(a.sigp) & isfinite(recon);      % planar-wave-significant (f,pos) cells
    if ~any(sel(:))
        warning('  %s: no significant (freq,position) cells -- nothing to test.', animalName); continue
    end
    obs = mean(recon(sel), 'omitnan');
    fprintf('  observed mean alignment over %d significant cells: %.1f deg\n', sum(sel(:)), obs);

    % =================================================================
    % NULL 1 (primary): permute which position's DIRECTION is paired with
    % which position's AXIS, within frequency. Both marginals are preserved.
    % =================================================================
    nullPos = nan(nPerm,1);
    for b = 1:nPerm
        v = nan(nFreq,nPos);
        for f = 1:nFreq
            pv = randperm(nPos);
            da = angle(exp(1i*(a.DIRp(f,:) - axis_ang(f,pv))));
            v(f,:) = rad2deg(min(abs(da), pi-abs(da)));
        end
        nullPos(b) = mean(v(sel), 'omitnan');
    end
    % one-sided: the claim is that alignment is BETTER (smaller) than chance
    p_pos = (1 + sum(nullPos <= obs)) / (1 + nPerm);

    % =================================================================
    % NULL 2 (supplementary): shuffle phases across electrode locations,
    % refit the direction, recompute the alignment. Same spatial shuffle the
    % PGD analysis uses; asks only whether the map being structured suffices.
    % =================================================================
    % The phase map of a given (freq, position) does not change between
    % permutations -- only which cell each phase lands in does. Build each map
    % ONCE and permute its values, rather than rebuilding it nPermPhase times.
    cells = find(sel);
    cache = cell(numel(cells),1);
    for q = 1:numel(cells)
        [f,p] = ind2sub(size(sel), cells(q));
        gp = build_grid_local(phi_all(:,f,p), coh_all(:,f,p), ch_row, ch_col, ...
            grid_rows, grid_cols, a.coh_sig(:,f));
        idx = find(~isnan(gp));
        cache{q} = struct('idx',idx,'vals',gp(idx),'f',f,'p',p);
    end
    nullPh = nan(nPermPhase,1);
    for b = 1:nPermPhase
        v = nan(numel(cells),1);
        for q = 1:numel(cells)
            c = cache{q};
            gs = nan(grid_rows,grid_cols);
            gs(c.idx) = c.vals(randperm(numel(c.vals)));
            da = angle(exp(1i*(pgd_dir_local(gs) - axis_ang(c.f, c.p))));
            v(q) = rad2deg(min(abs(da), pi-abs(da)));
        end
        nullPh(b) = mean(v, 'omitnan');
    end
    p_phase = (1 + sum(nullPh <= obs)) / (1 + nPermPhase);

    % ---- report ------------------------------------------------------
    fprintf('  NULL 1  position re-pairing : mean %.1f deg (95%% CI %.1f-%.1f)  p = %.4f  <-- report this\n', ...
        mean(nullPos), prctile(nullPos,2.5), prctile(nullPos,97.5), p_pos);
    fprintf('  NULL 2  spatial phase shuffle: mean %.1f deg (95%% CI %.1f-%.1f)  p = %.4f\n', ...
        mean(nullPh), prctile(nullPh,2.5), prctile(nullPh,97.5), p_phase);
    if p_pos < alpha
        verd = 'ORIGIN EFFECT: the wave radiates from the stimulus-driven patch';
    else
        verd = 'NO origin effect: alignment is what re-pairing positions gives by chance';
    end
    fprintf('  --> %s\n', verd);

    R(ia).animal = animalName;
    R(ia).freq = freq; R(ia).pos_keep = a.pos_keep;
    R(ia).recon = recon; R(ia).axis_ang = axis_ang; R(ia).sel = sel;
    R(ia).obs = obs; R(ia).nullPos = nullPos; R(ia).nullPh = nullPh;
    R(ia).p_pos = p_pos; R(ia).p_phase = p_phase; R(ia).verdict = verd;
    R(ia).nCells = sum(sel(:)); R(ia).recon_maxdiff = dmax;
end

%% --- Replication ----------------------------------------------------
fprintf('\n===== REPLICATION =====\n');
have = arrayfun(@(s) ~isempty(s.animal), R);
if sum(have) == numel(animals) && all([R(have).p_pos] < alpha)
    fprintf('  Origin effect significant in BOTH animals -- reconcile with the de-rotation result.\n');
elseif sum(have) == numel(animals)
    fprintf('  Not significant in both animals: the origin check is NEGATIVE, not untested.\n');
    fprintf('  The Limitations sentence calling it untested can be removed.\n');
else
    fprintf('  Only %d of %d animals produced a result -- no replication statement.\n', ...
        sum(have), numel(animals));
end

%% --- Figure ---------------------------------------------------------
idx = find(have); nA = numel(idx);
fg = new_fig(1200, 420*max(nA,1));
for j = 1:nA
    S2 = R(idx(j));
    ax = subplot(nA,2,(j-1)*2+1); hold(ax,'on');
    histogram(ax, S2.nullPos, 40, 'FaceColor',[.7 .7 .7],'EdgeColor','none', ...
        'Normalization','pdf','DisplayName','null: position re-pairing');
    histogram(ax, S2.nullPh, 30, 'FaceColor',[.9 .6 .2],'EdgeColor','none','FaceAlpha',.6, ...
        'Normalization','pdf','DisplayName','null: phase shuffle');
    xline(ax, S2.obs, 'r-','LineWidth',2,'DisplayName','observed');
    xline(ax, 45, 'k:','LineWidth',1.2,'DisplayName','45\circ (NOT the null)');
    xlabel(ax,'mean alignment, driven patch \rightarrow wave axis (deg)');
    ylabel(ax,'pdf');
    title(ax, {sprintf('%s -- wave-origin null (%d cells)', esc(S2.animal), S2.nCells), ...
               sprintf('observed %.1f\\circ, p = %.4f', S2.obs, S2.p_pos)},'FontSize',9);
    legend(ax,'Location','best','FontSize',7); grid(ax,'on');

    ax = subplot(nA,2,(j-1)*2+2);
    M = S2.recon; M(~S2.sel) = NaN;
    imagesc(ax, 1:numel(S2.pos_keep), S2.freq, M); set(ax,'YDir','normal');
    caxis(ax,[0 90]); colormap(ax, parula); cb = colorbar(ax); cb.Label.String = 'alignment (deg)';
    set(ax,'XTick',1:numel(S2.pos_keep),'XTickLabel',cellstr(num2str(S2.pos_keep(:))));
    xlabel(ax,'stimulus position'); ylabel(ax,'frequency (Hz)');
    title(ax, {'alignment per significant (freq, position) cell', ...
               'blue = radiates from the driven patch'},'FontSize',9);
end
sgtitle({'CHECK: does the planar wave radiate from the RF-driven patch? (a null for the origin statistic)', ...
         'The null holds both marginals fixed and permutes only which position''s direction pairs with which position''s axis.', ...
         '45\circ is NOT the chance level here: fitted directions cluster and driven patches are retinotopically arranged.'}, ...
         'FontSize', 9);
save_pdf(fg, fullfile(out_dir,'check_wave_origin_null.pdf'));

save(fullfile(res_dir,'check_wave_origin_null.mat'), 'R', 'nPerm','nPermPhase', ...
    'alpha','dv','MIN_CH','-v7.3');
fprintf('\nSaved figure to %s\n', fullfile(out_dir,'check_wave_origin_null.pdf'));

%% ===================== local functions ==============================

function grid = build_grid_local(phi_ch_in, coh_in, ch_row, ch_col, nR, nC, keep)
% IDENTICAL to build_grid in cortical_planar_wave_PGD.m. Kept local so this
% check cannot drift from the analysis it completes without the self-check
% above catching it.
phi_ch = phi_ch_in(:); rel = coh_in(:); keep = logical(keep(:));
ok = keep & isfinite(phi_ch) & isfinite(rel) & rel > 0;
grid = nan(nR,nC);
for ch = 1:numel(phi_ch)
    if ok(ch), grid(ch_row(ch),ch_col(ch)) = phi_ch(ch); end
end
end

function d = pgd_dir_local(grid)
% Propagation direction only -- the .dir field of pgd_metrics upstream.
[nR,nC] = size(grid);
gx = nan(nR,nC); gy = nan(nR,nC);
gx(:,1:nC-1) = angle(exp(1i*(grid(:,2:nC) - grid(:,1:nC-1))));
gy(1:nR-1,:) = angle(exp(1i*(grid(2:nR,:) - grid(1:nR-1,:))));
g = gx + 1i*gy; v = isfinite(g);
d = angle(mean(g(v)));
end

function s = esc(s), s = strrep(s, '_', '\_'); end

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
