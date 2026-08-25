% =====================================================================
% figures.m  —  paper figure set, rebuilt from saved results
%
% WHAT THIS IS
% ------------
% One script that re-renders every result figure listed in
%   writing/paper/V03/result_figs.docx
% in a single, consistent, Illustrator-editable house style, and writes
% them to  writing/figures/plots/*.pdf
%
% It is a RENDERING script, not an analysis script. Every number it draws is
% read from a .mat that the analysis pipelines already wrote. It recomputes no
% statistic, builds no permutation null, and writes nothing at all into the
% results tree - the only files it creates are the PDFs. The original analysis
% scripts are untouched and remain the source of truth.
%
% EVERYTHING HERE IS POOLED
% -------------------------
% Every panel in the set shows the two animals combined, or the monkey
% average. No figure draws a per-animal curve any more, which is why there is
% no per-animal colour: hermes and klecks are distinguished by which panel
% they are in (figures 7b and 7d), never by hue.
%
% FIGURES WRITTEN
% ---------------
%   fig03_erp_grandavg.pdf        RECTIFIED LFP and MUA evoked response,
%                                 hit against miss, channels and sessions
%                                 pooled; |hit-miss| significance shaded
%   fig04_tfr.pdf                 pooled hit-miss TFR and the pooled
%                                 fooof-corrected pre-stimulus spectrum
%   fig05_phase_measures.pdf      monkey-average coherence / correlation /
%                                 regression R^2 for MUA, LFP, RT, detection
%   fig06_hypotheses.pdf          H1-H4 overlay, monkey average
%   fig07a_phase_progression.pdf  preferred-phase systematicity vs frequency
%   fig07b_pgd.pdf                pooled z-PGD vs frequency, plus the
%                                 clearest example position in each animal
%   fig07c_derotation.pdf         frequency x speed de-rotation grids,
%                                 cortical and along the stimulus axis
%   fig07d_per_position_<animal>.pdf   per-position phase maps at the two
%                                 LFP-coherence peak frequencies, wrapped
%                                 over two rows per frequency
%
% Block header, DO_* toggle, figure variable and file name all carry the same
% letter, and the blocks appear in that order. Renaming a figure means
% changing all four together.
%
% ILLUSTRATOR-EDITABLE OUTPUT
% ---------------------------
% Everything is printed with -dpdf -painters, so the PDF is vector and the
% text stays live text (selectable, re-typable), not outlines.
%
% Two consequences that shaped the code:
%   1. NO TRANSPARENCY. MATLAB silently falls back to the OpenGL renderer
%      when a figure contains FaceAlpha/EdgeAlpha, which rasterises the whole
%      page. Every "translucent" fill here is instead an OPAQUE colour that
%      has been pre-blended toward white by tint(). Visually identical,
%      fully vector.
%   2. NO imagesc. imagesc embeds a bitmap that Illustrator can only move
%      and mask. Every heatmap is drawn by vec_heatmap() as one patch object
%      with one flat-shaded face per cell, and every colourbar by
%      vec_colorbar(). The grids here are small (at most 39x37), so the file
%      stays light and every cell is individually selectable.
%
% COLOURS - three rules, see palette() for the values and the reasoning
% ---------------------------------------------------------------------
%   R1  HUE NAMES THE MEASURED QUANTITY, and nothing else.
%       purple = MUA, teal = LFP, orange = reaction time, rose = detection.
%       A reader who learns the four hues from figure 5 carries them into
%       3, 4, 6 and 7a without relearning, which is the whole return on
%       having a palette at all.
%   R2  LIGHTNESS NAMES A LEVEL WITHIN THAT QUANTITY.
%       hit vs miss, and H1..H4, are levels of one measured thing rather
%       than four unrelated things, so they are rungs of that quantity's own
%       hue (PAL.pair.<key>, PAL.lad.<key>) and not four separate hues.
%       Lightness is also the channel that survives greyscale printing and
%       colour-vision deficiency.
%   R3  COLOUR NEVER ENCODES A STATISTIC.
%       Grey is the null, the threshold and the bias floor - anything grey is
%       chance level, never data. Significance is a band_tint() of the
%       curve's OWN hue, so the shading names the curve it belongs to rather
%       than adding a fifth meaning to teal.
%
% Separately, PURPLE -> PINK -> YELLOW (plasma, see cmap_ppy) carries every
% heatmap: the time-frequency map, the de-rotation grids and the array phase
% maps. It is used only for filled areas - no LINE in this set is yellow, and
% the four curve hues all sit outside the plasma ramp - so a yellow or pink
% AREA is always a colour scale and never a trace.
%
% SIGNIFICANCE, DRAWN THE SAME WAY EVERYWHERE
% -------------------------------------------
%   curves : dashed grey line = permutation threshold
%            pale band + solid bar along the top = significant, both in a
%            tint of the curve's own hue (band_tint, one constant everywhere)
%   maps   : WHITE outline around the significant region, padded by
%            sig_outline() so every loop CLOSES on the edge of the map - an
%            open arc leaves the reader unable to tell whether the region
%            stops there or runs off the plot
% Each panel's caption states which test produced it, because the tests
% differ (max-stat across frequency, cluster correction, sign-flip across
% sessions, pooled-z, BH-FDR across positions).
%
% FONTS AND LINE WIDTHS
% ---------------------
% Exactly TWO font sizes in the whole set, both set at the top: FONT_SIZE for
% everything inside a panel and FONT_SIZE_BIG for panel letters and column
% headers. Every data and threshold line inherits one width from
% DefaultLineLineWidth - no plot call states a size or a width of its own, so
% uniformity is structural rather than maintained by hand. Only the faint
% guides (zero lines, the stimulus-onset marker) set LW_GUIDE.
%
% LAYOUT
% ------
% Panels are placed by grid_pos(), whose margins and gaps are given in
% CENTIMETRES: these pages run from 9 to 24 cm tall and 15 to 21 cm wide, and
% a normalised gap that looks right on one is cramped on another. The margins
% already reserve room for the caption, the axis labels and the panel
% letters, so nothing is squeezed afterwards. Panels are lettered a, b, c ...
% in reading order by panel_label(), for the manuscript to refer to.
%
% USAGE
% -----
% Set the DO_* toggles below and run. Each figure block is independent, so
% a single figure can be rebuilt without touching the others.
% =====================================================================

clearvars; close all; clc

%% ─── Settings ────────────────────────────────────────────────────────
base     = '/mnt/hpc/projects/MWSampling/4Shivangi';
comb     = fullfile(base,'results_combined');
out_dir  = fullfile(base,'writing','figures','plots');
if ~exist(out_dir,'dir'), mkdir(out_dir); end

DO_FIG3  = 1;    % fig03  rectified evoked response, hit vs miss
DO_FIG4  = 1;    % fig04  TFR + pre-stimulus periodic spectrum
DO_FIG5  = 1;    % fig05  monkey-average phase measures
DO_FIG6  = 1;    % fig06  H1-H4 hypothesis overlay
DO_FIG7A = 1;    % fig07a phase-progression systematicity
DO_FIG7B = 1;    % fig07b pooled z-PGD + example position per animal
DO_FIG7C = 1;    % fig07c de-rotation grids
DO_FIG7D = 1;    % fig07d per-position phase maps, wrapped

erp_xlim  = [-0.10 0.30];   % s. Wide enough for the hermes LFP peak at ~136 ms
tfr_run   = 'unified_4to80Hz_2cyc';
pow_run   = 'prestim_2to80Hz';
dv_scan   = 'lfp';          % scanning analyses are LFP
cp        = 'cp10_till_100';

% Sessions excluded from the ERP. This is pool_both_animals.m's list, used for
% both signals, and it must stay that way: the finished pooled ERP result this
% figure reads was computed over exactly these sessions (klecks 21, hermes 28),
% so averaging the evoked waveforms over any other set would put the curves and
% the significance mask on different data.
erp_animals(1).name = 'hermes';
erp_animals(1).dir  = fullfile(base,'results_hermes');
erp_animals(1).drop = [1 3 21];
erp_animals(2).name = 'klecks';
erp_animals(2).dir  = fullfile(base,'results_klecks');
erp_animals(2).drop = [2 15 16 18];

% Canonical frequency axis. The coherence and correlation result files carry
% their own 'freq'; the regression files do not (their grid is the same one),
% so this is the fallback used whenever a file omits it and the lengths agree.
FREQ_AXIS = [];
fcan = fullfile(comb,'phase_coherence','complex',cp,'lfp','all_loc_difflev','monkey_avg_results.mat');
if isfile(fcan), Sc = load(fcan,'freq'); FREQ_AXIS = Sc.freq(:)'; end

% MATLAB's TeX interpreter mis-spaces a superscript in this release: 'R^{2}'
% prints the 2 detached and far to the right of the R. The Unicode superscript
% glyph sets correctly and survives the PDF export as live text, so R-squared
% is written with char(178) throughout. Built here rather than pasted as a
% literal so this line stays readable in any editor.
R2LAB = ['R' char(178)];

% Phase colour scales are ticked in pi, not in decimals. Unicode rather than
% TeX \pi for the same reason as R2LAB: TeX would be rendered at its own size
% and break the two-font-size rule.
PI_TICKS = {['-' char(960)], '0', char(960)};

% Colour family for every heatmap in the set: purple -> pink -> yellow.
% 'plasma' or 'magma'; see cmap_ppy().
MAP_FAMILY = 'magma';

% ─── Fonts ───────────────────────────────────────────────────────────
% TWO sizes in the whole set and nothing else: FONT_SIZE for everything a
% reader reads inside a panel (ticks, axis labels, titles, legends, the
% caption, colourbar labels), and FONT_SIZE_BIG for the panel letters and
% column headers only. Change either number here and every figure follows -
% no call below states a size of its own, they are inherited from the
% graphics defaults set immediately after.
FONT_NAME     = 'Helvetica';
FONT_SIZE     = 8;
FONT_SIZE_BIG = 10;

PAL = palette();
set(groot, 'DefaultAxesFontName',FONT_NAME, 'DefaultTextFontName',FONT_NAME, ...
           'DefaultLegendFontName',FONT_NAME, 'DefaultColorbarFontName',FONT_NAME, ...
           'DefaultAxesFontSize',FONT_SIZE, 'DefaultTextFontSize',FONT_SIZE, ...
           'DefaultLegendFontSize',FONT_SIZE, 'DefaultColorbarFontSize',FONT_SIZE, ...
           'DefaultAxesTitleFontSizeMultiplier',1, ...
           'DefaultAxesLabelFontSizeMultiplier',1, ...
           'DefaultAxesTitleFontWeight','normal', ...
           'DefaultAxesLineWidth',0.6, ...
           'DefaultAxesTickDir','out', 'DefaultAxesTickLength',[0 0], ...
           'DefaultAxesBox','off', 'DefaultAxesColor','none', ...
           'DefaultLineLineWidth',1.5, 'DefaultFigureColor','w');
% Every data and threshold line in the set inherits DefaultLineLineWidth, so
% no plot call below states a LineWidth of its own. The only exceptions are
% the faint guides - zero lines, the stimulus-onset marker - which are drawn
% at LW_GUIDE because they are not data.
LW_GUIDE = 0.6;


%% =====================================================================
%% Figure 3 — WRITES fig03_erp_grandavg.pdf
%% rectified LFP and MUA evoked response, hit against miss
%% =====================================================================
% Sources, both already computed - nothing here is re-run:
%   results_<animal>/<session>/ERP_<SIG>/ERP_real/norm_{hit,miss}_timelock.mat
%       the evoked waveforms. Each channel is RECTIFIED first, then averaged
%       within a session, then across the 49 sessions of both animals.
%   results_combined/group_ERP/<sig>/erp_hitmiss_pooled_both.mat
%       the finished pooled result. The shading is its RECTIFIED variant:
%       mean over channels of |hit - miss|, one-sided, corrected. That is the
%       tested quantity, so it is the only thing marked.

ERP_RECTIFY = 1;    % 0 = signed evoked response, the conventional ERP

if DO_FIG3
    fprintf('\n===== Figure 3: rectified evoked response =====\n');

    sigs = struct('key',{'lfp','mua'}, 'sub',{'ERP_LFP','ERP_MUA'}, 'lab',{'LFP','MUA'});

    f3   = new_fig(19, 9.0);
    % Margins chosen so the DRAWN BLOCK is centred, not the axes boxes: the
    % left margin also carries the y label, tick labels and the panel letter
    % (~1.7 cm of ink), so it has to exceed the right margin by that much for
    % the white space either side of the figure to come out equal. The top
    % margin clears the caption with room to spare.
    M3   = [2.9 1.3 1.2 3.1 3.0 0];
    lets = panel_letters(numel(sigs));

    for is = 1:numel(sigs)
        W  = load_erp_waveforms(erp_animals, sigs(is));
        P  = load_erp_pooled(sigs(is), comb, 'rect');
        ax = axes('Position', grid_pos(f3, 1, numel(sigs), is, M3)); hold(ax,'on');
        panel_label(ax, lets{is}, FONT_SIZE_BIG);

        if isempty(W)
            axis(ax,'off');
            text(0.5,0.5,'no data','Parent',ax,'HorizontalAlignment','center');
            continue
        end
        if ERP_RECTIFY
            hitv = W.hit_rect; missv = W.miss_rect;
            ylab = sprintf('%s mean |amplitude|', sigs(is).lab);
        else
            % reachable only when ERP_RECTIFY is set to 0 at the top; the
            % Code Analyzer cannot see that because the toggle is a constant
            hitv = W.hit; missv = W.miss;  %#ok<UNRCH>
            ylab = sprintf('%s (normalised)', sigs(is).lab);
        end

        w  = W.time >= erp_xlim(1) & W.time <= erp_xlim(2);
        yl = legend_headroom(pad_lim([hitv(w) missv(w)]), 0.34);

        % hue = the signal this panel measures; the two conditions are a
        % lightness pair inside it, and the significance band is a tint of
        % that same hue so it names the curves it belongs to
        hue = PAL.(sigs(is).key);
        pr  = PAL.pair.(sigs(is).key);
        if ~isempty(P)
            wm = P.time >= erp_xlim(1) & P.time <= erp_xlim(2);
            sig_band(ax, P.time(wm), P.mask(wm), yl, band_tint(hue));
        end
        plot(ax, erp_xlim, [0 0], '-', 'Color', PAL.grey_lt, 'LineWidth', LW_GUIDE);
        plot(ax, [0 0], yl,      '-', 'Color', PAL.grey_lt, 'LineWidth', LW_GUIDE);
        hh(1) = plot(ax, W.time(w), hitv(w),  '-', 'Color', pr(1,:));
        hh(2) = plot(ax, W.time(w), missv(w), '-', 'Color', pr(2,:));
        if ~isempty(P)
            sig_bar(ax, P.time(wm), P.mask(wm), yl, hue);
        end

        xlim(ax, erp_xlim); ylim(ax, yl);
        xlabel(ax,'Time from stimulus (s)'); ylabel(ax, ylab);
        % No session count on the panel. The ERP pools 49 sessions (klecks 21,
        % hermes 28) while the TFR and POW pool 48 (hermes 27) - one hermes
        % session has ERP data but no TFR file - and two different counts on
        % two figures of the same set invites the reader to think one is a
        % mistake. The counts are printed to the console on every run and
        % belong in the methods, stated once, with the reason.
        title(ax, sigs(is).lab, 'FontWeight','normal');
        lg = legend(hh, {'Hit','Miss'}, 'Location','northeast');
        style_legend(lg);
        place_legend(ax, lg, 'right');
        tidy(ax);
        fprintf('  %s: %d sessions (%s)%s\n', sigs(is).key, W.nS, W.split, ...
                ternary(isempty(P), ', no pooled result', ''));
    end
    supertitle(f3, {'Rectified evoked response, channels and sessions pooled, hit against miss', ...
        'Shaded = |hit - miss| significant in the finished pooled test (one-sided max-stat, p<0.05)'});
    save_fig(f3, out_dir, 'fig03_erp_grandavg');
end


%% =====================================================================
%% Figure 4 — WRITES fig04_tfr.pdf
%% hit-miss TFR and pre-stimulus periodic power
%% =====================================================================
% Source: results_combined/group_TFR/<run>/tfr_hitmiss_pooled_both.mat
%         results_combined/group_POW/<run>/pow_hitmiss_pooled_both.mat
%         both written by pool_both_animals.m and read here as stored.
%
% Everything in this figure is the two animals combined.
%
% TFR   unit of observation is the TRIAL (the permutation relabels trials
%       within a session), sessions are averaged, and the max-pixel
%       correction runs on that average. `mask` is the corrected map, drawn
%       as the outline.
% POW   unit of observation is the SESSION (fooof needs an averaged
%       spectrum), tested by sign-flip across the 48 sessions.
%
% The power panel shows the pooled fooof-corrected spectrum for each
% condition, with every frequency the sign-flip test marks as significant
% marked, over the whole range and with no band division. The marks refer to
% the hit-minus-miss difference, which is what the test was run on; the two
% conditions differ by far less than the spread of the spectrum itself, so
% the mark is what carries the result, not the visible gap between the
% curves.

if DO_FIG4
    fprintf('\n===== Figure 4: TFR and pre-stimulus power =====\n');

    T = load(fullfile(comb,'group_TFR',tfr_run,'tfr_hitmiss_pooled_both.mat'));
    P = load(fullfile(comb,'group_POW',pow_run,'pow_hitmiss_pooled_both.mat'));

    cmap_tfr = cmap_ppy(256, MAP_FAMILY);
    f4 = new_fig(19, 10.5);
    % generous horizontal gap: the left panel carries a colourbar and the
    % right one a y-label, and they must not meet
    % As figure 3, plus a wide column gap: the colourbar and its label live
    % in that gap, and panel b's y label after them.
    M4   = [2.9 1.4 1.2 3.1 4.6 0];
    lets = panel_letters(2);

    % ── left: the pooled TFR ──────────────────────────────────────────
    % The pooled log10(hit/miss) map is POSITIVE AT EVERY PIXEL: it runs from
    % +0.0027 to +0.0450, i.e. hits carry more power than misses everywhere in
    % the window, at every frequency. There is therefore no negative side to
    % represent, and a scale symmetric about zero would spend half the colour
    % range on values that do not occur - which is what flattened this map into
    % one shade. The limits run from 0 (no difference) to the 99th percentile,
    % so the whole purple-pink-yellow ramp is used by the data that exists:
    % purple = smallest hit advantage, yellow = largest.
    M    = double(T.real_avg);
    v    = M(~isnan(M));
    clim = [0 quantile(v, 0.99)];
    ax = axes('Position', grid_pos(f4, 1, 2, 1, M4)); hold(ax,'on');
    panel_label(ax, lets{1}, FONT_SIZE_BIG);
    vec_heatmap(ax, T.times, T.freqs, M, cmap_tfr, clim);
    sig_outline(ax, T.times, T.freqs, T.mask, [1 1 1]);
    plot(ax, [0 0], full_lims(T.freqs), '-', 'Color', PAL.ink, 'LineWidth', LW_GUIDE);
    xlim(ax, full_lims(T.times)); ylim(ax, full_lims(T.freqs));
    set(ax,'Layer','top','Box','on');
    xlabel(ax,'Time from stimulus (s)'); ylabel(ax,'Frequency (Hz)');
    title(ax, sprintf('Hit minus miss power, %d sessions', T.nSess), 'FontWeight','normal');
    tidy(ax);
    % the measure itself, not a shorthand: TFRdiff_hitmiss.m's statistic is
    % the log10 ratio of power between conditions
    % shrink = false: the column gap already reserves the strip's space, so
    % shrinking here as well would leave panel a narrower than panel b and
    % the pair visibly lopsided
    vec_colorbar(ax, cmap_tfr, clim, 'Log10(hit / miss)', false);
    fprintf('  TFR: %d sessions, %d/%d bins significant, range %+.4f..%+.4f (%.0f%% positive)\n', ...
            T.nSess, sum(T.mask(:)), numel(v), min(v), max(v), 100*mean(v>0));

    % ── right: the pooled fooof-corrected spectrum, hit against miss ──
    fr    = P.freqs(:)';
    dP    = P.stat_flat.d(:)';                % pooled hit-miss, periodic
    mk    = P.stat_flat.mask(:)';             % sign-flip, max-stat corrected
    hitS  = mean(P.FHa, 1, 'omitnan');
    missS = mean(P.FMa, 1, 'omitnan');

    ax = axes('Position', grid_pos(f4, 1, 2, 2, M4)); hold(ax,'on');
    panel_label(ax, lets{2}, FONT_SIZE_BIG);
    yl = legend_headroom(pad_lim([hitS missS]), 0.34);
    % the pre-stimulus spectrum is an LFP measurement, so it takes the LFP
    % hue and the LFP hit/miss pair
    sig_band(ax, fr, mk, yl, band_tint(PAL.lfp));
    hh(1) = plot(ax, fr, hitS,  '-', 'Color', PAL.pair.lfp(1,:));
    hh(2) = plot(ax, fr, missS, '-', 'Color', PAL.pair.lfp(2,:));
    sig_bar(ax, fr, mk, yl, PAL.lfp);
    xlim(ax,[min(fr) max(fr)]); ylim(ax, yl);
    xlabel(ax,'Frequency (Hz)'); ylabel(ax,'Periodic power (FOOOF)');
    title(ax, sprintf('Pre-stimulus periodic power, %d sessions', P.nS), 'FontWeight','normal');
    lg = legend(hh, {'Hit','Miss'}, 'Location','northeast');
    style_legend(lg);
    place_legend(ax, lg, 'right');
    tidy(ax);
    fprintf('  POW: %d sessions, hit-miss significant at %s Hz (mean diff %s)\n', ...
            P.nS, strjoin(compose('%.1f', fr(mk)), ', '), ...
            strjoin(compose('%+.4f', dP(mk)), ', '));

    supertitle(f4, {'Hit minus miss time-frequency response and pre-stimulus periodic power, both animals combined', ...
        'White outline = max-pixel corrected p<0.05; shaded = sign-flip across sessions, max-stat corrected'});
    save_fig(f4, out_dir, 'fig04_tfr');
end


%% =====================================================================
%% Figure 5 — WRITES fig05_phase_measures.pdf
%% monkey-average phase measures
%% =====================================================================
% Source: results_combined/{phase_coherence,phase_correlation,multi_lin_reg}/
%         complex/cp10_till_100/<dv>/.../monkey_avg_results.mat
% Same files, same variable names and same monkey-average curves that
% compare_all_measures.m plots; only the rendering differs.
%
% Rows are the dependent variables, columns the three ways of asking the
% same question of them. Each curve carries its own permutation threshold
% (max-statistic across frequency), drawn as the dashed grey line.
%
% H1 ('complex') throughout: the phase is pooled before the magnitude is
% taken, so nothing here inherits a Jensen advantage. Figure 6 is where the
% relaxed versions are compared.

if DO_FIG5
    fprintf('\n===== Figure 5: monkey-average phase measures =====\n');

    coh_root  = fullfile(comb,'phase_coherence',  'complex', cp);
    corr_root = fullfile(comb,'phase_correlation','complex', cp);
    reg_root  = fullfile(comb,'multi_lin_reg',    'complex', cp);

    rows = struct( ...
        'lab', {'LFP','MUA','Reaction time','Detection'}, ...
        'col', {PAL.lfp, PAL.mua, PAL.rt, PAL.det}, ...
        'coh_file', { fullfile(coh_root,'lfp','all_loc_difflev','monkey_avg_results.mat'), ...
                      fullfile(coh_root,'mua','all_loc_difflev','monkey_avg_results.mat'), ...
                      fullfile(coh_root,'RT', 'all_loc_difflev','monkey_avg_results.mat'), ...
                      fullfile(corr_root,'hit_miss_itc','all_loc_difflev','monkey_avg_results_itc.mat')}, ...
        'coh_var', {'coh_monkey_avg','coh_monkey_avg','coh_monkey_avg','itc_monkey_avg'}, ...
        'coh_thr', {'thresh_monkey_avg','thresh_monkey_avg','thresh_monkey_avg','thresh_monkey_avg_itc'}, ...
        'coh_lab', {'Coherence','Coherence','Coherence','ITC'}, ...
        'corr_file',{ fullfile(corr_root,'lfp','all_loc_difflev','monkey_avg_results.mat'), ...
                      fullfile(corr_root,'mua','all_loc_difflev','monkey_avg_results.mat'), ...
                      fullfile(corr_root,'RT', 'all_loc_difflev','monkey_avg_results.mat'), ...
                      fullfile(corr_root,'hit_miss','all_loc_difflev','monkey_avg_results_pos.mat')}, ...
        'corr_var',{'corr_monkey_avg','corr_monkey_avg','corr_monkey_avg','pos_monkey_avg'}, ...
        'corr_thr',{'thresh_monkey_avg','thresh_monkey_avg','thresh_monkey_avg','thresh_monkey_avg_pos'}, ...
        'corr_lab',{'Correlation','Correlation','Correlation','POS'}, ...
        'reg_dv',  {'LFP_ERP_ampl_all','MUA_ERP_ampl_all','RT','hit_miss'});

    col_head = {'Phase coherence','Phase correlation',['Regression ' R2LAB]};

    f5 = new_fig(19, 24);
    M5   = [2.1 1.3 0.5 2.9 2.0 1.7];
    lets = panel_letters(numel(rows)*3);
    for r = 1:numel(rows)

        % column 1 — coherence / ITC
        k  = (r-1)*3 + 1;
        ax = axes('Position', grid_pos(f5, numel(rows), 3, k, M5)); panel_label(ax, lets{k}, FONT_SIZE_BIG);
        [v,t,fq] = load_curve(rows(r).coh_file, rows(r).coh_var, rows(r).coh_thr, 'freq');
        curve_panel(ax, fq, v, t, rows(r).col, PAL, rows(r).coh_lab);
        ylabel(ax, sprintf('%s\n%s', rows(r).lab, rows(r).coh_lab));
        if r==1, col_header(ax, col_head{1}, FONT_SIZE_BIG); end

        % column 2 — circular-linear correlation / POS
        k  = (r-1)*3 + 2;
        ax = axes('Position', grid_pos(f5, numel(rows), 3, k, M5)); panel_label(ax, lets{k}, FONT_SIZE_BIG);
        [v,t,fq] = load_curve(rows(r).corr_file, rows(r).corr_var, rows(r).corr_thr, 'freq');
        curve_panel(ax, fq, v, t, rows(r).col, PAL, rows(r).corr_lab);
        if r==1, col_header(ax, col_head{2}, FONT_SIZE_BIG); end

        % column 3 — partial R^2 from the multiple regression
        k  = (r-1)*3 + 3;
        ax = axes('Position', grid_pos(f5, numel(rows), 3, k, M5)); panel_label(ax, lets{k}, FONT_SIZE_BIG);
        rf = fullfile(reg_root, rows(r).reg_dv, 'monkey_avg_results.mat');
        v = []; t = NaN; fq = [];
        if isfile(rf)
            S = load(rf);
            if isfield(S,'monkey_avg_obs') && isfield(S.monkey_avg_obs,'phase')
                v = S.monkey_avg_obs.phase;
            end
            if isfield(S,'thresh_monkey') && isfield(S.thresh_monkey,'phase')
                t = S.thresh_monkey.phase;
            end
            if isfield(S,'freq'), fq = S.freq; end
        end
        if isempty(fq) && numel(FREQ_AXIS) == numel(v), fq = FREQ_AXIS; end
        curve_panel(ax, fq, v, t, rows(r).col, PAL, ['Partial ' R2LAB]);
        if r==1, col_header(ax, col_head{3}, FONT_SIZE_BIG); end
    end
    supertitle(f5, {'Monkey average: how much of each measure the pre-stimulus phase explains', ...
        'Dashed grey = permutation threshold (max-stat across frequency, p<0.05); shading and bar = significant'});
    save_fig(f5, out_dir, 'fig05_phase_measures');
end


%% =====================================================================
%% Figure 6 — WRITES fig06_hypotheses.pdf
%% H1-H4 hypothesis overlay, monkey average
%% =====================================================================
% Source: the same monkey_avg_results.mat files as Figure 5, but read for
% all four aggregation levels:
%   H1 complex           phase pooled over everything, then magnitude
%   H2 abs_per_pos       magnitude per stimulus position, then averaged
%   H3 abs_per_pos_diff  magnitude per position x difficulty
%   H4 abs_per_chan      magnitude per channel
%
% WHY H1 KEPT DISAPPEARING, AND WHAT IS DONE ABOUT IT
% By Jensen's inequality mean(|x|) >= |mean(x)|, so H1 is by construction
% the SMALLEST of the four. Autoscaling then squashes it against the axis
% floor, and because the original overlay drew H1 first, H2-H4 were painted
% on top of it wherever the curves coincide - for correlation, H4 is
% identically equal to H1 (circ_corrcl is already non-negative per channel),
% so the H1 line sat exactly underneath H4 and was invisible.
%
% Rather than guess in advance which pair might coincide, the panel MEASURES
% it: any set of hypotheses whose curves are equal to within rounding is
% drawn by draw_curves() as a single interleaved dashed line that alternates
% through their colours, so every hypothesis sharing that curve is visible
% and it is obvious that they are the same curve rather than one hiding the
% other. Curves that do not coincide are drawn solid as usual, and H1 is
% drawn last so it is never overpainted by a near neighbour.
%
% Significance is the PER-HYPOTHESIS test (each curve against its own
% permutation threshold), drawn as four bars in a strip RESERVED above the
% data so a bar can never sit on a trace. The paired H_n - H_(n-1) test that
% removes the Jensen advantage lives in compare_hypotheses.m and is not
% repeated here.
%
% The legend is labelled H1..H4 only and lives in the top-right of panel l -
% one legend for the whole grid, showing that row's rungs. What the levels
% mean is in the caption, because four spelled-out labels do not fit a panel.

if DO_FIG6
    fprintf('\n===== Figure 6: H1-H4 overlay =====\n');

    hyp = struct( ...
        'key',   {'complex','abs_per_pos','abs_per_pos_diff','abs_per_chan'}, ...
        'lab',   {'H1','H2','H3','H4'}, ...
        'lw',    {1.5, 1.5, 1.5, 1.5});
    draw_order = [4 3 2 1];      % H1 last, so it stays on top

    % H1..H4 are a nested chain of pooling levels applied to ONE measured
    % quantity, not four unrelated quantities, so within a row they are four
    % rungs of that row's own hue rather than four different hues. Colour
    % then still says "this row measures MUA" while lightness says "this is
    % the pooled level"; the reader learns one ladder instead of relearning
    % four hues per row. H1 is the darkest, and because H1 is the smallest
    % curve by Jensen's inequality the curves rise up the panel as they get
    % lighter - the lightness order and the vertical order agree, which makes
    % the Jensen advantage a property of the picture rather than an assertion
    % in the caption.

    dvs = struct( ...
        'lab',      {'LFP','MUA','Reaction time','Detection'}, ...
        'key',      {'lfp','mua','rt','det'}, ...
        'coh_sub',  {'lfp','mua','RT','hit_miss_itc'}, ...
        'coh_file', {'monkey_avg_results.mat','monkey_avg_results.mat', ...
                     'monkey_avg_results.mat','monkey_avg_results_itc.mat'}, ...
        'coh_var',  {'coh_monkey_avg','coh_monkey_avg','coh_monkey_avg','itc_monkey_avg'}, ...
        'coh_thr',  {'thresh_monkey_avg','thresh_monkey_avg','thresh_monkey_avg','thresh_monkey_avg_itc'}, ...
        'corr_sub', {'lfp','mua','RT','hit_miss'}, ...
        'corr_file',{'monkey_avg_results.mat','monkey_avg_results.mat', ...
                     'monkey_avg_results.mat','monkey_avg_results_pos.mat'}, ...
        'corr_var', {'corr_monkey_avg','corr_monkey_avg','corr_monkey_avg','pos_monkey_avg'}, ...
        'corr_thr', {'thresh_monkey_avg','thresh_monkey_avg','thresh_monkey_avg','thresh_monkey_avg_pos'}, ...
        'reg_dv',   {'LFP_ERP_ampl_all','MUA_ERP_ampl_all','RT','hit_miss'});

    pipes = {'Phase coherence','Phase correlation',['Regression ' R2LAB]};
    nD = numel(dvs); nP = numel(pipes);

    f6 = new_fig(19, 24);
    M6   = [2.1 1.3 0.6 3.4 1.9 1.7];
    LEG_PANEL = 12;     % panel l
    lets = panel_letters(nD*nP);
    leg_h = gobjects(1,numel(hyp)); leg_done = false;

    for d = 1:nD
        for p = 1:nP
            k  = (d-1)*nP + p;
            ax = axes('Position', grid_pos(f6, nD, nP, k, M6)); hold(ax,'on');
            panel_label(ax, lets{k}, FONT_SIZE_BIG);

            lad = PAL.lad.(dvs(d).key);      % 4x3, dark to light
            for h = 1:numel(hyp), hyp(h).col = lad(h,:); end

            V = cell(1,numel(hyp)); TH = nan(1,numel(hyp)); fq = [];
            for h = 1:numel(hyp)
                switch p
                    case 1
                        f = fullfile(comb,'phase_coherence',hyp(h).key,cp, ...
                                     dvs(d).coh_sub,'all_loc_difflev',dvs(d).coh_file);
                        if strcmp(dvs(d).coh_sub,'hit_miss_itc')
                            f = fullfile(comb,'phase_correlation',hyp(h).key,cp, ...
                                         dvs(d).coh_sub,'all_loc_difflev',dvs(d).coh_file);
                        end
                        [V{h}, TH(h), fqh] = load_curve(f, dvs(d).coh_var, dvs(d).coh_thr, 'freq');
                    case 2
                        f = fullfile(comb,'phase_correlation',hyp(h).key,cp, ...
                                     dvs(d).corr_sub,'all_loc_difflev',dvs(d).corr_file);
                        [V{h}, TH(h), fqh] = load_curve(f, dvs(d).corr_var, dvs(d).corr_thr, 'freq');
                    case 3
                        f = fullfile(comb,'multi_lin_reg',hyp(h).key,cp,dvs(d).reg_dv, ...
                                     'monkey_avg_results.mat');
                        [V{h}, TH(h), fqh] = load_reg(f);
                end
                if isempty(fq), fq = fqh; end
                % A partial R^2 outside [0,1], or a non-finite threshold, is a
                % failed fit, not a small effect. The H3 x regression x hit/miss
                % cell is exactly that: the per-position x difficulty strata leave
                % almost no residual variance, so the ratio explodes to +/-1e10 and
                % Inf. Drop the curve and say so rather than let it set the axis.
                if ~isempty(V{h}) && (any(~isfinite(V{h})) || ...
                        (p == 3 && (max(V{h}) > 1 || min(V{h}) < 0)) || ...
                        (~isfinite(TH(h)) && p == 3))
                    fprintf('  dropped %s / %s / %s: degenerate values (%.3g to %.3g, thr %.3g)\n', ...
                            hyp(h).lab, pipes{p}, dvs(d).lab, min(V{h}), max(V{h}), TH(h));
                    V{h} = []; TH(h) = NaN;
                end
            end

            allv = [V{:}];
            if isempty(allv), axis(ax,'off'); continue; end
            nfq = max(cellfun(@numel, V));
            if isempty(fq) && numel(FREQ_AXIS) == nfq, fq = FREQ_AXIS; end
            if isempty(fq), fq = 1:nfq; end

            % Reserve a strip ABOVE the data for the four significance rows,
            % so a bar can never sit on a trace. yl_dat is where the curves
            % live; the strip runs from just above it to the top of the axes.
            yl_dat = pad_lim(allv);
            rng0   = yl_dat(2) - yl_dat(1);
            % The significance strip is anchored to the TOP OF THE AXIS, at
            % the same fraction of panel height in every panel, so the bars
            % line up across the whole grid and can be compared by eye. It
            % used to be anchored to the top of the DATA, which put it lower
            % in the one panel that carries extra headroom for the legend.
            % headroom then has to clear the strip as well as the curves.
            extra  = 0.50;
            yl     = [yl_dat(1), yl_dat(2) + extra*rng0];
            R      = yl(2) - yl(1);
            strip  = [yl(2) - 0.30*R, yl(2) - 0.02*R];

            % Neither a significance BAND nor the per-hypothesis threshold
            % LINES are drawn here. Four overlapping bands turn the panel into
            % mud, and four dotted thresholds add four more horizontal lines
            % competing with the curves - in panels j and k they sat right
            % where the significance strip does, so the two were hard to tell
            % apart. The four stacked bars along the top carry the same
            % information, one hypothesis per row, at the same height in every
            % panel.
            leg_h = draw_curves(ax, fq, V, {hyp.col}, [hyp.lw], draw_order, leg_h);
            % significance rows inside the reserved strip, one per hypothesis
            for h = 1:numel(hyp)
                if isempty(V{h}) || ~isfinite(TH(h)), continue; end
                sig_ticks(ax, fq, V{h} >= TH(h), strip, h, numel(hyp), hyp(h).col);
            end

            xlim(ax,[min(fq) max(fq)]); ylim(ax, yl);
            if d==nD, xlabel(ax,'Frequency (Hz)'); end
            if p==1,  ylabel(ax, sprintf('%s\nMagnitude', dvs(d).lab)); end
            if d==1,  col_header(ax, pipes{p}, FONT_SIZE_BIG); end
            tidy(ax);

            % ONE legend for the whole grid, drawn in NEUTRAL GREYS.
            %
            % A legend can only name lines that live in its own axes, so a
            % coloured key would necessarily show one row's hue - rose, here -
            % next to eleven panels that use a different one. What every row
            % actually shares is the LIGHTNESS ordering, not the hue, so the
            % key is four dummy greys at the ladder's own lightness steps: it
            % says "dark to light is H1 to H4" and claims no hue, which is
            % true of all twelve panels.
            if k == LEG_PANEL && ~leg_done && all(arrayfun(@(x) isgraphics(x), leg_h))
                gl = neutral_ladder(numel(hyp));
                dh = gobjects(1, numel(hyp));
                for h = 1:numel(hyp)
                    dh(h) = plot(ax, NaN, NaN, '-', 'Color', gl(h,:));
                end
                % Vertical, in the TOP-RIGHT corner. Measured, the four-row
                % key is 1.6 cm against a 3.55 cm panel - 45% of the height -
                % so it cannot sit below the significance strip (another 30%)
                % without squeezing this panel's curves into the bottom
                % quarter. It does not have to: the strip spans the full
                % width but the BARS in this panel only reach ~15 Hz, so the
                % key shares that band from the right and touches nothing.
                % Worth re-checking if this panel ever gains significance at
                % high frequencies - the key would then cover those bars.
                lg = legend(dh, {hyp.lab}, 'Location','northeast');
                style_legend(lg);
                place_legend(ax, lg, 'right', 0.99);
                leg_done = true;
            end
        end
    end
    % Single line. The hue/ladder explanation is a methods sentence, not a
    % figure caption - the picture reads without it, and the block comment
    % above records why the key is neutral.
    supertitle(f6, 'H1 pooled, H2 per position, H3 per position x difficulty, H4 per channel — monkey average');
    save_fig(f6, out_dir, 'fig06_hypotheses');
end


%% =====================================================================
%% Figure 7a — WRITES fig07a_phase_progression.pdf
%% preferred phase versus stimulus position
%% =====================================================================
% Source: results_combined/scanning/phase_progression/cp10_till_100/<dv>/
%         monkey_avg_results.mat, with the frequency axis taken from the
%         matching per-animal channel_avg_results.mat (the monkey-average
%         file does not store one).
%
% Systematicity is the circular-linear correlation between a channel's
% preferred phase and the stimulus POSITION: high means the preferred phase
% marches with where the stimulus is, which is what a scanning process would
% produce. Channels are averaged within animal, animals then averaged.
%
% All four dependent variables are shown, so the frequencies at which the
% phase tracks position can be compared against the frequencies at which
% phase predicts behaviour.

if DO_FIG7A
    fprintf('\n===== Figure 7a: phase progression =====\n');

    pp_dvs = {'lfp','mua','RT','hit_miss'};
    pp_lab = {'LFP','MUA','Reaction time','Detection'};
    pp_key = {'lfp','mua','rt','det'};
    % Each panel takes the hue of the quantity it measures, the same hue that
    % quantity has in Figures 3-6, so a reader who learned the four hues there
    % carries them straight into this figure.
    pp_col = cellfun(@(k) PAL.(k), pp_key, 'uni', 0);

    f7a = new_fig(19, 13.8);
    M7a  = [2.1 1.3 0.6 3.0 2.2 2.0];
    lets = panel_letters(numel(pp_dvs));
    for k = 1:numel(pp_dvs)
        ax = axes('Position', grid_pos(f7a, 2, 2, k, M7a)); hold(ax,'on');
        panel_label(ax, lets{k}, FONT_SIZE_BIG);
        f = fullfile(comb,'scanning','phase_progression',cp,pp_dvs{k},'monkey_avg_results.mat');
        if ~isfile(f), axis(ax,'off'); title(ax,[pp_lab{k} ' — no data']); continue; end
        S = load(f);

        % the frequency axis lives with the per-animal files
        fq = [];
        for ia = 1:numel(S.animals)
            fa = fullfile(base,['results_' S.animals{ia}],'scanning','phase_progression', ...
                          cp,pp_dvs{k},'channel_avg_results.mat');
            if isfile(fa), Sa = load(fa,'freq'); fq = Sa.freq(:)'; break; end
        end
        if isempty(fq), fq = 1:numel(S.R_monkey_avg); end

        v  = S.R_monkey_avg(:)';
        th = S.thresh_monkey;
        mk = v >= th;
        yl = legend_headroom(pad_lim([v th]), 0.30);

        sig_band(ax, fq, mk, yl, band_tint(pp_col{k}));
        hc(1) = plot(ax, [min(fq) max(fq)], [th th], '--', 'Color', PAL.grey);
        hc(2) = plot(ax, fq, v, '-', 'Color', pp_col{k});
        sig_bar(ax, fq, mk, yl, pp_col{k});

        xlim(ax,[min(fq) max(fq)]); ylim(ax, yl);
        xlabel(ax,'Frequency (Hz)'); ylabel(ax,'Systematicity (circular-linear r)');
        % no count in the title: the shading already shows where, and how
        % much, without asking the reader to hold a number in mind
        title(ax, pp_lab{k}, 'FontWeight','normal');
        if k==1
            % only the threshold needs naming: the single solid curve IS the
            % monkey average, which the y label and caption already say
            lg = legend(hc(1), {'Max-stat threshold'}, ...
                        'Location','northeast');
            style_legend(lg);
            place_legend(ax, lg, 'right');
        end
        tidy(ax);
    end
    supertitle(f7a, {'Does the preferred phase track the stimulus position?', ...
        'Monkey average; dashed grey = max-stat threshold across frequency, shading = significant (p<0.05)'});
    save_fig(f7a, out_dir, 'fig07a_phase_progression');
end


%% =====================================================================
%% Figure 7b — WRITES fig07b_pgd.pdf
%% planar-wave existence: pooled z-PGD, plus an example position per animal
%% =====================================================================
% Source: results_combined/scanning/planar_wave_existence/cp10_till_100/lfp/
%         planar_wave_existence.mat, written by cortical_planar_wave_PGD.m
%         results_<animal>/scanning/phase_progression/.../phase_progression.mat
%         for the raw phase of the example maps, which that file does not
%         duplicate (its rawphi/rawcoh are dropped before saving).
%
% Panel a. Phase-gradient directionality says how close the preferred-phase
% map across the 8x8 array is to a PLANE. Each animal's PGD is first
% standardised against its own shuffled-map null (zPGD), which is what makes
% the two comparable, and the pooled curve is the mean of those z-scores.
% zG is formed here by averaging the per-animal z curves the file already
% contains - the same line cortical_planar_wave_PGD.m runs before it saves
% sigG. No statistic is recomputed: sigG, the significance actually drawn, is
% read as stored.
%
% Panels b, c. The clearest single case the claim rests on, one per animal,
% so the pooled spectrum and the map it is made of are read together. The
% position is chosen from the stored per-position statistics (PGDp,
% peak_pos_fdr) - the highest-PGD position that survives BH-FDR - never by
% eye, and only at the peak frequency where something survives.

if DO_FIG7B
    fprintf('\n===== Figure 7b: pooled z-PGD =====\n');

    R = load(fullfile(comb,'scanning','planar_wave_existence',cp,dv_scan, ...
                      'planar_wave_existence.mat'));
    A     = R.results.A;
    valid = find(arrayfun(@(a) ~isempty(a.PGD), A));
    fr    = R.results.freq(:)';

    zG = zeros(numel(fr),1); zN = zeros(numel(fr), size(A(valid(1)).zPGD_null,2));
    for ia = valid
        zG = zG + A(ia).zPGD;
        zN = zN + A(ia).zPGD_null;
    end
    zG = (zG / numel(valid))';
    zN =  zN / numel(valid);
    zthr = quantile(zN, 0.95, 2)';        % per-frequency 95th percentile of the null
    sigG = R.results.sigG(:)';

    % Spectrum on top spanning the width, the two example maps beneath it,
    % so the claim and the clearest single case it rests on are read together
    % rather than on separate pages.
    f7b  = new_fig(15, 15);
    M7b  = [2.0 1.2 2.8 2.4 1.8 2.4];
    % b and c are SQUARE, so their boxes end up narrower than the grid cells
    % they sit in. Panel a is therefore spanned from b's left edge to c's
    % right edge, not across the raw cells - otherwise the top panel is wider
    % than the row beneath it and the three do not line up.
    pB   = square_pos(f7b, grid_pos(f7b, 2, 2, 3, M7b));
    pC   = square_pos(f7b, grid_pos(f7b, 2, 2, 4, M7b));
    posA = grid_pos(f7b, 2, 2, 1, M7b);
    lets = panel_letters(3);

    ax = axes('Position', [pB(1) posA(2) pC(1)+pC(3)-pB(1) posA(4)]);
    hold(ax,'on');
    panel_label(ax, lets{1}, FONT_SIZE_BIG);
    yl = pad_lim([zG 0]);

    % the scanning analyses are all LFP (dv_scan), so this takes the LFP hue
    sig_band(ax, fr, sigG, yl, band_tint(PAL.lfp));
    plot(ax, [min(fr) max(fr)], [0 0], '-', 'Color', PAL.grey_lt, 'LineWidth', LW_GUIDE);
    plot(ax, fr, zG, '-', 'Color', PAL.lfp);
    sig_bar(ax, fr, sigG, yl, PAL.lfp);

    xlim(ax,[min(fr) max(fr)]); ylim(ax, yl);
    xlabel(ax,'Frequency (Hz)');
    ylabel(ax,'Pooled z-PGD');
    title(ax, 'Pooled planar-wave evidence', 'FontWeight','normal');
    % no legend: the y label names the curve
    tidy(ax);

    % ── b, c: the clearest FDR-surviving position in each animal ──────
    % Chosen from the stored per-position statistics (PGDp, peak_pos_fdr),
    % not by eye: the highest-PGD position that survives BH-FDR at the peak
    % frequency. Only the peak where something survives is shown - at 13.3 Hz
    % nothing does in either animal, and the caption says so rather than the
    % figure carrying two empty boxes.
    grid_rows = 8; grid_cols = 8;
    ch_col = ceil((1:(grid_rows*grid_cols))' / grid_rows);
    ch_row = grid_rows - mod((1:(grid_rows*grid_cols))' - 1, grid_rows);
    cmap_ph = cmap_cyclic(256, MAP_FAMILY);

    % Which position to show, by animal. Empty = the highest-PGD position
    % that survives FDR; a number pins a specific one. A pinned position is
    % still CHECKED against peak_pos_fdr - if it does not survive, the panel
    % says so rather than showing it as though it did.
    EX_POS = struct('hermes', 9, 'klecks', []);

    nPkAll = numel(A(valid(1)).peak_freq);
    jbest  = 0;
    for j = 1:nPkAll
        if any(arrayfun(@(kk) any(A(valid(kk)).peak_pos_fdr(j,:)), 1:numel(valid)))
            jbest = j; break
        end
    end
    ex_freq = NaN; dropped = [];
    if jbest > 0
        ex_freq = A(valid(1)).peak_freq(jbest);
        dropped = A(valid(1)).peak_freq(setdiff(1:nPkAll, jbest));
    end

    for k = 1:numel(valid)
        ia = valid(k);
        if k==1, pk = pB; else, pk = pC; end
        ax = axes('Position', pk); hold(ax,'on');
        panel_label(ax, lets{1+k}, FONT_SIZE_BIG, 1.05);
        if jbest == 0, axis(ax,'off'); continue; end

        pf = fullfile(base,['results_' A(ia).animal],'scanning','phase_progression', ...
                      cp,dv_scan,'phase_progression.mat');
        if ~isfile(pf), axis(ax,'off'); continue; end
        S    = load(pf,'pref_phase','coh_mag','freq');
        pref = S.pref_phase(:,:,A(ia).pos_keep);
        cohm = S.coh_mag(:,:,A(ia).pos_keep);
        frq  = S.freq(:)';
        [~, fi] = min(abs(frq - ex_freq));

        pgd  = A(ia).PGDp(fi,:);
        pass = A(ia).peak_pos_fdr(jbest,:);
        want = [];
        if isfield(EX_POS, A(ia).animal), want = EX_POS.(A(ia).animal); end
        if ~isempty(want)
            ip = find(A(ia).pos_keep == want, 1);
            if isempty(ip)
                error('figures:exPos','%s: position %d is not in pos_keep', ...
                      A(ia).animal, want);
            end
            best = pgd(ip);
            if ~pass(ip)
                fprintf(['  NOTE %s: requested position %d does NOT survive ' ...
                         'FDR at %.1f Hz (PGD %.2f)\n'], A(ia).animal, want, ...
                        frq(fi), best);
            end
        else
            pgd(~pass) = -Inf;
            [best, ip] = max(pgd);
        end
        if ~isfinite(best)
            axis(ax,'off');
            text(0.5,0.5,sprintf('%s: none survives FDR', capitalise(A(ia).animal)), ...
                 'Parent',ax,'HorizontalAlignment','center');
            continue
        end

        gp = build_phase_grid(pref(:,fi,ip), cohm(:,fi,ip), ch_row, ch_col, ...
                              grid_rows, grid_cols, A(ia).coh_sig(:,fi));
        vec_heatmap(ax, 1:grid_cols, 1:grid_rows, gp, cmap_ph, [-pi pi]);
        set(ax,'YDir','reverse','XTick',[],'YTick',[],'Box','on');
        axis(ax,[0.5 grid_cols+0.5 0.5 grid_rows+0.5]);
        dd = A(ia).DIRp(fi,ip);
        if ~isnan(dd)
            quiver(ax, grid_cols/2, grid_rows/2, 2*cos(dd), 2*sin(dd), 0, ...
                   'Color', [0 0 0], 'LineWidth', 2.6, 'MaxHeadSize', 2);
        end
        % the frequency belongs on the panel, not only in the caption: these
        % maps are one frequency out of 35 and the reader should not have to
        % go looking for which
        title(ax, sprintf('%s, %.1f Hz\nPosition %d, PGD %.2f', ...
              capitalise(A(ia).animal), frq(fi), A(ia).pos_keep(ip), best), ...
              'FontWeight','normal');
        tidy(ax);
        if k==numel(valid)
            vec_colorbar(ax, cmap_ph, [-pi pi], 'Preferred phase', false, PI_TICKS);
        end
        fprintf('  example: %s @ %.1f Hz, position %d, PGD %.2f\n', ...
                A(ia).animal, frq(fi), A(ia).pos_keep(ip), best);
    end

    fprintf('  %d animals, significant at %s Hz\n', numel(valid), ...
            ternary(any(sigG), strjoin(compose('%.1f', fr(sigG)), ', '), '(none)'));
    cap = 'a: cluster-corrected significant';
    if jbest > 0
        cap = sprintf('%s.  b, c: PGD at an example position, %.1f Hz', cap, ex_freq);
    end
    % The frequencies where NO position survives FDR are no longer named in
    % the caption - they are still printed to the console below, and figure
    % 7d shows every position at both peaks, so the negative result is on
    % record; it just belongs in the text rather than in this caption.
    if ~isempty(dropped)
        fprintf('  no position survives FDR at %s Hz in either animal\n', ...
                strjoin(compose('%.1f', dropped), ', '));
    end
    supertitle(f7b, {'Is the phase map a plane, and what does one position look like?', cap});
    save_fig(f7b, out_dir, 'fig07b_pgd');
end


%% =====================================================================
%% Figure 7c — WRITES fig07c_derotation.pdf
%% frequency x speed de-rotation grids
%% =====================================================================
% Source: results_combined/scanning/planar_wave_derotation/.../planar_wave_derotation.mat
%         results_combined/scanning/stimulus_loc_wave/.../stimulus_loc_wave.mat
%
% A plane in the phase map says the phases are ORDERED in space; it does not
% say the pattern travels. The de-rotation test asks the propagation question
% directly: undo the delay a wave of speed v would produce, and see whether
% the map becomes more coherent than it was. The plotted quantity is that
% GAIN, coherence after de-rotation minus coherence at v = infinity, on a
% frequency x speed grid.
%
% Top row    ACROSS CORTEX, speeds in m/s along the array.
% Bottom row ALONG THE STIMULUS AXIS, speeds in deg/s in visual space -
%            i.e. a wave that follows the stimulus, not the cortical sheet.
% Columns are the two estimators, which differ in what they treat as the
% per-location quantity:
%   phase      resultant of the per-location preferred phases (NOT a coherence)
%   coherence  the complex coherence itself
% Both are shown because they can disagree, and a propagation claim should
% survive either.
%
% Pooled across animals with the pooled-z test; the outline is the corrected
% significant region.

if DO_FIG7C
    fprintf('\n===== Figure 7c: de-rotation grids =====\n');

    D = load(fullfile(comb,'scanning','planar_wave_derotation',cp,dv_scan, ...
                      'planar_wave_derotation.mat'));
    Sv = load(fullfile(comb,'scanning','stimulus_loc_wave',cp,dv_scan, ...
                      'stimulus_loc_wave.mat'));
    % Same purple-pink-yellow family as the TFR, with the centre pinned at
    % zero by the two-slope limits below: yellow = de-rotation helps, purple
    % = it hurts, pink = it changes nothing.
    cmap_gain = cmap_ppy(256, MAP_FAMILY);

    est   = D.results.ESTIMATORS;
    panels = {};
    for ie = 1:numel(est)
        C = D.results.C.(est{ie});
        panels{end+1} = struct('G', C.gain, 'sig', C.sig_pool, ...
            'f', C.fHz(:)', 'v', C.speeds(:)', 'vunit','m/s', ...
            'ttl', sprintf('Across cortex — %s', est{ie})); %#ok<SAGROW>
    end
    for ie = 1:numel(est)
        C = Sv.results.C.(est{ie}).visual;
        panels{end+1} = struct('G', C.gain, 'sig', C.sig_pool_G, ...
            'f', C.fHz(:)', 'v', C.speeds(:)', 'vunit','deg/s', ...
            'ttl', sprintf('Along stimulus axis — %s', est{ie})); %#ok<SAGROW>
    end

    f7c = new_fig(19, 16.8);
    M7c  = [2.0 1.4 2.2 3.0 3.4 2.2];
    lets = panel_letters(numel(panels));
    for k = 1:numel(panels)
        Pk = panels{k};
        ax = axes('Position', grid_pos(f7c, 2, 2, k, M7c)); hold(ax,'on');
        panel_label(ax, lets{k}, FONT_SIZE_BIG);
        g  = Pk.G(isfinite(Pk.G));
        lo = min(quantile(g, 0.02),  -eps);
        hi = max(quantile(g, 0.998),  eps);
        cl = [lo 0 hi];        % two-slope, pinned at "de-rotation changes nothing"
        % frequency on x, speed on y: the grid is stored [freq x speed], so
        % it is transposed to match, and so are the significance mask below
        % and the axis labels
        vec_heatmap(ax, Pk.f, log10(Pk.v), Pk.G', cmap_gain, cl);
        sig_outline(ax, Pk.f, log10(Pk.v), Pk.sig', [1 1 1]);
        xlim(ax, full_lims(Pk.f)); ylim(ax, full_lims(log10(Pk.v)));
        set(ax,'Layer','top','Box','on');
        log_yticks(ax, Pk.v);
        xlabel(ax,'Frequency (Hz)');
        ylabel(ax, sprintf('Speed (%s, log scale)', Pk.vunit));
        title(ax, Pk.ttl, 'FontWeight','normal');
        tidy(ax);
        % gain runs from ~-0.9 to ~+0.01, so plain decimals need four or five
        % places to say anything; powers of ten read at a glance
        vec_colorbar(ax, cmap_gain, cl, 'Gain', true, sci_labels(cl));
    end
    supertitle(f7c, {'De-rotation: does undoing a travelling delay make the map more coherent?', ...
        'Colour = gain over the no-travel (v = infinity) fit; white outline = pooled-z significant across animals'});
    save_fig(f7c, out_dir, 'fig07c_derotation');
end


%% =====================================================================
%% Figure 7d — WRITES fig07d_per_position_<animal>.pdf
%% per-position phase maps at the LFP-coherence peaks
%% =====================================================================
% Source: results_<animal>/scanning/phase_progression/.../phase_progression.mat
%         (pref_phase, coh_mag) combined with the channel mask, kept positions
%         and per-position statistics stored in planar_wave_existence.mat.
%         The raw phase arrays are deliberately not duplicated in that file,
%         so they are re-read here from phase_progression.mat.
%
% Why per position at all: the collapsed map vector-averages over positions,
% which cancels any wave whose direction depends on where the stimulus is.
% Each position is therefore tested on its own, at the two frequencies fixed
% a priori by the LFP coherence spectrum, and corrected across positions with
% Benjamini-Hochberg FDR (cluster correction across frequency is unavailable
% once the frequency is pinned).
%
% One page per animal: rows are the two peak frequencies, columns the
% stimulus positions. The arrow is the fitted propagation direction, heavy
% where that position's plane survives FDR.

if DO_FIG7D
    fprintf('\n===== Figure 7d: per-position phase maps =====\n');

    R = load(fullfile(comb,'scanning','planar_wave_existence',cp,dv_scan, ...
                      'planar_wave_existence.mat'));
    A     = R.results.A;
    valid = find(arrayfun(@(a) ~isempty(a.PGD), A));
    grid_rows = 8; grid_cols = 8;
    ch_col = ceil((1:(grid_rows*grid_cols))' / grid_rows);
    ch_row = grid_rows - mod((1:(grid_rows*grid_cols))' - 1, grid_rows);
    cmap_ph = cmap_cyclic(256, MAP_FAMILY);

    for k = 1:numel(valid)
        ia = valid(k);
        an = A(ia).animal;

        pf = fullfile(base,['results_' an],'scanning','phase_progression',cp,dv_scan, ...
                      'phase_progression.mat');
        if ~isfile(pf), fprintf('  %s: no phase_progression.mat, skipped\n', an); continue; end
        S = load(pf,'pref_phase','coh_mag','freq');
        % planar_wave_existence dropped under-sampled positions; index the
        % raw arrays with the SAME pos_keep so the columns line up
        pref = S.pref_phase(:,:,A(ia).pos_keep);
        cohm = S.coh_mag(:,:,A(ia).pos_keep);
        fr   = S.freq(:)';
        nPos = numel(A(ia).pos_keep);
        pk   = A(ia).peak_freq(:)';
        nPk  = numel(pk);

        % Wrap the positions over two rows per frequency ONLY when there are
        % enough of them to make the page an unusable shape otherwise. Hermes
        % has 15 positions, which unwrapped runs past 35 cm wide; klecks has
        % 9, which fits one row comfortably and reads better that way - the
        % positions are an ordered sequence, and a single row shows that
        % order directly instead of breaking it across a line.
        WRAP_ABOVE = 10;
        if nPos > WRAP_ABOVE
            nSub = 2;  nCol = ceil(nPos/nSub);
        else
            nSub = 1;  nCol = nPos;
        end
        nRow = nSub*nPk;
        M7d  = [1.6 0.8 2.6 3.3 0.35 1.4];      % l b r t hgap vgap
        fw   = max(15, 2.1*nCol + 3.8);
        % Height derived from the COLUMN WIDTH, not a fixed per-row figure:
        % these panels are square, so a row shorter than a column wide makes
        % square_axes shrink every map to the row height and the array maps
        % come out needlessly small. Sizing the row to the cell width lets
        % each map fill its cell.
        cellw = (fw - M7d(1) - M7d(3) - (nCol-1)*M7d(5)) / nCol;
        fd    = new_fig(fw, M7d(2) + M7d(4) + nRow*cellw + (nRow-1)*M7d(6));
        for j = 1:nPk
            [~, fi] = min(abs(fr - pk(j)));
            dirs = A(ia).DIRp(fi,:);
            for p = 1:nPos
                sub = ceil(p/nCol);                     % which sub-row
                col = p - (sub-1)*nCol;
                kk  = ((j-1)*nSub + sub - 1)*nCol + col;
                ax = axes('Position', grid_pos(fd, nRow, nCol, kk, M7d)); hold(ax,'on');
                gp = build_phase_grid(pref(:,fi,p), cohm(:,fi,p), ch_row, ch_col, ...
                                      grid_rows, grid_cols, A(ia).coh_sig(:,fi));
                square_axes(ax);
                vec_heatmap(ax, 1:grid_cols, 1:grid_rows, gp, cmap_ph, [-pi pi]);
                set(ax,'YDir','reverse','XTick',[],'YTick',[],'Box','on');
                axis(ax,[0.5 grid_cols+0.5 0.5 grid_rows+0.5]);
                passed = A(ia).peak_pos_fdr(j,p);
                if ~isnan(dirs(p))
                    if passed, acol = PAL.ink; alw = 2.0; else, acol = PAL.grey_mid; alw = 1.0; end
                    quiver(ax, grid_cols/2, grid_rows/2, ...
                           2*cos(dirs(p)), 2*sin(dirs(p)), 0, ...
                           'Color', acol, 'LineWidth', alw, 'MaxHeadSize', 2);
                end
                title(ax, sprintf('Position %d%s\nPGD %.2f', A(ia).pos_keep(p), ...
                      ternary(passed,' *',''), A(ia).PGDp(fi,p)), ...
                      'FontWeight', ternary(passed,'bold','normal'));
                if col==1
                    ylabel(ax, sprintf('%.1f Hz', fr(fi)), ...
                           'FontWeight','bold','Visible','on');
                end
                if p==1
                    % one letter per frequency block rather than per position:
                    % a letter on every near-identical array map would be
                    % noise, and the block is what the text refers to
                    rowlets = panel_letters(nPk);
                    panel_label(ax, rowlets{j}, FONT_SIZE_BIG, 1.05);
                end
                tidy(ax);
                if col==nCol && sub==1
                    vec_colorbar(ax, cmap_ph, [-pi pi], 'Preferred phase', false, PI_TICKS);
                end
            end
        end
        supertitle(fd, {sprintf('%s — preferred-phase map per stimulus position', capitalise(an)), ...
            '* and heavy arrow = plane survives BH-FDR q=0.05; blank cell = channel not significant'});
        save_fig(fd, out_dir, sprintf('fig07d_per_position_%s', an));
    end
end


fprintf('\nAll figures written to %s\n', out_dir);


%% =====================================================================
%% Local functions
%% =====================================================================

function PAL = palette()
% ONE HUE PER MEASURED QUANTITY, lightness for levels within it.
%
% The rule this replaces: four hues were shared out among four different
% dimensions at once, so purple meant miss, MUA, H1, the pooled z-PGD curve
% and the dark end of every heatmap, and teal meant hit, LFP, H2 AND
% "significant" - which in Figure 5 put a teal LFP curve inside teal
% significance shading, where the colour no longer said whether it was naming
% the quantity or the test result.
%
%   R1  hue names the measured quantity and nothing else
%       purple = MUA, teal = LFP, orange = reaction time, rose = detection
%   R2  lightness names a level within that quantity
%       hit vs miss, and H1..H4, are levels of one thing, not four things
%   R3  colour never encodes a statistic
%       grey stays the null/threshold; significance is a tint of the curve's
%       OWN hue, so the shading names the curve it belongs to
%
% The ladder rungs sit at equal CIELAB L* targets ACROSS hues (30/48/64/80,
% and 40/70 for the hit-miss pair), not at equal blend amounts. A fixed tint
% of teal (base L* 49) and of orange (base L* 61) would land on different
% lightnesses and the reader would have to relearn the ladder in every row;
% at equal L*, rung 2 of purple and rung 2 of orange are equally light and
% the ladder is learned once. Adjacent rungs are dE 16-39 in CIELAB and
% 15-35 under simulated deuteranopia - a lightness ramp is what survives
% colour-vision deficiency. Weakest step: teal rung 2->3 at dE 16.1.
%
% The heatmap ramps are deliberately outside these four hues (see cmap_ppy),
% which is what keeps the map family and the curve family apart.

% base hues - the quantity
PAL.mua = hx('#5C368C');
PAL.lfp = hx('#008287');
PAL.rt  = hx('#E6731A');
PAL.det = hx('#B83D70');

% four-rung ladders, dark to light: H1 H2 H3 H4
PAL.lad.mua = [hx('#583485'); hx('#8265A7'); hx('#A793C1'); hx('#CDC2DC')];
PAL.lad.lfp = [hx('#124E55'); hx('#017F84'); hx('#4EA8AC'); hx('#9ECFD1')];
PAL.lad.rt  = [hx('#653D24'); hx('#AF5C1E'); hx('#E8802F'); hx('#F3BB90')];
PAL.lad.det = [hx('#72304F'); hx('#BC4979'); hx('#D284A4'); hx('#E6BBCD')];

% hit/miss pairs, hit is the DARK rung because it is the larger response
PAL.pair.lfp = [hx('#09696F'); hx('#6DB7BA')];
PAL.pair.mua = [hx('#704E9A'); hx('#B5A4CB')];

% neutrals - never data
PAL.ink      = hx('#212129');   % structure, zero-crossings, arrows
PAL.grey     = hx('#7A7A80');   % threshold / null floor
PAL.grey_mid = hx('#9E9EA1');   % a direction that fails FDR
PAL.grey_lt  = hx('#CCCCD1');   % zero lines, guides

% aliases, so older call sites keep working
PAL.purple = PAL.mua;
PAL.teal   = PAL.lfp;
PAL.orange = PAL.rt;
PAL.rose   = PAL.det;
end


function c = hx(h)
% '#RRGGBB' -> [r g b] in 0..1. The palette is written in hex because that is
% the form it was specified and checked in.
c = double(sscanf(h(2:end), '%2x%2x%2x'))' / 255;
end


function c = band_tint(col)
% The ONE tint used for every significance band in the set. Having a single
% constant here is the point: the same shading weight in every figure means a
% shaded region always reads as the same thing, and a colour that looks
% different is then actually a different colour rather than a different tint
% of the same one.
c = tint(col, 0.86);
end


function c = tint(col, amount)
% Blend a colour toward white. This is how every "translucent" fill in this
% script is made: an OPAQUE pre-blended colour, because FaceAlpha would push
% the whole figure onto the OpenGL renderer and rasterise the PDF.
% amount = 0 is the colour itself, 1 is white.
c = col(:)' * (1-amount) + [1 1 1] * amount;
end


function f = new_fig(w_cm, h_cm)
f = figure('Visible','off','Color','w','Units','centimeters', ...
           'Position',[1 1 w_cm h_cm],'Renderer','painters', ...
           'InvertHardcopy','off');
end


function save_fig(fig, out_dir, name)
% Vector PDF with live text. -painters is what keeps the page editable in
% Illustrator; the paper size is matched to the figure so nothing is cropped.
%
% No layout fixing here any more. Every figure lays its panels out with
% grid_pos(), whose margins are given in CENTIMETRES and already reserve
% room for the caption, the axis labels and the panel letters. Squeezing
% the axes afterwards, as this used to do, fought with those margins and
% stranded anything positioned by hand.
p = get(fig,'Position');
set(fig,'PaperUnits','centimeters','PaperSize',p(3:4), ...
        'PaperPosition',[0 0 p(3:4)],'PaperPositionMode','manual');
align_ylabels(fig);
f = fullfile(out_dir,[name '.pdf']);
% '-painters', not '-vector': R2024b's Code Analyzer prefers '-vector', but
% that flag only exists from R2024a, and this script has to keep working on
% the older MATLAB the rest of the pipeline runs under. '-painters' produces
% the same vector PDF on every release.
print(fig, f, '-dpdf', '-painters');
fprintf('  saved %s\n', f);
close(fig);
end


function align_ylabels(fig)
% Centre every y label on its axis, and put all the labels in a COLUMN at the
% same distance from the plot.
%
% MATLAB places a y label just clear of whatever the tick labels happen to
% need, so a panel ticked "0.05" and one ticked "0.002" end up with their
% labels at different distances from the axis and the column reads as ragged.
% This finds the leftmost label in each column and moves the rest out to
% match it, and pins each one to the vertical middle of its own axis.
%
% Colourbar axes are skipped - their label is on the right by construction
% and belongs to the strip, not to a column of panels.
drawnow limitrate
ax = findall(fig,'Type','axes');
keep = false(size(ax));
for k = 1:numel(ax)
    if strcmp(get(ax(k),'YAxisLocation'),'right'), continue; end
    h = get(ax(k),'YLabel');
    if isempty(h) || isempty(get(h,'String')), continue; end
    keep(k) = true;
end
ax = ax(keep);
if isempty(ax), return; end

pos = zeros(numel(ax),4);
for k = 1:numel(ax), pos(k,:) = get(ax(k),'Position'); end
colkey = round(pos(:,1)*200);            % same left edge = same column

for c = unique(colkey)'
    g  = find(colkey == c);
    xf = nan(numel(g),1);
    for i = 1:numel(g)
        h = get(ax(g(i)),'YLabel');
        set(h,'Units','normalized');
        q = get(h,'Position');
        set(h,'Position',[q(1) 0.5 0]);                 % centred on the axis
        xf(i) = pos(g(i),1) + q(1)*pos(g(i),3);         % in figure units
    end
    xmin = min(xf);
    for i = 1:numel(g)
        h = get(ax(g(i)),'YLabel');
        set(h,'Position', [(xmin - pos(g(i),1))/pos(g(i),3), 0.5, 0]);
    end
end
end


function tidy(ax)
set(ax,'Box','off','TickDir','out','Layer','top','Color','none', ...
       'XColor',[0.25 0.25 0.28],'YColor',[0.25 0.25 0.28]);
% Kill the "x10^-3" exponent label. MATLAB draws it at 0.8x the axes font,
% which is a third font size in the figure however carefully everything else
% is pinned, and it is the only text on the page whose size is not ours to
% choose. Printing the tick values in full costs a couple of characters and
% removes the problem.
try
    ax.XAxis.Exponent = 0;
    ax.YAxis.Exponent = 0;
catch
end
end


function style_legend(lg, boxed)
% Legends do NOT honour DefaultLegendFontSize - they fall back to 0.9x the
% axes font - so the size is set here, on every legend, from the same
% default everything else inherits.
if nargin < 2, boxed = true; end
set(lg, 'FontSize', get(groot,'DefaultAxesFontSize'), ...
        'FontName', get(groot,'DefaultAxesFontName'));
if boxed
    set(lg, 'Box','on', 'EdgeColor','none', 'Color','w');
else
    set(lg, 'Box','off');
end
end


function pos = grid_pos(fig, nrow, ncol, k, mcm)
% Position of panel k in an nrow x ncol grid, with every margin and gap given
% in CENTIMETRES rather than normalised units.
%
% Centimetres because the figures in this set range from 7 cm to 22 cm tall
% and 19 to 34 cm wide: a normalised gap that looks right on one page is
% cramped on another, whereas "1.8 cm between columns" is 1.8 cm everywhere.
% The margins are what reserve room for the caption, the axis labels and the
% panel letters, so nothing has to be squeezed afterwards.
%
% mcm = [left bottom right top hgap vgap]. Panels are numbered in reading
% order, left to right then top to bottom.
fp = get(fig,'Position');
Wf = fp(3); Hf = fp(4);
L = mcm(1)/Wf; B = mcm(2)/Hf; R = mcm(3)/Wf; T = mcm(4)/Hf;
HG = mcm(5)/Wf; VG = mcm(6)/Hf;
w = (1 - L - R - (ncol-1)*HG) / ncol;
h = (1 - B - T - (nrow-1)*VG) / nrow;
r = ceil(k/ncol); c = k - (r-1)*ncol;
pos = [L + (c-1)*(w+HG), 1 - T - r*h - (r-1)*VG, w, h];
end


function square_axes(ax)
% Make the axes BOX square by trimming its Position, instead of calling
% axis('equal') and letting MATLAB shrink the drawn square inside a wider
% box. The difference matters twice over: panels whose boxes were shrunk by
% different amounts (say, one of them by a colourbar) end up drawing squares
% of DIFFERENT sizes, and normalised text - the panel letter - is positioned
% against the box, so it drifts away from a square that no longer fills it.
fp  = get(get(ax,'Parent'),'Position');
pos = get(ax,'Position');
w   = pos(3)*fp(3); h = pos(4)*fp(4);          % cm
side = min(w, h);
set(ax,'Position',[pos(1) + (w-side)/2/fp(3), pos(2) + (h-side)/2/fp(4), ...
                   side/fp(3), side/fp(4)]);
end


function panel_label(ax, txt, fs, up_cm)
% The a/b/c label above the panel's top-left corner, for the manuscript to
% refer to.
%
% The offset is in CENTIMETRES, converted to this axes' normalised units at
% the end. It used to be a fixed FRACTION of the axes width, which put the
% letter of a full-width panel far further left than the letter of a narrow
% panel beneath it - so the letters down a figure did not line up, and on a
% panel spanning two columns the letter drifted almost to the page edge.
%
% up_cm is how far above the axes top the letter sits: the default clears a
% one-line title, a TWO-line title needs about 1 cm.
if nargin < 4 || isempty(up_cm), up_cm = 0.55; end
LEFT_CM = 0.95;
fp  = get(get(ax,'Parent'),'Position');     % figure size, cm
pos = get(ax,'Position');
text(ax, -(LEFT_CM/fp(3))/pos(3), 1 + (up_cm/fp(4))/pos(4), txt, ...
     'Units','normalized', 'FontWeight','bold', ...
     'FontSize', fs, 'Color', [0.13 0.13 0.16], ...
     'HorizontalAlignment','left', 'VerticalAlignment','middle');
end


function pos = square_pos(fig, pos)
% The square, centred box that square_axes() would produce - returned WITHOUT
% needing an axes, so a spanning panel can be sized from the squares that
% will sit beneath it.
fp   = get(fig,'Position');
w    = pos(3)*fp(3); h = pos(4)*fp(4);
side = min(w, h);
pos  = [pos(1) + (w-side)/2/fp(3), pos(2) + (h-side)/2/fp(4), ...
        side/fp(3), side/fp(4)];
end


function letters = panel_letters(n)
a = 'abcdefghijklmnopqrstuvwxyz';
letters = cellstr(a(1:n)');
end


function out = capitalise(name)
out = name;
if ~isempty(out), out(1) = upper(out(1)); end
end


function supertitle(fig, txt, pad_cm)
% The page heading, centred on the PANEL BLOCK and held clear of the page top.
%
% sgtitle centres on the figure WIDTH, which is not the same thing: the left
% margin of these figures is wider than the right because it carries the y
% label, the tick labels and the panel letter, so a heading centred on the
% page sits noticeably left of the panels it describes. This measures where
% the axes actually are and centres on that. Colourbar axes are excluded -
% they are an annex to one panel, not part of the block.
%
% pad_cm is the clear space between the top of the page and the top of the
% heading; the gap between heading and panels is then set by that figure's
% top margin in grid_pos.
if nargin < 3 || isempty(pad_cm), pad_cm = 0.8; end
fp = get(fig,'Position');
fs = get(groot,'DefaultAxesFontSize');

ax = findall(fig,'Type','axes');
x0 = inf; x1 = -inf;
for k = 1:numel(ax)
    if strcmp(get(ax(k),'YAxisLocation'),'right'), continue; end
    q  = get(ax(k),'Position');
    x0 = min(x0, q(1));
    x1 = max(x1, q(1)+q(3));
end
if ~isfinite(x0), x0 = 0.05; x1 = 0.95; end

% Widest box that stays CENTRED ON THE BLOCK without running off the page.
% Sizing it to the block alone made long captions wrap into the panels on the
% figures with wide margins; sizing it to the page would re-centre it.
EDGE  = 0.5 / fp(3);
ctr   = (x0 + x1) / 2;
halfw = min(ctr - EDGE, 1 - EDGE - ctr);

% How many lines this will actually take, so the box is tall enough. A
% caption that wraps and a box sized for the unwrapped line count is how the
% text ends up sitting on the first row of panels.
cpl   = max(1, floor(2*halfw*fp(3) / (0.5*fs/72*2.54)));
lines = txt; if ~iscell(lines), lines = {lines}; end
nline = 0;
for k = 1:numel(lines), nline = nline + max(1, ceil(numel(lines{k})/cpl)); end

h = (nline * fs * 1.5 / 72 * 2.54) / fp(4);        % text block height, fig units
y = 1 - pad_cm/fp(4) - h;
annotation(fig, 'textbox', [ctr-halfw, y, 2*halfw, h], 'String', txt, ...
    'EdgeColor','none', 'FitBoxToText','off', 'Margin', 0, ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
    'FontName', get(groot,'DefaultAxesFontName'), 'FontSize', fs, ...
    'Color', [0.25 0.25 0.28]);
end


function col_header(ax, txt, fs)
% Column heading placed above the axes rather than as a title. MATLAB parks
% the axis exponent at the top-left INSIDE the axes box, which a normal title
% runs straight into; this sits clear above both, on the same line as the
% panel letter.
text(0.5, 1.13, txt, 'Parent',ax, 'Units','normalized', ...
     'HorizontalAlignment','center', 'FontWeight','bold', 'FontSize', fs, ...
     'Color',[0.13 0.13 0.16]);
end


function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end


function place_legend(ax, lg, side, topfrac)
% Park a legend in the headroom strip but BELOW the significance bar. The bar
% is drawn in the top 3% of the panel, and MATLAB's 'north...' inset is not
% always deep enough to clear it on a short panel, so the position is set
% explicitly: the legend's top edge lands at 88% of the axes height.
% side = 'left' (default) or 'right' picks which corner.
if nargin < 3 || isempty(side), side = 'left'; end
if nargin < 4 || isempty(topfrac), topfrac = 0.88; end
drawnow limitrate
ap = get(ax,'Position');
set(lg,'Units','normalized');
lp = get(lg,'Position');
if strcmp(side,'right')
    x = ap(1) + ap(3) - lp(3) - 0.008;
else
    x = ap(1) + 0.008;
end
set(lg,'Position',[x, ap(2)+topfrac*ap(4)-lp(4), lp(3), lp(4)]);
end


function g = neutral_ladder(n)
% n greys spanning the same lightness range as the coloured ladders, for a
% key that has to stand for every hue at once.
lo = 0.16; hi = 0.80;
g  = repmat(linspace(lo, hi, n)', 1, 3);
end


function yl = legend_headroom(yl, amount)
% Reserve a strip above the data for the legend. The legend is then placed at
% a 'north...' location, which is inside that strip and therefore over empty
% axes - the only way to guarantee it never lands on a trace whatever the y
% range works out to be, short of putting it outside the panel and losing
% plot width. A VERTICAL legend is taller than a horizontal one, so callers
% that stack their entries pass a larger amount.
if nargin < 2, amount = 0.26; end
yl(2) = yl(2) + amount*(yl(2)-yl(1));
end


function leg_h = draw_curves(ax, x, V, cols, lws, draw_order, leg_h)
% Draw a set of curves that may contain exact duplicates.
%
% The problem this solves: two hypotheses can produce the identical curve -
% for the circular-linear correlation, H4 equals H1 exactly, because
% circ_corrcl already returns a non-negative magnitude per channel, so
% averaging magnitudes and taking the magnitude of the average are the same
% arithmetic. Drawn as two solid lines, the second simply erases the first
% and the panel silently claims one hypothesis where it has two.
%
% So coincident curves are found by comparison rather than assumed, and each
% such group is drawn ONCE as an interleaved dashed line whose dashes cycle
% through the colours of its members. The result reads as a dashed
% two-colour (or three-, or four-) line: every member is visible, and the
% interleaving itself is the signal that these hypotheses share a curve.
% Curves with no partner are drawn solid.
%
% draw_order sets painting order (last drawn sits on top); a group inherits
% the position of its earliest-drawn member.
present = find(~cellfun(@isempty, V));
if isempty(present), return; end

% group indices of curves that agree to within rounding
grp = zeros(1, numel(V));
gid = 0;
for a = present
    if grp(a), continue; end
    gid = gid + 1; grp(a) = gid;
    for b = present
        if b <= a || grp(b), continue; end
        if curves_equal(V{a}, V{b}), grp(b) = gid; end
    end
end

% paint groups in draw_order, the group taking its LAST position in that
% order so a group containing H1 is still drawn on top
gids = unique(grp(present));
rank = nan(1, numel(gids));
for i = 1:numel(gids)
    mem = find(grp == gids(i));
    rank(i) = max(arrayfun(@(m) find(draw_order == m, 1), mem));
end
[~, ord] = sort(rank);

for i = ord
    mem = find(grp == gids(i));
    if isscalar(mem)
        leg_h(mem) = plot(ax, x, V{mem}, '-', 'Color', cols{mem}, ...
                          'LineWidth', lws(mem));
    else
        hs = plot_interleaved(ax, x, V{mem(1)}, cols(mem), max(lws(mem)));
        for k = 1:numel(mem), leg_h(mem(k)) = hs(k); end
    end
end
end


function tf = curves_equal(a, b)
% "The same curve", not "similar curves". The tolerance is there only to
% absorb floating-point rounding between two routes to the same arithmetic,
% so it is scaled to the values and stays far below any real difference.
a = a(:)'; b = b(:)';
if numel(a) ~= numel(b), tf = false; return; end
sc = max([abs(a) abs(b) eps]);
tf = max(abs(a - b)) <= 1e-9 * sc;
end


function hs = plot_interleaved(ax, x, y, cols, lw)
% One curve, drawn as a dashed line whose dashes cycle through several
% colours, so every series lying on it is visible. Each colour is one line
% object with the other colours' stretches blanked to NaN, which keeps it a
% single selectable object per series in the PDF and gives each a legend
% entry of its own.
k = numel(cols);
x = x(:)'; y = y(:)';
if numel(x) < 2
    hs = gobjects(1,k);
    for c = 1:k, hs(c) = plot(ax, x, y, '-', 'Color', cols{c}, 'LineWidth', lw); end
    return
end

% Densify first. Chunking the raw samples does not work: these curves have
% only ~35 points, so asking for ~50 dashes makes round() collapse the chunk
% edges onto each other (1 2 2 3 4 4 ...), most dashes become a single point,
% and a single point draws nothing - which silently left one colour of the
% pair invisible. Linear interpolation costs nothing in fidelity because
% plot() already renders the curve as straight segments between samples.
n     = max(400, 12*numel(x));
xd    = linspace(x(1), x(end), n);
yd    = interp1(x, y, xd);
nDash = 22 * k;
edges = round(linspace(1, n, nDash + 1));
hs    = gobjects(1, k);
for c = 1:k
    yy = nan(1, n);
    for j = c:k:nDash
        idx = edges(j):edges(j+1);
        yy(idx) = yd(idx);
    end
    hs(c) = plot(ax, xd, yy, '-', 'Color', cols{c}, 'LineWidth', lw);
end
end


function yl = pad_lim(v)
v = v(isfinite(v));
if isempty(v), yl = [0 1]; return; end
lo = min(v); hi = max(v);
if hi == lo, hi = lo + 1; end
pad = 0.08*(hi-lo);
yl = [lo-pad, hi+pad];
end


% ── Significance marking, identical in every panel ───────────────────

function sig_band(ax, x, mask, yl, col)
% Pale opaque band behind the curve over each run of significant samples.
mask = logical(mask(:)');
if ~any(mask) || numel(mask) ~= numel(x), return; end
[st, en] = runs(mask);
e = edges_from_centers(x);
for k = 1:numel(st)
    xa = e(st(k)); xb = e(en(k)+1);
    patch('Parent',ax,'XData',[xa xb xb xa],'YData',[yl(1) yl(1) yl(2) yl(2)], ...
          'FaceColor',col,'EdgeColor','none');
end
end


function sig_bar(ax, x, mask, yl, col)
% Solid bar along the top of the panel: readable even when the band behind
% a busy curve is not.
mask = logical(mask(:)');
if ~any(mask) || numel(mask) ~= numel(x), return; end
[st, en] = runs(mask);
e  = edges_from_centers(x);
h  = 0.030*(yl(2)-yl(1));
y0 = yl(2) - h;
for k = 1:numel(st)
    xa = e(st(k)); xb = e(en(k)+1);
    patch('Parent',ax,'XData',[xa xb xb xa],'YData',[y0 y0 y0+h y0+h], ...
          'FaceColor',col,'EdgeColor','none');
end
end


function sig_ticks(ax, x, mask, strip, row, nrow, col)
% One significance bar per hypothesis, stacked inside a strip the caller has
% reserved ABOVE the data. Passing the strip rather than the axis limits is
% what keeps the bars off the traces: the caller decides how much room the
% rows need and adds it to the y range before anything is drawn.
mask = logical(mask(:)');
if ~any(mask) || numel(mask) ~= numel(x), return; end
[st, en] = runs(mask);
e   = edges_from_centers(x);
band = (strip(2) - strip(1)) / nrow;
h    = 0.62 * band;
y0   = strip(2) - row*band + (band - h)/2;
for k = 1:numel(st)
    xa = e(st(k)); xb = e(en(k)+1);
    patch('Parent',ax,'XData',[xa xb xb xa],'YData',[y0 y0 y0+h y0+h], ...
          'FaceColor',col,'EdgeColor','none');
end
end


function sig_outline(ax, x, y, mask, col)
% Outline the significant region of a map, padded so every loop CLOSES.
%
% contour() traces the boundary only through the data it is given, so a
% region that runs into the edge of the map comes out as an open arc: the
% reader cannot tell whether it stops there or continues off the plot. The
% mask is therefore padded with a ring of "not significant" placed exactly on
% the outer CELL EDGES, which forces every contour to close along the border
% of the map. The axis limits are set to those same cell edges by the caller
% (full_lims), so the closing segments stay visible instead of being clipped.
if ~any(mask(:)), return; end
x = x(:)'; y = y(:)';
ex = edges_from_centers(x); ey = edges_from_centers(y);
xp = [ex(1) x ex(end)];
yp = [ey(1) y ey(end)];
Mp = zeros(numel(yp), numel(xp));
Mp(2:end-1, 2:end-1) = double(mask);
contour(ax, xp, yp, Mp, [0.5 0.5], 'LineColor', col, 'LineWidth', 1.2);
end


function lims = full_lims(v)
% Axis limits at the outer cell edges rather than the cell centres, so a map
% fills its axes and a boundary outline drawn on the edge is not clipped.
e = edges_from_centers(v);
lims = [e(1) e(end)];
end


function [st, en] = runs(mask)
d  = diff([false, mask, false]);
st = find(d ==  1);
en = find(d == -1) - 1;
end


function e = edges_from_centers(x)
% Bin edges from sample centres, valid for irregular and log-spaced axes.
x = x(:)';
if isempty(x), e = []; return; end
if isscalar(x), e = [x-0.5, x+0.5]; return; end
d = diff(x);
e = [x(1)-d(1)/2, x(1:end-1)+d/2, x(end)+d(end)/2];
end


function curve_panel(ax, fq, v, thr, col, PAL, ylab)
% The standard curve panel: threshold, significance band, significance bar,
% curve on top. Used for every monkey-average spectrum in Figures 5 and 7c.
hold(ax,'on');
if isempty(v)
    axis(ax,'off');
    text(0.5,0.5,'no data','Parent',ax,'HorizontalAlignment','center');
    return
end
v  = v(:)';
fq = fq(:)';
if numel(fq) ~= numel(v), fq = 1:numel(v); end
yl = pad_lim([v thr]);
if isfinite(thr)
    sig_band(ax, fq, v >= thr, yl, band_tint(col));
    plot(ax, [min(fq) max(fq)], [thr thr], '--', 'Color', PAL.grey);
end
plot(ax, fq, v, '-', 'Color', col);
if isfinite(thr)
    sig_bar(ax, fq, v >= thr, yl, col);
end
xlim(ax,[min(fq) max(fq)]); ylim(ax, yl);
xlabel(ax,'Frequency (Hz)'); ylabel(ax, ylab);
tidy(ax);
end


% ── Vector heatmaps and colourmaps ───────────────────────────────────

function h = vec_heatmap(ax, x, y, C, cmap, clim)
% A heatmap drawn as patches instead of imagesc, so the PDF stays vector.
% imagesc embeds a bitmap that Illustrator can only move and mask.
%
% One detail that dictates the shape of this function: MATLAB's painters
% renderer RASTERISES a patch whose colour comes from FaceVertexCData with
% FaceColor 'flat' - the obvious implementation produces exactly the bitmap
% it was meant to avoid. A patch with a scalar FaceColor is exported as true
% vector. So the cell values are quantised to NLEV colour levels and each
% level is drawn as ONE multi-face patch with a solid colour. That is
% vector, compact (at most NLEV objects per map), and in Illustrator each
% level is a single selectable object.
%
% C is [numel(y) x numel(x)]; NaN cells are simply not drawn, so a blank
% cell always means "no value here", never "value zero".
%
NLEV = 128;
x = x(:)'; y = y(:)';
xe = edges_from_centers(x);
ye = edges_from_centers(y);
nx = numel(x); ny = numel(y);

V = zeros(4*nx*ny, 2);
Fc = zeros(nx*ny, 4);
lev = zeros(nx*ny, 1);
k = 0;
for iy = 1:ny
    for ix = 1:nx
        val = C(iy,ix);
        if ~isfinite(val), continue; end
        k = k + 1;
        b = (k-1)*4;
        V(b+1,:) = [xe(ix)   ye(iy)  ];
        V(b+2,:) = [xe(ix+1) ye(iy)  ];
        V(b+3,:) = [xe(ix+1) ye(iy+1)];
        V(b+4,:) = [xe(ix)   ye(iy+1)];
        Fc(k,:)  = b + (1:4);
        lev(k)   = round(1 + clim_frac(val, clim)*(NLEV-1));
    end
end
if k == 0, h = []; return; end
V = V(1:4*k,:); Fc = Fc(1:k,:); lev = lev(1:k);

h = gobjects(0);
for L = unique(lev)'
    col = map_color(clim_value((L-1)/(NLEV-1), clim), cmap, clim);
    h(end+1) = patch('Parent',ax,'Faces',Fc(lev==L,:),'Vertices',V, ...
                     'FaceColor',col,'EdgeColor','none'); %#ok<AGROW>
end
% keep the axes colormap/limits in step, so any later colorbar agrees
colormap(ax, cmap);
set(ax,'CLim',[clim(1) clim(end)]);
end


function rgb = map_color(v, cmap, clim)
% clim is either [lo hi] (linear) or [lo centre hi] (TWO-SLOPE: lo maps to the
% bottom of the map, centre to its middle, hi to the top, each arm scaled
% independently).
%
% The two-slope form exists for the de-rotation gain, whose distribution is
% wildly asymmetric: it reaches -0.9 where de-rotation destroys the map but
% only +0.01 where it helps. On a symmetric scale the entire positive side -
% the only side that could support a travelling wave - collapses into the
% pale centre and becomes invisible. Pinning the centre at zero keeps the
% sign unambiguous while giving each arm its full colour range.
n = size(cmap,1);
t = clim_frac(v, clim);
i = 1 + t*(n-1);
i0 = floor(i); i1 = min(i0+1, n); w = i - i0;
rgb = (1-w)*cmap(i0,:) + w*cmap(i1,:);
end


function labs = sci_labels(clim)
% The three colourbar ticks written as m.m x 10^n, using UNICODE superscripts.
%
% Not TeX: '\times10^{-3}' is mis-spaced in this MATLAB release (the exponent
% detaches and drifts right, same bug as R^{2}), and it would also render the
% exponent at its own font size and break the two-size rule. The Unicode
% superscript digits and minus survive the PDF export as ordinary glyphs at
% the one font size - verified in the output.
tk   = linspace(0,1,3);
labs = arrayfun(@(t) sci_label(clim_value(t, clim)), tk, 'uni', 0);
end


function str = sci_label(v)
SUP = [char(8304) char(185) char(178) char(179) char(8308) ...
       char(8309) char(8310) char(8311) char(8312) char(8313)];   % 0..9
if v == 0 || ~isfinite(v), str = '0'; return; end
e = floor(log10(abs(v)));
m = v / 10^e;
if abs(m) >= 9.95, m = m/10; e = e + 1; end        % rounding can push to 10.0
ds = sprintf('%d', abs(e));
sup = SUP(ds - '0' + 1);
if e < 0, sup = [char(8315) sup]; end
str = sprintf('%.1f%c10%s', m, char(215), sup);
end


function t = clim_frac(v, clim)
% Where a value sits on the colour axis, in 0..1, under either clim form.
if numel(clim) == 3
    if v <= clim(2)
        t = 0.5 * (v - clim(1)) / max(clim(2)-clim(1), eps);
    else
        t = 0.5 + 0.5 * (v - clim(2)) / max(clim(3)-clim(2), eps);
    end
else
    t = (v - clim(1)) / max(clim(2)-clim(1), eps);
end
t = min(max(t,0),1);
end


function v = clim_value(t, clim)
% Inverse of clim_frac: the value at a given position on the colour axis.
if numel(clim) == 3
    if t <= 0.5
        v = clim(1) + 2*t*(clim(2)-clim(1));
    else
        v = clim(2) + 2*(t-0.5)*(clim(3)-clim(2));
    end
else
    v = clim(1) + t*(clim(2)-clim(1));
end
end


function vec_colorbar(ax, cmap, clim, lab, shrink, tlabs)
% A colourbar built from patches for the same reason as vec_heatmap: a
% MATLAB colorbar is exported as a bitmap strip. Placed just outside the
% right edge of ax, in figure-normalised coordinates.
p  = get(ax,'Position');
w  = 0.011;
gp = 0.008;
% shrink = false when the caller's layout already reserves a right margin;
% shrinking then would make this panel narrower than its siblings
if nargin < 5, shrink = true; end
fp0  = get(get(ax,'Parent'),'Position');
room = w + gp + 0.055 * (19/fp0(3));   % strip + gap + tick labels + the label
if shrink && p(3) > 2*room
    p(3) = p(3) - room;         % give the strip its own space, do not overlap
    set(ax,'Position',p);
end
cax = axes('Parent',get(ax,'Parent'), ...
           'Position',[p(1)+p(3)+gp, p(2)+0.12*p(4), w, 0.76*p(4)]);
hold(cax,'on');
n = 64;
yv = linspace(0,1,n+1);
for k = 1:n
    c = map_color(clim_value((k-0.5)/n, clim), cmap, clim);
    patch('Parent',cax,'XData',[0 1 1 0],'YData',[yv(k) yv(k) yv(k+1) yv(k+1)], ...
          'FaceColor',c,'EdgeColor','none');
end
set(cax,'XTick',[],'YAxisLocation','right','TickDir','out','TickLength',[0 0], ...
        'YLim',[0 1],'XLim',[0 1],'Box','on', ...
        'XColor',[0.25 0.25 0.28],'YColor',[0.25 0.25 0.28]);
tk = linspace(0,1,3);
if nargin >= 6 && ~isempty(tlabs)
    % caller-supplied labels, for scales whose natural units are not decimal
    % - a phase axis reads as -pi, 0, pi, never as -3.14, 0, 3.14
    labs = tlabs;
else
    labs = arrayfun(@(t) sprintf('%.3g', clim_value(t, clim)), tk, 'uni',0);
end
set(cax,'YTick',tk,'YTickLabel',labs);
% The label is placed just past the tick labels rather than left to MATLAB,
% which parks it far enough out that on these narrow strips it reads as
% belonging to the next panel.
%
% The width of the tick text is MEASURED, not estimated. A character-count
% estimate cannot work here: the same coefficient that puts "-8.6x10^-1"
% in the right place leaves a bare "pi" almost touching the strip, because
% the labels across this figure set range from one glyph to ten and include
% superscripts and multiplication signs of quite different widths. Measuring
% also gets the units right for free - Extent comes back normalised to the
% strip, which is exactly the space the label position is expressed in.
fp  = get(get(ax,'Parent'),'Position');       % cm
tmp = text(cax, 0, 0, labs{1}, 'Units','normalized', 'Visible','off', ...
           'FontName', get(cax,'FontName'), 'FontSize', get(cax,'FontSize'));
tickw = 0;
for k = 1:numel(labs)
    set(tmp,'String',labs{k});
    e = get(tmp,'Extent');
    tickw = max(tickw, e(3));
end
delete(tmp);
GAP_CM = 0.22;                       % tick offset + breathing room
gap    = GAP_CM / (w*fp(3));         % cm -> strip widths
set(get(cax,'YLabel'),'String',lab,'Rotation',270, ...
    'Units','normalized','Position',[1 + tickw + gap, 0.5, 0], ...
    'HorizontalAlignment','center','VerticalAlignment','middle');
end


function m = cmap_ppy(n, family)
% The purple-pink-yellow family used for every colour map in this figure set:
% matplotlib's PLASMA (deep purple -> purple -> pink/magenta -> orange ->
% yellow), the perceptually-uniform map common in the Fries-lab figures.
% Pass 'magma' for the darker sibling (near-black -> purple -> magenta ->
% orange -> pale yellow) if that is the one intended.
%
% USED ON SIGNED DATA. plasma is sequential, not diverging, so on its own it
% gives zero no special status. That is fixed by the caller, not here: every
% signed map in this script sets limits SYMMETRIC about zero (or, for the
% two-slope de-rotation grids, pins the centre at zero), which puts zero on
% plasma's midpoint. Zero therefore always reads as the same pink, negative
% runs into purple and positive into yellow. Do not use this map on signed
% data with asymmetric, un-centred limits - zero would land on an arbitrary
% colour and the sign would stop being readable.
if nargin < 2, family = 'plasma'; end
switch lower(family)
    case 'magma'
        anchors = [0.001462 0.000466 0.013866
                   0.078815 0.054184 0.211667
                   0.232077 0.059889 0.437695
                   0.390384 0.100379 0.501864
                   0.550287 0.161158 0.505719
                   0.716387 0.214982 0.475290
                   0.868793 0.287728 0.409303
                   0.967671 0.439703 0.359630
                   0.994738 0.624350 0.427397
                   0.996580 0.793927 0.545941
                   0.987053 0.991438 0.749504];
    otherwise   % plasma
        anchors = [0.050383 0.029803 0.527975
                   0.294279 0.009606 0.631017
                   0.417642 0.000564 0.658390
                   0.517933 0.021563 0.654109
                   0.610667 0.090204 0.619951
                   0.692840 0.165141 0.564522
                   0.764193 0.240396 0.502126
                   0.826588 0.315714 0.441316
                   0.881443 0.392529 0.383229
                   0.928329 0.472975 0.326067
                   0.965024 0.559118 0.268513
                   0.988260 0.652325 0.211364
                   0.994141 0.753137 0.161404
                   0.977995 0.861432 0.142808
                   0.940015 0.975158 0.131326];
end
m = interp_map(anchors, n);
end


function m = cmap_cyclic(n, family)
% Cyclic map for PHASE, built by MIRRORING the sequential purple-pink-yellow
% ramp: the scale runs dark at -pi, up to its bright end at phase 0, and back
% down to the identical dark at +pi.
%
% That construction is what makes -pi and +pi the same colour by definition
% rather than by careful choice of endpoints, and it makes the map visibly
% symmetric about zero phase, which is what was asked for.
%
% THE COST, stated plainly: a mirrored map is symmetric, so +phi and -phi are
% given the SAME colour. The sign of a phase difference can no longer be read
% off the colour - only its magnitude. In these panels the propagation
% direction is carried by the arrow and the PGD value printed above each map,
% so the figure does not lose the information, but do not use this map
% anywhere the sign of the phase has to be read from the colour alone.
if nargin < 2, family = 'plasma'; end
base = cmap_ppy(256, family);
u    = 1 - abs(linspace(-1, 1, n)');     % 0 at +/-pi, 1 at phase 0
m    = interp1((1:256)', base, 1 + u*255);
m    = min(max(m,0),1);
end


function m = interp_map(anchors, n)
xi = linspace(1, size(anchors,1), n)';
m  = [interp1((1:size(anchors,1))', anchors(:,1), xi), ...
      interp1((1:size(anchors,1))', anchors(:,2), xi), ...
      interp1((1:size(anchors,1))', anchors(:,3), xi)];
m  = min(max(m,0),1);
end


function log_yticks(ax, v)
% Ticks at round speeds on a log10 axis, labelled in the original units.
lo = min(v); hi = max(v);
cand = [1 2 3 5 10 20 30 50 100 200 300 500 1000];
cand = cand(cand >= lo & cand <= hi);
if numel(cand) < 2, cand = [lo hi]; end
set(ax,'YTick',log10(cand),'YTickLabel',arrayfun(@(c) sprintf('%g',c), cand,'uni',0));
end


function grid = build_phase_grid(phi, coh, ch_row, ch_col, nR, nC, keep)
% Same rule as build_grid in cortical_planar_wave_PGD.m: a channel's phase
% is placed on the array only if it is coherence-significant at this
% frequency AND has a defined phase and a positive resultant. Cells that
% fail are left NaN and are not drawn, so a blank cell means "not reliable
% here", never "phase zero".
phi = phi(:); coh = coh(:); keep = logical(keep(:));
ok  = keep & isfinite(phi) & isfinite(coh) & coh > 0;
grid = nan(nR, nC);
for ch = 1:numel(phi)
    if ok(ch), grid(ch_row(ch), ch_col(ch)) = phi(ch); end
end
end


% ── Loading saved results ────────────────────────────────────────────

function [val, thr, fq] = load_curve(filepath, val_field, thr_field, freq_field)
% Read one observed curve and its permutation threshold. A missing file is
% not an error: the panel is drawn empty and says so, which is what should
% happen when a pipeline has not finished rather than a silent blank.
val = []; thr = NaN; fq = [];
if ~isfile(filepath), return; end
S = load(filepath);
if isfield(S, val_field), val = S.(val_field)(:)'; end
if isfield(S, thr_field), thr = S.(thr_field); end
if isfield(S, freq_field), fq = S.(freq_field)(:)'; end
if numel(thr) > 1, thr = thr(1); end
end


function [val, thr, fq] = load_reg(filepath)
% The regression files nest the observed curve and its threshold one level
% deeper than the coherence/correlation files.
val = []; thr = NaN; fq = [];
if ~isfile(filepath), return; end
S = load(filepath);
if isfield(S,'monkey_avg_obs') && isfield(S.monkey_avg_obs,'phase')
    val = S.monkey_avg_obs.phase(:)';
end
if isfield(S,'thresh_monkey') && isfield(S.thresh_monkey,'phase')
    thr = S.thresh_monkey.phase;
end
if isfield(S,'freq'), fq = S.freq(:)'; end
end


function P = load_erp_pooled(sig, comb, variant)
% The FINISHED pooled hit-miss ERP result, read as stored. pool_both_animals.m
% computed it over 49 sessions with a 1000-permutation session-averaged null;
% this function does no arithmetic beyond unpacking the struct.
%
% variant picks which of the two stored tests to return, and they answer
% different questions, so the caller must choose deliberately:
%   'signed'  mean over channels of (hit - miss). Keeps the direction of the
%             effect. Two-sided, and it cancels across animals.
%   'rect'    mean over channels of |hit - miss|. Non-negative, so the animals
%             cannot cancel; the sign is gone and the test is ONE-sided
%             (lim_lo is NaN by construction).
P = [];
f = fullfile(comb,'group_ERP',sig.key,'erp_hitmiss_pooled_both.mat');
if ~isfile(f), return; end
S = load(f);
if ~isfield(S,'erp') || ~isfield(S.erp, variant), return; end
E = S.erp;
V = E.(variant);
P.variant    = variant;
P.time       = E.time(:)';
P.nS         = E.nS;
P.animal_of  = E.animal_of(:)';
P.per_animal = E.per_animal;
P.obs        = V.obs;
P.mean       = V.mean(:)';
P.lim_hi     = V.lim_hi;
P.lim_lo     = V.lim_lo;
P.mask       = logical(V.mask(:)');
end


function W = load_erp_waveforms(animals, sig)
% Hit and miss grand averages POOLED over both animals, from the per-session
% timelocks erpdiff_*.m already wrote. Channels are averaged WITHIN a session
% first, so a session with more electrodes does not count for more, then the
% sessions of both animals are averaged together - the same construction, and
% the same session list, as the pooled difference in the finished result file.
%
% Only the two conditions are read here. The difference and its significance
% come from that finished file; nothing is recomputed.
W = [];
hitm = []; missm = []; hitr = []; missr = []; tvec = []; n_an = zeros(1,numel(animals));

for ia = 1:numel(animals)
    an = animals(ia);
    sessions = dir(fullfile(an.dir, [an.name '*']));
    sessions = sessions([sessions.isdir]);
    names    = {sessions.name};
    if isempty(names), continue; end
    keep = setdiff(1:numel(names), an.drop);

    acf = fullfile(an.dir,'critical_time','all_channels.mat');
    if ~isfile(acf), fprintf('    %s: no all_channels.mat\n', an.name); continue; end
    AC = load(acf); all_channels = AC.all_channels(:);

    for k = 1:numel(keep)
        op = fullfile(an.dir, names{keep(k)}, sig.sub);
        fh = fullfile(op,'ERP_real','norm_hit_timelock.mat');
        fm = fullfile(op,'ERP_real','norm_miss_timelock.mat');
        if ~isfile(fh) || ~isfile(fm), continue; end
        H = load(fh); H = H.norm_hit_timelock;
        M = load(fm); M = M.norm_miss_timelock;
        kc = ismember(H.label, all_channels);
        if ~any(kc), continue; end
        if isempty(tvec), tvec = H.time(:)'; end
        if numel(H.time) ~= numel(tvec), continue; end
        hitm(end+1,:)  = mean(H.avg(kc,:), 1, 'omitnan'); %#ok<AGROW>
        missm(end+1,:) = mean(M.avg(kc,:), 1, 'omitnan'); %#ok<AGROW>
        % rectified: each channel's absolute value FIRST, then the channel
        % average, so rectification happens at the same stage as in the
        % pooled result's rect variant
        hitr(end+1,:)  = mean(abs(H.avg(kc,:)), 1, 'omitnan'); %#ok<AGROW>
        missr(end+1,:) = mean(abs(M.avg(kc,:)), 1, 'omitnan'); %#ok<AGROW>
        n_an(ia) = n_an(ia) + 1;
    end
end
if isempty(hitm), return; end
W.time  = tvec;
W.nS    = size(hitm,1);
W.hit       = mean(hitm,  1, 'omitnan');
W.miss      = mean(missm, 1, 'omitnan');
W.hit_rect  = mean(hitr,  1, 'omitnan');
W.miss_rect = mean(missr, 1, 'omitnan');
W.split = strjoin(arrayfun(@(i) sprintf('%s %d', animals(i).name, n_an(i)), ...
                           1:numel(animals), 'uni',0), ', ');
end
