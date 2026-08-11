% =====================================================================
% Planar traveling-wave EXISTENCE test (per animal + combined)
%
% Asks, per animal: IS there a planar traveling wave, and at WHICH
% frequencies?
%
% Method, per frequency:
%   - Build the preferred-phase map on the 8x8 array using only
%     coherence-SIGNIFICANT channels.
%   - PGD = |mean(grad phi)| / mean(|grad phi|): are the phase-gradient
%     arrows aligned (1 = planar wave, 0 = random)?
%   - Significance = shuffle phases across electrodes (null) + CLUSTER
%     permutation across frequency -> report a BAND, not bins.
%
% Three views of the same data:
%   COLLAPSED   positions combined (reliability-weighted circular mean);
%               a wave survives only if consistent across positions.
%   PER-POSITION wave fit separately per stimulus position, with a
%               wave-ORIGIN check (does it radiate from the RF-driven
%               patch?) and cross-position direction agreement.
%   CONSENSUS   one merged band = collapsed-sig AND positions agree in
%               direction AND significant in >= CONSENSUS_MIN_POS positions.
%
% Combine animals (never pool channels): REPLICATION (sig in both) +
% POOLED standardized-PGD evidence. Speed is a flagged secondary readout.
%
% ---------------------------------------------------------------------
% WHAT THIS SCRIPT DOES NOT ANSWER — does a significant band PROPAGATE?
%
% A significant PGD says the phase map is a plane. It does NOT say the
% plane is moving. PGD is scale-free and is evaluated one frequency at a
% time, so it cannot tell a traveling wave from a phase offset that is
% simply frozen in space:
%
%   a real wave      = one conduction speed v for every frequency
%                      -> k = 2*pi*f/v grows in proportion to f
%   a fixed offset   = one phase offset for every frequency
%                      -> k constant, so the IMPLIED v = 2*pi*f/k grows with f
%
% Both give a clean plane at every single frequency; they differ only in
% how the ramp scales ACROSS frequency. That question is settled in
% cortical_planar_wave_derotation.m, whose best-fit speed ridge is fitted
% to all electrodes at once. Do not try to settle it from the per-frequency
% v_f printed below: v_f is derived from a local finite-difference gradient
% whose bias floor (GMAG_null) is a large fraction of the signal, so it is
% both biased and frequency-dependent in its bias. v_f is reported here as
% a descriptive readout only, always alongside GMAG_null so the margin over
% the floor is visible.
% =====================================================================

clearvars; close all; clc

%% ─── Settings ────────────────────────────────────────────────────────
animals    = {'hermes','klecks'};
dv         = 'lfp';
base       = '/mnt/hpc/projects/MWSampling/4Shivangi';
grid_rows  = 8; grid_cols = 8;
SPACING_MM = 0.4;                 % electrode pitch 
nPerm      = 1000;
alpha      = 0.05;
CH_FILTER  = 'significant';       % 'significant' = keep only channels whose phase
                                  %   COHERENCE is significant at that freq (max-stat
                                  %   perm threshold from phase_coherence/complex);
                                  % 'all' = keep every channel with a defined phase.
COH_SIG_ALPHA = 0.05;             % per-channel coherence significance level
SPEED_OK   = [5 100];             % physiologically plausible cortical wave speed (cm/s)
MIN_CH     = 8;                   % min reliable channels to attempt a map
DO_PER_POSITION = true;           % also test the planar wave SEPARATELY per stimulus position
CONSENSUS  = true;                % merge collapsed + per-position into ONE "consensus"
                                  % band per animal (needs DO_PER_POSITION). A frequency
                                  % is a consensus wave if ALL of:
                                  %   (i)   collapsed PGD significant, AND
                                  %   (ii)  cross-position direction agreement above chance
                                  %         (Rayleigh p < alpha), AND
                                  %   (iii) individually significant in >= CONSENSUS_MIN_POS positions
CONSENSUS_MIN_POS = 0;            % min # positions individually significant (criterion iii)
rng(2025);

% RF-mapping summary files (to find which channels each stimulus position
% DRIVES, for the wave-origin check). Same files as traveling_wave_H2_H1.m.
rf_sessions = struct( ...
    'klecks','klecks_20170829_rfmapping_bar_1_channel_target_summary.txt', ...
    'hermes','hermes_20170829_rfmapping_bar_1_channel_target_summary.txt');

out_dir = fullfile(base,'Plots','scanning','planar_wave_existence','cp10_till_100', dv);
res_dir = fullfile(base,'results_combined','scanning','planar_wave_existence','cp10_till_100', dv);
if ~exist(out_dir,'dir'), mkdir(out_dir); end
if ~exist(res_dir,'dir'), mkdir(res_dir); end
ch_col = ceil((1:(grid_rows*grid_cols))' / grid_rows);
ch_row = grid_rows - mod((1:(grid_rows*grid_cols))' - 1, grid_rows);

A = struct();

%% ─── Per-animal existence test ───────────────────────────────────────
for ia = 1:numel(animals)
    animalName = animals{ia};
    pp = fullfile(base, ['results_' animalName], 'scanning', ...
        'phase_progression','cp10_till_100', dv, 'phase_progression.mat');
    if ~isfile(pp), warning('No data for %s — skipping.', animalName); continue; end
    S = load(pp, 'pref_phase','coh_mag','freq','positions');
    freq = S.freq(:); nFreq = numel(freq);
    [nCh,~,nPos] = size(S.pref_phase);
    fprintf('\n=== %s / %s : %d ch, %d freq, %d pos ===\n', animalName, upper(dv), nCh, nFreq, nPos);

    % Coherence-significance mask [nCh x nFreq]: keep only channels whose
    % phase is reliably present (significant) at each frequency. This
    % replaces any arbitrary quantile cut — we keep ALL significant channels.
    if strcmp(CH_FILTER,'significant')
        coh_root = fullfile(base, ['results_' animalName], 'phase_coherence', ...
            'complex','cp10_till_100', dv, 'all_loc_difflev');
        coh_sig = load_coh_sig_mask(coh_root, nCh, nFreq, COH_SIG_ALPHA);
        fprintf('  coherence-significant channels: %d–%d across freqs (median %d)\n', ...
            min(sum(coh_sig,1)), max(sum(coh_sig,1)), round(median(sum(coh_sig,1))));
    else
        coh_sig = true(nCh, nFreq);
    end

    PGD = nan(nFreq,1); PGD_thr = nan(nFreq,1); PGD_p = nan(nFreq,1);
    GMAG = nan(nFreq,1); GMAG_null = nan(nFreq,1);   % gradient magnitude + bias floor
    DIR  = nan(nFreq,1);                              % propagation direction
    fracpos = nan(nFreq,1);                           % fraction of positions individually organised
    nPGD = nan(nFreq,nPerm);

    for f = 1:nFreq
        % per-position PGD (robustness check)
        pgd_pos = nan(nPos,1);
        for p = 1:nPos
            gp = build_grid(S.pref_phase(:,f,p), S.coh_mag(:,f,p), ch_row, ch_col, grid_rows, grid_cols, coh_sig(:,f));
            if sum(~isnan(gp(:)))>=MIN_CH, mp = pgd_metrics(gp); pgd_pos(p)=mp.pgd; end
        end
        fracpos(f) = mean(pgd_pos > 0.4, 'omitnan');   % loose descriptive threshold

        % collapsed map (magnitude-weighted circular mean over positions)
        z = sum(S.coh_mag(:,f,:) .* exp(1i*S.pref_phase(:,f,:)), 3, 'omitnan');
        grid = build_grid_from_complex(z, ch_row, ch_col, grid_rows, grid_cols, coh_sig(:,f));
        if sum(~isnan(grid(:))) < MIN_CH, continue; end

        m = pgd_metrics(grid);
        PGD(f)=m.pgd; GMAG(f)=m.gmag; DIR(f)=m.dir;

        % permutation null: shuffle phases across electrode locations
        vals = grid(~isnan(grid)); idx = find(~isnan(grid));
        pn = nan(nPerm,1); gn = nan(nPerm,1);
        for b = 1:nPerm
            gsh = nan(grid_rows,grid_cols); gsh(idx) = vals(randperm(numel(vals)));
            mm = pgd_metrics(gsh); pn(b)=mm.pgd; gn(b)=mm.gmag;
        end
        PGD_thr(f)=quantile(pn,1-alpha); PGD_p(f)=mean(pn>=m.pgd);
        nPGD(f,:)=pn'; GMAG_null(f)=mean(gn);
    end

    % EXISTENCE: cluster-corrected PGD bands across frequency
    sig = cluster_correct(PGD, PGD_thr, nPGD, alpha);

    % standardized PGD (for combining animals)
    mu = mean(nPGD,2,'omitnan'); sd = std(nPGD,0,2,'omitnan');
    zPGD     = (PGD - mu) ./ max(sd,eps);
    zPGD_null = (nPGD - mu) ./ max(sd,eps);          % nFreq x nPerm

    % GRADIENT MAGNITUDE (how steep the phase tilts across the array):
    %   GMAG        COHERENT |mean(grad phi)| per electrode (rad/electrode) —
    %               the planar tilt; noise self-cancels so no de-bias needed.
    %   k_planar    wavenumber (rad/mm) = GMAG / spacing
    %   dphi_array  total phase change across the array span, in DEGREES —
    %               the most interpretable form (e.g. 30 deg = a twelfth of a
    %               cycle across the whole array = gentle tilt).
    array_span_mm = (grid_cols-1) * SPACING_MM;          % long-axis extent
    k_corr = GMAG / SPACING_MM;                          % rad/mm (coherent)
    dphi_array = rad2deg(k_corr * array_span_mm);        % deg across the array
    v_f = (2*pi*freq ./ max(k_corr,eps)) / 10;           % cm/s

    A(ia).animal=animalName; A(ia).freq=freq;
    A(ia).PGD=PGD; A(ia).PGD_thr=PGD_thr; A(ia).PGD_p=PGD_p;
    A(ia).sig=sig; A(ia).DIR=DIR; A(ia).fracpos=fracpos;
    A(ia).zPGD=zPGD; A(ia).zPGD_null=zPGD_null;
    A(ia).GMAG=GMAG; A(ia).GMAG_null=GMAG_null; A(ia).k_corr=k_corr;
    A(ia).dphi_array=dphi_array; A(ia).v_f=v_f;

    % per-band report
    fprintf('  --- %s planar-wave bands (cluster-corrected PGD) ---\n', animalName);
    runs = find_runs(sig);
    if isempty(runs), fprintf('    (none)\n'); end
    A(ia).bands = struct('f_lo',{},'f_hi',{},'PGD',{},'dir',{}, ...
        'grad_radmm',{},'dphi_array_deg',{},'v_lo',{},'v_hi',{},'v_med',{}, ...
        'speed_ok',{},'fracpos',{},'floor_ratio',{});
    for r = 1:numel(runs)
        ix = runs{r};
        vmed = median(v_f(ix),'omitnan');
        vlo  = v_f(ix(1)); vhi = v_f(ix(end));
        ok = vmed>=SPEED_OK(1) && vmed<=SPEED_OK(2);
        dmean = rad2deg(angle(mean(exp(1i*DIR(ix)),'omitnan')));
        kmean = mean(k_corr(ix),'omitnan');              % rad/mm
        dphi  = mean(dphi_array(ix),'omitnan');          % deg across array

        % MARGIN OVER THE BIAS FLOOR. GMAG is a vector mean over ~56 finite-
        % difference sites, so random scatter does NOT cancel completely —
        % GMAG_null is what a SHUFFLED map of the same phases still scores.
        % A ratio near 1 means the measured tilt (and every speed derived
        % from it) is mostly finite-array artifact. Read v_med with this.
        fratio = mean(GMAG(ix),'omitnan') / max(mean(GMAG_null(ix),'omitnan'),eps);

        fprintf(['    %.1f–%.1f Hz | mean PGD=%.2f | dir=%3.0f° | grad=%.3f rad/mm ' ...
                 '(%.0f° across array) | %.0f%% pos organised | speed %.0f–%.0f cm/s %s\n'], ...
            freq(ix(1)), freq(ix(end)), mean(PGD(ix)), dmean, kmean, dphi, ...
            100*mean(fracpos(ix),'omitnan'), min(vlo,vhi), max(vlo,vhi), ...
            ternary(ok,'(plausible)','(NON-PHYSICAL -> near-synchronous, not a real wave)'));
        fprintf('        gradient vs shuffle floor: GMAG/GMAG_null = %.2f%s\n', fratio, ...
            ternary(fratio<1.5, '  <- AT THE FLOOR: tilt and speed not interpretable', ''));

        A(ia).bands(end+1) = struct('f_lo',freq(ix(1)),'f_hi',freq(ix(end)), ...
            'PGD',mean(PGD(ix)),'dir',dmean,'grad_radmm',kmean,'dphi_array_deg',dphi, ...
            'v_lo',vlo,'v_hi',vhi,'v_med',vmed,'speed_ok',ok, ...
            'fracpos',mean(fracpos(ix),'omitnan'),'floor_ratio',fratio); %#ok<AGROW>
    end

    %% ── PER-POSITION mode (keeps the collapsed result above) ──────────
    % The collapse above can cancel a stimulus-evoked wave whose direction
    % differs per stimulus position. So also fit the planar wave SEPARATELY
    % per position, test direction AGREEMENT across positions, and check
    % whether the wave radiates from the RF-driven patch (origin check).
    if DO_PER_POSITION
        ch_covers = load_driven(rf_sessions.(animalName), base, animalName, nCh, nPos);
        PGDp = nan(nFreq,nPos); DIRp = nan(nFreq,nPos); sigp = false(nFreq,nPos);
        origin_align = nan(nFreq,nPos);     % |angle| between prop axis and source->centre (deg, 0..90)
        nrel_pos = nan(nFreq,nPos);         % # reliable channels per (freq,pos) = coverage
        for p = 1:nPos
            nP_pos = nan(nFreq,nPerm);
            for f = 1:nFreq
                gp = build_grid(S.pref_phase(:,f,p), S.coh_mag(:,f,p), ch_row, ch_col, grid_rows, grid_cols, coh_sig(:,f));
                nrel_pos(f,p) = sum(~isnan(gp(:)));
                if sum(~isnan(gp(:)))<MIN_CH, continue; end
                mp = pgd_metrics(gp); PGDp(f,p)=mp.pgd; DIRp(f,p)=mp.dir;
                vals=gp(~isnan(gp)); idx=find(~isnan(gp));
                pn=nan(nPerm,1);
                for b=1:nPerm
                    gs=nan(grid_rows,grid_cols); gs(idx)=vals(randperm(numel(vals)));
                    pn(b)=pgd_metrics(gs).pgd;
                end
                nP_pos(f,:)=pn';
                % wave-origin: axis of propagation vs source(driven)->array-centre axis
                drv = find(ch_covers(:,p));
                rc = ~isnan(gp);                                   % reliable cells used
                if ~isempty(drv) && any(rc(:))
                    src = [mean(ch_row(drv)), mean(ch_col(drv))];
                    [rr,cc]=find(rc); ctr=[mean(rr) mean(cc)];
                    outv = (ctr(2)-src(2)) + 1i*(ctr(1)-src(1));   % col + i*row, like m.dir
                    da = angle(exp(1i*(mp.dir - angle(outv))));
                    origin_align(f,p) = rad2deg(min(abs(da), pi-abs(da))); % fold to axis [0,90]
                end
            end
            thr_p = quantile(nP_pos,1-alpha,2);
            sigp(:,p) = cluster_correct(PGDp(:,p), thr_p, nP_pos, alpha);
        end
        % direction AGREEMENT across positions at each frequency (resultant R)
        dir_agree = nan(nFreq,1);
        for f=1:nFreq
            d = DIRp(f,~isnan(DIRp(f,:)));
            if numel(d)>=2, dir_agree(f)=abs(mean(exp(1i*d))); end
        end
        A(ia).PGDp=PGDp; A(ia).DIRp=DIRp; A(ia).sigp=sigp; A(ia).ch_covers=ch_covers;
        A(ia).dir_agree=dir_agree; A(ia).origin_align=origin_align;
        A(ia).nrel_pos=nrel_pos;
        A(ia).rawphi=S.pref_phase; A(ia).rawcoh=S.coh_mag;   % for the origin phase-map figure
        A(ia).coh_sig=coh_sig;

        % console
        nsig_any = sum(any(sigp,2));
        fprintf('  per-position: %d/%d freqs show a planar wave in >=1 position; ', nsig_any, nFreq);
        fprintf('mean cross-position direction agreement R=%.2f (at sig freqs)\n', ...
            mean(dir_agree(any(sigp,2)),'omitnan'));
        oa = origin_align(sigp);
        fprintf('  wave-origin: mean axis-alignment to driven->centre = %.0f deg ', mean(oa,'omitnan'));
        fprintf('(small => wave radiates from the RF-driven patch)\n');

        % coverage diagnostic: median reliable channels per position, and a
        % verdict for positions with no wave (too few channels vs. genuinely
        % well-covered but no wave).
        med_nrel = median(nrel_pos,1,'omitnan');     % per position
        for p = 1:nPos
            if any(sigp(:,p)), continue; end          % only the no-wave positions
            if med_nrel(p) < 2*MIN_CH
                verdict = sprintf('too few channels (median %d) -> cannot test', round(med_nrel(p)));
            else
                verdict = sprintf('well covered (median %d) -> genuinely no wave', round(med_nrel(p)));
            end
            fprintf('    pos %d: NO WAVE — %s\n', p, verdict);
        end

        % ── CONSENSUS wave: merge collapsed + per-position into one band ──
        % freq is a consensus wave if collapsed-significant AND positions
        % agree in direction (Rayleigh p<alpha) AND >=CONSENSUS_MIN_POS
        % positions are individually significant.
        if CONSENSUS
            npos_sig = sum(sigp,2);                   % nFreq x 1
            rayl_p   = nan(nFreq,1);
            for f=1:nFreq
                d = DIRp(f,~isnan(DIRp(f,:)));
                if numel(d)>=2, rayl_p(f) = rayleigh_p(d); end
            end
            consensus = sig & (rayl_p < alpha) & (npos_sig >= CONSENSUS_MIN_POS);
            A(ia).consensus=consensus; A(ia).rayl_p=rayl_p; A(ia).npos_sig=npos_sig;
            fprintf('  CONSENSUS planar-wave band(s) [collapsed & dir-agree & >=%d pos]: %s\n', ...
                CONSENSUS_MIN_POS, band_str(freq, consensus));
        end
    end
end
valid = find(arrayfun(@(s) ~isempty(s.animal), A));

%% ─── Combine animals ─────────────────────────────────────────────────
freqC = A(valid(1)).freq; nFreq = numel(freqC);

% 1) REPLICATION: significant in every animal
rep = true(nFreq,1);
for ia = valid, rep = rep & A(ia).sig; end

% 2) POOLED EVIDENCE: mean standardized PGD + cluster correction
zG = zeros(nFreq,1); zG_null = zeros(nFreq,nPerm);
for ia = valid
    zG = zG + A(ia).zPGD;
    zG_null = zG_null + A(ia).zPGD_null;
end
zG = zG/numel(valid); zG_null = zG_null/numel(valid);
zG_thr = quantile(zG_null,1-alpha,2);                 % per-freq threshold
sigG = cluster_correct(zG, zG_thr, zG_null, alpha);

% 3) CONSENSUS replicated across animals (if computed)
rep_cons = [];
if CONSENSUS && isfield(A,'consensus')
    rep_cons = true(nFreq,1);
    for ia = valid
        if isempty(A(ia).consensus), rep_cons = false(nFreq,1); break; end
        rep_cons = rep_cons & A(ia).consensus;
    end
end

fprintf('\n================ COMBINED ================\n');
fprintf('Replicated planar-wave freqs (sig in ALL animals): %s\n', band_str(freqC, rep));
fprintf('Pooled-evidence planar-wave band(s): %s\n', band_str(freqC, sigG));
if ~isempty(rep_cons)
    fprintf('CONSENSUS wave replicated in ALL animals: %s\n', band_str(freqC, rep_cons));
end

% Band summary across animals. These bands say the phase map is a PLANE.
% Whether that plane PROPAGATES is not answerable here — see the header
% and cortical_planar_wave_derotation.m. Speeds below are descriptive and
% only meaningful where GMAG/GMAG_null is comfortably above 1.
fprintf('\n--- SIGNIFICANT PGD BANDS (planar, propagation NOT tested here) ---\n');
fprintf('%-10s %-14s %-20s %10s %10s\n','animal','band','v across band','tilt (deg)','GMAG/null');
for k = 1:numel(valid)
    ia = valid(k);
    if ~isfield(A,'bands') || isempty(A(ia).bands)
        fprintf('%-10s (no significant band)\n', A(ia).animal); continue;
    end
    for r = 1:numel(A(ia).bands)
        B = A(ia).bands(r);
        fprintf('%-10s %5.1f–%-7.1f %7.0f -> %-10.0f %10.0f %10.2f%s\n', ...
            A(ia).animal, B.f_lo, B.f_hi, B.v_lo, B.v_hi, B.dphi_array_deg, ...
            B.floor_ratio, ternary(B.floor_ratio<1.5,'  <- at the floor',''));
    end
end
fprintf(['\n=> For propagation vs frozen offset, run cortical_planar_wave_derotation.m:\n' ...
         '   its best-fit speed ridge is fitted to all electrodes at once, so it does not\n' ...
         '   inherit the finite-difference bias floor that contaminates v above.\n']);

%% ─── Plots ───────────────────────────────────────────────────────────
cols = lines(numel(animals));

% Fig 1: per-animal PGD vs frequency + null + significant bands
f1 = figure('Position',[80 80 820 300*numel(valid)]);
for k=1:numel(valid)
    ia=valid(k); subplot(numel(valid),1,k); hold on;
    fr=A(ia).freq;
    shade_bands(fr, A(ia).sig, [0.2 0.5 0.9]);
    plot(fr, A(ia).PGD, '-', 'Color',cols(ia,:),'LineWidth',1.8,'DisplayName','PGD');
    plot(fr, A(ia).PGD_thr, '--', 'Color',[.4 .4 .4],'DisplayName','95% null');
    ylim([0 1]); ylabel('PGD'); title(sprintf('%s — planar wave where shaded',A(ia).animal));
    if k==numel(valid), xlabel('Frequency (Hz)'); end
    legend('Location','best'); grid on;
end
sgtitle('Planar-wave existence per animal (PGD vs frequency)');
saveas(f1, fullfile(out_dir,'pgd_existence_per_animal.pdf'));

% Fig 2: combined evidence
f2 = figure('Position',[80 80 820 420]); hold on;
shade_bands(freqC, sigG, [0.6 0.85 0.6]);
for k=1:numel(valid)
    ia=valid(k);
    plot(A(ia).freq, A(ia).zPGD, '-', 'Color',cols(ia,:),'LineWidth',1.2, ...
        'DisplayName',sprintf('%s z-PGD',A(ia).animal));
end
plot(freqC, zG, 'k-', 'LineWidth',2.2, 'DisplayName','pooled z-PGD');
% mark replicated frequencies
yl=ylim; plot(freqC(rep), repmat(yl(2)*0.95,sum(rep),1), 'kv','MarkerFaceColor','k', ...
    'DisplayName','replicated (sig in both)');
xlabel('Frequency (Hz)'); ylabel('standardized PGD (z vs own null)');
title('Combined planar-wave evidence (shaded=pooled sig band; ▼=replicated)');
legend('Location','best'); grid on;
saveas(f2, fullfile(out_dir,'pgd_existence_combined.pdf'));

% Fig 3: gradient magnitude (steepness of the phase tilt) per animal.
% Shown as total de-biased phase change across the array, in degrees:
% small = gentle tilt (slow/near-synchronous), large = steep tilt.
f3 = figure('Position',[80 80 820 300*numel(valid)]);
for k=1:numel(valid)
    ia=valid(k); subplot(numel(valid),1,k); hold on;
    fr=A(ia).freq;
    shade_bands(fr, A(ia).sig, [0.2 0.5 0.9]);
    plot(fr, A(ia).dphi_array, '-', 'Color',cols(ia,:),'LineWidth',1.8, ...
        'DisplayName','measured tilt');
    % The shuffle floor in the SAME units. The vector mean over ~56 sites
    % does not fully cancel, so this is what a structureless map scores.
    % Where the two curves meet, the measured tilt is uninterpretable.
    plot(fr, rad2deg(A(ia).GMAG_null*(grid_cols-1)), '--', 'Color',[.45 .45 .45], ...
        'LineWidth',1.2, 'DisplayName','shuffle floor');
    hl=yline(360,'k:','one full cycle'); set(hl,'YLimInclude','off'); % reference, don't rescale
    ylabel('\Delta\phi across array (deg)');
    title(sprintf('%s — phase tilt magnitude (shaded = planar-wave band)',A(ia).animal));
    if k==numel(valid), xlabel('Frequency (Hz)'); end
    legend('Location','northwest','FontSize',7); grid on;
end
sgtitle('Gradient magnitude: phase change across the array, against its shuffle floor');
saveas(f3, fullfile(out_dir,'gradient_magnitude_per_animal.pdf'));

%% ── Per-position figures ─────────────────────────────────────────────
if DO_PER_POSITION && isfield(A,'sigp')
    % Fig 4: per-position planar-wave significance (freq x position) + cross-
    % position direction agreement, per animal.
    f4 = figure('Position',[60 60 460*numel(valid) 640]);
    for k=1:numel(valid)
        ia=valid(k); fr=A(ia).freq; nPos=size(A(ia).sigp,2);
        % top: significance grid (PGD value, significant cells outlined)
        subplot(2,numel(valid),k);
        imagesc(1:nPos, fr, A(ia).PGDp); set(gca,'YDir','normal'); caxis([0 1]); hold on;
        [fi,pi_]=find(A(ia).sigp);
        plot(pi_, fr(fi), 'w.','MarkerSize',8);   % significant (freq,pos) cells
        xlabel('stimulus position'); ylabel('Frequency (Hz)');
        title(sprintf('%s: per-position PGD (• = sig)',A(ia).animal)); colorbar;
        % bottom: cross-position direction agreement vs freq
        subplot(2,numel(valid),numel(valid)+k);
        plot(fr, A(ia).dir_agree, '-', 'Color',cols(ia,:),'LineWidth',1.6); hold on;
        sg=any(A(ia).sigp,2);
        plot(fr(sg), A(ia).dir_agree(sg), 'o','Color',cols(ia,:),'MarkerFaceColor',cols(ia,:));
        ylim([0 1]); yline(0.5,'k:'); xlabel('Frequency (Hz)');
        ylabel('cross-position dir agreement R');
        title('R~1 = positions agree (combine OK); R~0 = disagree (don''t combine)');
    end
    sgtitle('Per-position planar wave + cross-position direction agreement');
    saveas(f4, fullfile(out_dir,'per_position_existence.pdf'));

    % Fig 5: wave-origin check — phase map for EVERY position (so coverage is
    % visible even where there is no wave), DRIVEN channels outlined, and the
    % propagation arrow where a wave exists. Title shows coverage + verdict.
    for k=1:numel(valid)
        ia=valid(k); nPos=size(A(ia).sigp,2); fr=A(ia).freq;
        nc=min(nPos,5); nr=ceil(nPos/nc);
        f5=figure('Position',[60 60 230*nc 250*nr]);
        med_nrel = median(A(ia).nrel_pos,1,'omitnan');
        for p=1:nPos
            subplot(nr,nc,p); hold on;
            sgf=find(A(ia).sigp(:,p));
            % pick the frequency to display: peak-PGD sig freq, else the
            % best-covered frequency (so no-wave maps still show coverage)
            if ~isempty(sgf)
                [~,bi]=max(A(ia).PGDp(sgf,p)); fpk=sgf(bi); haswave=true;
            else
                [~,fpk]=max(A(ia).nrel_pos(:,p)); haswave=false;
            end
            gp = build_grid(A(ia).rawphi(:,fpk,p), A(ia).rawcoh(:,fpk,p), ch_row, ch_col, grid_rows, grid_cols, A(ia).coh_sig(:,fpk));
            imagesc(1:grid_cols,1:grid_rows,gp,'AlphaData',~isnan(gp)); axis equal tight;
            colormap(gca,hsv); caxis([-pi pi]); set(gca,'YDir','reverse');
            axis([0.5 grid_cols+0.5 0.5 grid_rows+0.5]);
            drv=find(A(ia).ch_covers(:,p));
            plot(ch_col(drv), ch_row(drv),'ws','MarkerSize',9,'LineWidth',1.5); % driven channels
            if haswave
                d=A(ia).DIRp(fpk,p);
                quiver(grid_cols/2, grid_rows/2, 2*cos(d), 2*sin(d), 0,'k','LineWidth',2,'MaxHeadSize',2);
                ttl=sprintf('pos %d: WAVE %.0fHz\nalign=%.0f°, nrel=%d',p,fr(fpk),A(ia).origin_align(fpk,p),round(med_nrel(p)));
            elseif med_nrel(p) < 2*MIN_CH
                ttl=sprintf('pos %d: no wave\n(too few ch, nrel=%d)',p,round(med_nrel(p)));
            else
                ttl=sprintf('pos %d: no wave\n(covered, nrel=%d)',p,round(med_nrel(p)));
            end
            title(ttl,'FontSize',7); set(gca,'XTick',[],'YTick',[]);
        end
        sgtitle(sprintf('%s wave origin: phase map (HSV) + driven ch (□) + propagation (arrow); nrel = median reliable channels',A(ia).animal));
        saveas(f5, fullfile(out_dir, sprintf('wave_origin_%s.pdf',A(ia).animal)));
    end

    % Fig 6: coverage vs wave — # reliable channels per position, marking
    % which positions have a wave. Confirms whether "no wave" = low coverage.
    f6 = figure('Position',[60 60 460*numel(valid) 360]);
    for k=1:numel(valid)
        ia=valid(k); nPos=size(A(ia).sigp,2);
        med_nrel = median(A(ia).nrel_pos,1,'omitnan');
        haswave = any(A(ia).sigp,1);
        subplot(1,numel(valid),k); hold on;
        b=bar(1:nPos, med_nrel, 'FaceColor',[.7 .7 .7]);  %#ok<NASGU>
        bar(find(haswave), med_nrel(haswave), 'FaceColor',cols(ia,:));
        yline(MIN_CH,'r--','MIN\_CH'); yline(2*MIN_CH,'k:','"covered" cutoff');
        xlabel('stimulus position'); ylabel('median # reliable channels');
        title(sprintf('%s: coverage per position (colored = has wave)',A(ia).animal)); grid on;
    end
    sgtitle('Coverage per stimulus position — does "no wave" track low channel count?');
    saveas(f6, fullfile(out_dir,'coverage_per_position.pdf'));
end

% Fig 7: CONSENSUS wave — PGD with the three criteria and the merged band.
if CONSENSUS && isfield(A,'consensus')
    f7 = figure('Position',[60 60 820 320*numel(valid)]);
    for k=1:numel(valid)
        ia=valid(k); fr=A(ia).freq; subplot(numel(valid),1,k); hold on;
        shade_bands(fr, A(ia).sig,       [0.80 0.80 0.80]);   % collapsed-sig (grey)
        shade_bands(fr, A(ia).consensus, [0.20 0.60 0.20]);   % CONSENSUS (green)
        plot(fr, A(ia).PGD, '-', 'Color',cols(ia,:),'LineWidth',1.8,'DisplayName','PGD');
        plot(fr, A(ia).PGD_thr, '--','Color',[.4 .4 .4],'DisplayName','PGD 95% null');
        % criterion ticks near the top
        yc=[0.93 0.86 0.79];
        crit = {A(ia).sig, A(ia).rayl_p<alpha, A(ia).npos_sig>=CONSENSUS_MIN_POS};
        lab  = {'collapsed sig','dir agree','\geq pos'};
        for c=1:3
            fcc=fr(crit{c}); plot(fcc, repmat(yc(c),numel(fcc),1), 's', ...
                'Color',cols(ia,:),'MarkerFaceColor',cols(ia,:),'MarkerSize',3,'HandleVisibility','off');
            text(fr(end)*1.01, yc(c), lab{c},'FontSize',6,'Color',[.3 .3 .3]);
        end
        ylim([0 1]); ylabel('PGD'); xlim([0 max(fr)*1.12]);
        title(sprintf('%s — CONSENSUS band (green) = collapsed & dir-agree & \\geq%d pos',A(ia).animal,CONSENSUS_MIN_POS));
        if k==numel(valid), xlabel('Frequency (Hz)'); end
        grid on;
    end
    sgtitle('Consensus planar wave: merge of collapsed + per-position evidence');
    saveas(f7, fullfile(out_dir,'consensus_wave.pdf'));
end

Asave = A;   % drop the big raw phase arrays (already in phase_progression.mat)
if isfield(Asave,'rawphi'), Asave = rmfield(Asave,{'rawphi','rawcoh'}); end
results = struct('A',Asave,'rep',rep,'sigG',sigG,'freq',freqC,'animals',{animals}, ...
    'rep_cons',rep_cons,'CONSENSUS',CONSENSUS,'CONSENSUS_MIN_POS',CONSENSUS_MIN_POS, ...
    'dv',dv,'SPACING_MM',SPACING_MM,'SPEED_OK',SPEED_OK);
save(fullfile(res_dir,'planar_wave_existence.mat'),'results','-v7.3');
fprintf('\nSaved figures + results under %s\n', out_dir);

%% =====================================================================
%% Helpers
%% =====================================================================
function ch_covers = load_driven(rf_name, base, animalName, nCh, nPos)
% Which channels each stimulus position DRIVES (RF overlaps the stimulus),
% parsed from the RF-mapping summary (same logic as traveling_wave_H2_H1.m).
ch_covers = false(nCh, nPos);
rf_file = fullfile(base,'Plots','RF_Mapping',animalName,'loc_RF_map','gaussian_overlap', rf_name);
if ~isfile(rf_file), warning('RF file not found for %s — origin check skipped.', animalName); return; end
t = readtable(rf_file,'Delimiter','\t');
for ch = 1:min(nCh,height(t))
    s = t.Locations_Inside{ch};
    if strcmp(s,'none')||isempty(s), continue; end
    locs = str2double(strsplit(s,',')); locs = locs(~isnan(locs)&locs>=1&locs<=nPos);
    ch_covers(ch,locs) = true;
end
end

function [sig, freq] = load_coh_sig_mask(coh_root, nCh, nFreq, alpha)
% Per-channel coherence-significance mask [nCh x nFreq] from the
% phase_coherence/complex pipeline. Per channel, the max-statistic
% (freq-corrected) threshold is the (1-alpha) quantile of the per-perm max
% over frequency of |coh_perm_complex|; the channel is significant at freq f
% where |coh_complex(f)| exceeds that single threshold.
sig = false(nCh, nFreq); freq = [];
for ch = 1:nCh
    cfile = fullfile(coh_root, num2str(ch), 'coherence.mat');
    pfile = fullfile(coh_root, num2str(ch), 'coh_perm_complex.mat');
    if ~isfile(cfile) || ~isfile(pfile), continue; end
    Cc = load(cfile, 'coh', 'coh_complex');
    Pp = load(pfile, 'coh_perm_complex');
    if isfield(Cc,'coh'),             obs = Cc.coh(:)';
    elseif isfield(Cc,'coh_complex'), obs = abs(Cc.coh_complex(:))';
    else, continue; end
    perm_mag = abs(Pp.coh_perm_complex);
    if any(isnan(obs)) || any(isnan(perm_mag(:))), continue; end
    thr = quantile(max(perm_mag,[],2), 1-alpha);
    n = min(numel(obs), nFreq);
    sig(ch,1:n) = obs(1:n) >= thr;
end
end

function grid = build_grid(phi_ch_in, coh_in, ch_row, ch_col, nR, nC, keep)
% Place a channel's phase on the grid if it is in the KEEP set (e.g. the
% coherence-significant channels at this frequency) AND has a defined phase.
phi_ch = phi_ch_in(:); rel = coh_in(:); keep = logical(keep(:));
ok = keep & isfinite(phi_ch) & isfinite(rel) & rel>0;
grid = nan(nR,nC);
for ch = 1:numel(phi_ch)
    if ok(ch), grid(ch_row(ch),ch_col(ch)) = phi_ch(ch); end
end
end


function grid = build_grid_from_complex(z, ch_row, ch_col, nR, nC, keep)
grid = build_grid(angle(z), abs(z), ch_row, ch_col, nR, nC, keep);
end

function m = pgd_metrics(grid)
[gx,gy] = phase_grad(grid);
g = gx + 1i*gy; v = isfinite(g);
m.gtot = mean(abs(g(v)));                         % total local gradient (incl. noise)
m.gmag = abs(mean(g(v)));                         % COHERENT/planar gradient (rad/electrode):
                                                  % the vector-mean magnitude — random
                                                  % scatter self-cancels, so this is the
                                                  % true planar tilt (= PGD * gtot). No
                                                  % de-biasing needed.
m.pgd  = m.gmag / max(m.gtot, eps);              % directionality 0..1
m.dir  = angle(mean(g(v)));                       % propagation direction
end

function [gx,gy] = phase_grad(phi)
[nR,nC]=size(phi); gx=nan(nR,nC); gy=nan(nR,nC);
gx(:,1:nC-1)=angle(exp(1i*(phi(:,2:nC)-phi(:,1:nC-1))));
gy(1:nR-1,:)=angle(exp(1i*(phi(2:nR,:)-phi(1:nR-1,:))));
end

function sigmask = cluster_correct(obs, thr, nullmat, alpha)
nFreq=numel(obs); nPerm=size(nullmat,2); sigmask=false(nFreq,1);
supra = obs>=thr; stat=max(0,obs-thr); stat(~supra)=0;
runs=find_runs(supra); if isempty(runs), return; end
mass=cellfun(@(r) sum(stat(r)), runs);
mx=zeros(nPerm,1);
for p=1:nPerm
    o=nullmat(:,p); sp=o>=thr; st=max(0,o-thr); st(~sp)=0;
    rn=find_runs(sp);
    if ~isempty(rn), mx(p)=max(cellfun(@(r) sum(st(r)), rn)); end
end
cl=quantile(mx,1-alpha);
for c=1:numel(runs), if mass(c)>cl, sigmask(runs{c})=true; end, end
end

function runs = find_runs(mask)
mask=mask(:)'; runs={}; d=diff([0 double(mask) 0]);
s=find(d==1); e=find(d==-1)-1;
for i=1:numel(s), runs{end+1}=(s(i):e(i))'; end %#ok<AGROW>
end

function shade_bands(freq, sig, col)
runs=find_runs(sig);
for r=1:numel(runs)
    ix=runs{r}; xr=[freq(ix(1)) freq(ix(end))];
    h=patch([xr fliplr(xr)],[-1e6 -1e6 1e6 1e6],col, ...
        'FaceAlpha',0.25,'EdgeColor','none','HandleVisibility','off');
    set(h,'YLimInclude','off','XLimInclude','off');  % don't let shading rescale axes
end
end

function s = band_str(freq, sig)
runs=find_runs(sig);
if isempty(runs), s='(none)'; return; end
parts=cell(1,numel(runs));
for r=1:numel(runs), ix=runs{r}; parts{r}=sprintf('%.1f-%.1f Hz',freq(ix(1)),freq(ix(end))); end
s=strjoin(parts,', ');
end

function p = rayleigh_p(angles)
% Rayleigh test for non-uniformity of circular data (Zar approximation).
angles = angles(~isnan(angles)); n = numel(angles);
if n < 2, p = NaN; return; end
R = abs(mean(exp(1i*angles(:)))); z = n*R^2;
p = exp(-z) * (1 + (2*z - z^2)/(4*n) - (24*z - 132*z^2 + 76*z^3 - 9*z^4)/(288*n^2));
p = min(max(p,0),1);
end

function out = ternary(c,a,b), if c, out=a; else, out=b; end, end
