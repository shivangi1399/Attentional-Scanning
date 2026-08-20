% =====================================================================
% Planar traveling-wave EXISTENCE test (per animal + combined)
%
% Per animal: is there a planar traveling wave, and at which frequencies?
%
% Per frequency:
%   - build the preferred-phase map on the 8x8 array from coherence-
%     SIGNIFICANT channels only;
%   - PGD = |mean(grad phi)| / mean(|grad phi|): are the gradient arrows
%     aligned (1 = planar wave, 0 = random)?
%   - significance: shuffle phases across electrodes, then cluster-permute
%     across frequency -> report a BAND, not bins.
%
% Three views of the same data:
%   COLLAPSED     positions combined (reliability-weighted circular mean); a
%                 wave survives only if consistent across positions.
%   PER-POSITION  fit separately per stimulus position, with a wave-ORIGIN
%                 check (does it radiate from the RF-driven patch?) and
%                 cross-position direction agreement.
%   CONSENSUS     one merged band: collapsed-sig AND positions agree in
%                 direction AND sig in >= CONSENSUS_MIN_POS positions.
%
% Combining animals (channels never pooled): REPLICATION (sig in both) plus
% pooled standardised-PGD evidence.
%
% WHAT THIS DOES NOT ANSWER -- whether a significant band PROPAGATES. PGD says
% the phase map is a plane, not that the plane is moving; it is scale-free and
% evaluated one frequency at a time, so it cannot separate
%   a real wave    one speed v at every frequency -> k = 2*pi*f/v grows with f
%   a fixed offset one phase offset at every f    -> k constant, so the implied
%                                                    v = 2*pi*f/k grows with f
% Both give a clean plane at every single frequency and differ only in how the
% ramp scales ACROSS frequency, which cortical_planar_wave_derotation.m settles
% by fitting all electrodes at once.
%
% NO SPEED IS COMPUTED HERE, deliberately. A speed would have to come from
% v = 2*pi*f/k with k taken from GMAG, a vector mean of local finite
% differences whose bias floor (GMAG_null) is a large fraction of the signal.
% Any such speed is biased, and biased differently at each frequency, so it is
% not interpretable and was removed rather than left in as a tempting number.
% The same floor also qualifies the TILT (dphi_array_deg), which is therefore
% always reported alongside GMAG/GMAG_null so the margin over the floor is
% visible. Speeds belong to cortical_planar_wave_derotation.m.
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
MIN_CH     = 8;                   % min reliable channels to attempt a map
MIN_TRIALS_FRAC = 0.5;            % drop stimulus positions sampled with fewer than this
                                  % fraction of the animal's MEDIAN trials/position.
                                  % coh_mag is a resultant length, so its noise floor
                                  % scales as 1/sqrt(n): an under-sampled position gets an
                                  % INFLATED magnitude and therefore MORE weight in the
                                  % across-position circular mean, which is backwards.
                                  % Set to 0 to keep every position.
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
    S = load(pp, 'pref_phase','coh_mag','freq','positions','n_pos');
    freq = S.freq(:); nFreq = numel(freq);
    [nCh,~,nPos_all] = size(S.pref_phase);
    fprintf('\n=== %s / %s : %d ch, %d freq, %d pos ===\n', animalName, upper(dv), nCh, nFreq, nPos_all);

    % ── Drop under-sampled stimulus positions ─────────────────────────
    % n_pos is trials per (channel,position); dead channels sit at 0, so take
    % the max over channels as the position's trial count. pos_keep indexes
    % back into the ORIGINAL 1..nPos_all numbering (the RF file and the figure
    % labels both use those original indices).
    trials_pos = double(max(S.n_pos, [], 1));            % 1 x nPos_all
    if MIN_TRIALS_FRAC > 0
        keep_pos = trials_pos >= MIN_TRIALS_FRAC * median(trials_pos);
    else
        keep_pos = true(1, nPos_all);
    end
    pos_keep = find(keep_pos);
    if any(~keep_pos)
        for pd = find(~keep_pos)
            fprintf(['  DROP position %d (coord %g): %d trials vs median %d ' ...
                     '— under-sampled, its coh_mag weight would be noise-inflated\n'], ...
                pd, S.positions(pd), trials_pos(pd), round(median(trials_pos)));
        end
    end
    S.pref_phase = S.pref_phase(:,:,keep_pos);
    S.coh_mag    = S.coh_mag(:,:,keep_pos);
    S.positions  = S.positions(keep_pos);
    nPos = numel(pos_keep);
    if nPos < nPos_all
        fprintf('  keeping %d of %d positions\n', nPos, nPos_all);
    end

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
        % p = (1 + #{perm >= obs}) / (1 + nPerm), not mean(perm >= obs): under the
        % null the observed arrangement is exchangeable with the nPerm shuffled
        % ones, so it belongs in its own reference set. The naive form omits it,
        % is biased low, and can return exactly 0 — a certainty no finite number
        % of shuffles can support.
        PGD_thr(f)=quantile(pn,1-alpha); PGD_p(f)=(1+sum(pn>=m.pgd))/(1+nPerm);
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

    A(ia).animal=animalName; A(ia).freq=freq; A(ia).pos_keep=pos_keep;
    A(ia).PGD=PGD; A(ia).PGD_thr=PGD_thr; A(ia).PGD_p=PGD_p;
    A(ia).sig=sig; A(ia).DIR=DIR; A(ia).fracpos=fracpos;
    A(ia).zPGD=zPGD; A(ia).zPGD_null=zPGD_null;
    A(ia).GMAG=GMAG; A(ia).GMAG_null=GMAG_null; A(ia).k_corr=k_corr;
    A(ia).dphi_array=dphi_array;

    % per-band report
    fprintf('  --- %s planar-wave bands (cluster-corrected PGD) ---\n', animalName);
    runs = find_runs(sig);
    if isempty(runs), fprintf('    (none)\n'); end
    A(ia).bands = struct('f_lo',{},'f_hi',{},'PGD',{},'dir',{}, ...
        'grad_radmm',{},'dphi_array_deg',{},'fracpos',{},'floor_ratio',{});
    for r = 1:numel(runs)
        ix = runs{r};
        dmean = rad2deg(angle(mean(exp(1i*DIR(ix)),'omitnan')));
        kmean = mean(k_corr(ix),'omitnan');              % rad/mm
        dphi  = mean(dphi_array(ix),'omitnan');          % deg across array

        % MARGIN OVER THE BIAS FLOOR. GMAG is a vector mean over ~56 finite-
        % difference sites, so random scatter does NOT cancel completely —
        % GMAG_null is what a SHUFFLED map of the same phases still scores.
        % A ratio near 1 means the measured tilt is mostly finite-array
        % artifact. Read dphi_array_deg with this.
        fratio = mean(GMAG(ix),'omitnan') / max(mean(GMAG_null(ix),'omitnan'),eps);

        fprintf(['    %.1f–%.1f Hz | mean PGD=%.2f | dir=%3.0f° | grad=%.3f rad/mm ' ...
                 '(%.0f° across array) | %.0f%% pos organised\n'], ...
            freq(ix(1)), freq(ix(end)), mean(PGD(ix)), dmean, kmean, dphi, ...
            100*mean(fracpos(ix),'omitnan'));
        fprintf('        gradient vs shuffle floor: GMAG/GMAG_null = %.2f%s\n', fratio, ...
            ternary(fratio<1.5, '  <- AT THE FLOOR: tilt not interpretable', ''));

        A(ia).bands(end+1) = struct('f_lo',freq(ix(1)),'f_hi',freq(ix(end)), ...
            'PGD',mean(PGD(ix)),'dir',dmean,'grad_radmm',kmean,'dphi_array_deg',dphi, ...
            'fracpos',mean(fracpos(ix),'omitnan'),'floor_ratio',fratio); %#ok<AGROW>
    end

    %% ── PER-POSITION mode (keeps the collapsed result above) ──────────
    % The collapse above can cancel a stimulus-evoked wave whose direction
    % differs per stimulus position. So also fit the planar wave SEPARATELY
    % per position, test direction AGREEMENT across positions, and check
    % whether the wave radiates from the RF-driven patch (origin check).
    if DO_PER_POSITION
        % RF file lists ORIGINAL position indices, so parse with nPos_all and
        % then subset to the kept positions.
        ch_covers = load_driven(rf_sessions.(animalName), base, animalName, nCh, nPos_all);
        ch_covers = ch_covers(:, pos_keep);
        PGDp = nan(nFreq,nPos); DIRp = nan(nFreq,nPos); sigp = false(nFreq,nPos);
        origin_align = nan(nFreq,nPos);     % |angle| between prop axis and source->centre (deg, 0..90)
        nrel_pos = nan(nFreq,nPos);         % # reliable channels per (freq,pos) = coverage
        PGDp_p   = nan(nFreq,nPos);         % RAW per-(freq,pos) shuffle p, kept for the
                                            % fixed-frequency test at the coherence peaks
                                            % (sigp below is cluster-corrected over FREQ,
                                            % which a two-frequency test cannot use)
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
                nP_pos(f,:)=pn'; PGDp_p(f,p)=(1+sum(pn>=mp.pgd))/(1+nPerm);   % see PGD_p above
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
        A(ia).PGDp_p=PGDp_p;
        A(ia).dir_agree=dir_agree; A(ia).origin_align=origin_align;
        A(ia).nrel_pos=nrel_pos;
        A(ia).rawphi=S.pref_phase; A(ia).rawcoh=S.coh_mag;   % for the origin phase-map figure
        A(ia).coh_sig=coh_sig;
        A(ia).positions=S.positions; A(ia).trials_pos=trials_pos(pos_keep);

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
            fprintf('    pos %d: NO WAVE — %s\n', pos_keep(p), verdict);
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
% and cortical_planar_wave_derotation.m. The tilt below is descriptive and
% only meaningful where GMAG/GMAG_null is comfortably above 1.
fprintf('\n--- SIGNIFICANT PGD BANDS (planar, propagation NOT tested here) ---\n');
fprintf('%-10s %-14s %12s %12s\n','animal','band','tilt (deg)','GMAG/null');
for k = 1:numel(valid)
    ia = valid(k);
    if ~isfield(A,'bands') || isempty(A(ia).bands)
        fprintf('%-10s (no significant band)\n', A(ia).animal); continue;
    end
    for r = 1:numel(A(ia).bands)
        B = A(ia).bands(r);
        fprintf('%-10s %5.1f–%-7.1f %12.0f %12.2f%s\n', ...
            A(ia).animal, B.f_lo, B.f_hi, B.dphi_array_deg, ...
            B.floor_ratio, ternary(B.floor_ratio<1.5,'  <- at the floor',''));
    end
end
fprintf(['\n=> For propagation vs frozen offset, run cortical_planar_wave_derotation.m:\n' ...
         '   it fits all electrodes at once, so it does not inherit the finite-difference\n' ...
         '   bias floor that contaminates the local-gradient tilt above.\n']);

%% ─── Plots ───────────────────────────────────────────────────────────
% Figures are written with save_pdf (page = figure size), not saveas, which
% would crop them to a default letter page.
cols = lines(numel(animals));

% Fig 1: per-animal PGD vs frequency + null + significant bands
f1 = new_fig(820, 300*numel(valid));
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
save_pdf(f1, fullfile(out_dir,'pgd_existence_per_animal.pdf'));

% Fig 2: combined evidence
f2 = new_fig(820, 420); hold on;
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
save_pdf(f2, fullfile(out_dir,'pgd_existence_combined.pdf'));

% Fig 3: gradient magnitude (steepness of the phase tilt) per animal.
% Shown as total de-biased phase change across the array, in degrees:
% small = gentle tilt (slow/near-synchronous), large = steep tilt.
f3 = new_fig(820, 300*numel(valid));
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
save_pdf(f3, fullfile(out_dir,'gradient_magnitude_per_animal.pdf'));

%% ── Per-position figures ─────────────────────────────────────────────
if DO_PER_POSITION && isfield(A,'sigp')
    % Fig 4: per-position planar-wave significance (freq x position) + cross-
    % position direction agreement, per animal.
    f4 = new_fig(460*numel(valid), 640);
    for k=1:numel(valid)
        ia=valid(k); fr=A(ia).freq; nPos=size(A(ia).sigp,2);
        % top: significance grid (PGD value, significant cells outlined)
        subplot(2,numel(valid),k);
        imagesc(1:nPos, fr, A(ia).PGDp); set(gca,'YDir','normal'); caxis([0 1]); hold on;
        [fi,pi_]=find(A(ia).sigp);
        plot(pi_, fr(fi), 'w.','MarkerSize',8);   % significant (freq,pos) cells
        set(gca,'XTick',1:nPos,'XTickLabel',cellstr(num2str(A(ia).pos_keep(:))));  % original position indices
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
    save_pdf(f4, fullfile(out_dir,'per_position_existence.pdf'));

    % Fig 5: wave-origin check — phase map for EVERY position (so coverage is
    % visible even where there is no wave), DRIVEN channels outlined, and the
    % propagation arrow where a wave exists. Title shows coverage + verdict.
    for k=1:numel(valid)
        ia=valid(k); nPos=size(A(ia).sigp,2); fr=A(ia).freq;
        nc=min(nPos,5); nr=ceil(nPos/nc);
        cbw=95;                                  % right margin reserved for the colorbar
        f5=new_fig(260*nc+cbw, 285*nr);
        med_nrel = median(A(ia).nrel_pos,1,'omitnan');
        ax=gobjects(nPos,1);
        for p=1:nPos
            ax(p)=subplot(nr,nc,p); hold on;
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
                ttl=sprintf('pos %d: WAVE %.0fHz\nalign=%.0f°, nrel=%d',A(ia).pos_keep(p),fr(fpk),A(ia).origin_align(fpk,p),round(med_nrel(p)));
            elseif med_nrel(p) < 2*MIN_CH
                ttl=sprintf('pos %d: no wave\n(too few ch, nrel=%d)',A(ia).pos_keep(p),round(med_nrel(p)));
            else
                ttl=sprintf('pos %d: no wave\n(covered, nrel=%d)',A(ia).pos_keep(p),round(med_nrel(p)));
            end
            title(ttl,'FontSize',7); set(gca,'XTick',[],'YTick',[]);
        end
        % one shared phase colorbar — every panel uses the same HSV [-pi pi]
        % scale. Squeeze the panels left to free the right margin for it.
        fw=get(f5,'Position'); sc=1-cbw/fw(3);
        for p=1:nPos
            pp=get(ax(p),'Position'); set(ax(p),'Position',[pp(1)*sc pp(2) pp(3)*sc pp(4)]);
        end
        cax=axes('Parent',f5,'Position',[0 0 1 1],'Visible','off','Color','none');
        colormap(cax,hsv); caxis(cax,[-pi pi]);
        cb=colorbar(cax,'Position',[1-0.62*cbw/fw(3) 0.18 0.16*cbw/fw(3) 0.62]);
        set(cb,'Ticks',[-pi -pi/2 0 pi/2 pi],'TickLabels',{'-180','-90','0','90','180'}, ...
            'FontSize',8);
        set(get(cb,'Label'),'String','preferred phase (deg)','FontSize',8);
        sgtitle(sprintf('%s wave origin: phase map (HSV) + driven ch (□) + propagation (arrow); nrel = median reliable channels',A(ia).animal), ...
            'FontSize',10);
        save_pdf(f5, fullfile(out_dir, sprintf('wave_origin_%s.pdf',A(ia).animal)));
    end

    % Fig 6: coverage vs wave — # reliable channels per position, marking
    % which positions have a wave. Confirms whether "no wave" = low coverage.
    f6 = new_fig(460*numel(valid), 360);
    for k=1:numel(valid)
        ia=valid(k); nPos=size(A(ia).sigp,2);
        med_nrel = median(A(ia).nrel_pos,1,'omitnan');
        haswave = any(A(ia).sigp,1);
        subplot(1,numel(valid),k); hold on;
        b=bar(1:nPos, med_nrel, 'FaceColor',[.7 .7 .7]);  %#ok<NASGU>
        bar(find(haswave), med_nrel(haswave), 'FaceColor',cols(ia,:));
        yline(MIN_CH,'r--','MIN\_CH'); yline(2*MIN_CH,'k:','"covered" cutoff');
        set(gca,'XTick',1:nPos,'XTickLabel',cellstr(num2str(A(ia).pos_keep(:))));  % original position indices
        xlabel('stimulus position'); ylabel('median # reliable channels');
        title(sprintf('%s: coverage per position (colored = has wave)',A(ia).animal)); grid on;
    end
    sgtitle('Coverage per stimulus position — does "no wave" track low channel count?');
    save_pdf(f6, fullfile(out_dir,'coverage_per_position.pdf'));
end

% Fig 7: CONSENSUS wave — PGD with the three criteria and the merged band.
if CONSENSUS && isfield(A,'consensus')
    f7 = new_fig(820, 320*numel(valid));
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
    save_pdf(f7, fullfile(out_dir,'consensus_wave.pdf'));
end

%% ─── Fig 8: phase map at the LFP-coherence peak frequencies ──────────
% The two peaks of the monkey-average LFP coherence spectrum are the
% frequencies at which the LFP phase is most reliably locked to the stimulus.
% Show the phase map there, averaged ACROSS POSITIONS with the same
% reliability-weighted circular mean the collapsed test uses, plus the fitted
% propagation direction. Peaks are read from the spectrum rather than typed in,
% then snapped to the nearest bin of this animal's frequency axis.
coh_peak_file = fullfile(base,'results_combined','phase_coherence','abs_per_chan', ...
    'cp10_till_100','lfp','all_loc_difflev','monkey_avg_results.mat');
COH_PEAKS_HZ = [6.11 13.33];        % fallback if the spectrum file is missing
if isfile(coh_peak_file)
    Cm = load(coh_peak_file,'coh_monkey_avg','freq');
    pk = two_local_peaks(Cm.coh_monkey_avg(:), Cm.freq(:));
    if numel(pk)==2, COH_PEAKS_HZ = pk(:)'; end
else
    warning('Monkey-avg LFP coherence not found — using default peak frequencies.');
end
nPk = numel(COH_PEAKS_HZ);
fprintf('\n--- LFP coherence peaks: %s Hz ---\n', num2str(COH_PEAKS_HZ,'%.2f  '));

cbw8 = 95;
f8 = new_fig(300*nPk + cbw8, 330*numel(valid));
ax8 = gobjects(numel(valid)*nPk,1); q = 0;
for k = 1:numel(valid)
    ia = valid(k); fr = A(ia).freq;
    if isfield(A,'rawphi') && ~isempty(A(ia).rawphi)
        phi_all = A(ia).rawphi; coh_all = A(ia).rawcoh;
    else
        Sp = load(fullfile(base,['results_' A(ia).animal],'scanning','phase_progression', ...
            'cp10_till_100', dv,'phase_progression.mat'),'pref_phase','coh_mag');
        phi_all = Sp.pref_phase(:,:,A(ia).pos_keep);     % same position drop as above
        coh_all = Sp.coh_mag(:,:,A(ia).pos_keep);
    end
    A(ia).peak_freq = nan(1,nPk); A(ia).peak_pgd = nan(1,nPk);
    A(ia).peak_dir  = nan(1,nPk); A(ia).peak_nch = nan(1,nPk);
    for j = 1:nPk
        [~,fi] = min(abs(fr - COH_PEAKS_HZ(j)));
        z  = sum(coh_all(:,fi,:) .* exp(1i*phi_all(:,fi,:)), 3, 'omitnan');  % across positions
        gm = build_grid_from_complex(z, ch_row, ch_col, grid_rows, grid_cols, A(ia).coh_sig(:,fi));
        nch = sum(~isnan(gm(:)));
        q = q+1; ax8(q) = subplot(numel(valid), nPk, (k-1)*nPk + j); hold on;
        imagesc(1:grid_cols,1:grid_rows,gm,'AlphaData',~isnan(gm)); axis equal tight;
        colormap(gca,hsv); caxis([-pi pi]); set(gca,'YDir','reverse');
        axis([0.5 grid_cols+0.5 0.5 grid_rows+0.5]); set(gca,'XTick',[],'YTick',[]);
        if nch >= MIN_CH
            m = pgd_metrics(gm);
            quiver(grid_cols/2, grid_rows/2, 2*cos(m.dir), 2*sin(m.dir), 0, 'k', ...
                'LineWidth',2,'MaxHeadSize',2);
            title(sprintf('%s — %.2f Hz%s\nPGD=%.2f (p=%.3f), dir=%.0f°, n=%d ch', ...
                A(ia).animal, fr(fi), ternary(A(ia).sig(fi),' [sig band]',''), ...
                m.pgd, A(ia).PGD_p(fi), rad2deg(m.dir), nch),'FontSize',8);
            A(ia).peak_pgd(j) = m.pgd; A(ia).peak_dir(j) = rad2deg(m.dir);
            fprintf('  %-8s %6.2f Hz : PGD=%.2f (p=%.3f), dir=%4.0f deg, %d channels\n', ...
                A(ia).animal, fr(fi), m.pgd, A(ia).PGD_p(fi), rad2deg(m.dir), nch);
        else
            title(sprintf('%s — %.2f Hz\n(only %d reliable ch)',A(ia).animal,fr(fi),nch),'FontSize',8);
            fprintf('  %-8s %6.2f Hz : only %d reliable channels — no fit\n', A(ia).animal, fr(fi), nch);
        end
        A(ia).peak_freq(j) = fr(fi); A(ia).peak_nch(j) = nch;
    end
end
% shared phase colorbar (same HSV [-pi pi] scale in every panel)
fw8 = get(f8,'Position'); sc8 = 1 - cbw8/fw8(3);
for q = 1:numel(ax8)
    pp = get(ax8(q),'Position'); set(ax8(q),'Position',[pp(1)*sc8 pp(2) pp(3)*sc8 pp(4)]);
end
cax8 = axes('Parent',f8,'Position',[0 0 1 1],'Visible','off','Color','none');
colormap(cax8,hsv); caxis(cax8,[-pi pi]);
cb8 = colorbar(cax8,'Position',[1-0.62*cbw8/fw8(3) 0.18 0.16*cbw8/fw8(3) 0.62]);
set(cb8,'Ticks',[-pi -pi/2 0 pi/2 pi],'TickLabels',{'-180','-90','0','90','180'},'FontSize',8);
set(get(cb8,'Label'),'String','preferred phase (deg)','FontSize',8);
sgtitle('Phase map at the LFP-coherence peak frequencies (averaged across positions)','FontSize',10);
save_pdf(f8, fullfile(out_dir,'phase_map_coherence_peaks.pdf'));

%% ─── Fig 9/10: PER-POSITION planar wave AT the coherence peaks ───────
% The collapsed map above vector-averages over positions, which cancels any
% wave whose direction depends on the stimulus position. So test each position
% SEPARATELY, but at the two frequencies fixed a priori by the coherence
% spectrum rather than at each position's own best frequency.
%
% Why a different correction from sigp: sigp is cluster-corrected ACROSS
% FREQUENCY, which needs a contiguous frequency axis. Pinning two frequencies
% removes that, leaving nPos independent tests per frequency — so correct
% across POSITIONS with Benjamini-Hochberg FDR instead.
%
% Two questions, kept separate:
%   (1) does any position carry a plane?      per-position shuffle p -> FDR over positions
%   (2) do the positions AGREE on direction?  Rayleigh test on the nPos directions
% (1) can pass while (2) fails: individual planes that point every which way are
% not a stimulus-evoked traveling wave, and are exactly what the collapse hides.
if DO_PER_POSITION && isfield(A,'PGDp_p')
    FDR_Q     = 0.05;                       % FDR level across positions
    n_fam     = numel(valid)*nPk;           % Rayleigh tests in the family
    alpha_fam = alpha / n_fam;              % Bonferroni for the direction tests
    fprintf('\n=== PER-POSITION planar wave AT the LFP-coherence peaks ===\n');
    fprintf('    per-position: shuffle p, BH-FDR q=%.2f across positions\n', FDR_Q);
    fprintf('    direction agreement: Rayleigh, Bonferroni alpha=%.4f (%d tests)\n', alpha_fam, n_fam);

    f10 = new_fig(300*nPk + 40, 300*numel(valid));
    for k = 1:numel(valid)
        ia = valid(k); fr = A(ia).freq; nPos = size(A(ia).sigp,2);
        A(ia).peak_pos_p   = nan(nPk,nPos);  A(ia).peak_pos_fdr = false(nPk,nPos);
        A(ia).peak_rayl_R  = nan(1,nPk);     A(ia).peak_rayl_p  = nan(1,nPk);
        for j = 1:nPk
            [~,fi] = min(abs(fr - COH_PEAKS_HZ(j)));
            pv  = A(ia).PGDp_p(fi,:);                    % 1 x nPos raw shuffle p
            [fdr_sig, pcrit] = fdr_bh(pv, FDR_Q);
            dirs = A(ia).DIRp(fi,:);
            dv   = dirs(~isnan(dirs));
            Rres = abs(mean(exp(1i*dv)));  pray = rayleigh_p(dv);
            A(ia).peak_pos_p(j,:)=pv; A(ia).peak_pos_fdr(j,:)=fdr_sig;
            A(ia).peak_rayl_R(j)=Rres; A(ia).peak_rayl_p(j)=pray;

            % ── console: the per-position table at this frequency ──
            fprintf('\n  %s @ %.2f Hz  (%d positions)\n', A(ia).animal, fr(fi), nPos);
            fprintf('    pos : %s\n', sprintf('%6d', A(ia).pos_keep));
            fprintf('    PGD : %s\n', sprintf('%6.2f', A(ia).PGDp(fi,:)));
            fprintf('    p   : %s\n', sprintf('%6.3f', pv));
            marks = repmat({'.'},1,nPos); marks(fdr_sig) = {'*'};
            fprintf('    FDR : %s   (%d/%d pass, p_crit=%.4f)\n', ...
                sprintf('%6s', marks{:}), sum(fdr_sig), nPos, pcrit);
            fprintf('    FDR-significant positions: %s\n', ...
                ternary(any(fdr_sig), mat2str(A(ia).pos_keep(fdr_sig)), '(none)'));
            fprintf('    direction agreement: R=%.2f, Rayleigh p=%.4g %s\n', Rres, pray, ...
                ternary(pray<alpha_fam,'(survives Bonferroni)', ...
                ternary(pray<alpha,'(nominal only — FAILS Bonferroni)','(n.s.)')));

            % ── Fig 9: the per-position phase maps at this frequency ──
            nc = min(nPos,5); nr = ceil(nPos/nc);
            cbw9 = 95; f9 = new_fig(260*nc + cbw9, 285*nr); ax9 = gobjects(nPos,1);
            for p = 1:nPos
                ax9(p) = subplot(nr,nc,p); hold on;
                gp = build_grid(A(ia).rawphi(:,fi,p), A(ia).rawcoh(:,fi,p), ...
                    ch_row, ch_col, grid_rows, grid_cols, A(ia).coh_sig(:,fi));
                imagesc(1:grid_cols,1:grid_rows,gp,'AlphaData',~isnan(gp)); axis equal tight;
                colormap(gca,hsv); caxis([-pi pi]); set(gca,'YDir','reverse');
                axis([0.5 grid_cols+0.5 0.5 grid_rows+0.5]); set(gca,'XTick',[],'YTick',[]);
                drv = find(A(ia).ch_covers(:,p));
                plot(ch_col(drv), ch_row(drv),'ws','MarkerSize',9,'LineWidth',1.5);
                if ~isnan(dirs(p))
                    % heavy black arrow = survives FDR; thin grey = does not
                    if fdr_sig(p), acol='k'; alw=2.4; else, acol=[.45 .45 .45]; alw=1.2; end
                    quiver(grid_cols/2, grid_rows/2, 2*cos(dirs(p)), 2*sin(dirs(p)), 0, ...
                        'Color',acol,'LineWidth',alw,'MaxHeadSize',2);
                end
                title(sprintf('pos %d%s\nPGD=%.2f, p=%.3f', A(ia).pos_keep(p), ...
                    ternary(fdr_sig(p),' *',''), A(ia).PGDp(fi,p), pv(p)),'FontSize',7);
            end
            fw9 = get(f9,'Position'); sc9 = 1 - cbw9/fw9(3);
            for p = 1:nPos
                pp = get(ax9(p),'Position'); set(ax9(p),'Position',[pp(1)*sc9 pp(2) pp(3)*sc9 pp(4)]);
            end
            cax9 = axes('Parent',f9,'Position',[0 0 1 1],'Visible','off','Color','none');
            colormap(cax9,hsv); caxis(cax9,[-pi pi]);
            cb9 = colorbar(cax9,'Position',[1-0.62*cbw9/fw9(3) 0.18 0.16*cbw9/fw9(3) 0.62]);
            set(cb9,'Ticks',[-pi -pi/2 0 pi/2 pi],'TickLabels',{'-180','-90','0','90','180'},'FontSize',8);
            set(get(cb9,'Label'),'String','preferred phase (deg)','FontSize',8);
            sgtitle(sprintf(['%s @ %.2f Hz — per-position phase map (* = plane survives FDR q=%.2f); ' ...
                'arrow = propagation, \\square = RF-driven channels'], ...
                A(ia).animal, fr(fi), FDR_Q),'FontSize',10);
            save_pdf(f9, fullfile(out_dir, sprintf('per_position_peak_%s_%.0fHz.pdf', ...
                A(ia).animal, fr(fi))));

            % ── Fig 10 panel: the direction-agreement (Rayleigh) test ──
            figure(f10); subplot(numel(valid), nPk, (k-1)*nPk + j); hold on;
            th = linspace(0,2*pi,181);
            plot(cos(th), sin(th), '-', 'Color',[.82 .82 .82]);
            plot([-1 1],[0 0],':','Color',[.85 .85 .85]); plot([0 0],[-1 1],':','Color',[.85 .85 .85]);
            for p = 1:nPos
                if isnan(dirs(p)), continue; end
                if fdr_sig(p), lc = cols(ia,:); lw = 2.0; else, lc = [.62 .62 .62]; lw = 1.0; end
                plot([0 cos(dirs(p))],[0 sin(dirs(p))],'-','Color',lc,'LineWidth',lw);
            end
            % resultant: length R, the statistic the Rayleigh test evaluates
            if ~isempty(dv)
                ang = angle(mean(exp(1i*dv)));
                quiver(0,0,Rres*cos(ang),Rres*sin(ang),0,'k','LineWidth',2.6,'MaxHeadSize',1.2);
            end
            axis equal; axis([-1.15 1.15 -1.15 1.15]); set(gca,'XTick',[],'YTick',[]); box on;
            title(sprintf('%s @ %.2f Hz\nR=%.2f, Rayleigh p=%.3g%s', A(ia).animal, fr(fi), ...
                Rres, pray, ternary(pray<alpha_fam,' *','')),'FontSize',9);
        end
    end
    figure(f10);
    sgtitle(['Do the positions agree on direction? Unit vector per position (colored = plane ' ...
             'survives FDR); black = resultant. * = Rayleigh survives Bonferroni'],'FontSize',10);
    save_pdf(f10, fullfile(out_dir,'peak_direction_agreement.pdf'));
end

Asave = A;   % drop the big raw phase arrays (already in phase_progression.mat)
if isfield(Asave,'rawphi'), Asave = rmfield(Asave,{'rawphi','rawcoh'}); end
results = struct('A',Asave,'rep',rep,'sigG',sigG,'freq',freqC,'animals',{animals}, ...
    'rep_cons',rep_cons,'CONSENSUS',CONSENSUS,'CONSENSUS_MIN_POS',CONSENSUS_MIN_POS, ...
    'dv',dv,'SPACING_MM',SPACING_MM, ...
    'COH_PEAKS_HZ',COH_PEAKS_HZ,'MIN_TRIALS_FRAC',MIN_TRIALS_FRAC);
save(fullfile(res_dir,'planar_wave_existence.mat'),'results','-v7.3');
fprintf('\nSaved figures + results under %s\n', out_dir);

%% =====================================================================
%% Helpers
%% =====================================================================
function [sig, pcrit] = fdr_bh(p, q)
% Benjamini-Hochberg FDR across the positions tested at ONE frequency.
% Sort the nPos p-values ascending, find the largest rank i with p(i) <= i/m*q,
% and declare everything at or below that p-value significant. Controls the
% expected PROPORTION of false discoveries at q (not the family-wise error), and
% assumes the tests are independent or positively dependent — reasonable here:
% the positions are separate trial sets on a shared electrode layout.
p = p(:)'; ok = isfinite(p); sig = false(size(p)); pcrit = NaN;
pv = sort(p(ok)); m = numel(pv);
if m == 0, return; end
below = find(pv <= (1:m)/m * q, 1, 'last');
if isempty(below), return; end
pcrit = pv(below);
sig(ok & p <= pcrit) = true;
end

function pk = two_local_peaks(c, f)
% Frequencies of the two highest interior local maxima of the spectrum c(f),
% returned in ascending frequency order.
c = c(:); f = f(:); n = numel(c); pk = [];
if n < 3, return; end
loc = find(c(2:n-1) > c(1:n-2) & c(2:n-1) >= c(3:n)) + 1;
if isempty(loc), return; end
[~,o] = sort(c(loc),'descend'); loc = loc(o);
pk = sort(f(loc(1:min(2,numel(loc)))));
end

function f = new_fig(w, h)
% Figure with a requested pixel canvas, kept on-screen so nothing is clipped.
ss = get(0,'ScreenSize');
w = min(w, ss(3)-80); h = min(h, ss(4)-120);
f = figure('Units','pixels','Position',[40 40 w h],'Color','w');
end

function save_pdf(fig, fname)
% Print to a PDF page exactly the size of the figure. saveas(...,'.pdf')
% instead uses a default letter page and crops anything that overflows.
drawnow;                                   % let sgtitle/legends settle first
set(fig,'Color','w','InvertHardcopy','off');
set(fig,'Units','inches'); p = get(fig,'Position');
set(fig,'PaperUnits','inches','PaperSize',[p(3) p(4)], ...
        'PaperPosition',[0 0 p(3) p(4)],'PaperPositionMode','manual');
set(fig,'Units','pixels');
try
    exportgraphics(fig, fname, 'ContentType','vector','BackgroundColor','white');
catch
    print(fig, fname, '-dpdf', '-painters', '-r300');   % pre-R2020a
end
end

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
