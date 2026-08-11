% =====================================================================
% Wave-TYPE classification from saved phase-progression data
%
% Phase can be spatially organised in more ways than a plane. Per frequency,
% which pattern does the preferred-phase map look like?
%
%   PLANE / travelling   phase increases linearly across the array, gradient
%                        vectors all point the same way. Metric: phase-gradient
%                        directionality PGD = |mean(grad)|/mean(|grad|) (->1).
%   RADIAL / source-sink phase increases (source) or decreases (sink) from a
%                        centre, gradient points radially. Metric: net
%                        divergence (mean Laplacian). + source, - sink.
%   ROTATING / spiral    phase winds around a centre, gradient circulates.
%                        Metrics: net curl, and the number of phase
%                        singularities (2x2 loops whose wrapped phase
%                        differences sum to +/-2*pi).
%   SYNCHRONOUS / none   little spatial gradient at all (mean|grad| not above
%                        its null) -> no wave.
%
% Unlike the plane fit, this works on WRAPPED phase using local neighbour
% differences and does not 2-D unwrap: global unwrapping would destroy the
% singularities a spiral is defined by.
%
% Significance: phases shuffled across electrode locations and every metric
% recomputed (same null logic as the speed script). A type is present at a
% frequency if its metric beats the 95th-percentile null.
%
% Outputs (Plots/scanning/wave_type/cp10_till_100/<dv>/):
%   wave_type_vs_freq.pdf   metric vs frequency per animal, with nulls
%   wave_type_strip.pdf     dominant significant type per frequency
%   phase_maps_focus.pdf    phase grid + gradient quiver at FOCUS_HZ
% =====================================================================
clearvars; close all; clc
addpath /mnt/hpc/projects/MWSampling/4Shivangi/software_folder/CircStat2012a

%% ─── Settings ────────────────────────────────────────────────────────
animals    = {'hermes','klecks'};
dv         = 'lfp';
base       = '/mnt/hpc/projects/MWSampling/4Shivangi';
grid_rows  = 8; grid_cols = 8;
FOCUS_HZ   = [6 20 45];      % frequencies to draw phase maps for
nPerm      = 1000;
alpha      = 0.05;
RELIABLE_Q = 0.5;            % drop bottom quantile of coh_mag per freq (NaN out)
rng(2025);

out_dir = fullfile(base,'Plots','scanning','wave_type','cp10_till_100', dv);
if ~exist(out_dir,'dir'), mkdir(out_dir); end
ch_col = ceil((1:(grid_rows*grid_cols))' / grid_rows);
ch_row = grid_rows - mod((1:(grid_rows*grid_cols))' - 1, grid_rows);

types = {'planar','radial','rotational'};
A = struct();

%% ─── Per animal ──────────────────────────────────────────────────────
for ia = 1:numel(animals)
    animalName = animals{ia};
    pp = fullfile(base, ['results_' animalName], 'scanning', ...
        'phase_progression','cp10_till_100', dv, 'phase_progression.mat');
    if ~isfile(pp), warning('No data for %s — skipping.', animalName); continue; end
    S = load(pp, 'pref_phase','coh_mag','freq','positions');
    freq = S.freq(:); nFreq = numel(freq);
    [nCh,~,nPos] = size(S.pref_phase);

    fprintf('\n=== %s / %s : %d ch, %d freq, %d pos ===\n', animalName, upper(dv), nCh, nFreq, nPos);

    PGD   = nan(nFreq,1);  PGD_thr  = nan(nFreq,1);  PGD_p  = nan(nFreq,1);
    DIV   = nan(nFreq,1);  DIV_thr  = nan(nFreq,1);  DIV_p  = nan(nFreq,1);
    CURL  = nan(nFreq,1);  CURL_thr = nan(nFreq,1);  CURL_p = nan(nFreq,1);
    NSING = nan(nFreq,1);  NS_thr   = nan(nFreq,1);  NS_p   = nan(nFreq,1);
    GSTR  = nan(nFreq,1);  GSTR_thr = nan(nFreq,1);  GSTR_p = nan(nFreq,1);
    nPGD = nan(nFreq,nPerm); nDIV = nan(nFreq,nPerm); nCURL = nan(nFreq,nPerm); % per-perm nulls (cluster test)
    phi_grids = cell(nFreq,1);

    for f = 1:nFreq
        % collapse positions: magnitude-weighted circular mean per channel
        z = sum(S.coh_mag(:,f,:) .* exp(1i*S.pref_phase(:,f,:)), 3, 'omitnan'); % nCh x 1
        phi_ch = angle(z);
        rel    = abs(z);                       % reliability weight
        good   = isfinite(phi_ch) & isfinite(rel) & rel>0;
        % drop least-reliable channels
        if any(good)
            cutoff = quantile(rel(good), RELIABLE_Q);
            phi_ch(rel < cutoff) = NaN;
        end
        grid = nan(grid_rows, grid_cols);
        for ch = 1:nCh
            if ~isnan(phi_ch(ch)), grid(ch_row(ch), ch_col(ch)) = phi_ch(ch); end
        end
        phi_grids{f} = grid;
        if sum(~isnan(grid(:))) < 8, continue; end

        m = wave_metrics(grid);
        PGD(f)=m.pgd; DIV(f)=m.div; CURL(f)=m.curl; NSING(f)=m.nsing; GSTR(f)=m.gstr;

        % permutation null: shuffle phases across valid locations
        vals = grid(~isnan(grid));
        idx  = find(~isnan(grid));
        np = nan(nPerm,5);
        for b = 1:nPerm
            gsh = nan(grid_rows, grid_cols);
            gsh(idx) = vals(randperm(numel(vals)));
            mm = wave_metrics(gsh);
            np(b,:) = [mm.pgd mm.div mm.curl mm.nsing mm.gstr];
        end
        PGD_thr(f)=quantile(np(:,1),1-alpha);  PGD_p(f)=mean(np(:,1)>=m.pgd);
        DIV_thr(f)=quantile(abs(np(:,2)),1-alpha); DIV_p(f)=mean(abs(np(:,2))>=abs(m.div));
        CURL_thr(f)=quantile(abs(np(:,3)),1-alpha); CURL_p(f)=mean(abs(np(:,3))>=abs(m.curl));
        NS_thr(f)=quantile(np(:,4),1-alpha);   NS_p(f)=mean(np(:,4)>=m.nsing);
        GSTR_thr(f)=quantile(np(:,5),1-alpha); GSTR_p(f)=mean(np(:,5)>=m.gstr);
        nPGD(f,:)=np(:,1)'; nDIV(f,:)=abs(np(:,2))'; nCURL(f,:)=abs(np(:,3))';
    end

    % Significance per wave TYPE, cluster-corrected across frequency, so a
    % band of contiguous significant frequencies survives and isolated chance
    % pokes do not. Detection uses the DIRECTIONAL metrics (PGD/div/curl) only,
    % never gated on gradient magnitude (GSTR): a slow planar wave has a small
    % gradient but a consistent direction, so gating on magnitude would discard
    % it. GSTR only tells slow from fast, and is reported separately.
    sig_planar = cluster_correct(PGD,       PGD_thr,  nPGD,  alpha);
    sig_radial = cluster_correct(abs(DIV),  DIV_thr,  nDIV,  alpha);
    sig_rot    = cluster_correct(abs(CURL), CURL_thr, nCURL, alpha);
    sigmat = [sig_planar, sig_radial, sig_rot];

    % dominant type per frequency: among cluster-significant types, the one
    % with the largest relative excess over its own null threshold.
    rel = [ (PGD-PGD_thr)./max(PGD_thr,eps), ...
            (abs(DIV)-DIV_thr)./max(DIV_thr,eps), ...
            (abs(CURL)-CURL_thr)./max(CURL_thr,eps) ];
    dom = zeros(nFreq,1);                  % 0=none,1=planar,2=radial,3=rotational
    for f=1:nFreq
        if ~any(sigmat(f,:)), continue; end
        r = rel(f,:); r(~sigmat(f,:)) = -inf;
        [~,dom(f)] = max(r);
    end

    A(ia).animal=animalName; A(ia).freq=freq;
    A(ia).PGD=PGD; A(ia).PGD_thr=PGD_thr; A(ia).PGD_p=PGD_p;
    A(ia).DIV=DIV; A(ia).DIV_thr=DIV_thr; A(ia).DIV_p=DIV_p;
    A(ia).CURL=CURL; A(ia).CURL_thr=CURL_thr; A(ia).CURL_p=CURL_p;
    A(ia).NSING=NSING; A(ia).NS_thr=NS_thr; A(ia).NS_p=NS_p;
    A(ia).GSTR=GSTR; A(ia).GSTR_p=GSTR_p;
    A(ia).dom=dom; A(ia).phi_grids={phi_grids};
    A(ia).sig=sigmat;   % [nFreq x 3] cluster-significant planar/radial/rot

    % console summary: cluster-significant bands per type
    fprintf('  --- %s significant wave bands (cluster-corrected) ---\n', animalName);
    any_band = false;
    for t=1:3
        runs = find_runs(sigmat(:,t));
        for rr=1:numel(runs)
            ix = runs{rr}; any_band = true;
            meanPGD = mean(PGD(ix)); meanGSTR = mean(GSTR(ix));
            fprintf('    %-11s %5.1f–%-5.1f Hz  (mean PGD=%.2f, mean tilt=%.3f rad/el, %s)\n', ...
                types{t}, freq(ix(1)), freq(ix(end)), meanPGD, meanGSTR, ...
                ternary(meanGSTR<0.15,'gentle tilt','steep tilt'));
        end
    end
    if ~any_band, fprintf('    (none)\n'); end
end
valid = find(arrayfun(@(s) ~isempty(s.animal), A));

%% ─── Fig 1: metric vs frequency (planar / radial / rotational) ───────
cols = lines(numel(animals));
f1 = figure('Position',[60 60 460*numel(valid) 720]);
rowlab = {'PGD (planar)','|net div| (radial)','|net curl| (rotational)','# singularities'};
for k=1:numel(valid)
    ia=valid(k);
    obs={A(ia).PGD, abs(A(ia).DIV), abs(A(ia).CURL), A(ia).NSING};
    thr={A(ia).PGD_thr, A(ia).DIV_thr, A(ia).CURL_thr, A(ia).NS_thr};
    for rr=1:4
        subplot(4,numel(valid),(rr-1)*numel(valid)+k); hold on;
        plot(A(ia).freq, obs{rr}, '-', 'Color',cols(ia,:),'LineWidth',1.6);
        plot(A(ia).freq, thr{rr}, '--', 'Color',[.4 .4 .4]);
        if rr==1, title(A(ia).animal); end
        if k==1, ylabel(rowlab{rr},'FontSize',8); end
        if rr==4, xlabel('Frequency (Hz)'); end
        grid on;
    end
end
sgtitle(sprintf('Wave-type metrics vs frequency — %s (solid=obs, dashed=95%% null)', upper(dv)));
saveas(f1, fullfile(out_dir,'wave_type_vs_freq.pdf'));

%% ─── Fig 2: dominant-type strip per frequency ───────────────────────
f2 = figure('Position',[60 60 820 120*numel(valid)+80]);
cmap = [0.85 0.85 0.85; 0.10 0.45 0.80; 0.85 0.35 0.10; 0.20 0.65 0.25]; % none,planar,radial,rot
for k=1:numel(valid)
    ia=valid(k);
    subplot(numel(valid),1,k);
    imagesc(A(ia).freq, 1, A(ia).dom'); colormap(cmap); caxis([0 3]);
    set(gca,'YTick',[]); xlabel('Frequency (Hz)'); title(A(ia).animal);
end
cb=colorbar('Ticks',[0.4 1.1 1.9 2.6],'TickLabels',{'none','planar','radial','rotational'}, ...
    'Position',[0.92 0.3 0.02 0.4]);
sgtitle('Dominant significant wave type per frequency');
saveas(f2, fullfile(out_dir,'wave_type_strip.pdf'));

%% ─── Fig 3: phase maps + gradient quiver at FOCUS_HZ ────────────────
[X,Y] = meshgrid(1:grid_cols, 1:grid_rows);
f3 = figure('Position',[60 60 320*numel(FOCUS_HZ) 300*numel(valid)]);
for k=1:numel(valid)
    ia=valid(k); pg=A(ia).phi_grids{1};
    for j=1:numel(FOCUS_HZ)
        [~,fi]=min(abs(A(ia).freq - FOCUS_HZ(j)));
        grid = pg{fi};
        subplot(numel(valid),numel(FOCUS_HZ),(k-1)*numel(FOCUS_HZ)+j);
        imagesc(1:grid_cols,1:grid_rows,grid,'AlphaData',~isnan(grid)); axis equal tight;
        colormap(gca,hsv); caxis([-pi pi]); hold on;
        [gx,gy]=phase_grad(grid);
        quiver(X,Y,gx,gy,0.7,'k','LineWidth',0.6);
        ns = count_singularities(grid);
        ttype = '';
        if A(ia).dom(fi)>0, ttype = types{A(ia).dom(fi)}; else, ttype='none'; end
        title(sprintf('%s %.0fHz: %s (ns=%d)', A(ia).animal, A(ia).freq(fi), ttype, ns), 'FontSize',8);
        set(gca,'XTick',[],'YTick',[]);
    end
end
sgtitle('Preferred-phase map (HSV) + gradient field at focus frequencies');
saveas(f3, fullfile(out_dir,'phase_maps_focus.pdf'));

fprintf('\nSaved wave-type figures under %s\n', out_dir);

%% =====================================================================
%% Cluster / utility helpers
%% =====================================================================
function sigmask = cluster_correct(obs, thr, nullmat, alpha)
% Cluster-based permutation across frequency. Cluster-forming threshold =
% each frequency's own 95th-pctile null (thr). Observed contiguous runs of
% obs>=thr are kept only if their mass (summed excess) beats the null
% distribution of the MAX cluster mass. FWER-controlled, band-tuned.
nFreq = numel(obs); nPerm = size(nullmat,2);
sigmask = false(nFreq,1);
supra = obs >= thr;
stat  = max(0, obs - thr); stat(~supra) = 0;
runs = find_runs(supra);
if isempty(runs), return; end
mass_obs = cellfun(@(r) sum(stat(r)), runs);
maxmass = zeros(nPerm,1);
for p = 1:nPerm
    o  = nullmat(:,p);
    sp = o >= thr;
    st = max(0, o - thr); st(~sp) = 0;
    rn = find_runs(sp);
    if isempty(rn), maxmass(p) = 0;
    else, maxmass(p) = max(cellfun(@(r) sum(st(r)), rn)); end
end
clthr = quantile(maxmass, 1-alpha);
for c = 1:numel(runs)
    if mass_obs(c) > clthr, sigmask(runs{c}) = true; end
end
end

function runs = find_runs(mask)
mask = mask(:)'; runs = {};
d = diff([0, double(mask), 0]);
starts = find(d == 1); ends = find(d == -1) - 1;
for i = 1:numel(starts), runs{end+1} = (starts(i):ends(i))'; end %#ok<AGROW>
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

%% =====================================================================
%% Helpers
%% =====================================================================
function m = wave_metrics(grid)
% Phase-gradient metrics on a WRAPPED phase grid (no global unwrap).
[gx,gy] = phase_grad(grid);
g = gx + 1i*gy;                       % complex gradient (col + i*row)
v = isfinite(g);
m.gstr = mean(abs(g(v)));             % overall gradient strength
m.pgd  = abs(mean(g(v))) / max(mean(abs(g(v))), eps);   % directionality 0..1
% divergence = d(gx)/dx + d(gy)/dy ; curl = d(gy)/dx - d(gx)/dy
[nR,nC]=size(grid);
dgx_dx=nan(nR,nC); dgy_dy=nan(nR,nC); dgy_dx=nan(nR,nC); dgx_dy=nan(nR,nC);
dgx_dx(:,2:nC-1)=(gx(:,3:nC)-gx(:,1:nC-2))/2;
dgy_dy(2:nR-1,:)=(gy(3:nR,:)-gy(1:nR-2,:))/2;
dgy_dx(:,2:nC-1)=(gy(:,3:nC)-gy(:,1:nC-2))/2;
dgx_dy(2:nR-1,:)=(gx(3:nR,:)-gx(1:nR-2,:))/2;
divg = dgx_dx + dgy_dy;  curlg = dgy_dx - dgx_dy;
m.div  = mean(divg(isfinite(divg)));   % signed net divergence
m.curl = mean(curlg(isfinite(curlg))); % signed net curl
m.nsing = count_singularities(grid);
end

function [gx,gy] = phase_grad(phi)
% Forward wrapped phase differences (rad per electrode). gx=col, gy=row.
[nR,nC]=size(phi);
gx=nan(nR,nC); gy=nan(nR,nC);
gx(:,1:nC-1)=angle(exp(1i*(phi(:,2:nC)-phi(:,1:nC-1))));
gy(1:nR-1,:)=angle(exp(1i*(phi(2:nR,:)-phi(1:nR-1,:))));
end

function ns = count_singularities(phi)
% Count 2x2 plaquettes whose wrapped phase circulation is +/- 2*pi
% (phase singularities = spiral centres / topological charges).
[nR,nC]=size(phi); ns=0;
for r=1:nR-1
  for c=1:nC-1
    loop=[phi(r,c) phi(r,c+1) phi(r+1,c+1) phi(r+1,c)];
    if any(isnan(loop)), continue; end
    d=angle(exp(1i*diff([loop loop(1)])));
    if abs(abs(sum(d))-2*pi) < 0.6, ns=ns+1; end
  end
end
end

