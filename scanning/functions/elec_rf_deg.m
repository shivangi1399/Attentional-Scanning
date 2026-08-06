function [xy_deg, valid_gauss] = elec_rf_deg(base, animal, nCh, rf_date, ppd, screen_xy)
% Shared helper — extracted from stimulus_loc_traveling_wave.m so that
% erp_latency_wave.m (and anything else needing electrode RF centres) uses the
% SAME coordinate-frame conversion. Getting this frame wrong is off by
% (840, 525) px, so it must exist in exactly one place.
% Per-channel electrode RF centre in FIXATION-CENTRED DEGREES:
%   xy_deg [nCh x 2] = (x, y) with (0,0) = fovea — the SAME frame the stimulus
%   positions are converted to, so the two may be differenced safely.
%
% Source: the Gaussian-fit centre written by RF_Mapping/mapping_lfp.m into
%   Plots/RF_Mapping/<animal>/loc_RF_map/gaussian_overlap/
%       <animal>_<date>_..._channel_target_summary.txt
% RF_Center_X/Y there are SCREEN pixels — chan_loc_mua.m adds centerX/centerY
% when writing them — so the screen centre is SUBTRACTED here before /ppd.
% The stimulus values (trialinfo col 16/17) are already fixation-centred and
% therefore need no such offset; mixing the two frames without this
% subtraction is off by (840, 525) px. Rows with no fit hold '-' -> NaN.
% Channel index is 1:1 with the phase-progression array (64 rows).
%
% valid_gauss(ch) is TRUE only where Status == 'Valid_Gaussian', i.e. a real
% 2-D Gaussian fit to that channel's RF map. Status 'Extrapolated' means the
% fit FAILED and the centre was filled in from the median of the channel's
% 8-connected array neighbours, or from a plane fit across the array
% (RF_Mapping/mapping_lfp.m:513) — inferred under a retinotopic-smoothness
% assumption, not measured. 'No_Data' rows are NaN and never valid.
xy_deg = nan(nCh,2); valid_gauss = false(nCh,1);
odir = fullfile(base,'Plots','RF_Mapping',animal,'loc_RF_map','gaussian_overlap');
rf_file = fullfile(odir, sprintf('%s_%s_rfmapping_bar_1_channel_target_summary.txt', animal, rf_date));
if ~isfile(rf_file)                                   % fall back to any bar_* for that date
    d = dir(fullfile(odir, sprintf('%s_%s_*channel_target_summary.txt', animal, rf_date)));
    if isempty(d), warning('RF summary for %s (%s) not found — no electrode centres.', animal, rf_date); return; end
    rf_file = fullfile(d(1).folder, d(1).name);
end
t = readtable(rf_file, 'Delimiter','\t');
x = t.RF_Center_X; y = t.RF_Center_Y;                 % may import as cell if '-' present
if iscell(x), x = str2double(x); end
if iscell(y), y = str2double(y); end
st = t.Status;                                        % char/cell/string/categorical
if ~iscell(st), st = cellstr(string(st)); end
cx = screen_xy(1)/2; cy = screen_xy(2)/2;             % fixation = screen centre (840, 525)
n = min(nCh, height(t));
xy_deg(1:n,1) = (x(1:n) - cx) / ppd;                  % NaN passes through
xy_deg(1:n,2) = (y(1:n) - cy) / ppd;
valid_gauss(1:n) = strcmp(strtrim(st(1:n)), 'Valid_Gaussian') & ...
                   isfinite(xy_deg(1:n,1)) & isfinite(xy_deg(1:n,2));
end
