function fun_powspec_session(cfg)
% fun_powspec_session
% -------------------
% One slurm job = one session.
%
% Pre-stimulus power spectrum for hits and misses separately, with the 1/f
% background removed by fooof.
%
% Trials are combined as a GEOMETRIC mean: log10 is taken on single trials and
% then averaged, rather than averaging power and logging afterwards. This
% matches fun_tfr_perm_session.m, and is the standard choice for power, which
% is close to log-normal. The two differ by the Jensen gap, which depends on
% trial-to-trial variance - for this data the arithmetic mean shrinks the
% hit-miss difference roughly tenfold, so the choice is not cosmetic.
%
% fooof is fit on that geometric-mean spectrum, so the session is the unit of
% observation for the group test.
%
% input (struct):
%   cfg.inputfile   path to zlfptrials.mat
%   cfg.outputfile  path of the .mat to write
%   cfg.toilim      pre-stimulus window, e.g. [-1 0]
%   cfg.foi         frequencies (Hz)
%   cfg.min_trials  skip the session below this many hits or misses
%   cfg.run_name    stored for traceability
%
% output (saved to cfg.outputfile):
%   flat_hit, flat_miss     chan x freq, log10 periodic power
%   pow_hit, pow_miss       chan x freq, log10 power (geometric mean)
%   off_hit, off_miss       chan x 1, aperiodic offset
%   exp_hit, exp_miss       chan x 1, aperiodic exponent
%   label, freq, n_hit, n_miss, skipped, reason, settings

%% load

tmp  = load(cfg.inputfile);
fn   = fieldnames(tmp);
data = tmp.(fn{1});
clear tmp

%% trial selection

hitIdx  = find(data.trialinfo(:,20) == 1);
missIdx = find(data.trialinfo(:,20) == 5);

if length(hitIdx) < cfg.min_trials || length(missIdx) < cfg.min_trials
    skipped = true; %#ok<NASGU>
    reason  = sprintf('only %d hits / %d misses', length(hitIdx), length(missIdx));
    warning('%s: %s - skipped', cfg.inputfile, reason)
    save(cfg.outputfile, 'skipped', 'reason')
    return
end

t_start = max(cellfun(@(t) t(1), data.time));
if cfg.toilim(1) < t_start
    skipped = true; %#ok<NASGU>
    reason  = sprintf('epoch starts at %g s, window needs %g s', t_start, cfg.toilim(1));
    warning('%s: %s - skipped', cfg.inputfile, reason)
    save(cfg.outputfile, 'skipped', 'reason')
    return
end

cfgr        = [];
cfgr.toilim = cfg.toilim;
dat         = ft_redefinetrial(cfgr, data);
clear data

%% spectra per condition

[pow_hit,  ap_hit,  freq, label] = local_spec(dat, hitIdx,  cfg.foi);
[pow_miss, ap_miss]              = local_spec(dat, missIdx, cfg.foi);
clear dat

% flattened (periodic) spectra: everything above the fitted background
flat_hit  = pow_hit  - ap_hit;
flat_miss = pow_miss - ap_miss;

% aperiodic parameters: log10(A) = offset - exponent*log10(f)
[off_hit,  exp_hit ] = local_aperiodic(ap_hit,  freq);
[off_miss, exp_miss] = local_aperiodic(ap_miss, freq);

%% save

n_hit    = length(hitIdx);
n_miss   = length(missIdx);
skipped  = false;
reason   = '';
settings = cfg;
save(cfg.outputfile, 'flat_hit','flat_miss','pow_hit','pow_miss', ...
    'off_hit','off_miss','exp_hit','exp_miss', ...
    'label','freq','n_hit','n_miss','skipped','reason','settings')

end


function [pow_log, ap_log, freq, label] = local_spec(dat, trials, foi)
% geometric-mean log power spectrum and the fooof aperiodic fit of it

cfg            = [];
cfg.method     = 'mtmfft';
cfg.taper      = 'hann';
cfg.output     = 'pow';
cfg.keeptrials = 'yes';           % single trials, so the log comes first
cfg.foi        = foi;
cfg.pad        = 'nextpow2';
cfg.trials     = trials';
p = ft_freqanalysis(cfg, dat);

P = p.powspctrm;                  % rpt x chan x freq
P(P <= 0) = NaN;
pow_log = squeeze(nanmean(log10(P), 1));   % chan x freq, geometric mean

freq  = p.freq;
label = p.label;

ap_log = log10(local_fooof_aperiodic(10.^pow_log, freq));
end


function ap = local_fooof_aperiodic(pow, freq)
% Aperiodic component of a power spectrum, using the same brainstorm fooof
% routine and the same defaults that ft_freqanalysis uses for
% cfg.output = 'fooof_aperiodic'. Called directly here because the spectrum
% we want to fit is the geometric mean, which ft_freqanalysis cannot produce.

ft_hastoolbox('brainstorm', 1);

TF(:,1,:) = pow;
F = freq; F(F==0) = [];

opts_bst = getfield(process_fooof('GetDescription'), 'options'); %#ok<GFLD>
opt                     = struct();
opt.freq_range          = F([1 end]);
opt.peak_width_limits   = opts_bst.peakwidth.Value{1};
opt.max_peaks           = opts_bst.maxpeaks.Value{1};
opt.min_peak_height     = opts_bst.minpeakheight.Value{1}/10;   % dB -> B
opt.aperiodic_mode      = opts_bst.apermode.Value;
opt.peak_threshold      = 2;
opt.return_spectrum     = 1;
opt.border_threshold    = 1;
opt.power_line          = '50';
opt.peak_type           = opts_bst.peaktype.Value;
opt.proximity_threshold = opts_bst.proxthresh.Value{1};
opt.guess_weight        = opts_bst.guessweight.Value;
opt.thresh_after        = true;
opt.sort_type           = opts_bst.sorttype.Value;
opt.sort_param          = opts_bst.sortparam.Value;
opt.sort_bands          = opts_bst.sortbands.Value;

hasOptimTools = exist('fmincon','file') > 0;
[fs, fg] = process_fooof('FOOOF_matlab', TF, freq, opt, hasOptimTools);

apf = cat(1, fg.ap_fit);
ap  = nan(size(pow));
for k = 1:size(apf,1)
    ap(k,:) = interp1(fs, apf(k,:), freq, 'linear', nan);
end
end


function [offset, expo] = local_aperiodic(logA, freq)
% offset and exponent of the aperiodic component, from
%   log10(A) = offset - exponent*log10(f)

x = log10(freq(:));
X = [ones(length(x),1) -x];
offset = nan(size(logA,1),1);
expo   = nan(size(logA,1),1);
for ichan = 1:size(logA,1)
    y  = logA(ichan,:)';
    ok = isfinite(y) & isfinite(x);
    if sum(ok) < 3, continue, end
    b = X(ok,:) \ y(ok);
    offset(ichan) = b(1);
    expo(ichan)   = b(2);
end
end
