function fun_powspec_session(cfg)
% fun_powspec_session
% -------------------
% One slurm job = one session.
%
% Pre-stimulus power spectrum for hits and misses separately, with the 1/f
% background removed by fooof (as in pow_freq.m). No time axis.
%
% For each condition it returns
%   the flattened (periodic) spectrum   log10(pow) - log10(aperiodic fit)
%   the aperiodic offset and exponent   from log10(A) = offset - exponent*log10(f)
%
% fooof is fit on the condition-average spectrum, not on single trials, so
% the session is the unit of observation for the group test.
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
%   pow_hit, pow_miss       chan x freq, log10 raw power
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

% pre-stimulus window
t_start = max(cellfun(@(t) t(1),   data.time));
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

[pow_hit,  frac_hit,  freq, label] = local_spec(dat, hitIdx,  cfg.foi);
[pow_miss, frac_miss]              = local_spec(dat, missIdx, cfg.foi);
clear dat

% flattened (periodic) spectra, as in pow_freq.m
flat_hit  = log10(pow_hit)  - log10(frac_hit);
flat_miss = log10(pow_miss) - log10(frac_miss);

% aperiodic parameters: log10(A) = offset - exponent*log10(f)
[off_hit,  exp_hit ] = local_aperiodic(frac_hit,  freq);
[off_miss, exp_miss] = local_aperiodic(frac_miss, freq);

pow_hit  = log10(pow_hit);
pow_miss = log10(pow_miss);

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


function [pow, frac, freq, label] = local_spec(dat, trials, foi)
% trial-averaged power spectrum and its fooof aperiodic fit

cfg            = [];
cfg.method     = 'mtmfft';
cfg.taper      = 'hann';
cfg.keeptrials = 'no';
cfg.foi        = foi;
cfg.pad        = 'nextpow2';
cfg.trials     = trials';

cfg.output = 'pow';
p          = ft_freqanalysis(cfg, dat);

cfg.output = 'fooof_aperiodic';
a          = ft_freqanalysis(cfg, dat);

pow   = p.powspctrm;
frac  = a.powspctrm;
freq  = p.freq;
label = p.label;
end


function [offset, expo] = local_aperiodic(frac, freq)
% offset and exponent of the aperiodic component.
% Fitted rather than read out of the fooof struct so this does not depend on
% which FieldTrip version stores fooofparams.

x = log10(freq(:));
X = [ones(length(x),1) -x];              % log10(A) = offset - exponent*log10(f)
offset = nan(size(frac,1),1);
expo   = nan(size(frac,1),1);
for ichan = 1:size(frac,1)
    y = log10(frac(ichan,:))';
    ok = isfinite(y);
    if sum(ok) < 3, continue, end
    b = X(ok,:) \ y(ok);
    offset(ichan) = b(1);
    expo(ichan)   = b(2);
end
end
