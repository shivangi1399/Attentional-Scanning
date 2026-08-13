function fun_tfr_perm_session(cfg)
% fun_tfr_perm_session
% --------------------
% One slurm job = one session.
%
% Computes the single-trial TFR of a session once, then the real hit-miss
% difference map and cfg.nperm permuted difference maps. Hit/miss labels are
% relabelled within this session only, keeping the session's own hit and miss
% counts.
%
% input (struct):
%   cfg.inputfile   path to zlfptrials.mat
%   cfg.outputfile  path of the .mat to write
%   cfg.foi         frequencies (Hz)
%   cfg.ncycles     sliding window length, in cycles
%   cfg.toi         time points (s)
%   cfg.nperm       number of within-session relabelings
%   cfg.min_trials  skip the session below this many hits or misses
%   cfg.seed        rng seed, must differ between sessions
%
% output (saved to cfg.outputfile):
%   diff_real   chan x freq x time, log10(hit) - log10(miss)
%   diff_perm   nperm x chan x freq x time, the same under relabelling
%   label, freq, time, n_hit, n_miss
%
% No file is written if the session is skipped.

rng(cfg.seed)

%% load

tmp  = load(cfg.inputfile);
fn   = fieldnames(tmp);
data = tmp.(fn{1});
clear tmp

%% trial selection

hitIdx  = find(data.trialinfo(:,20) == 1);
missIdx = find(data.trialinfo(:,20) == 5);

if length(hitIdx) < cfg.min_trials || length(missIdx) < cfg.min_trials
    warning('%s: %d hits, %d misses - skipped', cfg.inputfile, ...
        length(hitIdx), length(missIdx))
    return
end

% the sliding window has to fit inside the epoch, otherwise mtmconvol
% returns NaN and those NaNs spread through the permutation arithmetic
t_start = max(cellfun(@(t) t(1),   data.time));
t_end   = min(cellfun(@(t) t(end), data.time));
maxwin  = max(cfg.ncycles ./ cfg.foi);
if cfg.toi(1) - maxwin/2 < t_start || cfg.toi(end) + maxwin/2 > t_end
    warning('%s: epoch %g to %g s too short for toi %g to %g s with a %g s window - skipped', ...
        cfg.inputfile, t_start, t_end, cfg.toi(1), cfg.toi(end), maxwin)
    return
end

%% single-trial TFR

cfgf            = [];
cfgf.output     = 'pow';
cfgf.method     = 'mtmconvol';
cfgf.taper      = 'hanning';
cfgf.foi        = cfg.foi;
cfgf.t_ftimwin  = cfg.ncycles ./ cfg.foi;
cfgf.toi        = cfg.toi;
cfgf.pad        = 'nextpow2';
cfgf.keeptrials = 'yes';
cfgf.trials     = [hitIdx; missIdx];      % hits first, then misses
freqpow         = ft_freqanalysis(cfgf, data);
clear data

label = freqpow.label;
freq  = freqpow.freq;
time  = freqpow.time;

% log10 on single trials: the condition average becomes a geometric mean, so
% the contrast is a log power ratio and the skewed power distribution is made
% roughly symmetric
pow = single(freqpow.powspctrm);
clear freqpow
pow(pow <= 0) = NaN;
pow = log10(pow);

nTr = size(pow,1);
nH  = length(hitIdx);
nM  = nTr - nH;
nC  = length(label);
nF  = length(freq);
nT  = length(time);

P = reshape(pow, nTr, []);                % trials x voxels
clear pow
if any(isnan(P(:)))
    warning('%s: NaNs in the TFR, the difference maps will contain NaNs', cfg.inputfile)
end

%% real difference

w         = zeros(1, nTr, 'single');
w(1:nH)   =  1/nH;
w(nH+1:end) = -1/nM;
diff_real = reshape(w * P, nC, nF, nT);

%% permuted differences
% written as a weight matrix times the trial x voxel matrix, so all
% permutations are one matrix product rather than nperm passes over the data

W = zeros(cfg.nperm, nTr, 'single');
for p = 1:cfg.nperm
    r = randperm(nTr);
    W(p, r(1:nH))     =  1/nH;
    W(p, r(nH+1:end)) = -1/nM;
end
diff_perm = reshape(W * P, cfg.nperm, nC, nF, nT);
clear P W

%% save

n_hit  = nH;
n_miss = nM;
save(cfg.outputfile, 'diff_real','diff_perm','label','freq','time', ...
    'n_hit','n_miss','-v7.3')

end
