function [perm_cond_max,perm_cond_min] = timelock_analysis(cfg_fun)

load(cfg_fun.inputfile)

iperm  = cfg_fun.iperm;

cd(cfg_fun.output_paths),
load('trial_perm_ind.mat')

cd(fullfile(cfg_fun.output_paths,'ERP_real')),
load('bsl_avg')

% hits
cfg = [];
cfg.trials = trial_perm_ind.rand_matrix(iperm,1:length(trial_perm_ind.hit));
lfp_hits = ft_timelockanalysis(cfg,zlfpTrials);
norm_hit = lfp_hits.avg-bsl_avg;

% misses
cfg = [];
cfg.trials = trial_perm_ind.rand_matrix(iperm,length(trial_perm_ind.hit)+1:end);
lfp_misses = ft_timelockanalysis(cfg,zlfpTrials);
norm_miss = lfp_misses.avg-bsl_avg;

% max min
perm_cond_diff = [];
perm_cond_diff = norm_hit-norm_miss;

perm_cond_max = max(perm_cond_diff,[],2);
perm_cond_min = min(perm_cond_diff,[],2);

cd(cfg_fun.outputfile)
save perm_cond_max perm_cond_max
save perm_cond_min perm_cond_min



