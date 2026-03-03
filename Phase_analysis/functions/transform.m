function [phase_ar] = transform(cfg_diff)

ichan  = cfg_diff.ichan;
iter_n = cfg_diff.iter_n;
toi    = cfg_diff.toi;

% cd('/mnt/hpx/projects/MWSampling/4Shivangi/results/klecks_20170804_attentional-sampling_1/Phase_analysis/hit_miss/50/time_res_0.04_c')
% load('freqpow_ogs.mat')

% phase

for iter=1:iter_n
    cd(cfg_diff.inputfile)
    cd(num2str(iter))
    load('freqpow.mat')
    %ph(iter,:,:) = squeeze(angle(freqpow.fourierspctrm(:,ichan,:,toi)));
    phf(iter,:,:) = squeeze(freqpow.fourierspctrm(:,ichan,:,toi));
    clear freqpow
end

cd(cfg_diff.inputfile)
cd(num2str(iter))
load('freqpow.mat')

ph_vec= squeeze(mean(phf,1)./abs(mean(phf,1)));
amp_vec=squeeze(mean(abs(phf),1));
Pf_avg = amp_vec.*(ph_vec);   
phase_ar = angle(Pf_avg);  % the phase we decided to use
%phase_og = squeeze(angle(freqpow_ogs.fourierspctrm(:,ichan,:,toi)));

% diff between AR phase and original signal phase

%diff= angdiff(phase_ar,phase_og);

% save
transf = [];
transf.trialinfo = freqpow.trialinfo;
transf.label = freqpow.label;
transf.ar_phase = phase_ar;
transf.ar_transform = Pf_avg;
%phase.ogs_phase = phase_og;
% phase.diff = diff;

cd(cfg_diff.outputfile)
ESIsave transf transf







