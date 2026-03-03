function [ph,diff] = diff_ar_ogs(cfg_diff)

ichan  = cfg_diff.ichan;
iter_n = cfg_diff.iter_n;
trl    = cfg_diff.trl;
toi    = cfg_diff.toi;

cd('/mnt/hpx/projects/MWSampling/4Shivangi/results/klecks_20170804_attentional-sampling_1/Phase_analysis/hit_miss/50/time_res_0.04_c')
load('freqpow_ogs.mat')
    
    
% phase

f= 12%length(cfg.foi);

for iter=1:iter_n
    cd(cfg_diff.inputfile)
    cd(num2str(iter))
    load('freqpow.mat')
    P{iter} = angle(freqpow.fourierspctrm(trl,ichan,:,toi));
    ph(iter,:) = reshape(P{iter},1,f);
    Pf{iter} = freqpow.fourierspctrm(trl,ichan,:,toi);
    phf(iter,:) = reshape(Pf{iter},1,f);
end


ph_vec= mean(phf,1)./abs(mean(phf,1));
amp_vec=mean(abs(phf),1);
Pf_avg = amp_vec.*(ph_vec);   % the phase we decided to use
ph(iter_n+1,:) = angle(Pf_avg);
ph(iter_n+2,:) = angle(freqpow_ogs.fourierspctrm(trl,ichan,:,toi));

% diff between AR phase and original signal phase

diff= angdiff(ph(iter_n+1,:),ph(iter_n+2,:));  


% cd(cfg_diff.outputfile)
% save diff diff
% save ph ph







