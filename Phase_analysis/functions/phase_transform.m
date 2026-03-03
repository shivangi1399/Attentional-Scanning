function [phase_ar] = phase_transform(cfg)

ichan  = cfg.ichan;
iter_n = cfg.iter_n;
toi    = cfg.toi;

% phase

for iter=1:iter_n
    cd(cfg.inputfile)
    cd(num2str(iter))
    load('freqpow.mat')
    phf(iter,:,:) = squeeze(freqpow.fourierspctrm(:,ichan,:,toi));
    clear freqpow
end

cd(cfg.inputfile)
cd(num2str(iter))
load('freqpow.mat')

ph_vec= squeeze(mean(phf,1)./abs(mean(phf,1)));
amp_vec=squeeze(mean(abs(phf),1));
Pf_avg = amp_vec.*(ph_vec);   % the phase we decided to use
phase_ar = angle(Pf_avg);

% save
phase = [];
phase.trialinfo = freqpow.trialinfo;
phase.label = freqpow.label;
phase.freq = freqpow.freq;
phase.iter_phase = angle(phf);
phase.ph_vec = ph_vec;
phase.amp_vec = amp_vec;
phase.Pf_avg = Pf_avg;
phase.ar_phase = phase_ar;
%phase.ogs_phase = phase_og;
% phase.diff = diff;

cd(cfg.outputfile)
ESIsave phase phase







