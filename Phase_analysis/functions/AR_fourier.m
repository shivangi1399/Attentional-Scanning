function [datar] = AR_fourier(cfg_fun)

load(cfg_fun.inputfile)

exms           = cfg_fun.exms; % Extra sample (in ms) at which the phases will be extrapolated
arord          = cfg_fun.arord; % AR model order-it should be nfsample/20
cycles         = cfg_fun.cycles;
t              = cfg_fun.toi;
f              = cfg_fun.foi;
fsample        = cfg_fun.fsample; % Sampling rate of original zdata
nfsample       = cfg_fun.nfsample; % Sampling rate after downsampling

fexdat = [];

dsrat = fsample/nfsample; % order of downsampling. It is assumed that fsample is
        % an integer multiple of nfsample
 
%% AR and extrapolation

for ichan= 1:length(zdat.label)
    
    data = cellfun(@(x) x(ichan,:),zdat.trial,'UniformOutput',false);
    
    % Downsampling
    
    if  nfsample<fsample
        
        dsfiltord = 8; % order of butterworth filter applied before downsampling
        lpcuttoff = 150; % frequency in Hz
        
        lpdat = cellfun(@(x) subfunc_butterlp(x,lpcuttoff,dsfiltord,fsample), ...
            data,'UniformOutput',false); % Lowpass filter
        
        dat = cellfun(@(x) downsample(x,dsrat),lpdat,'UniformOutput',false); % Downsample
        
        clear dfiltord lpcuttoff data
    else
        
        dat = data;
    end
    
    % AR and extrapolation
    
    [West,Aest,Cest,~,th] = subfunc_ARfit(dat,arord);
    
    disp(['Extrapolating for channel ' num2str(ichan)])
    
    exdat(ichan,:) = cellfun(@(x) ...
        subfunc_ARextrap(x,exms,dsrat,Aest,Cest), ...
        dat,'UniformOutput',false);
    
    for itrial = 1:length(zdat.trial)
        fexdat{1,itrial}(ichan,:) = exdat{ichan,itrial};
    end
end

%% get zdata into fieldtrip format

datar = [];
datar = zdat;
extrapfinaltime = max(zdat.time{1,1})+exms/zdat.fsample;
datar.time  = cellfun(@(x) [x [max(zdat.time{1,1})+1/zdat.fsample:1/zdat.fsample:extrapfinaltime]],zdat.time,'uniformoutput',false);  % create new time axis MAKE more general ..
datar.trial = [];
datar.trial = fexdat;

% % ds
% datar = [];
% datar = zdat;
% sampr = (zdat.fsample)/dsrat;
% time = downsample(zdat.time{1,1},dsrat);
% extrapfinaltime = max(time)+(exms/dsrat)/sampr;
% datar.time  = cellfun(@(x) [downsample(x,dsrat) [max(time)+1/sampr:1/sampr:extrapfinaltime]],zdat.time,'uniformoutput',false);  % create new time axis MAKE more general ..
% datar.trial = [];
% datar.trial = fexdat;

%% fourier 

cfg              = [];
cfg.output       = 'fourier';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hann';
cfg.keeptrials   = 'yes';
cfg.foi          = f;
cfg.t_ftimwin    = cycles./cfg.foi;
cfg.toi          = t;
freqpow          = ft_freqanalysis(cfg, datar); %how to do this for downsampled data?

cd(cfg_fun.outputfile),
save('datar','datar')
save('freqpow','freqpow')




