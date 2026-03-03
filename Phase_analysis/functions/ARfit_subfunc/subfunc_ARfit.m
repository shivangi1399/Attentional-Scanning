function [W,A,C,sbc,fpe,th] = subfunc_ARfit(datf,arord)

% Data should ideally go in a time x channel x trial form. If data is 
% concatenated noise variance will be overestimated

datf = vertcat(datf{:});
datf = datf';

datf = reshape(datf,[size(datf,1) 1 size(datf,2)]);

[W,A,C,sbc,fpe,th] = arfit(datf,arord,arord,'sbc');

% cfg = [];
% cfg.order = arord;
% cfg.demean = 'no';
% 
% fsample = 1000;
% 
% data.trial = {dat};
% data.label = {'1'};
% data.time  = {0:(1/fsample):((numel(dat)-1)/fsample)};
% data.fsample = fsample;
% data.trialdef = [1 numel(dat) 0];
% 
% [mvardat] = ft_mvaranalysis(cfg, data);
% 
% arcoeff=squeeze(mvardat.coeffs)';
% noisecoeff=mvardat.noisecov;


