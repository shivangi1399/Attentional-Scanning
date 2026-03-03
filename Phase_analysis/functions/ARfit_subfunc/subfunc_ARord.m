function [A,C,sbc,fpe] = subfunc_ARord(dat,arordmin,arordmax)

% Data should ideally go in a time x channel x trial form. If data is 
% concatenated noise variance will be overestimated

dat = vertcat(dat{:});
dat = dat';

dat = reshape(dat,[size(dat,1) 1 size(dat,2)]);

[~,A,C,sbc,fpe,~] = arfit(dat,arordmin,arordmax,'sbc');

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


