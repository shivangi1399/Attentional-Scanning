function [fldat,h] = subfunc_butter(dat,ord,cfreq,smo,fsample)

smo = [(cfreq-smo)' (cfreq+smo)'];

[z,p,k] = butter(ord,2*smo/fsample,'bandpass');
[sos,g] = zp2sos(z,p,k);

h = dfilt.df2sos(sos,g);

fldat  = filter(h,dat);

