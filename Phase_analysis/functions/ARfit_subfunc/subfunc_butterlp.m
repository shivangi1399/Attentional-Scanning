function [fdat] = subfunc_butterlp(dat,hfreq,ord,fsample)

[z,p,k] = butter(ord,2*hfreq/fsample,'low');
[sos,g] = zp2sos(z,p,k);

h = dfilt.df2sos(sos,g);

fdat = filter(h,dat);