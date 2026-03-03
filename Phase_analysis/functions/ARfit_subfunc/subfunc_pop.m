function [met] = subfunc_pop(inp)

% See VanRullen 2016 Frontiers in Neuroscience

v1 = inp.v1;
v2 = inp.v2;

clear inp

met = 0;

itc1 = abs(subfunc_ITC(v1));
itc2 = abs(subfunc_ITC(v2));
itca = abs(subfunc_ITC([v1 ; v2]));

met = (itc1.*itc2)-itca.^2;

met = squeeze(met);