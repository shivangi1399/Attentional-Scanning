function [met] = subfunc_pos(inp)

% See VanRullen 2016 Frontiers in Neuroscience

v1 = inp.v1;
v2 = inp.v2;

clear inp

met = 0;

% itc1 = abs(subfunc_ITC(v1));
% itc2 = abs(subfunc_ITC(v2));
% itca = abs(subfunc_ITC([v1 ; v2]));
% 
% met = itc1+itc2-2.*itca;
% 
% met = squeeze(met);

itc1 = abs(subfunc_ITC(v1))./size(v1,1);
itc2 = abs(subfunc_ITC(v2))./size(v2,1);
itca = abs(subfunc_ITC([v1 ; v2]))./size([v1;v2],1);

met = itc1+itc2-2.*itca;

met = squeeze(met);