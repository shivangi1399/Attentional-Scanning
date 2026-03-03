function map = back_project(data, directions, method)
% BACK_PROJECT Returns a back-projected MAP of the
%              intersected DATA.
% MAP = back_project(DATA, DIRECTIONS, <METHOD>);
% MAP is a square array (DATA length x DATA length)
%   of the intersected (averaged) DATA.
%   obs. it uses array coordinates (ij). To see (plot) it
%     use: [axis ij] or plot a flipped MAP [flipud(MAP)].
% DATA is an array of N lines, representing the vectors
%   (projections) to be intersected.
% DIRECTIONS are the N stimulus directions for the N DAT vectors,
%         in DEGREES (zero at left, counterclockwise).
% METHOD <optional>: 0 = ARITHMETIC (default);
%                1 = GEOMETRIC(*);
%                2 = PRODUCT.
% (*)negative number sigs are removed before rooting and then reassigned.
%
% mario fiorani, last modified November 2012 (first version, 2001)

if (size(data,1) ~= numel(directions)), error('dimension mismatch'), end
if (nargin == 2), method = 0; end
sz = size(data,2);
pad = ceil((sqrt(2)*sz-sz)/2);
pad_dat(1:sz + 2*pad) = nan;
cj = repmat(-(sz-1)/2:1:(sz-1)/2, sz, 1);
ci = cj';
map(1:sz,1:sz) = ~~method;
for n = 1:length(directions)
    t = directions(n)*pi/180;
    pad_dat(pad + 1:pad + sz) = data(n,:);
    rcj = round(cj*cos(t) - ci*sin(t) + ceil(pad + sz/2));
    dat_array = pad_dat(rcj);
    if method
        map = map.*dat_array; 
    else
        map = map + dat_array; 
    end
end
switch method
    case 0 
        map = map./n;
    case 1
        s = sign(map); 
        in = isnan(map); 
        map(in) = 1;
        map = abs(map).^(1/n).*s; 
        map(in) = nan;
end