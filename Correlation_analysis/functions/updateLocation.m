function updateLocation(direction)
currentLocIdx = currentLocIdx + direction;
if currentLocIdx < 1
    currentLocIdx = length(targ_loc);
elseif currentLocIdx > length(targ_loc)
    currentLocIdx = 1;
end

set(hImg,'CData',sig_grids{currentLocIdx});
title(ax, sprintf('Significant Channels - Loc %d Dlev %d-%d', ...
    targ_loc(currentLocIdx), min_difflev(currentLocIdx), max_difflev(currentLocIdx)));
drawnow
end