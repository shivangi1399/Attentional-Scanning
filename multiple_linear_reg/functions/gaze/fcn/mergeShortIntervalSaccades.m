function a = mergeShortIntervalSaccades(saccades,maxISI)
%MERGESHORTINTERVALSACCADES
%   mergeShortIntervalSaccades(saccades) merge with default maxISI
%   mergeShortIntervalSaccades(saccades, maxISI)
    fprintf('\n[mergeShortIntervalSaccades]\n')
    assert(isa(saccades,'Saccade'),'merge failed, input 1 is not Saccade');
    assert(all(arrayfun(@(x) ~isempty(x.parent), saccades)),'not all saccades have parent');
    assert(numel(unique(arrayfun(@(x) x.parent, saccades)))==1, 'saccades have no commmon parent');
    while true
        t1 = arrayfun(@(x) x.time(1), saccades);
        t2 = arrayfun(@(x) x.time(end), saccades);
        isi = t1(2:end) - t2(1:end-1);
        if nargin == 1
            ind = find(isi<0.02);
        else
            ind = find(isi<maxISI);
        end
        
        if ~isempty(ind)
            ind2reserve = union(ind,setdiff(1:numel(saccades),union(ind,ind+1)));
            sac1Array = saccades(ind);
            sac2Array = saccades(ind+1);
            if nargin == 1
                mergeFunc = @(sac1,sac2) mergeTwoSaccades(sac1,sac2);
            else
                mergeFunc = @(sac1,sac2) mergeTwoSaccades(sac1,sac2,maxISI);
            end
            tempSaccades = arrayfun(mergeFunc,sac1Array,sac2Array);
            saccades(ind) = tempSaccades;

            saccades = saccades(ind2reserve);
            fprintf('\n%d saccades remained, continue...',numel(saccades));
        else
            break
        end
                
    end
    fprintf('\n%d saccades remained',numel(saccades));
    a = saccades;

end