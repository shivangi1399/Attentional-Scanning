function sacArray = cell2sac(sacCellArray)
% convert saccade cell array into sacArray
    nCell = cellfun(@(x) numel(x), sacCellArray);
    n = sum(nCell);
    sacArray(n,1) = Saccade();
    counter = 1;
    for iCell = 1 : nCell
        try
            sacArray(counter:counter+nCell(iCell)-1) = sacCellArray{iCell};
        catch
            warning('Cell %d conversion failed, leave empty Saccade instead',ii);
        end
            counter = counter + nCell(iCell);
    end
end