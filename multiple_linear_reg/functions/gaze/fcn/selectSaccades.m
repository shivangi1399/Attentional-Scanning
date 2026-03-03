function selectedSaccades = selectSaccades(saccades,n,cfg)
%SELECTSACCADES Select Saccade array
%   selectSaccade(saccades,n,cfg) returns maxium n saccades from saccades.
%   criteria is specifed in cfg struct, avalible options are:
%       cfg.starttime_range = [minStarttim,maxStarttime]
%       cfg.endtime_range = [minEndtime, maxEndtime]
%       cfg.magnitude_range = [minSize, maxSize]
%       cfg.direction_range = [minDir, maxDir]
%       cfg.duration_range = [minDuration, maxDuration]
%       cfg.peakVel_range = [minPeakVel, maxPeakVel]
%       cfg.startEcc_range = [minStartEcc, maxStartEcc]

    fprintf('\n[select_saccades]\n')
    fprintf('check input...\n')
    assert(isa(saccades,'Saccade'), 'Not a Saccade array!')  
    assert(isa(cfg,'struct'),'cfg not a struct!')
    
    allowedFieldNames = {
        'starttime_range'
        'endtime_range'
        'magnitude_range'
        'direction_range'
        'duration_range'
        'peakVel_range'
        'startEcc_range'
        };   
    assert(all(ismember(fieldnames(cfg),allowedFieldNames)),'cfg has invalid field(s)')
    for ii = 1 : numel(allowedFieldNames)
        try
            eval(sprintf('%s=cfg.%s;',allowedFieldNames{ii},allowedFieldNames{ii}));
            eval(sprintf('assert(isa(%s,''double''))',allowedFieldNames{ii}));
            eval(sprintf('assert(numel(%s)==2)',allowedFieldNames{ii}));
        catch
            eval(sprintf('%s=[];',allowedFieldNames{ii}))
        end
    end
    
    % apply filters    
    flag_selected = ones(numel(saccades),1);
    
    if ~isempty(starttime_range)
        fprintf('apply starttime filter [%f,%f]\n',starttime_range(1), starttime_range(2));
        fcn_isInRange = @(x) x.time(1)>=starttime_range(1) && x.time(1)<=starttime_range(2);
        flag_selected = flag_selected .* arrayfun(fcn_isInRange,saccades);
    end
    
    if ~isempty(endtime_range)
        fprintf('apply endtime filter [%f,%f]\n',endtime_range(1), endtime_range(2));
        fcn_isInRange = @(x) x.time(end)>=endtime_range(1) && x.time(end)<=endtime_range(2);
        flag_selected = flag_selected .* arrayfun(fcn_isInRange,saccades);
    end
    
    if ~isempty(magnitude_range)
        fprintf('apply size filter [%f,%f]\n',magnitude_range(1), magnitude_range(2));
        fcn_isInRange = @(x) x.magnitude>=magnitude_range(1) && x.magnitude<=magnitude_range(2);
        flag_selected = flag_selected .* arrayfun(fcn_isInRange,saccades);
    end
    
    if ~isempty(direction_range)
        fprintf('apply direction filter [%f,%f]\n',direction_range(1), direction_range(2));
        fcn_isInRange = @(x) x.direction>=direction_range(1) && x.direction<=direction_range(2);
        flag_selected = flag_selected .* arrayfun(fcn_isInRange,saccades);
    end
    
    if ~isempty(duration_range)
        fprintf('apply duration filter [%f,%f]\n',duration_range(1), duration_range(2));
        fcn_isInRange = @(x) x.duration>=duration_range(1) && x.duration<=duration_range(2);
        flag_selected = flag_selected .* arrayfun(fcn_isInRange,saccades);
    end
    
    if ~isempty(peakVel_range)
        fprintf('apply peakVel filter [%f,%f]\n',peakVel_range(1), peakVel_range(2));
        fcn_isInRange = @(x) x.peakVel>=peakVel_range(1) && x.peakVel<=peakVel_range(2);
        flag_selected = flag_selected .* arrayfun(fcn_isInRange,saccades);
    end
    
    if ~isempty(startEcc_range)
        fprintf('apply startEcc filter [%f,%f]\n',startEcc_range(1), startEcc_range(2));
        fcn_getStartEcc = @(x) abs(complex(x.position(1,1),x.position(1,2)))./x.pixPerDeg;     
        fcn_isInRange = @(x) fcn_getStartEcc(x)>=startEcc_range(1) && fcn_getStartEcc(x)<=startEcc_range(2);
        flag_selected = flag_selected .* arrayfun(fcn_isInRange,saccades);
    end
    
    
    nSaccadesRemained = sum(flag_selected);
    
    if nSaccadesRemained == 0
        fprintf('no saccades remained, empty array will be returned\n')
        selectedSaccades = Saccade.empty();        
    elseif nSaccadesRemained <= n    
        fprintf('%.0d saccades remained and will be returned\n',nSaccadesRemained)
        selectedSaccades = saccades(logical(flag_selected));
    else
        fprintf('more saccades (%.0d) than requsted (%.0d) remained...\n',nSaccadesRemained,n);
        fprintf('> Press [v] to proceed for visual selection\n')
        fprintf('> Press [a] to return all remained saccades\n')
        while (true)
            option = input('> ','s');
            switch option
                case 'v'
                    fprintf('proceed to visual selection')
                    selectedSaccades = selectSaccades_visual(saccades(logical(flag_selected)),n);
                    break
                case 'a'
                    fprintf('ok, all remained saccades will be returned')
                    selectedSaccades = saccades(logical(flag_selected));
                    break
                otherwise
                    warning('(press [v] or [a])');
                    continue
            end
        end
    end   
end