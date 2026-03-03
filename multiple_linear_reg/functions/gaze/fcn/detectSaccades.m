function saccades = detectSaccades(gazeEpoch,cfg)
%GET_SACCADES   Detect saccades in GazeEpoch instance
%   saccades = get_saccades(gazeEpoch,cfg) returns a Saccade array as
%   children of the gazeEpoch.
%
%   cfg structure must have the following fields
%   method-> 'Engbert2003'
%   params-> {param1, param2, ...}
%    
%   `Engbert2003` parames:
%   {lambda}
   
    fprintf('[detect_saccades]\n')
    fprintf('check input...\n')
    assert(isa(gazeEpoch,'GazeEpoch'), 'Not a GazeEpoch instance!')
    
    if nargin == 1
        cfg.method = 'Engbert2003';
        cfg.params = {6};
    end
    
    assert(isa(cfg,'struct'),'cfg not a struct!')
    try
        assert(all(isfield(cfg, {...
            'method'
            'params'})));
        assert(ismember(cfg.method,{'Engbert2003'}))
        assert(isa(cfg.params,'cell'));
        switch cfg.method
            case 'Engbert2003'
                assert(length(cfg.params)==1);
                assert(isa(cfg.params{1},'double'));
                fprintf('use method `%s` with ', cfg.method)
                fprintf('lambda = %.1f\n',cfg.params{1})
        end             
    catch
        error('Invalid cfg!');
    end
        
    switch cfg.method
        case 'Engbert2003'
            lambda = cfg.params{1};
            fcn_is_outside_ellipse = @(a,b,x0,y0,x,y) (((x-x0)/a).^2 + ((y-y0)/b).^2) > 1;
            
            v = gazeEpoch.velocity;
            v_median = median(v,1,'omitnan');
            v_sqr = v .^2;
            v_sqr_median = median(v_sqr,1,'omitnan');            
            v_sigma = sqrt(v_sqr_median - v_median.^2);
            eta = lambda *v_sigma;
            
            flag_during_saccade = fcn_is_outside_ellipse(eta(1),eta(2),v_median(1),v_median(2),...
                v(:,1),v(:,2));            
            saccade_start = find(diff(flag_during_saccade)==1)+1;
            saccade_end   = find(diff(flag_during_saccade)==-1);                       
    end
    
    try
        if saccade_start(1)<= saccade_end(1) % start-end-start-end-...
            nSaccades = length(saccade_end);
            saccade_start = saccade_start(1:nSaccades);                             
        else % end-start-end-...
            nSaccades = length(saccade_end)-1;
            saccade_start = saccade_start(1:nSaccades);
            saccade_end = saccade_end(2:end);
        end
    catch
        nSaccades = 0;
    end
    
    % init Saccade array
    fprintf('%d saccades detected',nSaccades)
    
    clear saccades
    if nSaccades == 0
        saccades = Saccade.empty;
    else
        saccades(nSaccades,1) = Saccade();
        % set parent gazeEpoch
        arrayfun(@(x) x.setParent(gazeEpoch),...
            saccades);
        % fetch data from parent
        arrayfun(@(x,j,k) x.fetchDataFromParent(j,k),...
            saccades,saccade_start, saccade_end);
        % log
            arrayfun(@(x) x.logInfo(cfg),saccades),...        
    end
    
    fprintf('\n')
end
