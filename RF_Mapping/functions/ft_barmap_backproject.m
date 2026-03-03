function [allMaps, bestLatency] = ft_barmap_backproject(data, directionCol,latencies)
% ft_barmap_backproject Calculate response maps to moving bars 
% 
% This function applies Mario Fiorani's back_project to a Fieldtrip data set.
%
%   [allMaps, bestLatency] = ft_barmap_backproject(data, directionCol, {latencies})
%
% INPUT
% -----
%   data : Fieldtrip data struct similar to ft_datatype_raw.
%          Bar sweep direction must be provided in a column of the 
%          data.trialinfo matrix (degrees, 0=rightward, 180=leftward,
%          counter-clockwise). Duration of all trials must be identical, 
%          with the trial start being bar onset and trial end the offset.
%   directionCol : index of column in data.trialinfo that specifies the bar
%                  direction for each trial
%   latencies : optional, vector of response latencies (s) to be tested,
%               default=0.030:0.002:0.070
%               For each value the data is shifted backwards and the
%               response map calculated. The latency with the maximum
%               overall response will be selected.
%
% OUTPUT
% -----
%    allMaps : [nSamples x nSamples x nChannel] matrix of back-projected
%              responses. The square axes of the map correspond to the time
%              axis of the trials.
%    bestLatency : vector of maximum response latencies for each channel
%
%
%
% See also back_project, ft_datatype_raw

ft_defaults

if nargin == 2
    latencies = 0.040:0.005:0.070;    
end

directions = unique(data.trialinfo(:, directionCol));

avgCfg = [];
tl = cell(1, length(directions));
for iDirection = 1:length(directions)
    
    % select trials based on data.trialinfo
    currentDirection = directions(iDirection);
    selector = data.trialinfo(:,directionCol) == currentDirection;
    avgCfg.trials = find(selector)';
    
    % average across trials
    tl{iDirection} = ft_timelockanalysis(avgCfg, data);
end

% created back-projected maps for all latencies
allMaps = zeros(length(tl{1}.time), length(tl{1}.time), length(data.label));
bestLatency = zeros(1, length(data.label));

for iChannel = 1:length(data.label)
    
    % create [nDirections x nSamples] array for back_project function
    projections = cellfun(@(x) x.avg(iChannel,:), tl, 'unif', false);
    projections = cat(1, projections{:});
    
    map = zeros(size(projections, 2));
    for iLatency = 1:length(latencies)
        
        % shift data by latency
        shiftSamples = round(latencies(iLatency)*data.fsample);
        shiftedProjections = circshift(projections, -shiftSamples, 2);
        shiftedProjections(end-shiftSamples:end) = 0;
        
        % calculate back-projected 2D map
        newMap = back_project(shiftedProjections, directions);
        
        % use new map if maximum response is higher than previous
        if max(newMap(:)) > max(map(:))
            map = newMap;
            bestLatency(iChannel) = latencies(iLatency);
        end
    end
  
    allMaps(:,:,iChannel) = map;
end
