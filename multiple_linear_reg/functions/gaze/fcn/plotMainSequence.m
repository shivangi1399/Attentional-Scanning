function [h] = plotMainSequence(saccades, varargin)
    allowedOptions = {
        'peakVelocity'
    };

    if nargin == 1; varargin = {'peakVelocity'};end
    if ~isempty(varargin)
        assert(all(ismember(varargin,allowedOptions)),...
            'Not valid option');
    end
    
    flag_peakVelocity = any(cellfun(@(x) strcmp(x,'peakVelocity'),varargin));

    
    assert(isa(saccades,'Saccade'));
        
    flag_valid = arrayfun(@(x) ~isempty(x.time), saccades);  

    saccades  = saccades(flag_valid);   
    magnitude = arrayfun(@(x) x.magnitude, saccades);
    duration  = arrayfun(@(x) x.duration, saccades);
    peakVel   = arrayfun(@(x) abs(complex(x.peakVel(1),x.peakVel(2))), saccades);   
    %plot magnitude vs peakvelocity
    if flag_peakVelocity
        h = plot(magnitude,peakVel,'.k',...
            'markersize',6)
        %set(gca,'YScale','log');
        %set(gca,'XScale','log');
        title('Main Sequence');
        xlabel('Magnitude [dva]');
        ylabel('Peak Velocity [dva/sec]');
    end
   
end