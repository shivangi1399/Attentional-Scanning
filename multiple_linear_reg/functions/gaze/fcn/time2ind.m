function ind = time2ind(time,timeSeries)
    distance = abs(timeSeries - time);
    [~,ind] = min(distance);
end