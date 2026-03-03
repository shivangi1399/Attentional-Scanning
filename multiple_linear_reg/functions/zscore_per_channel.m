function data_z = zscore_per_channel(data, active_ch)
% data: chan × trials
data_z = data;
for ichan = active_ch'
    x = data(ichan,:);
    mu = mean(x,'omitnan');
    sd = std(x,0,'omitnan');
    if sd>0
        data_z(ichan,:) = (x-mu)./sd;
    end
end
end