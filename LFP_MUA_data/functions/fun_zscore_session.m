% fun_zscore_session
% -------------------
% function to zscore session over all trials
% this step should be done befrore averaging data from different sessions

% input: data in fieldtrip format
% output: zscored data in fieldtrip format

% steps:
% ------
% 1. concatenate all trials
% 2. compute SD and mean over ALL trials of the session
% 3. Use this SD and mean to z-score every single trial


function zdata = fun_zscore_session(data)

% 1. concatenate all trials
SD_chan  = nan(length(data.label),1); 
Avg_chan = nan(length(data.label),1); 

for ichan = 1:length(data.label)

    temp = [];
    temp = data.trial{1,1}(ichan,:);
    for itrial = 2:length(data.trial)
        temp = [temp data.trial{1,itrial}(ichan,:)];
    end

    % 2. compute SD and mean per channel over ALL trials of the session

    SD_chan(ichan,1)  = std(temp,0,2,"omitnan");
    Avg_chan(ichan,1) = mean(temp,2,'omitnan');

end

clear temp

% 3. Use this SD and mean to z-score every single trial
% zdata = (data-avg)/sd
z_trials = cell(1,length(data.trial));
z_trials = cellfun(@(x) (x-Avg_chan)./SD_chan, data.trial,'uniformoutput',0);

zdata = data;

clear data


zdata.trial = [];
zdata.trial = z_trials; % zscored data in fieldtrip format





% visual inspection
% figure
% for itrial = 1:length(data.trial)
%     
%     plot(zdata.time{1,itrial},zdata.trial{1,itrial}(30,:))
%     title(num2str(data.trialinfo(itrial,20))) % behavioral outcome
%     
%     pause
%     
%     clf
% end


