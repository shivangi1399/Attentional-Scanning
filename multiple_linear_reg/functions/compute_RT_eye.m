function out = compute_RT_eye(cfg) % based on Zhang2025

ESIload(cfg.eye_file)   % loads eyedata
eyeTrials = eyeData;
nTrl = numel(eyeTrials.trial);
RT   = nan(nTrl,1);

% visual inspection
% figure
% for itrial = 1:size(eyeTrials.trialinfo,1)
%     plot(eyeTrials.time{1,itrial}, eyeTrials.trial{1,itrial}) % plots x and y for both eyes
%     title(sprintf('Trial %d', itrial))
%     t = eyeTrials.time{1,itrial};
%     xticks(min(t):0.1:max(t))
%     pause
%     clf
% end


RT = NaN(nTrl,1);  % preallocate all RTs as NaN

for iTrl = 1:nTrl
    fprintf('Processing trial %d\n', iTrl);

    % keep hits only
    if eyeTrials.trialinfo(iTrl,15) ~= 1 || ...
       ~(eyeTrials.trialinfo(iTrl,20)==1 || eyeTrials.trialinfo(iTrl,20)==5)
        continue
    end

    screenX = 1680;
    screenY = 1050;
    ppd     = eyeTrials.trialinfo(iTrl,30);

    eyeTrace = eyeTrials.trial{iTrl}(2:3,:)' * [screenX/20, 0; 0, screenY/20]; % converts volts -> pixel
    eyeTime  = eyeTrials.time{iTrl}';

    gazeEpoch = GazeEpoch(eyeTrace, eyeTime, ppd); % converts pixels -> degrees

    % Detect saccades
    saccades = detectSaccades(gazeEpoch, cfg.cfg_detectSac);

    if ~isempty(saccades)
        % Merge short-interval saccades only if they exist
        saccades = mergeShortIntervalSaccades(saccades);

        % Select target saccade
        tarSac = selectSaccades(saccades, 1, cfg.cfg_selSac);
        % reaction time - calculated as the duration b/w when saccade is above and below the EK threshold
        if ~isempty(tarSac)
            RT(iTrl) = tarSac(1).time(1); % this gives the first time point when the saccade was initiated
            
        else
            % no target saccade after selection
            RT(iTrl) = NaN;
            warning('No target saccade selected in trial %d', iTrl);
        end
    else
        % no saccades detected at all
        RT(iTrl) = NaN;
        warning('No saccades detected in trial %d', iTrl);
    end
end


out.RT        = RT;
out.trialinfo = eyeTrials.trialinfo;

end
