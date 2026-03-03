function amp_z = zscore_amp_per_freq(amp, active_ch)
% amp: trials × freq × chan
amp_z = amp;
for ichan = active_ch'
    for ifreq = 1:size(amp,2)
        x = squeeze(amp(:,ifreq,ichan));
        mu = mean(x,'omitnan');
        sd = std(x,0,'omitnan');
        if sd>0
            amp_z(:,ifreq,ichan) = (x-mu)./sd;
        end
    end
end
end


