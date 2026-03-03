function perm_mua_lfp_coherence(cfg)
    cfg_freq = [];
    cfg_freq.output     = 'fourier';
    cfg_freq.method     = 'mtmfft';
    cfg_freq.foilim     = cfg.foilim;
    cfg_freq.tapsmofrq  = cfg.tapsmofrq;
    cfg_freq.keeptrials = 'yes';
    cfg_freq.channel    = cfg.channel;
    cfg_freq.channelcmb = cfg.channelcmb;

    freq_shuff = ft_freqanalysis(cfg_freq, cfg.data);

    cfg_coh = [];
    cfg_coh.method = 'coh';
    cfg_coh.channelcmb = cfg.channelcmb;

    fd_shuff = ft_connectivityanalysis(cfg_coh, freq_shuff);

    save(cfg.outputfile, 'fd_shuff');
end