function [corr_perm] = circlin_regression_perm(cfg_fun)

ichan = cfg_fun.ichan;
permut_n = cfg_fun.permut_n;
perm_indices_hits = cfg_fun.perm_ind_hits;
perm_indices_miss = cfg_fun.perm_ind_miss;
hits2use = cfg_fun.hits2use;
miss2use = cfg_fun.miss2use;

cd(cfg_fun.infile)
load('ph_all_sess.mat')

ERP_hits = []; ERP_hits = ph_comb.ERP_ampl_all(hits2use,:);
DE_hits = []; DE_hits = ph_comb.ERP_trialinfo(hits2use,18);
ERP_miss = []; ERP_miss = ph_comb.ERP_ampl_all(miss2use,:);
DE_miss = []; DE_miss = ph_comb.ERP_trialinfo(miss2use,18);

nFreq = size(ph_comb.phase_all, 2);
corr_perm = nan(permut_n, nFreq);
pvalue_perm = nan(permut_n, nFreq);

for perm = 1:permut_n
    
    % hits reg
    tbl = [];
    tbl = table(DE_hits,ERP_hits(:,ichan));
    model_hits  = fitlm(tbl);
    
    % misses reg
    tbl = [];
    tbl = table(DE_miss,ERP_miss(:,ichan));
    model_miss  = fitlm(tbl);
    
    % perm
    temp = [];
    temp = table2array(model_hits.Residuals(:,1));
    hits2use = find(~isnan(temp));
    ampl_resdl = temp(hits2use);
    ampl_resdl_perm = ampl_resdl(perm_indices_hits{perm}); % same shuffle for all channels
    
    temp = [];
    temp = table2array(model_miss.Residuals(:,1));
    miss2use = find(~isnan(temp));
    temp_m = temp(miss2use);
    ampl_resdl_perm = [ampl_resdl_perm; temp_m(perm_indices_miss{perm})]; % same shuffle for all channels
    
    % correlation
    for ifreq = 1:nFreq
        
        phase2use = [ph_comb.phase_all(hits2use,ifreq,ichan); ph_comb.phase_all(miss2use,ifreq,ichan)];
        
        if~isempty(ampl_resdl)
            [corr_perm(perm,ifreq),pvalue_perm(perm,ifreq)] = circ_corrcl(phase2use, ampl_resdl_perm);
        end
    end
end

cd(cfg_fun.outfile)
if ~exist(num2str(ichan), 'dir')
    mkdir(num2str(ichan))
end
cd(num2str(ichan))

ESIsave corr_perm corr_perm
ESIsave pvalue_perm pvalue_perm
ESIsave perm_indices_hits perm_indices_hits
ESIsave perm_indices_miss perm_indices_miss
end
