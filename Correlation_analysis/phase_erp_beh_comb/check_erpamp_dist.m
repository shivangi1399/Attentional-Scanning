clear all
close all
clc

%% Dependencies
addpath /opt/fieldtrip_github/
ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/
clc

%% Create data paths
datafolder = '/mnt/hpc/projects/MWSampling/4Shivangi/results_klecks';
cd(datafolder)
animalName = 'klecks';
temp = dir;
session_names = [];
ii = 0;
for i = 1:length(temp)
    if strfind(temp(i).name,animalName)
        ii = ii+1;
        session_names{ii,1} = temp(i).name;
    end
end

%% Load phase amp data

phase_folder = fullfile(datafolder, 'phase_coherence', 'cp10_till_100');
cd(phase_folder);
load('ph_all_sess.mat')

lfp_ampl = ph_comb.LFP_ERP_ampl_all;   % [nTrials x nChannels]
mua_ampl = ph_comb.MUA_ERP_ampl_all;   % [nTrials x nChannels]
nChannels = size(lfp_ampl,2);

% Colors
color_data    = [0 128 128]/255; % dark turquoise
color_outlier = [128 0 0]/255;   % maroon
line_color    = [0 0 0];         % black

% Output folder
outfolder = '/mnt/hpc/projects/MWSampling/4Shivangi/Plots/coherence/klecks/cp10_till_100';
if ~exist(outfolder,'dir'), mkdir(outfolder); end

%% Boxplots per channel
fig1 = figure('Name','Amplitude Boxplots - All Sessions','Position',[100 100 2000 800]);

% LFP boxplot
subplot(1,2,1);
boxplot(lfp_ampl, 'Labels', 1:nChannels, 'Symbol','', 'Colors','k'); % hide default outliers
hold on;
h = findobj(gca,'Tag','Box');
for j = 1:length(h)
    patch(get(h(j),'XData'), get(h(j),'YData'), color_data,'FaceAlpha',0.3,'EdgeColor','k','LineStyle','-');
end
for ch = 1:nChannels
    data_ch = lfp_ampl(:,ch);
    data_ch = data_ch(~isnan(data_ch));
    if isempty(data_ch), continue; end
    p1 = prctile(data_ch,1);      % lower percentile
    p99 = prctile(data_ch,99);    % upper percentile
    out_idx = find(data_ch < p1 | data_ch > p99);
    scatter(ch*ones(size(out_idx)), data_ch(out_idx), 50, color_outlier, '+', 'LineWidth',1.5);
    m   = mean(data_ch); med = median(data_ch); mo  = mode(data_ch);
    xj = [ch-0.4, ch+0.4];
    plot(xj, [m m],'k:','LineWidth',2); plot(xj, [med med],'k-.','LineWidth',2); plot(xj, [mo mo],'k--','LineWidth',2);
    plot([ch ch],[p1 p99],'k:','LineWidth',1.5);
end
title('LFP Amplitude - All Sessions','FontSize',14); xlabel('Channel','FontSize',12); ylabel('Amplitude','FontSize',12);
legend({'Box','Mean','Median','Mode','1st-99th percentile outliers'},'Location','northeast'); 
set(gca,'FontSize',10); xtickangle(45); grid on; hold off;

% MUA boxplot
subplot(1,2,2);
boxplot(mua_ampl, 'Labels', 1:nChannels, 'Symbol','', 'Colors','k'); % hide default outliers
hold on;
h = findobj(gca,'Tag','Box');
for j = 1:length(h)
    patch(get(h(j),'XData'), get(h(j),'YData'), color_data,'FaceAlpha',0.3,'EdgeColor','k','LineStyle','-');
end
for ch = 1:nChannels
    data_ch = mua_ampl(:,ch);
    data_ch = data_ch(~isnan(data_ch));
    if isempty(data_ch), continue; end
    p1 = prctile(data_ch,1);      % lower percentile
    p99 = prctile(data_ch,99);    % upper percentile
    out_idx = find(data_ch < p1 | data_ch > p99);
    scatter(ch*ones(size(out_idx)), data_ch(out_idx), 50, color_outlier, '+', 'LineWidth',1.5);
    m   = mean(data_ch); med = median(data_ch); mo  = mode(data_ch);
    xj = [ch-0.4, ch+0.4];
    plot(xj, [m m],'k:','LineWidth',2); plot(xj, [med med],'k-.','LineWidth',2); plot(xj, [mo mo],'k--','LineWidth',2);
    plot([ch ch],[p1 p99],'k:','LineWidth',1.5);
end
title('MUA Amplitude - All Sessions','FontSize',14); xlabel('Channel','FontSize',12); ylabel('Amplitude','FontSize',12);
legend({'Box','Mean','Median','Mode','1st-99th percentile outliers'},'Location','northeast'); 
set(gca,'FontSize',10); xtickangle(45); grid on; hold off;

% Save boxplot figure
saveas(fig1, fullfile(outfolder,'Amplitude_Boxplots_AllSessions.png'));

%% Histograms
fig2 = figure('Name','Amplitude Histograms - All Sessions','Position',[100 100 1600 700]);

% LFP histogram
subplot(1,2,1);
histogram(lfp_ampl(:),100,'FaceColor',color_data,'EdgeColor','none'); hold on;
data_ch = lfp_ampl(:); data_ch = data_ch(~isnan(data_ch));
p1 = prctile(data_ch,1);      
p99 = prctile(data_ch,99);    
m   = mean(data_ch); med = median(data_ch); mo  = mode(data_ch);
h1 = xline(m,'k:','LineWidth',2,'DisplayName','Mean'); 
h2 = xline(med,'k-.','LineWidth',2,'DisplayName','Median'); 
h3 = xline(mo,'k--','LineWidth',2,'DisplayName','Mode');
h4 = xline(p1,'Color',color_outlier,'LineWidth',2,'LineStyle','--','DisplayName','1st percentile'); 
h5 = xline(p99,'Color',color_outlier,'LineWidth',2,'LineStyle','--','DisplayName','99th percentile');
legend([h1,h2,h3,h4,h5],'Location','northeast');
title('LFP Histogram - All Sessions','FontSize',14); xlabel('Amplitude','FontSize',12); ylabel('Count','FontSize',12); grid on; hold off;

% MUA histogram
subplot(1,2,2);
histogram(mua_ampl(:),100,'FaceColor',color_data,'EdgeColor','none'); hold on;
data_ch = mua_ampl(:); data_ch = data_ch(~isnan(data_ch));
p1 = prctile(data_ch,1);      
p99 = prctile(data_ch,99);    
m   = mean(data_ch); med = median(data_ch); mo  = mode(data_ch);
h1 = xline(m,'k:','LineWidth',2,'DisplayName','Mean'); 
h2 = xline(med,'k-.','LineWidth',2,'DisplayName','Median'); 
h3 = xline(mo,'k--','LineWidth',2,'DisplayName','Mode');
h4 = xline(p1,'Color',color_outlier,'LineWidth',2,'LineStyle','--','DisplayName','1st percentile'); 
h5 = xline(p99,'Color',color_outlier,'LineWidth',2,'LineStyle','--','DisplayName','99th percentile');
legend([h1,h2,h3,h4,h5],'Location','northeast');
title('MUA Histogram - All Sessions','FontSize',14); xlabel('Amplitude','FontSize',12); ylabel('Count','FontSize',12); grid on; hold off;

% Save histogram figure
saveas(fig2, fullfile(outfolder,'Amplitude_Histograms_AllSessions.png'));
