% Description:
% -------------
% Goal: combine phases from all sessions that have been estimated before
% ----
% if concatenate_phases:
% Load the pre-computed spectra per trial and session and save in a common
% structure that also includes the trialinfo
% ----
% if concatenate_lfp_erp_amp:
% estimate erp amplitude and save in a common structure with the phase
% ----
% if concatenate_mua_erp_amp:
% estimate mua amplitude and save in a common structure with the phase
% ----
% if concatenate_mua_baseline:
% estimate mua baseline amplitude and save in a common structure with the phase
% ----
% if concatenate_pupil_baseline:
% estimate pupil amplitude and save in a common structure with the phase
% ----
% if concatenate_RT:
% load estimated reaction time and save in a common structure with the phase

clearvars
close all
clc

%% Dependecies

addpath /opt/fieldtrip_github/
ft_defaults
addpath /opt/ESIsoftware/matlab/tdt_preprocessing/
addpath /mnt/hpc/opt/ESIsoftware/matlab/esi-nbf
addpath /opt/ESIsoftware/matlab/slurmfun/
addpath /mnt/hpc/projects/MWSampling/4Shivangi/
clc

%% logicals

concatenate_phases = 1;
concatenate_lfp_erp_amp = 1;
concatenate_mua_erp_amp = 1;
concatenate_mua_baseline = 1;
concatenate_pupil_baseline = 1;
concatenate_RT = 1;

%% Create data paths

datafolder   = '/mnt/hpc/projects/MWSampling/4Shivangi/results_hermes';

cd(datafolder),
animalName = 'hermes';
temp = dir;
session_names = [];
ii = 0;
for i = 1:length(temp)
    if strfind(temp(i).name,animalName)
        ii = ii+1;
        session_names{ii,1} = temp(i).name;
    end
end

session_paths_files_lfp = [];
session_paths_files_lfp = cellfun(@(x) fullfile(datafolder,x, 'clean_lfp.mat'), session_names, 'uniform',0);

session_paths_files_mua = [];
session_paths_files_mua = cellfun(@(x) fullfile(datafolder,x, 'clean_mua.mat'), session_names, 'uniform',0);

session_paths_files_pup = [];
session_paths_files_pup = cellfun(@(x) fullfile(datafolder,x, 'pupData.mat'), session_names, 'uniform',0);

session_paths_files_RT = [];
session_paths_files_RT = cellfun(@(x) fullfile(datafolder,x, 'RT_sess.mat'), session_names, 'uniform',0);

phase_paths = cellfun(@(x) fullfile(datafolder, x,'Phase_analysis/hit_miss'),session_names, 'uniform',0);
info_folder = '/mnt/hpc/projects/MWSampling/4Shivangi/results_hermes/multi_lin_reg/cp10_till_100';

%% load channel per session matrix

% cd(datafolder)
% load('V4_lfp_data.mat')
% 
% all_channels = lfpTrials.label(:);
% num_channels = length(all_channels);
% num_sessions = length(phase_paths);
% 
% channels_sessions = NaN(num_channels, num_sessions);
% 
% for isess = 1:num_sessions
%     ESIload(session_paths_files_lfp{isess});
% 
%     for ichan = 1:num_channels
%         if any(strcmp(all_channels{ichan}, clean_data.label))
%             channels_sessions(ichan, isess) = 1;
%         end
%     end
% end
% 
% cd(datafolder)
% save channels_sessions channels_sessions

%% load critical time per channel matrix

cd('/mnt/hpc/projects/MWSampling/4Shivangi/results_hermes/critical_time')
load('CriticalTime.mat')
load('all_channels.mat')

cd(datafolder)
load('V4_lfp_data.mat')
full_channels = lfpTrials.label(:);

CriticalTime_mat = NaN(length(full_channels), 1);

for i = 1:length(all_channels)
    idx = find(strcmp(full_channels, all_channels{i}));
    if ~isempty(idx)
        CriticalTime_mat(idx) = CriticalTime(i);
    end
end

%% combine phases across session along with erp amplitude

cd(datafolder)
load('channels_sessions.mat')
chan_orig = 1:64;
cut_here = 0.010;

if concatenate_phases

    ph_all = []; amp_all = []; trlinfo_all = [];
    session_name_all = {};
    LFP_ERP_ampl_all = []; LFP_ERP_trialinfo = [];
    MUA_ERP_ampl_all = []; MUA_ERP_trialinfo = [];
    MUA_baseline = []; MUA_baseline_trialinfo  = [];
    pup_baseline = []; pup_baseline_trialinfo  = [];
    RT = []; RT_trialinfo = [];

    for isess = 1:length(phase_paths)
        isess

        %% load phase estimations

        % loop through sessions and concatenate the phases at cp over all sessions

        cd(fullfile(phase_paths{isess}, '100iter_cut@cp_m10/phase'));

        % list of available channels
        dr = dir;
        dr = dr(3:end);
        chan_nums = sort(str2double({dr.name}));

        % Preallocate for all channels
        cd(fullfile(phase_paths{isess}, '100iter_cut@cp_m10/phase', num2str(chan_nums(1))));
        load('phase.mat');
        ph = nan([size(phase.ar_phase), length(chan_orig)]); % size: (samples × trials × 64)
        amp = nan([size(phase.amp_vec), length(chan_orig)]);

        % Loop only over active channels
        active_ch = find(channels_sessions(:, isess) == 1);
        for idx = 1:length(active_ch)
            ichan = active_ch(idx);
            cd(fullfile(phase_paths{isess}, '100iter_cut@cp_m10/phase', num2str(chan_nums(idx))));
            load('phase.mat');
            ph(:, :, ichan) = phase.ar_phase;
            amp(:, :, ichan) = phase.amp_vec;
        end

        amp = zscore_amp_per_freq(amp, active_ch); % Zscore amplitude per freq × chan × session

        ph_all = [ph_all; ph];
        amp_all = [amp_all; amp];
        trlinfo_all = [trlinfo_all; phase.trialinfo];
        session_name_all = [session_name_all; repmat({session_names{isess}}, size(phase.trialinfo, 1), 1)];

        %% ERP amplitude lfp

        if concatenate_lfp_erp_amp

            ESIload(session_paths_files_lfp{isess}) % load broadband data for ERP-amplitude estimation

            [~, ia_phase, ib_clean] = intersect(phase.trialinfo(:,14), clean_data.trialinfo(:,14));%checks the trial numbers
            nTrials = length(ia_phase);

            ERP_ampl = nan(length(chan_orig), nTrials);

            active_ch = find(channels_sessions(:, isess) == 1);
            for ichan = 1:length(active_ch')
                ct = CriticalTime_mat(active_ch(ichan))-cut_here;
                if isnan(ct), continue; end

                cfg = [];
                cfg.channel = ichan;
                cfg.latency = [ct 0.1];
                data_chan = ft_selectdata(cfg, clean_data);

                % ERP amplitude (RMS)
                ampl = cellfun(@(x) sqrt(mean((x(:) - mean(x(:))).^2)), data_chan.trial);

                % Align trials
                [~, ia_chan, ib_chan] = intersect(phase.trialinfo(:,14), data_chan.trialinfo(:,14));

                % Store
                ch_data = nan(1, nTrials);
                ch_data(ia_chan) = ampl(ib_chan);
                ERP_ampl(active_ch(ichan),:) = ch_data;
            end

            ERP_ampl = zscore_per_channel(ERP_ampl, active_ch);

            LFP_ERP_ampl_all  = [LFP_ERP_ampl_all; ERP_ampl'];
            LFP_ERP_trialinfo = [LFP_ERP_trialinfo; clean_data.trialinfo(ib_clean,:)];
        end

        %% Baseline mua

        if concatenate_mua_baseline

            ESIload(session_paths_files_mua{isess})

            [~, ia_phase, ib_clean] = intersect(phase.trialinfo(:,14), clean_mua.trialinfo(:,14));%checks the trial numbers
            nTrials = length(ia_phase);

            ERP_ampl = nan(length(chan_orig), nTrials);

            active_ch = find(channels_sessions(:, isess) == 1);
            for ichan = 1:length(active_ch')
                ct = CriticalTime_mat(active_ch(ichan))-cut_here;
                if isnan(ct), continue; end

                cfg = [];
                cfg.channel = ichan;
                cfg.latency = [-0.8 0];
                data_chan = ft_selectdata(cfg, clean_mua);

                % ERP amplitude
                ampl = cellfun(@(x) sum(x(:)), data_chan.trial);

                % Align trials
                [~, ia_chan, ib_chan] = intersect(phase.trialinfo(:,14), data_chan.trialinfo(:,14));

                % Store
                ch_data = nan(1, nTrials);
                ch_data(ia_chan) = ampl(ib_chan);
                ERP_ampl(active_ch(ichan),:) = ch_data;
            end

            ERP_ampl = zscore_per_channel(ERP_ampl, active_ch);

            MUA_baseline  = [MUA_baseline; ERP_ampl'];
            MUA_baseline_trialinfo = [MUA_baseline_trialinfo; clean_mua.trialinfo(ib_clean,:)];
        end

        %% ERP amplitude mua

        if concatenate_mua_erp_amp

            ESIload(session_paths_files_mua{isess})

            [~, ia_phase, ib_clean] = intersect(phase.trialinfo(:,14), clean_mua.trialinfo(:,14));%checks the trial numbers
            nTrials = length(ia_phase);

            ERP_ampl = nan(length(chan_orig), nTrials);

            active_ch = find(channels_sessions(:, isess) == 1);
            for ichan = 1:length(active_ch')
                ct = CriticalTime_mat(active_ch(ichan))-cut_here;
                if isnan(ct), continue; end

                cfg = [];
                cfg.channel = ichan;
                cfg.latency = [ct 0.1];
                data_chan = ft_selectdata(cfg, clean_mua);

                % ERP amplitude
                ampl = cellfun(@(x) sum(x(:)), data_chan.trial);

                % Align trials
                [~, ia_chan, ib_chan] = intersect(phase.trialinfo(:,14), data_chan.trialinfo(:,14));

                % Store
                ch_data = nan(1, nTrials);
                ch_data(ia_chan) = ampl(ib_chan);
                ERP_ampl(active_ch(ichan),:) = ch_data;
            end

            ERP_ampl = zscore_per_channel(ERP_ampl, active_ch);

            MUA_ERP_ampl_all  = [MUA_ERP_ampl_all; ERP_ampl'];
            MUA_ERP_trialinfo = [MUA_ERP_trialinfo; clean_mua.trialinfo(ib_clean,:)];
        end

        %% Baseline pupil data

        if concatenate_pupil_baseline

            ESIload(session_paths_files_pup{isess})

            [~, ia_phase, ib_clean] = intersect(phase.trialinfo(:,14), pupData.trialinfo(:,11)); %checks the trial numbers
            nTrials = length(ia_phase);

            pup_ampl = nan(length(chan_orig), nTrials);

            active_ch = find(channels_sessions(:, isess) == 1);
            for ichan = 1:length(active_ch')
                ct = CriticalTime_mat(active_ch(ichan))-cut_here;
                if isnan(ct), continue; end

                % zscoring
                cfg = [];
                cfg.latency = [-1 0];
                data_chanr = ft_selectdata(cfg, pupData);

                t1 = cellfun(@(x) x(1,:), data_chanr.trial, 'UniformOutput', false); t1 = [t1{:}];
                t2 = cellfun(@(x) x(2,:), data_chanr.trial, 'UniformOutput', false); t2 = [t2{:}];

                m1 = mean(t1,'omitnan');  s1 = std(t1,0,'omitnan'); s1(s1==0) = NaN;
                m2 = mean(t2,'omitnan');  s2 = std(t2,0,'omitnan'); s2(s2==0) = NaN;

                zdata = data_chanr;
                zdata.trial = cellfun(@(x) ...
                    [(x(1,:)-m1)./s1; (x(2,:)-m2)./s2], ...
                    data_chanr.trial, 'uni', 0);

                % Mean pupil amplitude
                cfg = [];
                cfg.latency = [-0.8 0];
                data_chan = ft_selectdata(cfg, zdata);

                ampl = cellfun(@(x) mean(abs(x(1:2,:)), 'all'), data_chan.trial); %zscore within sess and avg 
                % them for 200ms before stimonset, the average is done of the abs of the signal

                % Align trials
                [~, ia_chan, ib_chan] = intersect(phase.trialinfo(:,14), data_chan.trialinfo(:,11));

                % Store
                ch_data = nan(1, nTrials);
                ch_data(ia_chan) = ampl(ib_chan);
                pup_ampl(active_ch(ichan),:) = ch_data;
            end

            pup_ampl = zscore_per_channel(pup_ampl, 1);

            pup_baseline  = [pup_baseline; pup_ampl'];
            pup_baseline_trialinfo = [pup_baseline_trialinfo; pupData.trialinfo(ib_clean,:)];
        end

        %% Reaction time data

        if concatenate_RT

            load(session_paths_files_RT{isess})

            [~, ia_phase, ib_RT] = intersect(phase.trialinfo(:,14), RT_sess.trialinfo(:,14)); %checks the trial numbers
            nTrials = length(ia_phase);

            RT_ampl = nan(length(chan_orig), nTrials);

            % Extract RT values
            RT_vals = RT_sess.RT(:);
            RT_vals = (RT_vals - mean(RT_vals,'omitnan')) ./ std(RT_vals,0,'omitnan'); %zscore

            % Align trials
            RT_aligned = nan(1, nTrials);
            RT_aligned(ia_phase) = RT_vals(ib_RT);

            % Store
            active_ch = find(channels_sessions(:, isess) == 1);

            for ichan = 1:numel(active_ch)
                ct = CriticalTime_mat(active_ch(ichan)) - cut_here;
                if isnan(ct), continue; end

                RT_ampl(active_ch(ichan), :) = RT_aligned;
            end

            RT  = [RT; RT_ampl'];
            RT_trialinfo = [RT_trialinfo; RT_sess.trialinfo(ib_RT,:)];

        end

    end

end

ph_comb =[];
ph_comb.phase_all = ph_all;
ph_comb.amp_all = amp_all;
ph_comb.trialinfo = trlinfo_all;
ph_comb.session_names = session_name_all;
ph_comb.dimord    = 'trlx freq x chan';

if concatenate_lfp_erp_amp
    ph_comb.LFP_ERP_ampl_all  = LFP_ERP_ampl_all;
    ph_comb.LFP_ERP_trialinfo = LFP_ERP_trialinfo;
end

if concatenate_mua_erp_amp
    ph_comb.MUA_ERP_ampl_all  = MUA_ERP_ampl_all;
    ph_comb.MUA_ERP_trialinfo = MUA_ERP_trialinfo;
end

if concatenate_mua_baseline
    ph_comb.MUA_baseline   = MUA_baseline;
    ph_comb.MUA_baseline_trialinfo = MUA_baseline_trialinfo;
end

if concatenate_pupil_baseline
    ph_comb.pup_baseline   = pup_baseline;
    ph_comb.pup_baseline_trialinfo = pup_baseline_trialinfo;
end

if concatenate_RT
    ph_comb.RT   = RT;
    ph_comb.RT_trialinfo = RT_trialinfo;
end

if ~isfolder(info_folder)
    mkdir(fullfile(info_folder))
end

cd(info_folder)
save('ph_all_sess','ph_comb')
clear ph_all_sess ph_comb



