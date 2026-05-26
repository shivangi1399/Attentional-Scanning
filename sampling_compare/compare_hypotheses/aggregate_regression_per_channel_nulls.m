% =====================================================================
% One-time aggregator: consolidates per-channel × per-perm regression
% null shards into a single per-channel mat so the per-channel paired
% H2−H1 regression test in compare_hypotheses_per_chan.m can do a
% single load per channel instead of 1000.
%
% INPUT shards (one file per perm per channel, written by
% regress_perm_R.m / regress_perm_R_pos.m):
%   results_<animal>/multi_lin_reg/<hyp>/cp10_till_100/
%       perm_R[_pos]/<dv_key>/<ch>/perm_NNNN.mat
%   each containing struct `results` with:
%     .null_R_phase  (1×nFreq)  R² for phase predictors
%     H1: .null_b_sin, .null_b_cos      (nFreq×1)   full-model β
%     H2: .null_b_sin_pos, .null_b_cos_pos (nPos×nFreq) per-position β
%
% OUTPUT (one per channel):
%   <ch>/per_channel_null.mat with:
%     null_R2_phase     (nPerm × nFreq)  partial R² null (same as
%                       .null_R_phase from each shard, just stacked)
%     null_R_phase_mag  (nPerm × nFreq)  R_phase = |complex β| null,
%                       computed per the hypothesis recipe:
%                         H1: sqrt(b_sin² + b_cos²)
%                         H2: mean_pos sqrt(b_sin_pos² + b_cos_pos²)
%
% After this runs, compare_hypotheses_per_chan.m can compute the
% paired test for reg_R2 and reg_Rphase per channel with one file
% per channel per hypothesis.
% =====================================================================

clearvars; close all; clc;

%% Settings
animals = {'hermes','klecks'};
hypotheses = struct( ...
    'label',    {'H1','H2'}, ...
    'folder',   {'complex','abs_per_pos'}, ...
    'perm_sub', {'perm_R','perm_R_pos'});
dvs    = {'MUA_ERP_ampl_all','LFP_ERP_ampl_all','RT','hit_miss'};
nCh    = 64;
base   = '/mnt/hpc/projects/MWSampling/4Shivangi';
overwrite = false;   % set true to re-aggregate even when output exists

total_files_written = 0;

for a = 1:numel(animals)
    animalName = animals{a};
    fprintf('\n=== %s ===\n', animalName);
    for h = 1:numel(hypotheses)
        H = hypotheses(h);
        for d = 1:numel(dvs)
            dv = dvs{d};
            dv_dir = fullfile(base, ['results_' animalName], 'multi_lin_reg', ...
                H.folder, 'cp10_till_100', H.perm_sub, dv);
            if ~exist(dv_dir, 'dir')
                fprintf('  SKIP %s/%s/%s — %s not found\n', animalName, H.label, dv, dv_dir);
                continue
            end
            fprintf('  %s / %s / %s\n', animalName, H.label, dv);

            for ch = 1:nCh
                ch_dir = fullfile(dv_dir, num2str(ch));
                if ~exist(ch_dir, 'dir'), continue; end
                out_file = fullfile(ch_dir, 'per_channel_null.mat');
                if ~overwrite && isfile(out_file)
                    continue
                end

                perm_files = dir(fullfile(ch_dir, 'perm_*.mat'));
                if isempty(perm_files), continue; end
                nPerm = numel(perm_files);

                % Probe first shard to learn nFreq (and confirm fields)
                tmp = load(fullfile(ch_dir, perm_files(1).name), 'results');
                if ~isfield(tmp, 'results') || ...
                   ~isfield(tmp.results, 'null_R_phase')
                    fprintf('    ch%d: probe missing null_R_phase, skip\n', ch);
                    continue
                end
                nFreq = numel(tmp.results.null_R_phase);

                null_R2_phase    = nan(nPerm, nFreq);
                null_R_phase_mag = nan(nPerm, nFreq);

                for p = 1:nPerm
                    pf = fullfile(ch_dir, perm_files(p).name);
                    try
                        S = load(pf, 'results');
                        r = S.results;
                    catch err
                        warning('    %s ch%d %s: load failed (%s)', ...
                            dv, ch, perm_files(p).name, err.message);
                        continue
                    end

                    % R² for phase — present in both hypotheses
                    if isfield(r, 'null_R_phase')
                        null_R2_phase(p, :) = r.null_R_phase(:)';
                    end

                    % R_phase magnitude — recipe depends on hypothesis
                    switch H.label
                        case 'H1'
                            if isfield(r,'null_b_sin') && isfield(r,'null_b_cos')
                                bs = r.null_b_sin(:)';   % 1 × nFreq
                                bc = r.null_b_cos(:)';
                                null_R_phase_mag(p, :) = sqrt(bs.^2 + bc.^2);
                            end
                        case 'H2'
                            if isfield(r,'null_b_sin_pos') && isfield(r,'null_b_cos_pos')
                                bs = r.null_b_sin_pos;    % nPos × nFreq
                                bc = r.null_b_cos_pos;
                                mag_per_pos = sqrt(bs.^2 + bc.^2);
                                null_R_phase_mag(p, :) = ...
                                    mean(mag_per_pos, 1, 'omitnan');
                            end
                    end
                end

                save(out_file, 'null_R2_phase', 'null_R_phase_mag', '-v7.3');
                total_files_written = total_files_written + 1;
                if mod(ch, 8) == 0 || ch == 1 || ch == nCh
                    fprintf('    ch%d done (%d perms)\n', ch, nPerm);
                end
            end
        end
    end
end

fprintf('\nDone. %d per-channel aggregates written.\n', total_files_written);
