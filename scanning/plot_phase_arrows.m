% =====================================================================
% Per-position phase-vector arrows on an 8×8 channel grid
%
% For each animal, at a set of fixed frequencies, plots an 8×8 grid
% (matching the electrode layout) where each subplot shows the
% per-position c_p vectors as colored arrows from the origin.
%
% Arrow = complex coherence c_p = mean(y · exp(i·φ)) for that position.
%   angle  → preferred phase at that position
%   length → coherence magnitude (coupling strength)
%
% Positions are color-coded peripheral (blue) → foveal (red).
% Channels significant in the H2−H1 test (per-channel, max-stat
% corrected) get a bold black border; non-significant channels are
% drawn faded so you can compare the spatial structure.
%
% Outputs one multi-page PDF per animal at:
%   Plots/scanning/phase_progression/cp10_till_100/lfp/
%       <animal>_phase_arrows.pdf
% =====================================================================
clearvars; close all; clc

%% Settings
animals       = {'klecks','hermes'};
dv_list       = {'lfp','mua'};
dv_key_map    = struct('lfp','LFP_ERP_ampl_all','mua','MUA_ERP_ampl_all');
target_freqs  = [5 10 15 20];      % Hz — one 8×8 grid page per frequency
alpha         = 0.05;
nCh           = 64;
grid_rows     = 8;
grid_cols     = 8;
base          = '/mnt/hpc/projects/MWSampling/4Shivangi';

for dv_idx = 1:numel(dv_list)
dv = dv_list{dv_idx};
for a = 1:numel(animals)
    animalName = animals{a};
    fprintf('\n=== %s / %s ===\n', animalName, upper(dv));

    %% Load phase-progression data (has per-position c_p via pref_phase + coh_mag)
    pp_file = fullfile(base, ['results_' animalName], 'scanning', ...
        'phase_progression','cp10_till_100', dv, 'phase_progression.mat');
    if ~isfile(pp_file)
        warning('No phase progression data for %s/%s.', animalName, dv);
        continue
    end
    S = load(pp_file, 'pref_phase','coh_mag','R_obs','R_null','freq','positions');
    % pref_phase: nCh × nFreq × nPos (angle of c_p)
    % coh_mag:    nCh × nFreq × nPos (|c_p|)
    freq      = S.freq;
    positions = S.positions;
    nPos      = numel(positions);
    nFreq     = numel(freq);

    %% Per-channel H2−H1 significance (regression Reg R² paired test)
    % This marks channels where positions significantly disagree on
    % preferred phase in the regression framework (max-stat corrected
    % within each channel). Uses the per_channel_null.mat files from
    % the aggregate_regression_per_channel_nulls pipeline.
    sig_mask = false(nCh, nFreq);
    dv_key = dv_key_map.(dv);

    H1_obs_file = fullfile(base, ['results_' animalName], ...
        'multi_lin_reg','complex','cp10_till_100', ...
        'multi_regression_channelwise_R2.mat');
    H2_obs_file = fullfile(base, ['results_' animalName], ...
        'multi_lin_reg','abs_per_pos','cp10_till_100', ...
        'multi_regression_channelwise_R2_abs_per_pos.mat');

    if isfile(H1_obs_file) && isfile(H2_obs_file)
        S1 = load(H1_obs_file, 'reg_results');
        S2 = load(H2_obs_file, 'reg_results');
        obs_H1 = S1.reg_results.(dv_key).R2_phase;   % nCh × nFreq
        obs_H2 = S2.reg_results.(dv_key).R2_phase;
        obs_diff = obs_H2 - obs_H1;

        h1_root = fullfile(base, ['results_' animalName], ...
            'multi_lin_reg','complex','cp10_till_100','perm_R', dv_key);
        h2_root = fullfile(base, ['results_' animalName], ...
            'multi_lin_reg','abs_per_pos','cp10_till_100','perm_R_pos', dv_key);

        for ch = 1:nCh
            h1f = fullfile(h1_root, num2str(ch), 'per_channel_null.mat');
            h2f = fullfile(h2_root, num2str(ch), 'per_channel_null.mat');
            if ~isfile(h1f) || ~isfile(h2f), continue; end
            P1 = load(h1f, 'null_R2_phase');
            P2 = load(h2f, 'null_R2_phase');
            if ~isfield(P1,'null_R2_phase') || ~isfield(P2,'null_R2_phase'), continue; end
            p1 = P1.null_R2_phase; p2 = P2.null_R2_phase;
            n_pair = min(size(p1,1), size(p2,1));
            if n_pair < 1 || size(p1,2) ~= nFreq || size(p2,2) ~= nFreq, continue; end
            diff_perm = p2(1:n_pair,:) - p1(1:n_pair,:);
            thr = quantile(max(diff_perm, [], 2), 1 - alpha);
            if isfinite(thr)
                sig_mask(ch,:) = obs_diff(ch,:) >= thr;
            end
        end
        fprintf('  H2-H1 sig: %d channels have >= 1 sig freq\n', sum(any(sig_mask,2)));
    else
        warning('H2-H1 regression results not found for %s — no significance overlay.', animalName);
    end

    %% Position colormap: peripheral (blue) → foveal (red)
    pos_cmap = [linspace(0.15, 0.85, nPos)', ...
                linspace(0.30, 0.15, nPos)', ...
                linspace(0.80, 0.20, nPos)'];

    %% Output setup
    plot_dir = fullfile(base, 'Plots','scanning','phase_progression', ...
        'cp10_till_100', dv);
    if ~exist(plot_dir,'dir'), mkdir(plot_dir); end
    out_pdf = fullfile(plot_dir, sprintf('%s_phase_arrows.pdf', animalName));
    out_ps  = strrep(out_pdf, '.pdf', '.ps');
    if isfile(out_pdf), delete(out_pdf); end
    if isfile(out_ps),  delete(out_ps);  end

    %% One page per target frequency
    for tf = 1:numel(target_freqs)
        [~, fi] = min(abs(freq - target_freqs(tf)));
        f_hz = freq(fi);

        fig = figure('Visible','off','Units','centimeters', ...
            'Position',[1 1 32 32]);
        set(fig,'PaperUnits','centimeters','PaperSize',[32 32], ...
            'PaperPosition',[0 0 32 32]);

        % Global max magnitude for consistent arrow scaling across channels
        mag_at_f = squeeze(S.coh_mag(:, fi, :));   % nCh × nPos
        rmax = max(mag_at_f(:), [], 'omitnan');
        if ~isfinite(rmax) || rmax == 0, rmax = 1; end

        for ch = 1:nCh
            % Column-wise bottom-to-top (matches loc_chan_pairing.m)
            col = ceil(ch / grid_rows);
            row = grid_rows - mod(ch - 1, grid_rows);
            % Subplot index: row 1 = top of figure
            subplot_idx = (row-1)*grid_cols + col;
            ax = subplot(grid_rows, grid_cols, subplot_idx);
            hold(ax, 'on');

            phi = squeeze(S.pref_phase(ch, fi, :));   % nPos × 1
            mag = squeeze(S.coh_mag(ch, fi, :));       % nPos × 1

            if all(isnan(phi))
                text(ax, 0, 0, 'no data', 'FontSize', 5, ...
                    'HorizontalAlignment','center');
                set(ax, 'XLim',[-1 1],'YLim',[-1 1]);
                axis(ax, 'square'); axis(ax, 'off');
                title(ax, sprintf('ch%d', ch), 'FontSize', 6);
                continue
            end

            % Draw unit circle for reference
            theta_c = linspace(0, 2*pi, 100);
            plot(ax, rmax*cos(theta_c), rmax*sin(theta_c), ...
                'Color',[0.85 0.85 0.85], 'LineWidth', 0.3);
            plot(ax, [-rmax rmax], [0 0], 'Color',[0.9 0.9 0.9], 'LineWidth', 0.2);
            plot(ax, [0 0], [-rmax rmax], 'Color',[0.9 0.9 0.9], 'LineWidth', 0.2);

            % Determine visual style based on significance at this freq
            is_sig = sig_mask(ch, fi);
            if is_sig
                arrow_alpha = 1.0;
                lw = 1.5;
            else
                arrow_alpha = 0.3;
                lw = 0.8;
            end

            % Draw arrows from origin to c_p for each position
            for p = 1:nPos
                if isnan(phi(p)) || isnan(mag(p)), continue; end
                x_end = mag(p) * cos(phi(p));
                y_end = mag(p) * sin(phi(p));
                col_p = pos_cmap(p,:);
                if ~is_sig
                    col_p = 0.5 + 0.5*col_p;  % wash out
                end
                quiver(ax, 0, 0, x_end, y_end, 0, ...
                    'Color', [col_p, arrow_alpha], ...
                    'LineWidth', lw, ...
                    'MaxHeadSize', 0.6);
            end

            lim = rmax * 1.15;
            set(ax, 'XLim',[-lim lim], 'YLim',[-lim lim]);
            axis(ax, 'square');
            set(ax, 'XTick',[],'YTick',[]);

            n_sig_f = sum(sig_mask(ch,:));
            if is_sig
                title(ax, sprintf('ch%d *', ch), ...
                    'FontSize', 6, 'FontWeight', 'bold', 'Color', [0.1 0.5 0.2]);
            else
                title(ax, sprintf('ch%d', ch), 'FontSize', 6, 'Color', [0.5 0.5 0.5]);
            end

            % Bold border for significant channels
            if is_sig
                set(ax, 'Box','on', 'LineWidth', 1.5, ...
                    'XColor',[0.1 0.5 0.2], 'YColor',[0.1 0.5 0.2]);
            else
                set(ax, 'Box','on', 'LineWidth', 0.3, ...
                    'XColor',[0.8 0.8 0.8], 'YColor',[0.8 0.8 0.8]);
            end
        end

        % Position colorbar legend
        cb_ax = axes('Position',[0.92 0.05 0.015 0.25],'Visible','off');
        hold(cb_ax,'on');
        for p = 1:nPos
            y_p = (p-1)/(nPos-1);
            plot(cb_ax, 0.5, y_p, 's', 'MarkerSize', 8, ...
                'MarkerFaceColor', pos_cmap(p,:), 'MarkerEdgeColor','none');
            text(cb_ax, 1.2, y_p, sprintf('p%g', positions(p)), ...
                'FontSize', 5, 'VerticalAlignment','middle');
        end
        set(cb_ax, 'XLim',[0 3],'YLim',[-0.1 1.1],'Visible','off');

        sgtitle(sprintf('%s — %s — %.1f Hz :: per-position phase vectors\n(blue=peripheral, red=foveal; bold border = sig at this freq)', ...
            animalName, upper(dv), f_hz), ...
            'FontSize', 10, 'FontWeight', 'bold');

        print(fig, '-dpsc', '-append', '-painters', out_ps);
        close(fig);
        fprintf('  Page %d: %.1f Hz\n', tf, f_hz);
    end

    %% Convert PS → PDF
    [s, msg] = system(sprintf('ps2pdf "%s" "%s"', out_ps, out_pdf));
    if s == 0 && isfile(out_pdf)
        delete(out_ps);
        fprintf('Saved %s\n', out_pdf);
    else
        warning('ps2pdf failed (%d): %s', s, strtrim(msg));
    end
end
end % dv loop

fprintf('\nDone.\n');
