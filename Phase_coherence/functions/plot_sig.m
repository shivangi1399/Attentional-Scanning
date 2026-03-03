function plot_sig(freq, values, limit, xlabel_str, ylabel_str)
    % plot_sigcoh: plot values vs frequency with shaded significance regions
    %
    % Usage:
    %   plot_sigcoh(freq, values, limit)
    %   plot_sigcoh(freq, values, limit, 'X Label', 'Y Label')

    if nargin < 4 || isempty(xlabel_str)
        xlabel_str = 'Frequency (Hz)';
    end
    if nargin < 5 || isempty(ylabel_str)
        ylabel_str = 'Correlation';
    end

    % Identify significant frequency ranges
    sig_freq = (values >= limit);
    sigonset = find(conv(sig_freq, [1 -1]) == 1);
    sigoffset = find(conv(sig_freq, [1 -1]) == -1) - 1;
    durplot = [freq(sigonset); freq(sigoffset)];

    % Plot values
    plot(freq, values, 'k'); hold on;

    % Get current y-axis limits
    ylims = ylim;

    % Add shaded patches spanning full y-axis range
    for i = 1:size(durplot,2)
        v = [durplot(1,i) ylims(1); ...
             durplot(2,i) ylims(1); ...
             durplot(2,i) ylims(2); ...
             durplot(1,i) ylims(2)];
        patch('Faces', [1 2 3 4], 'Vertices', v, ...
              'FaceColor', 'blue', 'FaceAlpha', 0.2, ...
              'EdgeColor', 'none');
    end

    % Red significance threshold line
    yline(limit, 'r--', 'LineWidth', 1.2);

    % Reset ylim in case patch modified it
    ylim(ylims);

    xlabel(xlabel_str); 
    ylabel(ylabel_str);
end
