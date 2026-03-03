function plot_sigfreq(freq, corr_vals, limit)
sig_freq = (corr_vals >= limit);
sigonset = find(conv(sig_freq, [1 -1]) == 1);
sigoffset = find(conv(sig_freq, [1 -1]) == -1) - 1;
durplot = [freq(sigonset); freq(sigoffset)];

plot(freq, corr_vals, 'k'); hold on;
for i = 1:size(durplot,2)
    v = [durplot(1,i) 0; durplot(2,i) 0; durplot(2,i) 0.2; durplot(1,i) 0.2];
    patch('Faces', [1 2 3 4], 'Vertices', v, ...
        'FaceColor', 'blue', 'FaceAlpha', 0.2, ...
        'EdgeColor', 'blue', 'EdgeAlpha', 0.2);
end
yline(limit, 'r--', 'LineWidth', 1.2);
ylim([0 0.2]); xlabel('Frequency (Hz)'); ylabel('Correlation');
end
