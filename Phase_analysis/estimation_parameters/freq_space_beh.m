clear all
close all
clc

%% Frequency space behavior: comparison of different log-spaced methods

% Standard logspace
logsp1 = logspace(log10(2), log10(80), 40);

% Power-law warped linspace in log space
N = 30;
log_start = log10(2);
log_end = log10(80);
s2 = linspace(0, 1, N).^1.5;  % Bias toward low frequencies
log_vals2 = log_start + (log_end - log_start) * s2;
logsp2 = 10.^log_vals2;

% Exponential stretch in log space
alpha = 2;  % Controls curvature; alpha > 1 favors low freqs
s3 = linspace(0, 1, N);
log_vals3 = log10(2) + (log10(80) - log10(2)) * (exp(alpha * s3) - 1) / (exp(alpha) - 1);
logsp3 = 10.^log_vals3;

% Hyperbolic Tangent Mapping
s4 = linspace(-1, 1, N);
warp = (tanh(2 * s4) + 1) / 2;  % Range scaled to [0, 1]
log_vals4 = log10(2) + warp * (log10(80) - log10(2));
logsp_tanh = 10.^log_vals4;

% Plotting
figure;
hold on;
grid on;

plot(logsp1, '-o', 'DisplayName', 'logspace (40 pts)', 'LineWidth', 1.5);
plot(logsp2, '-s', 'DisplayName', 'power-law (30 pts)', 'LineWidth', 1.5);
plot(logsp3, '-^', 'DisplayName', 'exponential map (30 pts)', 'LineWidth', 1.5);
plot(logsp_tanh, '-d', 'DisplayName', 'tanh mapping (30 pts)', 'LineWidth', 1.5);

xlabel('Index');
ylabel('Frequency (Hz)');
title('Comparison of log-like frequency spacings');
legend('Location', 'northwest');

% stick with the standard logspace one