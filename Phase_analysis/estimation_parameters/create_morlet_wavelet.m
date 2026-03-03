function wavelet = create_morlet_wavelet(frequency, sampling_rate, num_points)
    % Create a complex Morlet wavelet for a specific frequency.
    %
    % Args:
    %     frequency (scalar): Central frequency of the wavelet in Hz.
    %     sampling_rate (scalar): Sampling rate in Hz.
    %     num_points (scalar): Number of points in the wavelet.
    %
    % Returns:
    %     wavelet (vector): Complex Morlet wavelet.
    
    % Time vector centered at zero
    t = (-num_points/2 : num_points/2-1) / sampling_rate;
    
    % Define standard deviation for Gaussian taper
    sigma_t = 1 / (2 * pi * frequency);
    
    % Create Gaussian window and sinusoidal carrier
    gaussian_window = exp(-t.^2 / (2 * sigma_t^2));
    sinusoidal_carrier = exp(2i * pi * frequency * t);
    
    % Create the wavelet
    wavelet = gaussian_window .* sinusoidal_carrier;
end