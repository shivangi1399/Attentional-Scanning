% File: calculate_itc.m
function itc = calculate_itc(signal, frequencies, sampling_rate)
    % Calculate the Inter-Trial Coherence (ITC) for given signal trials.
    % 
    % Args:
    %     signal (matrix): Trials x Time matrix of signals.
    %     frequencies (vector): Vector of frequencies for wavelet analysis.
    %     sampling_rate (scalar): Sampling rate of the signal in Hz.
    % 
    % Returns:
    %     itc (matrix): Frequencies x Time ITC matrix.

    [num_trials, num_timepoints] = size(signal);
    num_frequencies = length(frequencies);
    
    % Initialize the ITC matrix
    itc = zeros(num_frequencies, num_timepoints);
    
    for f_idx = 1:num_frequencies
        % Get the current frequency
        freq = frequencies(f_idx);
        
        % Create the Morlet wavelet
        wavelet = create_morlet_wavelet(freq, sampling_rate, num_timepoints);
        
        % Convolve each trial with the wavelet
        convolution_results = zeros(num_trials, num_timepoints);
        for trial = 1:num_trials
            convolution_results(trial, :) = conv(signal(trial, :), wavelet, 'same');
        end
        
        % Extract the phase information
        phase_data = angle(convolution_results);
        
        % Compute ITC
        complex_vectors = exp(1i * phase_data); % Convert phase to unit-length complex vectors
        itc(f_idx, :) = abs(mean(complex_vectors, 1)); % Compute ITC magnitude across trials
    end
end


