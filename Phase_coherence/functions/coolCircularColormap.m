function cmap = coolCircularColormap(n)
    % coolCircularColormap Generate a smooth cyclic colormap in cool tones
    %   cmap = coolCircularColormap(n) returns an n×3 RGB array
    
    if nargin < 1, n = 256; end
    theta = linspace(0, 2*pi, n)';   % angle
    
    % Define RGB channels with cosine functions
    R = 0.4 + 0.6*cos(theta);
    G = 0.6 + 0.4*cos(theta - 2*pi/3);
    B = 0.6 + 0.4*cos(theta - 4*pi/3);
    
    % Bias toward cooler colors (reduce red contribution)
    R = R * 0.6;
    
    % Normalize to [0,1] and clip to avoid out-of-range values
    cmap = [R G B];
    cmap = max(min(cmap, 1), 0);
end
