function [coh, phase_spec] = phase_coherence(phase, erp_amp)
% phase_coherence computes phase coherence between phase and ERP amplitude
%
% INPUTS:
%   phase    - trials x frequencies matrix (in radians)
%   erp_amp  - trials x 1 vector of ERP amplitude
%
% OUTPUTS:
%   coh       - 1 x frequencies, coherence magnitude
%   phase_spec- 1 x frequencies, coherence phase

nFreq = size(phase, 2);

coh = nan(1, nFreq);
phase_spec = nan(1, nFreq);

for f = 1:nFreq
    % Scale unit vectors by ERP amplitude
    vec = exp(1i * phase(:,f)) .* erp_amp(:);
    
    % Complex average
    cavg = mean(vec);
    
    % Coherence magnitude and phase
    coh(f) = abs(cavg);
    phase_spec(f) = angle(cavg);
end
end
