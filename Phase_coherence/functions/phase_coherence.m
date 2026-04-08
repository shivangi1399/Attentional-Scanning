function [coh, phase_spec, coh_complex] = phase_coherence(phase, erp_amp)
% phase_coherence computes phase coherence between phase and ERP amplitude
%
% INPUTS:
%   phase       - trials x frequencies matrix (in radians)
%   erp_amp     - trials x 1 vector of ERP amplitude
%
% OUTPUTS:
%   coh         - 1 x frequencies, coherence magnitude
%   phase_spec  - 1 x frequencies, coherence phase
%   coh_complex - 1 x frequencies, complex coherence (magnitude + phase in one)

nFreq = size(phase, 2);

coh_complex = nan(1, nFreq);

for f = 1:nFreq
    % Scale unit vectors by ERP amplitude
    vec = exp(1i * phase(:,f)) .* erp_amp(:);

    % Complex average — preserve full complex value for downstream averaging
    coh_complex(f) = mean(vec);
end

% Decompose only when needed
coh       = abs(coh_complex);
phase_spec = angle(coh_complex);
end
