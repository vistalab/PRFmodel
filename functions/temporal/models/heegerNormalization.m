function output = trf_heegerNormalization(params, stim, t)

% INPUTS:
% params have the following fields
%       tau1  - time constant for linear impulse response
%       sigma - semi-saturation constant
%       alpha - weight on feedback (0 = no feedback)
%       n     - output non-linearity exponent

% OUTPUTS:

%% initialize response

R = zeros(1, length(stim)); % Normalized response
G = zeros(size(R));         % Feedback signal
F = zeros(size(R));         % Multiplicative feedback
K = 1;                      % Determines maximum responses

%% make impulse response function

irf     = gammaPDF(t, params.tau1, 2);

%% compute normalization response

for k = 1 : size(stim, 1)
    L       = convCut(stim(k, :), irf, length(stim));
    R(1) = max(L(1),0)^params.n;
    G(1) = params.alpha * R(1);
    F(1) = sqrt(K-G(1))/ params.sigma;
    
    for ii = 2:length(t)
        R(ii) = max(L(ii) * F(ii-1),0)^params.n;
        G(ii) = (1-params.alpha) * G(ii-1) + params.alpha * R(ii);
        G(ii) = min(G(ii), K);
        F(ii) = sqrt(K-G(ii-1)) / params.sigma;
    end
    output(k, :) = R./max(R);
end

end