function C = spm_drift(N, K)
% Implementation of the R low frequency drift calculation
% Belgium neuRoSim
% 
% Inputs:
%  N - The number of sample points over time
%  K - number of cosine basis functions
%
% Outputs
%   C - A matrix with N time points and K columns of the noise basis functions
%       C is scaled pretty arbitrarily by sqrt(2/N)*10
%
% The returned basis functions are NOT signal dependent.  So when used with a
% BOLD signal the basis terms must be both combined and scaled.
%
% GL
%
% See also
%   pmNoise


n = 0:(N - 1);
C = zeros(N,K);
C(:, 1) = 1/sqrt(N);
for k = 2:K
    C(:, k) = sqrt(2/N) * 10 * cos(pi * (2 * n + 1) * (k - 1)/(2 * N));
end

C = C*diag([1 0.5 1]);

end

