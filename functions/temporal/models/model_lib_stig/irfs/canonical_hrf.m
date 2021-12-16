function hrf = canonical_hrf(TR, P)
% Derives the canonical difference of gammas HRF implemented in SPM
% 
% INPUTS
%   TR: fMRI repetition time (s)
%   P: parameters of the response function in units of seconds (optional)
%	  P(1) = delay of response relative to onset (default = 6)
%	  P(2) = delay of undershoot relative to onset (default = 16)
%	  P(3) = length of kernel (default = 32)
% 
% Example of generating the custom HRF used in TemporalChannels code:
%   hrf = spm_hrf(1, [5 14 28]);
% 
% Adapted from vistasoft (https://github.com/vistalab/vistasoft) from SPM
% AS 2/2017

% sampling rate parameter
sr = 16; dt = TR / sr; fratio = 6;

% define default HRF parameters
p = [6 16 32]; if nargin > 1; p(1:length(P)) = P; end

% derive positive gamma function
x1 = 0:p(3) / dt; h1 = p(1); f1 = zeros(size(x1)); Q1  = find(x1 > 0);
f1(Q1) = exp((h1 - 1) .* log(x1(Q1)) + h1 .* log(dt) - dt .* x1(Q1) - gammaln(h1));

% derive negative gamma function
x2 = 0:p(3) / dt; h2 = p(2); f2 = zeros(size(x2)); Q2 = find(x2 > 0);
f2(Q2) = exp((h2 - 1) .* log(x2(Q2)) + h2 .* log(dt) - dt .* x2(Q2) - gammaln(h2));

% find difference of gamma functions and normalize
hrf = f1 - f2 / fratio;
hrf = hrf((0:(p(3) / TR)) * sr + 1);
hrf = hrf' / sum(hrf);

end
