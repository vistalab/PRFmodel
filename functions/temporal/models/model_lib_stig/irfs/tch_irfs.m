function irf = tch_irfs(channel, tau, n1, n2, kappa, fs)
% Derive the sustained and transient impulse response functions based on
% the formulation outlined by Watson (1986).
% 
% INPUT
%   1) channel: 'S' (sustained) or 'T' (transient)
%   2) tau: time constant of excitatory mechanism (ms)
%   3) n1: number of stages in excitatory mechanism
%   4) n2: number of stages in inhibitory mechanism
%   5) kappa: ratio of time constants for primary/secondary filters
%   6) fs: sampling rate of IRF (Hz)
% 
% OUTPUT
%   irf: either the sustained or transient IRF (sampled at fs Hz)
% 
% AS 2/2017

channel = upper(channel(1));
if sum(strcmp(channel, {'S' 'T'})) ~= 1
    error('Channel input must be "S" or "T"');
end
if nargin < 2; tau = 4.93; end
if nargin < 3; n1 = 9; end
if nargin < 4; n2 = 10; end
if nargin < 5; kappa = 1.33; end
if nargin < 6; fs = 1000; end
tau2 = kappa * tau;

% generate filters
% f1 = []; f2 = [];
f1 = zeros(1,length(0:tau * 1000));
f2=f1;
for t = 0:tau * 1000 % time in ms
    % excitatory filter
    f1(t + 1) = ((tau * gamma(n1)) ^ -1) * ...
        ((t / tau) ^ (n1 - 1)) * exp(-t / tau);
    % inhibitory filter
    f2(t + 1) = ((tau2 * gamma(n2)) ^ -1) * ...
        ((t / tau2) ^ (n2 - 1)) * exp(-t / tau2);
end

% derive sustained and transient IRFs
irfS = f1;
irfT = f1 - f2;

irfS = resample(irfS, 1, 1000 / fs)';
irfT = resample(irfT, 1, 1000 / fs)';

% % normalize sum of IRFs

irfT = irfT * (max(irfS) / max(irfT));

switch channel
    case 'S'
        irf = irfS;
    case 'T'
        irf = irfT;
end


% % normalize max of IRFs
% irfT = irfT * (max(irfS) / max(irfT));
% 
% % output sustained or transient IRF
% switch channel
%     case 'S'
%         irf = resample(irfS, 1, 1000 / fs)';
%     case 'T'
%         irf = resample(irfT, 1, 1000 / fs)';
% end


% clip values from end that are close to zero for efficiency
zclip = length(irf); zthresh = max(abs(irf)) / 1000;
while abs(irf(zclip)) < zthresh
    zclip = zclip - 1;
end
irf = irf(1:zclip);

end

