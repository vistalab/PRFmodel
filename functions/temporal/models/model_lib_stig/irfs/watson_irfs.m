function irf = watson_irfs(channel, fs)
% Derive the sustained and transient impulse response functions proposed
% by Watson (1986).
% 
% INPUT
%   1) channel: either 'S' (sustained) or 'T' (transient)
%   2) fs: sampling rate of IRF (Hz)
% 
% OUTPUT
%   irf: either the sustained or transient IRF (sampled at fs Hz)
% 
% AS 2/2017

channel = upper(channel(1));

% filter parameters from Watson
k = 1.33;    % ratio of time constants for different response stages
t1 = 4.94;   % time constant of excitatory mechanism
t2 = k * t1; % time constant of inhibitory mechanism
n1 = 9;      % number of stages in exctatory mechanism
n2 = 10;     % number of stages in inhibitory mechanism

% generate filters
time = 0:149; % time in ms
for t = time
    % excitatory filter
    f1(t + 1) = ((t1 * factorial(n1 - 1)) ^ -1) * ((t / t1) ^ (n1 - 1)) * exp(-t / t1);
    % inhibitory filter
    f2(t + 1) = ((t2 * factorial(n2 - 1)) ^ -1) * ((t / t2) ^ (n2 - 1)) * exp(-t / t2);
end

% derive sustained and transient IRFs
irfS = f1;
irfT = f1 - f2;
% normalize max of S and T IRFs
irfT = irfT * (max(irfS) / max(irfT));

% output sustained or transient IRF
if strcmp(channel,'S')
    irf = resample(irfS, 1, 1000 / fs);
elseif strcmp(channel,'T')
    irf = resample(irfT, 1, 1000 / fs);
else
    error('unexpected input aguement');
end

end
