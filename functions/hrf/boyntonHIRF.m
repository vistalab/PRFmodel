function [hirf, t, parms] = boyntonHIRF(t,parms)
%
%    [hirf, t, parms] = boyntonHIRF(t,parms)
%
% Examples:
%    [hirf,t] = boyntonHIRF(t,parms);
%    [hirf,t,parms] = boyntonHIRF(t);
%
%Author:   Wandell
%Purpose:
%   Compute the Boynton et al. HIRF function.  The return values include
% both the hirf and the values of time and the parameters used in the
% computation.
%

disp('Boynton HIRF')

if ~exist('parms','var')
    parms.delay = 2.05;
    parms.n = 3;
    parms.tau = 1.08;
end

tau = parms.tau;
n = parms.n;
delay = parms.delay;

hirf = (t/tau).^(n-1).*exp(-(t/tau))/(tau*(factorial(n-1)));

% Now, account for the early delay.
dt = t(2) - t(1);
t = [0:dt:(delay - dt),t + delay];
nZeros = length(t) - length(hirf);
earlyZeros = zeros(1,nZeros);
hirf = [earlyZeros, hirf];

if length(hirf) ~= length(t)
    warning('Bad hirf,t');
end

return;
