function [hirf, parms] = fristonHIRF(t,parms)
% Illustrate Friston-Worsley HIRF. 
%
%    [hirf, parms] = fristonHIRF(t,parms)
%
% Inputs:
%   t:  Temporal samples
%   parms:  Parameters for the Friston function.  Default is from
%           Friston-Worsley
%
%Brief description:
%
%  Illustrate Friston-Worsley HIRF.  I think this is defined in a
%  Worsely chapter in that Toga book.  Or at least it is defined
%  online, or hopefully in one of the readings for Psych 204.  Figure
%  it out.  Write down the answer here.
%   
% Example:
%    t = 0:0.1:15
%    [hirf, params] = fristonHIRF(t);
%    plot(t,hirf); 
%    xlabel('Time (sec)'); ylabel('Relative amp'); grid on;
%
% See also
%   boyntonHIRF
%

%%
if nargin < 2
    % Use the default parameters
    a(1) = 6; a(2) = 12;
    b(1:2) = 0.9;
    c = 0.35;
else
    a = parms.a; b = parms.b; c = parms.c;
end

for ii = 1:2
    d(ii) = a(ii)*b(ii);
end

% So we can return the parameters
parms.a = a;
parms.b = b;
parms.c = c;

% This is their parameterized function
hirf = (t/d(1)).^a(1).*exp(-(t - d(1))/b(1)) - c*(t/d(2)).^a(2).*exp(-(t-d(2))/b(2));

end

