function rate = dvdt(t, v, in_flow)
% Calculates rate of change in volume as a function of volume
% 
% Adapted from code by JG (gru.stanford.edu/svn/matlab/balloonmodel.m)
% AS 9/2017

global which_tau tauN tauP tauMTT alpha;

% get in_flow if it is passed by runge_kutta
if (iscell(in_flow))
    in_flow = in_flow{1};
end

if (which_tau == 1)
    rate = ((1 / tauMTT) * (in_flow - v ^ (1 / alpha))) / (1 + (tauP / tauMTT));
else
    rate = ((1 / tauMTT) * (in_flow - v ^ (1 / alpha))) / (1 + (tauN / tauMTT));
end

if rate < 0
    which_tau = 2;
else
    which_tau = 1;
end

end
