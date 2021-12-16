function rate = dqdt(t, q, v, in_flow)
% Calculates rate of change in normalized deoxyhemoglobin as a function of 
% deoxyhemoglobin, normalized volume, and oxygen extraction.
% 
% Adapted from code by JG (gru.stanford.edu/svn/matlab/balloonmodel.m)
% AS 9/2017

global E0  tauMTT
if (iscell(v))
    in_flow = v{2};
    v = v{1};
end

% relate oxygen extraction to rate of flow
qtE = 1 - (1 - E0) .^ (1 ./ in_flow);
vt = flow_out(v, t, in_flow);
qt = (in_flow * (qtE / E0) - vt * (q / v));
rate = (1 / tauMTT) * qt;

end
