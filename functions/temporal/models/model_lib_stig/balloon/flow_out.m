function flow = flow_out(v, t, in_flow)
% Calculates the output flow as a function of volume.
% 
% Adapted from code by JG (gru.stanford.edu/svn/matlab/balloonmodel.m)
% AS 9/2017

global which_tau tauN tauP tau alpha;

vt = dvdt(t, v, in_flow);
if (which_tau == 1)
    flow = v ^ (1 / alpha) + tauP * vt;
    tau = tauP;
else
    flow = v ^ (1 / alpha) + tauN * vt;
    tau = tauN;
end

end