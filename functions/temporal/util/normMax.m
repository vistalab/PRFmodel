function out_ts = normMax(in_ts)
% DESCRIPTION -------------------------------------------------------------
% function out_ts - normMax(in_ts)
% This file normalizes time course in each electrode to its maximum
% response. This function assumes that each time course is base-line
% reduced, i.e. the average pre-stimulus (first 200ms) response is 0.
%
% INPUTS ------------------------------------------------------------------
% in_ts is in the form of time courses x n electrodes
%
% OUTPUTS -----------------------------------------------------------------
% out_ts: the same shape as in_ts, except that all values within a time
% course are less than 1.

%% USEFUL FUNCTION

normMax = @(x) x./max(x);

%% DERIVED PARAMETERS

tLth  = size(in_ts, 1);
nElec = size(in_ts, 2);

%% COMPUTE NORMALIZING TO THE MAX

if tLth ~= 1 & nElec ~= 1
    for iElec = 1 : nElec
        out_ts(:, iElec) = normMax(in_ts(:, iElec));
    end
else
    out_ts = normMax(in_ts);
end

end