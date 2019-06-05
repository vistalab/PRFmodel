%% MAIN SCRIPT: for a single synthetic BOLD signal generation and testing
% 
% This script is a wrapper that generates A SINGLE synthetic BOLD signal and 
% estimates calculated with different pRF implementations. It is intended to be
% used as a learning and testing tool. Once it is mastered and tested with a
% single BOLD series that this is what we want to do, we can create multiple
% variations in s_main_table.m, analyze them with different implementations and
% test them.
% 
% The unit element of this software is a prfModel class. 
                pm = prfModel_basic;
                pm.plot('both')  % Shows a BOLD signal with and without noise
% If any element is changed, the result will be changed:
                pm.TR = 1.82;
                pm.Noise{3}.params.amplitude = 0.9;  % Respiratory noise
                pm.plot('both')
%                 
% Any changes made to any of the "editable" parameters, generates a new signal
% The editable parameters can be seen with the following function:
                pm.defaultsTable
