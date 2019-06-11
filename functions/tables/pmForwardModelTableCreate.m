function synthDT = pmForwardModelTableCreate()
% Creates default format for a PRF model table
%
% Syntax
%  synthDT = forwardModelTableCreate()
%
% Brief description
%   This function defines the default parameters (defaults) required
%   to perform a forward calculation
% 
% Inputs
%  N/A
% 
% Outputs
%  Matlab table
% 
% 

%%
% p = inputParser;


%%

HRF      = cell2table(     {"friston",  20, [6,12],  [0.9,0.9],   0.35}, ...
           'VariableNames',{'Type','duration','params_a','params_b','params_c'});
Stimulus = cell2table(     {"103"    , true   , 20      , 20}, ...
           'VariableNames',{'ExpName','Binary','fovHorz','fovVert'});
RF       = cell2table(     {0 ,    0,      0,         1,      1}, ...
           'VariableNames',{'x0','y0','theta','sigMajor','sigMinor'});
Noise    = cell2table(     {"white", 0.5}, ...
           'VariableNames',{'Type' ,'white_k'});

% Create the main table with one row and default values
pm = prfModel_basic;
synthDT = cell2table(      {"1" , 'basic',  1, HRF , Stimulus , RF , Noise , pm}, ...
           'VariableNames',{'ID','Type','TR','HRF','Stimulus','RF','Noise','pm'});


end

