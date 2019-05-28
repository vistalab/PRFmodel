function synthDT = forwardModelTableCreate()
% Creates a table with all parameters (defaults) required to perform a forward
% calculation
% 
%  Inputs: none
% 
%  Outputs: Matlab table
% 
% 


HRF      = cell2table(     {"Friston",  [6,12],  [0.9,0.9],   0.35}, ...
           'VariableNames',{'Type','Friston_a','Friston_b','Friston_c'});
Stimulus = cell2table(     {"103"    , true   , 100      , 100}, ...
           'VariableNames',{'ExpName','Binary','fovHorz','fovVert'});
RF       = cell2table(     {0 ,    0,      0,         1,      1}, ...
           'VariableNames',{'x0','y0','theta','sigMajor','sigMinor'});
Noise    = cell2table(     {"white", 0.5}, ...
           'VariableNames',{'Type' ,'white_k'});

% Create the main table with one row and default values
synthDT = cell2table(      {"1" , HRF , Stimulus , RF , Noise , prfModel()}, ...
           'VariableNames',{'ID','HRF','Stimulus','RF','Noise','pm'});

synthDT = forwardModelCalculate(synthDT);

end

