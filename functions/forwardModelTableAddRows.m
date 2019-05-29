function synthDT = forwardModelTableAddRows(synthDT, variableName, values)
% Adds rows combining the new values with the existing ones. 
% 
%  Inputs: variable name and value(s)
% 
%  Outputs: Matlab table
% 
% 

vars = strsplit(variableName,'.');

% TODO: add checks to values to validate that the value type corresponds to the
% variable that we want to expand
tmp = synthDT;
for ii=1:length(values)
    if length(vars) == 2
        tmp.(vars{1}).(vars{2}) = repmat(values(ii), [height(tmp),1]);
    else
        tmp.(vars{1})           = repmat(values(ii), [height(tmp),1]);
    end
    synthDT = [synthDT; tmp];
end
       
       
       
       
% Calculate the models       
synthDT = forwardModelCalculate(synthDT);

end

