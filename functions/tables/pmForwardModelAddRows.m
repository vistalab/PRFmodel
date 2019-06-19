function synthDT = pmForwardModelAddRows(synthDT, variableName, values)
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
    if length(vars) == 1
        tmp.(vars{1})           = repmat(values(ii), [height(tmp),1]);
    end
    if length(vars) == 2
        tmp.(vars{1}).(vars{2}) = repmat(values(ii), [height(tmp),1]);
    end
    if length(vars) == 3
        tmp.(vars{1}).(vars{2}).(vars{3}) = repmat(values(ii), [height(tmp),1]);
    end
    if length(vars) >= 4
        error('Only implemented until 3 levels of nesting in the table')
    end
    synthDT = [synthDT; tmp];
end

end

