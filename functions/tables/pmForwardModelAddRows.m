function synthDT = pmForwardModelAddRows(synthDT, variableName, values)
% Adds rows combining the new values with the existing ones.
%
%  Inputs: variable name and value(s)
%
%  Outputs: Matlab table
%
%

% We already have the defaults, which are the first value
if ~isa(values, 'table')
    values = values(2:end);
end

% TODO: add checks to values to validate that the value type corresponds to the
% variable that we want to expand
tmp = synthDT;

switch variableName
    case 'HRF'
        tmp.HRF          = repmat(values, [height(tmp),1]);
        % Solving the concatenation problem for Type
        tmp.HRF.Type     = categorical(cellstr(tmp.HRF.Type));
        synthDT.HRF.Type = categorical(synthDT.HRF.Type);
        % Solving the concatenation problem for Duration
        tmp.HRF.Duration     = double(tmp.HRF.Duration);
        synthDT.HRF.Duration = double(synthDT.HRF.Duration);
        
        % Now concatenate the original and the recently created one.
        synthDT          = [synthDT; tmp];
    case 'Noise'
        warning('Noise not implemented yet.')
    otherwise
        vars = strsplit(variableName,'.');
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
end

