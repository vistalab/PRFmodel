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

% Obtain a pm to have defaults
pm = prfModel;

% TODO: add checks to values to validate that the value type corresponds to the
% variable that we want to expand
tmp = synthDT;
switch variableName
    case {'HRF'}
        for ii=1:length(values)
            % Add rows with the combinations of parameters we want to check
            % Then check that the whole thing should be complete
            complete = pmParamsCompletenessCheck(values(ii), ...
                table2struct(pm.defaultsTable.HRF));
            % Check if only some of the fields in params where set
            complete.params = pmParamsCompletenessCheck(complete.params, ...
                pm.defaultsTable.HRF.params);
            % Convert it to table
            fieldValuesTable = struct2table(complete,'AsArray',true);
            % Take the copy of the existing table, and make everything equal
            % except the one thing we are changing
            tmp.(variableName)          = repmat(fieldValuesTable, [height(tmp),1]);
            
            % Now concatenate the original and the recently created one.
            synthDT          = [synthDT; tmp];
        end
    case {'Noise'}
        for ii=1:length(values)
            % Add rows with the combinations of parameters we want to check
            % Then check that the whole thing should be complete
            % If we are providing a specific type of noise voxel, obtain
            % such defaults
            if isfield(values,'voxel')
                completeNoise = pmParamsCompletenessCheck(values(ii), ...
                                            table2struct(pm.Noise.defaultsGet('voxel',values(ii).voxel)));
            else
                completeNoise = pmParamsCompletenessCheck(values(ii), ...
                                            table2struct(pm.defaultsTable.Noise));
            end
            complete = pmParamsCompletenessCheck(values(ii), ...
                table2struct(pm.defaultsTable.Noise));
            
            % Convert it to table
            fieldValuesTable = struct2table(complete,'AsArray',true);
            % Take the copy of the existing table, and make everything equal
            % except the one thing we are changing
            tmp.(variableName)          = repmat(fieldValuesTable, [height(tmp),1]);
            
            % Now concatenate the original and the recently created one.
            synthDT          = [synthDT; tmp];
        end
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

