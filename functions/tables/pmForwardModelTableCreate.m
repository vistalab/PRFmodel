function synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS)
% Creates default format for a PRF model table
%
% Syntax
%  synthDT = pmForwardModelTableCreate()
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


% Create the main table with one row and default values
pm              = prfModel;
synthDT         = pm.defaultsTable;
synthDT.pm      = pm;


% We don't want to have the defaults in the table if we did not specify them. 
% First row is going to be the the first values of all the options. 
% This means that we will run the loop twice, first populating the defaults, 
% second, creating the additionnal rows. 
fieldsToCombine = fieldnames(COMBINE_PARAMETERS);
REDUCED_COMBINE_PARAMETERS = COMBINE_PARAMETERS;
for ii=1:length(fieldsToCombine)
    % Construct fieldname
    fieldName   = fieldsToCombine{ii};
    fieldValues = COMBINE_PARAMETERS.(fieldName);
    switch fieldName
        case 'HRF'
            for nh = 1:length(fieldValues)
                % Add rows with the combinations of parameters we want to check
                % Check if only some of the fields in params where set
                % Check the struct params first
                fieldValues(nh).params = pmParamsCompletenessCheck(fieldValues(nh).params, ...
                                                   pm.defaultsTable.HRF.params);
                % Then check that the whole thing should be complete
                fieldValues(nh)        = pmParamsCompletenessCheck(fieldValues(nh), ...
                                            table2struct(pm.defaultsTable.HRF));
                % Convert it to table
                fieldValuesTable = struct2table(fieldValues(nh),'AsArray',true);
                % Change the default if provided
                synthDT.(fieldName) = fieldValuesTable;
            end
        case 'Noise'
        otherwise
            if ~isstruct(fieldValues)
                % Change the default if provided
                synthDT.(fieldName) = fieldValues(1);
                % If it is only one value, add it to defaults and delete it
                if length(fieldValues)==1
                    REDUCED_COMBINE_PARAMETERS = rmfield(REDUCED_COMBINE_PARAMETERS, (fieldName));
                end
            else
                subFieldsToCombine = fieldnames(fieldValues);
                for ii=1:length(subFieldsToCombine)
                    % Construct fieldname
                    subFieldName  = subFieldsToCombine{ii};
                    fieldValues   = getfield(COMBINE_PARAMETERS,fieldName,subFieldName);
                    if ~isstruct(fieldValues)
                        % Change the default if provided
                        synthDT.(fieldName).(subFieldsToCombine{ii}) = fieldValues(1);
                        % If it is only one value, add it to defaults and delete it
                        if length(fieldValues)==1
                            tmpStruct = REDUCED_COMBINE_PARAMETERS.(fieldName);
                            tmpStruct = rmfield(tmpStruct, subFieldsToCombine{ii});
                            REDUCED_COMBINE_PARAMETERS.(fieldName) = tmpStruct;
                        end
                    else
                        error('Only two levels of nesting implemented');
                    end
                    
                end
            end
    end
end



% Extract parameters and keed adding rows.
% BEWARE: THIS GROWS VERY FAST: each line multiplyes the rows of the
% previous one, accumulatively
% The behavior is different for HRF and Noise parameters, a set of parameters
% will generate an HRF or Noise, not their independent components
fieldsToCombine = fieldnames(REDUCED_COMBINE_PARAMETERS);
for ii=1:length(fieldsToCombine)
    % Construct fieldname
    fieldName   = fieldsToCombine{ii};
    fieldValues = REDUCED_COMBINE_PARAMETERS.(fieldName);
    switch fieldName
        case 'HRF'
            for nh = 1:length(fieldValues)
                % Add rows with the combinations of parameters we want to check
                % Check if only some of the fields in params where set
                % Check the struct params first
                fieldValues(nh).params = pmParamsCompletenessCheck(fieldValues(nh).params, ...
                                                   pm.defaultsTable.HRF.params);
                % Then check that the whole thing should be complete
                fieldValues(nh)        = pmParamsCompletenessCheck(fieldValues(nh), ...
                                            table2struct(pm.defaultsTable.HRF));
                % Convert it to table
                fieldValuesTable = struct2table(fieldValues(nh),'AsArray',true);
                % Add all possible combinations based on this options
                synthDT          = pmForwardModelAddRows(synthDT, fieldName,fieldValuesTable);
            end
        case 'Noise'
        otherwise
            if ~isstruct(fieldValues)
                % Add rows with the combinations of parameters we want to check
                synthDT = pmForwardModelAddRows(synthDT, fieldName,fieldValues);
            else
                subFieldsToCombine = fieldnames(fieldValues);
                for ii=1:length(subFieldsToCombine)
                    % Construct fieldname
                    subFieldName  = subFieldsToCombine{ii};
                    fieldValues   = getfield(REDUCED_COMBINE_PARAMETERS,fieldName,subFieldName);
                    if ~isstruct(fieldValues)
                        % Add rows with the combinations of parameters we want to check
                        subFieldName  = [fieldName '.' subFieldsToCombine{ii}];
                        synthDT = pmForwardModelAddRows(synthDT, subFieldName,fieldValues);
                    else
                        error('Only two levels of nesting implemented');
                    end
                    
                end
            end
    end
end

