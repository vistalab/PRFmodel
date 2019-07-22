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
            % We are only interested in the first array if there are more than one
            nh          = 1;
            completeHRF = pmParamsCompletenessCheck(fieldValues(nh), ...
                                            table2struct(pm.defaultsTable.HRF));
            % Do the same with params
            completeHRF.params = pmParamsCompletenessCheck(completeHRF.params, ...
                                                    pm.defaultsTable.HRF.params);
            % Convert it to table
            fieldValuesTable = struct2table(completeHRF,'AsArray',true);
            % Change the default
            synthDT.(fieldName) = fieldValuesTable;
            % If there is only only one HRF, then delete the variable
            if length(fieldValues)==1
                REDUCED_COMBINE_PARAMETERS = rmfield(REDUCED_COMBINE_PARAMETERS,'HRF');
            end
        % case 'Noise'
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
                    if strcmp(subFieldName,'sigmaMinor') && strcmp(fieldValues,'same')
                        fieldValues = getfield(COMBINE_PARAMETERS,fieldName,'sigmaMajor');
                    end
                    if ~isstruct(fieldValues)
                        % Change the default if provided
                        if strcmp(fieldName,'Noise')
                            
                        else
                            synthDT.(fieldName).(subFieldsToCombine{ii}) = fieldValues(1);
                        end
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
            for nh = 2:length(fieldValues)  % the first one is the default
                % Add rows with the combinations of parameters we want to check
                % Then check that the whole thing should be complete
                completeHRF        = pmParamsCompletenessCheck(fieldValues(nh), ...
                                            table2struct(pm.defaultsTable.HRF));
                % Check if only some of the fields in params where set
                % Check the struct params first
                completeHRF.params = pmParamsCompletenessCheck(completeHRF.params, ...
                                                   pm.defaultsTable.HRF.params);
                
                % Convert it to table
                fieldValuesTable = struct2table(completeHRF,'AsArray',true);
                % Add all possible combinations based on this options
                synthDT          = pmForwardModelAddRows(synthDT, fieldName, fieldValuesTable);
            end
        % case 'Noise'
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
                        % Check if sigmaMinor is 'same'. If it is do nothing
                        % here, but everytime sigmaMajor is changed, change
                        % minor too. 
                        if strcmp(subFieldName,'sigmaMinor') && strcmp(fieldValues,'same')
                            continue
                        elseif strcmp(subFieldName,'sigmaMajor') && strcmp(REDUCED_COMBINE_PARAMETERS.RF.sigmaMinor,'same')
                            % Add rows with the combinations of parameters we want to check
                            subFieldName          = [fieldName '.' subFieldsToCombine{ii}];
                            synthDT               = pmForwardModelAddRows(synthDT, subFieldName,fieldValues);
                            synthDT.RF.sigmaMinor = synthDT.RF.sigmaMajor;
                        else
                            % Add rows with the combinations of parameters we want to check
                            subFieldName  = [fieldName '.' subFieldsToCombine{ii}];
                            synthDT = pmForwardModelAddRows(synthDT, subFieldName,fieldValues);
                        end
                        
                    else
                        error('Only two levels of nesting implemented');
                    end
                    
                end
            end
    end
end

