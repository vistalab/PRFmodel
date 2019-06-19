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
pm         = prfModel;
synthDT    = pm.defaultsTable;
synthDT.pm = pm;

% Extract parameters and keed adding rows.
% BEWARE: THIS GROWS VERY FAST: each line multiplyes the rows of the
% previous one, accumulatively

fieldsToCombine = fieldnames(COMBINE_PARAMETERS);
for ii=1:length(fieldsToCombine)
    % Construct fieldname
    fieldName  = fieldsToCombine{ii};
    fieldValues=getfield(COMBINE_PARAMETERS,fieldName);
    if ~isstruct(fieldValues)
        % Add rows with the combinations of parameters we want to check
        synthDT = pmForwardModelAddRows(synthDT, fieldName,fieldValues);
    else
        subFieldsToCombine = fieldnames(fieldValues);
        for ii=1:length(subFieldsToCombine)
            % Construct fieldname
            subFieldName  = subFieldsToCombine{ii};
            fieldValues   = getfield(COMBINE_PARAMETERS,fieldName,subFieldName);
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

