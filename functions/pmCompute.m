function pmEstimates = pmCompute(inputTable, prfImplementation)
% Select and apply a PRF model to estimate model parameters
% 
% Inputs:
%   inputTable  - Data table including stimulus parameters and BOLD time series
%   prfImplementation - String defining the model
% 
% Outputs: 
%   pmEstimates: Data table of the estimated pRF model parameters
% 
% Key/val parameters (Optional)
%   N/A
%
% GLU Vistalab 05.2019
%
% See also
%     pmXXX

%  TODO: Make all the outputs the same so that we can use the same function to
%        compare them to the synthetic data. 
% 

% Examples:
%{
  pmCompute...
%}

%%
p = inputParser;
p.addRequired('dataTable',@istable);
p.addRequired('prfImplementation',@ischar);

%% Choose the analysis case

prfImplementation = mrvParamFormat(prfImplementation);

switch prfImplementation
    case {'analyzeprf'}
        % Let's define the format for the estimates so that this is
        % the same for all the methods.
        pmEstimates = table();
        
        % TODO: use parfor if the number of rows is larger than XX
        for ii=1:height(inputTable)
            
            % Obtain the required values for this pRF model
            pm       = inputTable.pm(ii);
            stimulus = stimValuesRead(pm.stimulus.values);
            data     = pm.BOLD.predictedWithNoise;
            TR       = pm.TR;
            options  = struct('seedmode',[0 1],'display','off');
            
            % Calculate PRF
            results  = analyzePRF({stimulus},{data},TR, options);
            
            % Add a new row of results
            pmEstimates = [pmEstimates; struct2table(results,'AsArray',true)];
            
        end
    case {'afni'}
        disp('NYI');
    case {'popeye'}
        disp('NYI');
    otherwise
        error('Method %s not implemented yet.', prfImplementation)
end


end

