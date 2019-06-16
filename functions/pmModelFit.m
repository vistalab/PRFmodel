function [pmEstimates, results] = pmModelFit(input, prfImplementation)
% Select and apply a PRF model to estimate model parameters
% 
% Syntax:
%    [pmEstimates, results] = pmModelFit(inputTable, prfImplementation)
%  
% Brief description:
%    Takes a table whose rows describe the BOLD time series in a
%    single voxel (along with its metadata) and estimates the pRF
%    model parameters for some implementation.
%
% Inputs:
%   inputTable  - Data table including stimulus parameters and BOLD
%                 time series. Each row of the data table is another
%                 voxel
%   prfImplementation - String defining the model
%
% Outputs: 
%   pmEstimates: Table format of the pRF model parameters in results
%   results:     The struct from analyzePRF
%
% Key/val parameters (Optional)
%   N/A
%
% GLU Vistalab 05.2019
%
% See also:
%     pmXXX
%

%  TODO: Make all the outputs the same so that we can use the same function to
%        compare them to the synthetic data. 
% 

% Examples:
%{
  pmCompute...
%}

%%
p = inputParser;
p.addRequired('input');
p.addRequired('prfImplementation',@ischar);

%% Choose the analysis case

prfImplementation = mrvParamFormat(prfImplementation);

switch prfImplementation
    case {'analyzeprf'}
        % TODO: Let's define the format for the estimates so that this is
        %       the same for all the methods.
        pmEstimates = table();
        
        % Check if the pm-s come in a table or alone
        if ~istable(input)
            temp = table();
            temp.pm = input;
            input = temp;
        end
        
        % Go line by line and compute required values for each pRF model
        for ii=1:height(input)
            % TODO: use parfor if the number of rows if the table is larger than XX
            pm       = input.pm(ii);
            stimulus = double(pm.Stimulus.getStimValues);
            data     = pm.BOLDnoise;
            TR       = pm.TR;
            options  = struct('seedmode',[0 1],...
                              'display','off',...
                              'dosave','modelpred');
            
            % Calculate PRF
            results  = analyzePRF({stimulus}, {data}, TR, options);
            
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

