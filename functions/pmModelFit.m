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
        % Let's define the format for the estimates so that this is
        % the same for all the methods.
        pmEstimates = table();
        
        if istable(input)
            % TODO: if is not a table, create a table of 1 row and use
            % the same loop.
            %
            % TODO: use parfor if the number of rows if the table is
            % larger than XX
            for ii=1:height(input)

                % Obtain the required values for this pRF model
                pm       = input.pm(ii);
                stimulus = pm.Stimulus.getStimValues;
                data     = pm.BOLDnoise;
                TR       = pm.TR;
                options  = struct('seedmode',[0 1],...
                    'display','off',...
                    'dosave','modelpred');

                % Calculate PRF
                results  = analyzePRF({stimulus}, {data}, TR, options);
                
                % TODO: make "results" the same format for everybody
                 
                % Add a new row of results
                pmEstimates = [pmEstimates; struct2table(results,'AsArray',true)];
            end
        else
            % Obtain the required values for this pRF model
            pm       = input;
            stimulus = pm.Stimulus.getStimValues;
            data     = pm.BOLDnoise;
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

