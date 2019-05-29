function estimates = calculateFit(dt, prfImplementation)
% Calculates the PRF and returns a table with estimations 
% 
%  Inputs: variable name and value(s)
% 
%  Outputs: Matlab table
% 
% 

switch lower(prfImplementation)
    case {'analyzeprf'}
        estimates = table();
        % TODO: use parfor
        for ii=1:height(dt)
            % Obtain all the required values
            pm       = dt.pm(ii);
            stimulus = stimValuesRead(pm.stimulus.values);
            data     = pm.BOLD.predictedWithNoise;
            TR       = pm.TR;
            options  = struct('seedmode',[0 1],'display','off');
            % Calculate PRF
            results  = analyzePRF(stimulus,data,TR, options);
            % Add a new row of results
            estimates = [estimates; struct2table(results,'AsArray',true)];
        end
    otherwise
        error('Method %s not implemented yet.', prfImplementation)
end


end

