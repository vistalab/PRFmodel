classdef pmHRF_friston <  pmHRF
    % This is "friston" implementation of the Hemodynamic Response Function
    % Friston et al (1994)
    %
    % Syntax:
    %      hrf = friston();
    %
    % Inputs:
    %
    % Outputs:
    %
    % Optional key/value pairs:
    %
    % Description
    %
    % See also
    %    boyntonHIRF
    
    % Examples
    %{
       
    %}
    
    
    % TODO: copy somewhere above
    % Inputs:
    %   t:  Temporal samples
    %   parms:  Parameters for the Friston function.  Default is from
    %           Friston-Worsley
    %
    %Brief description:
    %
    %   
    % Example:
    %    t = 0:0.1:15
    %    [hirf, params] = fristonHIRF(t);
    %    plot(t,hirf); 
    %    xlabel('Time (sec)'); ylabel('Relative amp'); grid on;
    %
    
    %
    
    
    
    properties (GetAccess=public, SetAccess=public)
        params;
    end
    
    properties (Dependent = true)
        values;
    end
    
    
    %%
    methods
        % Constructor
        function hrf = pmHRF_friston(pm,varargin)
            hrf.PM = pm;
            % Default parameters
            a  = [6  ,  12];
            b  = [0.9, 0.9];
            c  = 0.35;
            params.a = a; params.b = b; params.c = c;
            hrf.params = params;
        end
        
        
        function values = get.values(hrf)
            a = hrf.params.a;
            b = hrf.params.b;
            c = hrf.params.c;
            
            t = hrf.tSteps;
            % Calculate d
            for ii = 1:2
                d(ii) = a(ii)*b(ii);
            end
            % Calculate actual values
            values = (t/d(1)).^a(1)   .* exp(-(t - d(1))/b(1)) ...
                   - c*(t/d(2)).^a(2) .* exp(-(t-d(2))/b(2));
            
        end
    end
    
end