classdef pmHRF_boynton <  pmHRF
    % This is "boynton" implementation of the Hemodynamic Response Function
    % Boynton et al (1996)
    %
    % Syntax:
    %      hrf = pmHRF_boynton;
    %      hrf.plot
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
    %
    
    % Examples
    %{
       
    %}
   
    
    
    %%
    properties (GetAccess=public, SetAccess=public)
        params;
    end
    
    properties (Dependent = true)
        values;
    end
    
    
    %%
    methods
        % Constructor
        function hrf = pmHRF_boynton
            hrf.Type = 'boynton';
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
            values = (t/d(1)).^a(1) .* exp(-(t - d(1))/b(1)) ...
                   - c*(t/d(2)).^a(2) .* exp(-(t-d(2))/b(2));
            
        end
    end
    
end



