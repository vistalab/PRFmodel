classdef HRF <  handle
    % This is a superclass for HRF-s. Every particular instance of this class
    % will have different parameters, so it will be a children class. For
    % example, "Friston" or "Boynton"
    %
    % Syntax:
    %      hrf = HRF();
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
    
    properties (GetAccess=public, SetAccess=public)
                
        Duration  = 20;   % Seconds
        TR        =  1;   % Seconds
        Type      = 'Friston';
        
    end
    
    properties (Dependent = true)
        tSteps;
        values;
    end
    
    
    %%
    methods
        function tSteps = get.tSteps(hrf)
             tSteps     = 0:(hrf.TR):hrf.Duration;
          end

    end
    
end



