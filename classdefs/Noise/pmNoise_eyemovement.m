classdef pmNoise_eyemovement < pmNoise
    % This is the white noise implementation of a Noise class. 
    % 
    %
    % Syntax:
    %      whiteNoise = white;
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
    
    properties (GetAccess=public, SetAccess=private)
        params;  
    end
    
    properties (Dependent = true)
        values;
    end
    
    
    %%
    methods
        % Constructor
        function noise = pmNoise_eyemovement
             noise.params.k      =  0.5; % Value between 0 (no noise) and 1 (proportional to mean signal)
        end
        % Methods available to this class and his childrens (friston, boynton... classes)
        function values = get.values(noise)
            n      = noise.params.k * mean(noise.baselineBOLD);
            values = n * randn(size(noise.baselineBOLD));
        end

    end
    
end


