classdef pmNoise_white < pmNoise
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
    
    properties (GetAccess=public, SetAccess=public)
        params;
        values;
    end
    
    properties (Dependent = true)
        % values;
    end
    
    
    %%
    methods
        % Constructor
        function noise = pmNoise_white(pm)
            noise.PM            = pm;
            noise.Type          = 'white';
            noise.params.k      =  0.1; % Value between 0 (no noise) and 1 (proportional to amplitude of signal)
        end
        
        % Methods available to this class 
        
        % COMPUTE
        function compute(noise)
            n            = noise.params.k * (max(noise.PM.BOLD) - min(noise.PM.BOLD));
            noise.values = n * randn([1,noise.PM.timePointsN]);
        end
        
    end
    
end


