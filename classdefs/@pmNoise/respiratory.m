classdef pmNoise_respiratory < pmNoise
    % This is the respiratory noise implementation of a Noise class. 
    %   Implements a sinusoidal noise with 0.3 frequency rate and amplitude
    %   relative to the mean BOLD signal
    %
    % Syntax:
    %      noise = pmNoise_respiratory;
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
    
    properties
        params;
        values
    end
    
    % properties (Dependent = true)
    %    values;
    % end
    
    
    %%
    methods
        % Constructor
        function noise = pmNoise_respiratory(pm)
            noise.PM                    =  pm;
            noise.Type                  =  'respiratory';
            noise.params.frequency      =  0.3;  % 0.3 Hz : 18 breaths/min
            noise.params.amplitude      =  0.1;  % 0-1 proportion over BOLD signal amplitude
            
        end
        % Methods available to this class 
                % COMPUTE
        function compute(noise)
            signal = noise.PM.BOLD;
            % Frequency [Hz]
            fNoise = noise.params.frequency;
            % Amplitude
            aNoise = noise.params.amplitude * (max(noise.PM.BOLD) - min(noise.PM.BOLD)); 
            t      = noise.PM.timePointsSeries;
            % Calculate the noise
            noise.values = aNoise*sin(2*pi.*t.*fNoise);
            
        end

    end
    
end


