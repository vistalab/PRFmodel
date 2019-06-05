classdef pmNoise_cardiac < pmNoise
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
    
    properties (GetAccess=public, SetAccess=private)
        params;  
    end
    
    properties (Dependent = true)
        values;
    end
    
    
    %%
    methods
        % Constructor
        function noise = pmNoise_cardiac
             noise.Type                  =  'cardiac';  
             noise.params.frequency      =  1.25;  % 1.25 Hz:75 beats/min
             noise.params.amplitude      =  0.50;  % 0-1 proportion over mean BOLD signal
             
        end
        % Methods available to this class and his childrens (friston, boynton... classes)
        function values = get.values(noise)
            signal = noise.PM.BOLD;
            fNoise = noise.params.frequency;                % Frequency [Hz]
            aNoise = noise.params.amplitude * mean(signal); % Amplitude
            t      = noise.PM.timePointsSeries;
            % Calculate the noise
            values = aNoise*sin(2*pi.*t.*fNoise);
        end

    end
    
end

