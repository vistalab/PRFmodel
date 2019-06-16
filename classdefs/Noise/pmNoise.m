classdef pmNoise <  matlab.mixin.SetGet
    % This is a superclass for Noise-s. Every particular instance of this class
    % will have different parameters, so it will be a children class. For
    % example:
    %   - White noise (white)
    %
    %   - Eye motion jitter
    %     eyeMotionJitter = 1;  % Deg
    %
    %   - Motion related (translation and rotation)
    %
    %   - Cardiac Related
    %
    %   - Respiration related
    %
    %   - Low frequency physiological fluctuations
    %
    %   - Draining veins
    %
    %   - Low frequency drifts
    %     (slow head displacements, scanner related (e.g. heating...)
    %
    %   - Hardware related instabilities
    %
    % Syntax:
    %      noise = Noise;
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
        PM;           % prfModel that has some of the variables we need, such as TR
        Type;         % 'White', ...
    end
    
    properties(Dependent= true, SetAccess = private, GetAccess = public)
        TR;            % Seconds, it will be read from the parent class pm
    end
    
    
    
    %%
    methods
        % Constructor
        function noise = pmNoise
        end
        % Methods available to this class and his childrens (friston, boynton... classes)
        function TR = get.TR(noise)
            TR = noise.PM.TR;
        end
        
        % Plot it
        function plot(noise)
            % I think we don't want this to be stored in the object.
            % Calculate it and return every time we need it.
            mrvNewGraphWin('Noise to be added to signal');
            plot(noise.PM.timePointsSeries, noise.values);
            grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');
        end
    end
    
end


