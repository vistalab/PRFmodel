classdef pmRF <  matlab.mixin.SetGet
    % This is a superclass for RF-s. Every particular instance of this class
    % will have different parameters, so it will be a children class.
    %
    % Syntax:
    %      rf = RF();
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
        PM;         % prfModel that has some of the variables we need, such as TR
        Center;     % Deg
        Theta;      % Radians
        sigmaMajor; % Deg
        sigmaMinor; % Deg
    end
    
    properties (Dependent = true, Access = public)
        values;
    end
    properties(Dependent= true, SetAccess = private, GetAccess = public)
        TR;            % Seconds, it will be read from the parent class pm
    end
    
    
    
    %%
    methods
        % Constructor
        function rf = pmRF
            %% The receptive field parameters
            rf.Center     = [0 0];    % Deg
            rf.Theta      = 0;        % Radians
            rf.sigmaMajor = 1;        % deg
            rf.sigmaMinor = 1;        % deg
        end
        function v = get.TR(rf)
            v      = rf.PM.TR;
        end
        % Methods available to this class and childrens, if any
        function values = get.values(rf)
            XY = rf.PM.Stimulus.XY;
            values = rfGaussian2d(XY{1}, XY{2}, ...
                rf.sigmaMajor,rf.sigmaMinor,rf.Theta, ...
                rf.Center(1),rf.Center(2));
        end
        % Plot it
        function plot(rf)
            % I think we don't want this to be stored in the object.
            % Calculate it and return every time we need it.
            mrvNewGraphWin('Receptive Field');
            mesh(rf.values);
            grid on; xlabel('x'); ylabel('y');
        end
    end
    
end



