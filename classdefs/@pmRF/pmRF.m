classdef pmRF <   matlab.mixin.SetGet & matlab.mixin.Copyable
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
        Centerx0;   % Deg
        Centery0;   % Deg
        Theta;      % Radians
        sigmaMajor; % Deg
        sigmaMinor; % Deg
        % values;
    end
    properties (SetAccess = private, GetAccess = public)
         values;    % Result. Only can be changes from within this func.
    end
    properties(Dependent= true, SetAccess = private, GetAccess = public)
        TR;            % Seconds, it will be read from the parent class pm
    end
    
    
    
    %%
    methods (Static)
        function d = defaultsGet
            d.Centerx0   = 0;        % Deg
            d.Centery0   = 0;        % Deg
            d.Theta      = 0;        % Radians
            d.sigmaMajor = 1;        % deg
            d.sigmaMinor = 1;        % deg
            % Convert to table and return
            d = struct2table(d,'AsArray',true);
        end
    end
    methods
        % Constructor
        function rf = pmRF(pm, varargin)
            % Obtain defaults table. If a parameters is not passed, it will use
            % the default one defined in the static function
            d = rf.defaultsGet;
            % Read the inputs
            varargin = mrvParamFormat(varargin);
            p = inputParser;
            p.addRequired ('pm'        ,              @(x)(isa(x,'prfModel')));
            p.addParameter('centerx0'  ,d.Centerx0  , @isnumeric);
            p.addParameter('centery0'  ,d.Centery0  , @isnumeric);
            p.addParameter('theta'     ,d.Theta     , @isnumeric);
            p.addParameter('sigmamajor',d.sigmaMajor, @isnumeric);
            p.addParameter('sigmaminor',d.sigmaMinor, @isnumeric);
            p.parse(pm,varargin{:});
            
            % Initialize the pm model and hrf model parameters
            rf.PM           = pm;
            % The receptive field parameters
            rf.Centerx0     = p.Results.centerx0;
            rf.Centery0     = p.Results.centery0;
            rf.Theta        = p.Results.theta;
            rf.sigmaMajor   = p.Results.sigmamajor;
            rf.sigmaMinor   = p.Results.sigmaminor;
        end
        
        function v = get.TR(rf)
            v      = rf.PM.TR;
        end
        
        % Methods available to this class and childrens, if any
        function compute(rf)
            % Compute stimulus just in case
            rf.PM.Stimulus.compute;
            % Obtain grid XY
            XY = rf.PM.Stimulus.XY;
            % Calculate values
            rf.values = rfGaussian2d(XY{1}, XY{2}, ...
                        rf.sigmaMajor,rf.sigmaMinor,rf.Theta, ...
                        rf.Centerx0,rf.Centery0);
        end
                
        
        
        % Plot it
        function plot(rf)
            % Compute before plotting
            rf.compute
            % Plot it
            mrvNewGraphWin('Receptive Field');
            mesh(rf.values);
            grid on; xlabel('x'); ylabel('y');
        end
    end
    
end



