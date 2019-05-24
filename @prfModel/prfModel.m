classdef prfModel 
    % Initialize a prf model object for the forward time series
    % calculation
    %
    % Syntax:
    %      pm = prfModel(varargin);
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
       pm = prfModel;
       pm.rfCompute;
    %}
    
    properties (GetAccess=public, SetAccess=public)
        
        TR;        % Repetition time of the MR acquisition
        
        HRF;       % Struct of hemodynamic response function variables
        
        stimulus;  % Struct of stimulus properties, including images
        
        RF;        % Struct of receptive field properties
        
        noise;     % Struct of BOLD noise properties
        
    end
    
    %%
    methods
        
        % Construct
        function pm = prfModel(varargin)
            
            varargin = mrvParamFormat(varargin);
            
            p = inputParser;
            p.addParameter('tr',2,@isnumeric);              % Seconds
            p.addParameter('binarystimulus',[],@isnumeric); % Binary stimulus
            p.addParameter('fieldofview',[],@isnumeric);    % Binary stimulus
            p.addParameter('hrfduration',20,@isnumeric);    % Seconds
            
            p.parse(varargin{:});
            
            %% MR parameters
            pm.TR  = p.Results.tr;
            
            
            %% Initialize time steps and HRF
            pm.HRF.duration          = p.Results.hrfduration;
            pm.HRF.tSteps = 0:(pm.TR):pm.HRF.duration;   % For now, always to 20 sec
            pm.HRF.values = pmHRF('friston','time steps',pm.HRF.tSteps);
            
            %% The sequence of stimulus images
            
            % We will probably attach the stimulus movie (time series of images)
            % and write a method that converts the movie to the binary stimulus
            % for us.
            pm.stimulus.fieldofview = p.Results.fieldofview;
            pm.stimulus.binary      = p.Results.binarystimulus;
            
            % Spatial field of view
            
            % This will be a derived quantity
            spatialSample = pm.stimulus.fieldofview/size(pm.stimulus.binary,2);
            
            x = (spatialSample:spatialSample:pm.stimulus.fieldofview);
            x = x - mean(x);
            % mrvNewGraphWin; plot(x,x); grid on;
            
            % Set the spatial sampling parameters
            y = x;
            [X,Y] = meshgrid(x,y);
            pm.stimulus.X = X;
            pm.stimulus.Y = Y;
            
            %% The receptive field parameters
            
            pm.RF.center = [0 0];    % Deg
            pm.RF.theta =   0;       % Radians
            pm.RF.sigmaMajor = 1;    % deg
            pm.RF.sigmaMinor = 1;    % deg
            
        end
    end
    
end
