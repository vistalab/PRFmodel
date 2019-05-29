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
        
        BOLD;      % Struct of the predicted time series and params
        
    end
    
    %%
    methods
        
        % Construct
        function pm = prfModel(varargin)
            
            varargin = mrvParamFormat(varargin);
            
            p = inputParser;
            p.addParameter('tr',2,@isnumeric);                  % Seconds
            p.addParameter('values', [],@isnumeric);            % actual values, matrix
            p.addParameter('fieldofviewHorz',[],@isnumeric);    % Deg
            p.addParameter('fieldofviewVert',[],@isnumeric);    % Deg
            p.addParameter('hrfduration'    ,20,@isnumeric);    % Seconds
            
            p.parse(varargin{:});
            
            %% MR parameters
            pm.TR         = p.Results.tr;
            
               
            %% Initialize time steps and HRF
            pm.HRF.duration  = p.Results.hrfduration;
            pm.HRF.tSteps    = 0:(pm.TR):pm.HRF.duration;   % For now, always to 20 sec
            pm.HRF.modelName = 'friston';
            pm.HRF.Friston_a = [];
            pm.HRF.Friston_b = [];
            pm.HRF.Friston_c = [];
            pm               = pm.HRFget;
            
            %% The sequence of stimulus images
            
            % We will probably attach the stimulus movie (time series of images)
            % and write a method that converts the movie to the binary stimulus
            % for us.
            pm.stimulus.fieldofviewHorz = p.Results.fieldofviewHorz;
            pm.stimulus.fieldofviewVert = p.Results.fieldofviewVert;
            pm.stimulus.values          = p.Results.values;
            
            % Obtain the X and Y values
            pm                          = pm.spatialSampleCompute;
            
            %% The receptive field parameters
            
            pm.RF.center     = [0 0];    % Deg
            pm.RF.theta      = 0;        % Radians
            pm.RF.sigmaMajor = 1;        % deg
            pm.RF.sigmaMinor = 1;        % deg
            
            
            %% The predicted BOLD signal
            pm.BOLD.timeSeries          = [];
            pm.BOLD.tSamples            = [];
            pm.BOLD.predicted           = [];
            
            %% Noise
            pm.noise.Type               = "white";
            pm.noise.white_k            = 0.5;
            pm.BOLD.predictedWithNoise  = [];
            
            
        end
    end
    
end
