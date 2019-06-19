classdef prfModel_CSS < prfModel
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
        TR;                % Repetition time of the MR acquisition     
    end
    
    properties (Access = private)
        HRF;      % Class: Hemodynamic Response Function
        Stimulus; % Class: Stimulus.
        RF;       % Class: Receptive field. 
        Noise;    % Class: Noise
        BOLD;     % Class: Predicted synthetic time series
    end
    
    %%
    methods        
        % Constructs 
        function pm = prfModel
            pm.TR       = 1;
            pm.HRF      = HRF;
            pm.Stimulus = Stimulus;
            pm.RF       = RF;
            pm.Noise    = Noise;
            pm.BOLD     = BOLD;
        end
        
        function set.TR(pm, tr)
            pm.HRF.TR       = tr;
            pm.Stimulus.TR  = tr;
            pm.BOLD.TR      = tr;
        end
        
        function TR = get.TR(pm)
            TR = pm.TR;
        end
        
    end
    
end

%{

            

            
            
            %% The predicted BOLD signal
            pm.BOLD.timeSeries          = [];
            pm.BOLD.tSamples            = [];
            pm.BOLD.predicted           = [];
            
            %% Noise
            pm.noise.Type               = "white";
            pm.noise.white_k            = 0.5;
            pm.BOLD.predictedWithNoise  = [];
            
            
            %}