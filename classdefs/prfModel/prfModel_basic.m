classdef prfModel_basic < prfModel
    % Initialize a prf model object for the forward time series
    % calculation
    %
    % Syntax:
    %      pm = prfModel_basic(varargin);
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
       pm = prfModel_basic;  % It is a special case of a prfModel, for the basic
                             % implementation, others are CSS...
    %}
    
    %     properties (GetAccess=public, SetAccess=public)
    %
    %     end
    %
    
    properties (SetAccess = private)
        % This is the synthetic BOLD signal with several types of noise
        BOLD; 
    end
    
    
    methods
         % Constructs
         function pm = prfModel_basic
             pm.BOLD = [];
         end
         
    end
       
    methods (Access = public)
        function computeBOLD(pm,varargin)
            % For this linear model we take the inner product of the
            % receptive field with the stimulus contrast.  Then we
            % convolve that value with the HRF.
            %
            
            varargin = mrvParamFormat(varargin);
            p = inputParser;
            p.addRequired('pm',@(x)(isa(x,'prfModel')));
            p.addParameter('randomseed',1000,@(x)(round(x) == x && x > 0));
            p.parse(pm,varargin{:});
            
            % Load example stimulus
            stimValues = pm.Stimulus.getStimValues;
            
            % Initialize timeSeries, it is the signal prior to convolution
            [r,c,t] = size(stimValues);
            spaceStim = reshape(stimValues,r*c,t);

            % I had didactic code here that GL used at first. But that
            % took 9 secs.  This code takes 0.1 sec.  So harder to
            % read, but let's use it.
            timeSeries  = pm.RF.values(:)' * spaceStim;
            
            % Convolution between the timeSeries and the HRF
            pm.BOLD = conv(timeSeries, pm.HRF.values, 'same');
        end
        
    end
    
end
