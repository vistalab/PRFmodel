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
        BOLD            ; % This is the synthetic BOLD signal after the added noise(s)
    end
    
    
    methods
         % Constructs
         function pm = prfModel_basic
             pm.BOLD = [];
         end
         % Methods
    end
       
    methods (Access = public)
        function computeBOLD(pm)
            % Load example stimulus
            stimValues = pm.Stimulus.getStimValues;
            
            % Initialize timeSeries, it is the signal prior to convolution
            timeSeries = zeros(1, pm.timePointsN);
            for tt = 1:pm.timePointsN
                % This is called the hadamard product.  It is the pointwise
                % multiplication of the RF with the stimulus.  The hadProduct is
                % the same size as the stimulus
                hadProduct = stimValues(:,:,tt) .* pm.RF.values;
                % Now, we add up all of the hadProduct values
                timeSeries(tt) = sum(hadProduct(:));
            end
            
            % Convolution between the timeSeries and the HRF
            pm.BOLD = conv(timeSeries, pm.HRF.values, 'same');
        end
        
    end
    
end
