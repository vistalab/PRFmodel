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
        Type;
    end
    
    
    methods
         % Constructs
         function pm = prfModel_basic
             pm.BOLD = [];
             pm.Type = 'basic';
         end
         
    end
       
    methods (Access = public)
 
        
    end
    
end








%{
plot(timeSeries)
plot(pm.BOLD)

100*(   (max(pm.BOLD)-min(pm.BOLD))  /mean(pm.BOLD))
%}













