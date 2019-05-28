function pm = getHRF(pm)
% Create the HRF from different models
%
%
%
%
%

% Examples:
%{
  [HRF, tSteps] = pmHRF('friston');
  mrvNewGraphWin;
  plot(tSteps, HRF, '--'); 
  grid on; xlabel('secs');
%}
%{
  [HRF, tSteps] = pmHRF('friston','time steps',(0:2:20));
  mrvNewGraphWin;
  plot(tSteps, HRF, '--'); 
  grid on; xlabel('secs');
%}

%% delete this, now is a method in a class
% %% Interpret input arguments
% varargin = mrvParamFormat(varargin);
% 
% p = inputParser;
% p.addRequired('hrfname',@ischar);     % Which HRF model?
% p.addParameter('timesteps',(0:0.1:20),@isvector);  % Time steps in seconds
% p.addParameter('params',[],@isstruct);          % Structure contaning parameters for HRF model
% 
% p.parse(hrfName,varargin{:});
% 
% tSteps  = p.Results.timesteps;
% hrfName = p.Results.hrfname;
% params  = p.Results.params;

%%

switch pm.HRF.modelName
    case 'friston'
        [HRF, params] = fristonHIRF(pm.HRF.tSteps);
      
    otherwise
        error('Unknown HRF type %s\n',hrfName);
        
end

pm.HRF.params = params;
pm.HRF.values = HRF;

end
