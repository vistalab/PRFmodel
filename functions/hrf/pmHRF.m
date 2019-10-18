function [HRF, tSteps, params] = getHRF(hrfName, varargin)
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

%% Interpret input arguments
varargin = mrvParamFormat(varargin);

p = inputParser;
p.addRequired('hrfname',@ischar);     % Which HRF model?
p.addParameter('timesteps',(0:0.1:20),@isvector);  % Time steps in seconds
p.addParameter('params',[],@isstruct);          % Structure contaning parameters for HRF model

p.parse(hrfName,varargin{:});

tSteps  = p.Results.timesteps;
hrfName = p.Results.hrfname;
params  = p.Results.params;

%%

switch hrfName
    case 'friston'
        
        if isempty(params), [HRF,params] = fristonHIRF(tSteps);
        else, fristonHIRF(tSteps,params);
        end
        

    otherwise
        error('Unknown HRF type %s\n',hrfName);
        
end


end


figure
for HRFtype=pm.HRF.Types
    pm.HRF.Type = HRFtype{:};
    pm.HRF.plot; hold on;
end
legend(strrep(pm.HRF.Types,'_','\_'))
title('Different HRF implementations')