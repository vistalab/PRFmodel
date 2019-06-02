function pm = HRFget(pm)
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

switch lower(pm.HRF.modelName)
    case 'friston'
        if (isempty(pm.HRF.Friston_a) || isempty(pm.HRF.Friston_b) || isempty(pm.HRF.Friston_c))
            [HRF, params] = fristonHIRF(pm.HRF.tSteps);
        else
            params.a = pm.HRF.Friston_a;
            params.b = pm.HRF.Friston_b;
            params.c = pm.HRF.Friston_c;
            HRF      = fristonHIRF(pm.HRF.tSteps, params);
        end
    otherwise
        error('Unknown HRF type %s\n',hrfName);
        
end

pm.HRF.params = params;
pm.HRF.values = HRF;



%% MORE HRF stuff
    %{
    This is one from the default at the Winawer lab in analyzePRF
    testHIRF = getcanonicalhrf(TR,TR);
    mrvNewGraphWin; plot(testHIRF);
    set(gca,'xlim',[0 20]);
    grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');

    %% Another random one

    TR  = 1;    % Imagine we want the HRF at every TR
    tSteps = 0:TR/4:20;
    HRF = fristonHIRF(tSteps,params);
    mrvNewGraphWin; plot(tSteps,HRF)
    grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');

    %}



end
