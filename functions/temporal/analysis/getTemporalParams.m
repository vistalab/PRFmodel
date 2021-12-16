function params = getTemporalParams(params)
%% Get temporal parameters (getTemporalParams)

% Grabs defualt temporal params for each model (params.analysis.temporalModel) 

% The function first checks if there is a Constant file (user-defined temporal parameter values).
% if the function does not exist, defualt temporal parameters are loaded

% input:
% params.analysis.temporalModel: '1ch-dcts', '2ch-exp-sig','1ch-glm'

% Return values:
% fields: temporal parameter names
% param:  temporal parameters
% fs:     sampling rate (ms)
% tr:     TR (sec)
% num_channels: number of channels

% [ISK] note: st_getTemporalAttributes is the previous function name

%%

% check if Constant file is there. If there is a constant file load params
% from the Constants file
if exist('Constants') == 2
    c = Constants.getTemporalParams.temporalParams;
    for i = 1:length(c)
        if strcmp(c{i}.type, params.analysis.temporalModel)
            idx = i;
        end
    end
    temporal_param = c{idx}.prm;
    fs             = c{idx}.fs;
    num_channels   = c{idx}.num_channels;
    fields         = c{idx}.fields;
    tr = Constants.getTemporalParams.tr; % seconds
    
else % load defualt temporal params if Constant file is not there   
    % defualt TR is 1 sec
    tr = 1; % sec
    switch params.analysis.temporalModel
        case {'1ch-dcts'} % load 1ch-dcts (DN) model params
            num_channels = 1;
            fs = 1000;
            fields = ["tau1", "weight", "tau2", "nn", "delay", "shift", "scale"];
            temporal_param = [0.05 0 0.1 2 0.1 0 1];
        case {'2ch-exp-sig'}  % load 2ch-exp-sig (2ch) model params
            num_channels = 2;
            fs = 1000;
            fields = ["tau_s", "tau_ae", "Lp", "Kp", "Kn", "weight","shift"];
            temporal_param = [4.93 10000 0.1 3 3 0.5 0];
        case {'1ch-glm'} % load 1ch-glm (linear) model params
            num_channels = 1;
            fs = 1000;
            fields = ["shift", "scale"];
            temporal_param = [0 1];
    end
      
end

% pass values to params
params.analysis.temporal.model        = params.analysis.temporalModel;
params.analysis.temporal.fields       = fields;
params.analysis.temporal.fs  = fs;
params.analysis.temporal.param        = temporal_param;
params.analysis.temporal.num_channels = num_channels;
params.analysis.temporal.tr           = tr;

end