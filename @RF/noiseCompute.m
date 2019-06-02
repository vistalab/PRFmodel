function pm = noiseCompute(pm)
% Add noise to the predicted signal
%
%
% See also


% Eye motion jitter
% eyeMotionJitter = 1;  % Deg

% Motion related (translation and rotation)

% Cardiac Related

% Respiration related

% Low frequency physiological fluctuations

% Draining veins

% Low frequency drifts
% (slow head displacements, scanner related (e.g. heating...)

% Hardware related instabilities


%% Check
if isempty(pm.BOLD.predicted)
    error("There are no values in BOLD.predicted where to add noise")
end


%% Calculate
switch lower(pm.noise.Type)
    case 'white'
        if isempty(pm.noise.white_k)
            scaleNoise = 0.5;  % Multiplies the mean signal value
        else
            scaleNoise = pm.noise.white_k;
        end
        pm.BOLD.predictedWithNoise = pmNoiseWhite(pm.BOLD.predicted, scaleNoise);
    otherwise
        error('Unknown noise type %s\n', pm.noise.Type);
end



end