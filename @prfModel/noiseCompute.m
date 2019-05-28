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





%% Calculate
scaleNoise = 0.5;  % Multiplies the mean signal value
pm.BOLD.predictedWithNoise = pmNoiseWhite(pm.BOLD.predicted, scaleNoise);


end