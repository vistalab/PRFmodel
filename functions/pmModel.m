% Read cached stimuli
%
% It was created by this script:
%   s_pmStimulusInterface
%
% mrvNewGraphWin('check');
% imagesc(stim(:,:,66)); colormap gray; axis equal tight off;
%
% Example:
%{
  pm = prfModel;                       % Initialize an instance of a prf model for one voxel
  pm.set('thisParameter',thisValue);   % Adjust the parameters as you like
  pm.set ....                          % Other parameters
  pm.predictTimeSeries;
  thisTS = pm.get('time series');
  pm.plot('time series');
%}
%% Clear up.

% close all;

%%

% Produced by s_pmStimulusInterface
thisStimulus = load(fullfile(pmRootPath,'data','Exp103_onlyMask_Downsampled_Resized.mat'));

pm = prfModel('TR', 1,...
    'binary stimulus', thisStimulus.stim,...
    'field of view', 20);

pm = pm.rfCompute;

%% Make the binarized mask only for the forward calculation (Linear or CSS style)


%% pm.stimulusBinarize

% This might be
%
%    pm.stimulusCreate('aperture type');
%    pm.stimulusBinarize;
%    pm.plot('stimulus movie');

% Make decision regarding binarization, color use, etc. and add them to s_pmStimulusInterface.m
nstim   = stim - min(stim(:));
nstim   = nstim ./ max(nstim(:));
pm.stimulus.binary = imbinarize(nstim,.5);

% pm.plot('binary stimulus movie');
% mrvNewGraphWin; plot(tSteps,HRF)
% grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');

%% pm.rfCompute
RF = pm.rfCompute;

% Gaussian2d(X,Y,sigmaMajor,sigmaMinor,theta, x0,y0);
% pm.plot('receptive field');

% mesh(x,y,RF)


%% pm.timeSeriesCompute

% BOLD signal
% Predicted prf response
pm.timeSeries = zeros(1,tSamples);
for tt = 1:tSamples
    % This is called the hadamard product.  It is the pointwise
    % multiplication of the RF with the stimulus.  The hadProduct is
    % the same size as the stimulus
    hadProduct = pm.stimulus.binary(:,:,tt) .* pm.RF;
    
    % Now, we add up all of the hadProduct values
    pm.timeSeries(tt) = sum(hadProduct(:));
end

% Here is the time series.  This is the signal prior to convolution
% with the hemodynamic response function
% mrvNewGraphWin; plot(1:tSamples,timeSeries);

%% pm.plot('time series');

predictedTS = conv(timeSeries,HRF,'same');
mrvNewGraphWin('predicted'); plot(predictedTS)
grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');

%% Set the parameters and build the second RF


% mesh(x,y,RF)

%% Apply different noise models
scaleNoise = 0.5;  % Multiplies the mean signal value
predictedTS = pmNoiseWhite(predictedTS, scaleNoise);

% Eye motion jitter
eyeMotionJitter = 1;  % Deg

% Motion related (translation and rotation)

% Cardiac Related

% Respiration related

% Low frequency physiological fluctuations

% Draining veins

% Low frequency drifts
% (slow head displacements, scanner related (e.g. heating...)

% Hardware related instabilities

hold on; 
plot(predictedTS)
grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');

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


