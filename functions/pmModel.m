% Read cached stimuli
%
% It was created by this script:
%   s_pmStimulusInterface

%% Clear up.
clear all; close all; clc; 

%% Make the binarized mask only for the forward calculation (Linear or CSS style)

% Produced by s_pmStimulusInterface
load(fullfile(pmRootPath,'data','Exp103_onlyMask_Downsampled_Resized.mat'))

%%
TR  = 1;    % Imagine we want the HRF at every TR

% Make decision regarding binarization, color use, etc. and add them to s_pmStimulusInterface.m
nstim   = stim - min(stim(:));
nstim   = nstim ./ max(nstim(:));
pm.stim = imbinarize(nstim,.5);
tSamples = size(stim,3);

% mrvNewGraphWin('check');
% imagesc(stim(:,:,66)); colormap gray; axis equal tight off;

%%  Get an HRF
tSteps = 0:TR:20;
[HRF,params] = fristonHIRF(tSteps);
% mrvNewGraphWin; plot(tSteps,HRF)
% grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');

%%
% Define variables that will be looped afterwards
pm.fieldofview = 20;  % Deg

% This will be a derived quantity
spatialSample = 20/size(pm.stim,2);

x = (spatialSample:spatialSample:pm.fieldofview);
x = x - mean(x);
% mrvNewGraphWin; plot(x,x); grid on;

%% Set the parameters and build one RF

y = x;
[X,Y] = meshgrid(x,y);

% Center
x0 = 2;         % Deg
y0 = 3;         % Deg

% Spread
sigma = 1;      % Deg

% Orientation
theta = pi/4;   % Radians
sigmaMajor = 4;
sigmaMinor = 4;

RF = rfGaussian2d(X,Y,sigmaMajor,sigmaMinor,theta, x0,y0);
% mesh(x,y,RF)

% BOLD signal
% Predicted prf response
timeSeries = zeros(1,tSamples);
for tt = 1:tSamples
    % This is called the hadamard product.  It is the pointwise
    % multiplication of the RF with the stimulus.  The hadProduct is
    % the same size as the stimulus
    hadProduct = pm.stim(:,:,tt) .* RF;
    
    % Now, we add up all of the hadProduct values
    timeSeries(tt) = sum(hadProduct(:));
end

% Here is the time series.  This is the signal prior to convolution
% with the hemodynamic response function
% mrvNewGraphWin; plot(1:tSamples,timeSeries);

predictedTS = conv(timeSeries,HRF,'same');
mrvNewGraphWin('predicted'); plot(predictedTS)
grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');

%% Set the parameters and build the second RF

y = x;
[X,Y] = meshgrid(x,y);

% Center
x0 = 2;         % Deg
y0 = 3;         % Deg

% Orientation
theta = pi/4;   % Radians
sigmaMajor = 4;
sigmaMinor = 4;

RF = rfGaussian2d(X,Y,sigmaMajor,sigmaMinor,theta, x0,y0);
% mesh(x,y,RF)

% BOLD signal
% Predicted prf response
timeSeries = zeros(1,tSamples);
for tt = 1:tSamples
    % This is called the hadamard product.  It is the pointwise
    % multiplication of the RF with the stimulus.  The hadProduct is
    % the same size as the stimulus
    hadProduct = pm.stim(:,:,tt) .* RF;
    
    % Now, we add up all of the hadProduct values
    timeSeries(tt) = sum(hadProduct(:));
end

predictedTS = conv(timeSeries,HRF,'same');


% Apply different noise models
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


