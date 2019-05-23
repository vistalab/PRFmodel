% Read cached stimuli
clear all; close all; clc; 
load(fullfile(pmRootPath,'data','Exp103_onlyMask_Downsampled_Resized.mat'))

% Working with mask only, make sure the stim is binarized
% Make decision regarding binarization, color use, etc. and add them to s_pmStimulusInterface.m
nstim = stim - min(stim(:));
nstim = nstim ./ max(nstim(:));
stim = imbinarize(nstim,.5);
mrvNewGraphWin; imagesc(stim(:,:,66)); colormap gray; axis equal tight off;




% Define variables that will be looped afterwards
fieldRange = 20;  % Deg
sampleRate = 0.2; % Deg
x = [-fieldRange:sampleRate:fieldRange];
y = x;
[X,Y] = meshgrid(x,y);
sigma = 5;  % Deg
rf = rfGaussian2d(X,Y,sigma);



%% Function definition

% PRF model
g = exp(  -((x-x0)^2+(y-y0)^2)/2*sigma^2     );

% Stimulus aperture
% s = (x,y,t), binary, one per time point
% Read one example:
load(fullfile(pmRootPath,'data','examples','stimEx.mat'));

% Predicted prf response
r = s .* g;
