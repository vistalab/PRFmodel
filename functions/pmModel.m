
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
