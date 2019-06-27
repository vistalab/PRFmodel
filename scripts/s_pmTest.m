%% TEST SCRIPT: A single synthetic BOLD signal generation and testing
%
% This script is a wrapper that generates A synthetic BOLD signal and
% estimates calculated with different pRF implementations. It is intended
% to be used as a learning and testing tool. Once it is mastered and tested
% with a single BOLD series that this is what we want to do, we can create
% multiple variations in s_main_table.m, analyze them with different
% implementations and test them.
%
% The unit element of this software is a prfModel class.
%       pm = prfModel;

%% Create default values:
pm = prfModel;

% TODO: if not strictly necessary, change Center to an array
pm.RF.Centerx0    = 0;
pm.RF.Centery0    = 0;
pm.RF.Theta       = 0;        % Degrees, x-axis = 0, positive y-axis 90
pm.RF.sigmaMajor  = 1;      % Degrees, x-axis = 0, positive y-axis 90
pm.RF.sigmaMinor  = 1;      % Degrees, x-axis = 0, positive y-axis 90

pm.TR               = 2;

% Compute
pm.compute; 

% Fit model
options = struct('seedmode', [0 1], 'display' , 'iter');
results = pmModelFit(pm, 'analyzePRF', 'options',options);
mrvNewGraphWin('data and fit'); 
plot(1:pm.TR:pm.TR*pm.timePointsN,results.testdata(1,:));hold on;
plot(1:pm.TR:pm.TR*pm.timePointsN,results.modelpred(1,:));
legend({'testdata','modelpred'})
text(1, max(results.testdata(1,:)) , sprintf('R^2:%2.2f',results.R2(1)))


%{
% How to transform the results parameters into a center and size all
% in degrees
%
rfThetaRadians = deg2rad(pm.RF.Theta);
rfCenter  = pm.RF.Center;
rfSizeDeg = [pm.RF.sigmaMajor, pm.RF.sigmaMinor];
mrvNewGraphWin;
ellipsePlot(rfCenter,rfSizeDeg,rfThetaRadians);
%}

% Computes noise free and noise BOLD signals.
% NOTE: the call to compute calls the compute functions for the rest of
%       subclasses. In this case, after the RF parameter changes, the call to
%       pm.RF.compute only would necessary if we wanted to visualize the RF with
%       pm.RF.plot
pm.compute; 


% Visualize them
pm.plot('with noise');
            
pm.Noise{1}.params.noise2signal = 0;
pm.BOLDmeanValue = 1000;
pm.compute;
pm.plot('with noise');

% Visualize them
pm.plot('no noise');

% Visualize both on the same graph
pm.plot('both');

% And compute it with, for example, analyzePRF
results = pmModelFit(pm, 'analyzePRF');


%{
% How to transform the results parameters into a center and size all
% in degrees
%
thetaRadians = deg2rad(results.ang);
[cX,cY] = pol2cart(theta,results.ecc);
[cX,cY] = pol2cart(thetaRadians,results.ecc);
rfSizeDeg = pm.RF.sigmaMajor;
mrvNewGraphWin;
ellipsePlot(rfSizeDeg,rfSizeDeg,thetaRadians);
%}
% Evaluation: TODO
% TODO: make sense of the analyzePRF output and relate to the input

%%
% We can make changes to the parameters to check both the signal and the results.
pm.TR = 1.82;

% We need to compute again the BOLD signal after any change in params
pm.compute;
pm.plot('both')

% CHECK/VISUALIZE/CHANGE THE INDIVIDUAL COMPONENTS.
% The logical order is:

%% 1./ STIMULUS:
% Visualize
pm.Stimulus.plot

% pm.Stimulus.toVideo (creates a video and saves it in folder local)
% Change a parameter
pm.Stimulus.Binary   = false; % TODO: it will generate the same thing
pm.Stimulus.barWidth =     2; % TODO: it will generate the same thing

% Compute new stimulus based on the new parameters
pm.Stimulus.compute; % If exists the file, read it, otherwise, create it and save it.

% Visualize it
pm.Stimulus.plot

% See how the predicted synthetic BOLD signal changed
pm.compute;
pm.plot('with noise');

%% 2./ RF: Receptive Field
% Visualize
pm.RF.plot
% Change a couple of parameters
pm.RF.sigmaMajor = 3;
pm.RF.Theta      = pi/2;
% Compute 
pm.RF.compute;
% Visualize
pm.RF.plot
% See the new predicted synthetic BOLD signal
% We need to compute it first.
pm.compute;
pm.plot('with noise');

%{
     radSpacing = 0.01;
     a = 1; b = 3; theta = pi/4;
     [x,y] = ellipsePoints(a,b,theta,radSpacing);
     x = [x,x(1)]; y = [y, y(1)];   % Close it up for plotting
     plot(x,y,'-'); axis equal
 
%}

%% 3./ HRF: Hemodynamic Response Function
% Visualize
pm.HRF.plot
% Change a parameter
pm.HRF.params.c = 0.7;
% Compute 
pm.HRF.compute;
% Visualize 
pm.HRF.plot
% Change HRF from Friston to Canonical
pm.HRF = pmHRF(pm, 'Type','boynton');
% Compute 
pm.HRF.compute;
% Visualize 
pm.HRF.plot

%% 4./ Noise
% Visualize the default noise models
pm.Noise
% Visualize respiratory
pm.Noise{3}.plot
% Change a parameter
pm.Noise{3}.params.amplitude = 0.9;
% Compute the changes
pm.Noise{3}.compute;
% Visualize it 
pm.Noise{3}.plot
% Compute and visualize the new synthetic BOLD with noise
pm.compute;
pm.plot('with noise');


%% Some PRF tests: use demo data from Winawer's tutorials
addpath(genpath('~/winawerlab/analyzePRF'))
setup
load('exampledataset.mat');
data = data{1};
data = data(1,:);
stimulus = stimulus{1};
stimulus = stimulus(:,:,1:2:300);
data = double(data);
stimulus = double(stimulus);
size(data)
size(stimulus)

TR       = 2;
options  = struct('seedmode',[0 1],...
    'display','iter');
% Calculate PRF
results  = analyzePRF({stimulus}, {data}, TR, options);
mrvNewGraphWin('data and fit'); 
plot(1:TR:300,results.testdata(1,:));hold on;
plot(1:TR:300,results.modelpred(1,:))
legend({'testdata','modelpred'})
text(1, max(results.testdata(1,:)) , sprintf('R^2:%2.2f',results.R2(1)))


%% END