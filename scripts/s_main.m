% Use this script to create the different calls to the generators

% Create table with all the different options and then store in the last columns
% the pm with the generated synthetic

% Parameters affecting the simulations
% Use the same classification as the class structs. 
synthDT = forwardModelTableCreate()

       
       
% Now keep adding rows with all the combinations we want to check       
       
       

%% Initialize the basic model with defaults
pm = prfModel();

%% Stimulus
% Load stimulus produced by s_pmStimulusInterface
s = load(fullfile(pmRootPath,'data','Exp103_onlyMask_Downsampled_Resized.mat'));
% Normalize to 0>1 and binarize 
nstim   = s.stim - min(s.stim(:));
nstim   = nstim ./ max(nstim(:));
% Add it to the pm we just created. 
pm.stimulus.binary = imbinarize(nstim,.5);
% We can add it, or create it here
%    pm.stimulusCreate('aperture type');
%    pm.stimulusBinarize;
%    pm.plot('stimulus movie');

% Add the field of view information
%   This should be pm.set() and the set() function should make sure that all the
%   places where fieldofview is used, is updated. 
pm.stimulus.fieldofviewHorz = 20;
pm.stimulus.fieldofviewVert = 20;

%% Receptive field (RF)
pm = pm.rfCompute;
% Visualize the receptive field (RF)
% pm.plot('receptive field')

%% Synthetic time series
% This function performs the Hadamard product between the RF and the stimuli,
% and then the convolution with the hrf signal. 

% Use default 20s Friston HRF. Otherwise change it now. 
% pm = pm.getHRF('Boynton');  % For example
pm = pm.timeSeriesCompute;
% Visualize the predicted time series
% mrvNewGraphWin; plot(tSteps,HRF)
% grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');
mrvNewGraphWin('predicted'); plot(pm.BOLD.predicted)
grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');



%% Apply different noise models
pm = pm.noiseCompute;

hold on; 
plot(pm.BOLD.predictedWithNoise)
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


