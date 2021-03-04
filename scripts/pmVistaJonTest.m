% I would like to have the files and code needed to call pmModelFit. And to make
% it fast, I would like the data files to be small (say, 10 voxels). And to be
% able to test multiple scan functionality, I would like to have predicted BOLD
% time series for the three different stimuli. Doesn't have to be the three
% stimuli from the HCP. It could just be three sets of bar



%% SELECT SYNTH PARAMS
COMBINE_PARAMETERS.TR                    = [1.5];
COMBINE_PARAMETERS.RF                    = struct();
COMBINE_PARAMETERS.RF.Centerx0           = [0,3]; 
COMBINE_PARAMETERS.RF.Centery0           = [0,3];
COMBINE_PARAMETERS.RF.Theta              = [0]; %, deg2rad(45)];
COMBINE_PARAMETERS.RF.sigmaMajor         = [1,2];% [1,2];
COMBINE_PARAMETERS.RF.sigmaMinor         = "same";
COMBINE_PARAMETERS.Stimulus.durationSecs = 300;
HRF(1).Type = 'vista_twogammas';
COMBINE_PARAMETERS.HRF           = HRF;
Noise(1).seed                    = 'none'; % 'random'
COMBINE_PARAMETERS.Noise         = Noise;
    
%% SYNTHESIZE
synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS, 'repeats', 1);
pmForwardModelCalculate(synthDT, 'useparallel',false,...
                                 'writefiles',true,...
                                 'outputdir',fullfile(pmRootPath,'local'),...
                                 'subjectname','vista2test'); 

%% CREATE INPUT
niftiBOLDfile  = fullfile(pmRootPath,'local','vista2test.nii.gz');
jsonSynthFile  = fullfile(pmRootPath,'local','vista2test.json');
stimNiftiFname = fullfile(pmRootPath,'local','vista2test_Stim.nii.gz');

input          = {niftiBOLDfile, jsonSynthFile, stimNiftiFname};

%% VISTA OPTIONS
allOptions       = struct();
allOptions.vista = struct();
allOptions       = pmParamsCompletenessCheck(allOptions, allOptions);
% Edit vista options, otherwise uses defaults
%{
    options.vista            = allOptions.vista;
    options.vista.model      = 'one gaussian';
    options.vista.grid       = false;  % if true, returns gFit
    options.vista.wSearch    = 'coarse to fine and hrf'; 
    options.vista.detrend    = 0;
    options.vista.keepAllPoints            = true; 
    options.vista.numberStimulusGridPoints =  50;  
%}    


%% SOLVE IT 
results = pmModelFit(input,'vistasoft','options',options);