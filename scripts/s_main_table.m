%% MAIN SCRIPT (it works again, but I still need to do many other things)
% 
% This script is a wrapper that generates ARRAYS OF synthetic BOLD signal and 
% estimates calculated with different pRF implementations. It is intended to be
% used as a data analysis tool to combine a selected combinations of
% parameters and compare performamces of different tools. 
% 
% 
%
% STEPS:
% 1.- Create a table with parameters required to generate the synthetic BOLD
%     signal(s) we want to use. It can be one or many. 
%     We populate a table with all the parameters)
%   1.1.- Stimulus
%   1.2.- RF  (Receptive field)
%   1.3.- HRF (Hemodynamic Response Function)
%   1.4.- Noise
% 
% 2.- Calculate the forward model BOLD time series per every possible
%     combination of parameters provided in STEP 1. 
%     This creates a synthetic BOLD signal TEST dataset, which is matlab table. 
% 
% 3.- Save or export the TEST dataset to other formats. For example: nifti for AFNI. 
% 
% 4.- Estimate parameters for a given model. For example: analyzePRF, or AFNI. 
%     The result will be formated into a table, with the same format for all
%     different models. 
%     We will have a table per each set of results, that we will store as well. 
% 
% 5.- Compare synthetic BOLD data with the different result. 
%   5.1.- Create stats  
%   5.2.- Create visualizations
% 
% 
%
% To make this work:
%    0. add PRFmodel to the path
%    1. to test analyzePRF add winawerlab/analyzePRF to the path
%    2. (TODO: vistasoft implementation, maybe we'll add just the functions we want)
%    3. (TODO: test an AFNI wrapper, and requirement will be to have it installed)
%    4. 
% 
% 
% % See also:
%     s_pmTest.m

%% STEP 1: create a table with parameters that will generate synthetic BOLD signal
% The values that are not set-up below will use defaults. 
% To see default parameters, execute:
% clear classes
pm       = prfModel;
DEFAULTS = pm.defaultsTable
% 
% Each parameter below can be unique or an array of values: 
%     - If only one value is provided per every parameter, it will substitute
%       the defaults, and generate only one synthetic BOLD signal. 
%     - If one or more arrays of values are provided, it will generate all
%       possible combinations between parameters. TODO: provide more control to this...
% 
% How to include parameters:
%     - Visualize default values
%     - Create struct with the same organization as the table. 
%     - (DO NOT ADD DEFAULTS, they will be in the first row of the synthDT table)

% COMBINE_PARAMETERS.TR            = [2];      % Don't add 1
% COMBINE_PARAMETERS.Type          = {'CSS'};  % Don't add basic
COMBINE_PARAMETERS.RF.Centerx0   = [1,2];    % etc ...
COMBINE_PARAMETERS.RF.Centery0   = [1];
COMBINE_PARAMETERS.RF.sigmaMinor = [1];
COMBINE_PARAMETERS.RF.sigmaMajor = [1,2];
% Some parameters go in groups, because params are specific to the choice
% If we change only some, the rest will remain defaults. 
% HRF
HRF(1).Type                      = 'friston';
HRF(1).Duration                  = 18;
HRF(1).params.c                  = 0.5;

HRF(2).Type                      = 'canonical';
HRF(2).params.stimDuration       = 2;

COMBINE_PARAMETERS.HRF           = HRF;
% Noise
% Noise is an array of noises. 
% For simplicity, we are going to make it fixed: i.e. always the same noise 
% models will be present. As with the other cases, just specify what we want to
% change, and the rest will be maintained with the default
% Noise(1).Type                      = 'friston';
% Noise(1).params.noise2                  = 0.5;
% % Another one with the three default ones
% HRF(2).Type                      = 'canonical';
% HRF(2).params.stimDuration       = 2;


% TODO: add the possibility to repeat the same line N times (to check for random effects)

% 
% The entries in these different tables will be used to generate
% synthetic BOLD timeseries.  These will be analyzed and the estimates
% will be compared with the parameters set here.
% NOTE: row(1)from synthDT will always be the defaults, for testing purposes
synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS);



%% Compute the pm-s using the parameters in each row 
% To see the big list of cases, just type
%   synthDT
% The table, which might be quite large, will print out
       
% Once all the desired combinations of parameters have been created, 
% compute the forward models       
synthDT = pmForwardModelCalculate(synthDT);

% 
% The last column is a pm, for example
%{
thisPM = synthDT.pm(11);
thisPM.HRF.plot;
thisPM.Stimulus.plot;
thisPM.plot('both');

thisPM = synthDT.pm(111);
thisPM.plot('both');

%}

%% Visualize examples of the generated time series (UPDATE)
%{
    mrvNewGraphWin('predicted with noise');
    pm = synthDT.pm(11);
    pm = synthDT.pm(synthDT.RF.x0==1 & synthDT.RF.y0==2 & ...
                    synthDT.RF.sigMajor==3 & synthDT.RF.sigMinor==1);
    % Multiple
    pm = synthDT.pm(synthDT.RF.x0==1 &  ...
                    synthDT.RF.sigMajor==3 & synthDT.RF.sigMinor==1);

    if numel(pm)==1
        plot(pm.BOLD.tSamples, pm.BOLD.predictedWithNoise);
        grid on;xlabel('Time (sec)');ylabel('Relative amplitude');
    else
        error('Filter pm to be length one or assign to a specific row')
    end
%}

%% Convert the synthetic BOLD time series to other formats

% To nifti

% To csv (for Python implementations?)

% To JSON



%% Run different PRF analysis models

% In the first version run one model per row. 
% We will change the stimuli and most of the parameters, even the algorithm
% config, so it makes sense looping all the rows now. 

% analyzePRF
analyzePRF_estimates = pmModelFit(synthDT, 'analyzePRF');

%% Compare models
% Method 1



% Method 2

%% Plot synthetic data versus predicted
% Select the 11th row parameters and results. 
% We can use filtering over the parameters as well
ii = 11;

% Select if we want to save the plots as images for offline validation/publishing
savePlot     = false;
FileName = fullfile(pmRootPath,'data','figures',['analyzePRF_Voxel-' num2str(ii) '.svg']);  % .svg for publishing
plotSynthvsPredicted(synthDT.pm(ii), analyzePRF_estimates(ii,:),'savePlot',false,'plotFileName',FileName);
                            
%%























