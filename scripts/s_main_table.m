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
%     s_main_table.m

%% STEP 1: create a table with parameters that will generate synthetic BOLD signal
% The values that are not set-up below will use defaults. To see default parameters, execute:
%      forwardModelTableCreate()   % to see a table
%              or 
%      prfModel()                  % to see a struct
% 
% Each parameter below can be unique or an array of values: 
%     - If only one value is provided per every parameter, it will substitute
%       the defaults, and generate only one synthetic BOLD signal. 
%     - If one or more arrays of values are provided, it will generate all
%       possible combinations between parameters. TODO: provide more control to this...
% 
% How to include parameters: 
%     - TO DO
% 
% The entries in these different tables will be used to generate
% synthetic BOLD timeseries.  These will be analyzed and the estimates
% will be compared with the parameters set here.

synthDT = pmForwardModelTableCreate();
% TODO: make this into a function
            % Generate a default seed table with default values
            synthDT = pmForwardModelTableCreate();

            % Add rows with the combinations of parameters we want to check
            % BEWARE: THIS GROWS VERY FAST: each line multiplyes the rows of the
            % previous one, accumulatively
            synthDT = pmForwardModelAddRows(synthDT, 'RF.x0',[1,2,3]);
            synthDT = pmForwardModelAddRows(synthDT, 'RF.y0',[1,2,3]);
            synthDT = pmForwardModelAddRows(synthDT, 'RF.sigMajor',[2,3]);
            synthDT = pmForwardModelAddRows(synthDT, 'RF.sigMinor',[2,3]);

%% Create tables with different pRF parameters 
       
% Once all the desired combinations of parameters have been created, 
% compute the forward models       
synthDT = pmForwardModelCalculate(synthDT);

% To see the big list of cases, just type
%   synthDT
% The table, which might be quite large, will print out
% 
% The last column is a pm, for example
%
%    thisPM = synthDT.pm(11);
%    mrvNewGraphWin;  plot(thisPM.HRF.tSteps,thisPM.HRF.values)
%    grid on;

%% Visualize examples of the generated time series
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























