%% MAIN SCRIPT
%
% Generate the parameter table, generate synthetic bold signal, run
% PRF model, and compare/visualize results with the synthetically
% created ones.
%
% Add to path:
%    1. winawerlab/analyzePRF
%    2. 

%% Include the stimulus generation here

%% Create tables with different pRF parameters 
%
% The entries in these different tables will be used to generate
% synthetic BOLD timeseries.  These will be analyzed and the estimates
% will be compared with the parameters set here.

% Generate a default seed table with default values
synthDT = forwardModelTableCreate();
             
% Add rows with the combinations of parameters we want to check
% BEWARE: THIS GROWS VERY FAST: each line multiplyes the rows of the
% previous one, accumulatively

synthDT = forwardModelTableAddRows(synthDT, 'RF.x0',[1,2,3]);
synthDT = forwardModelTableAddRows(synthDT, 'RF.y0',[1,2,3]);
synthDT = forwardModelTableAddRows(synthDT, 'RF.sigMajor',[2,3]);
synthDT = forwardModelTableAddRows(synthDT, 'RF.sigMinor',[2,3]);

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
analyzePRF_estimates = pmCompute(synthDT, 'analyzePRF');

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























