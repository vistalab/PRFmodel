%% MAIN SCRIPT
% Generate the parameter table, generate synthetic bold signal, run PRF model, 
%           and compare/visualize results with the synthetically created ones. 
% Add to path:
%    1. winawerlab/analyzePRF
%    2. 

%% Create tables with different parameters and generate synthetic BOLD timeseries

% Generate a default seed table with default values
synthDT = forwardModelTableCreate();
% Visualize the default
%{
    mrvNewGraphWin('predicted with noise');
    plot(synthDT.pm.BOLD.tSamples, synthDT.pm.BOLD.predictedWithNoise);
    grid on;xlabel('Time (sec)');ylabel('Relative amplitude');
%}
              
% Add rows with the combinations of parameters we want to check
% BEWARE: THIS GROWS VERY FAST: each line multiplyes the rows of the previous
%         one, accumulatively
synthDT = forwardModelTableAddRows(synthDT, 'RF.x0',[1,2,3]);
synthDT = forwardModelTableAddRows(synthDT, 'RF.y0',[1,2,3]);
synthDT = forwardModelTableAddRows(synthDT, 'RF.sigMajor',[2,3]);
synthDT = forwardModelTableAddRows(synthDT, 'RF.sigMinor',[2,3]);
% Visualize some examples
%{
    mrvNewGraphWin('predicted with noise');
    pm = synthDT.pm(11);
    pm = synthDT.pm(synthDT.RF.x0==1 & synthDT.RF.y0==2 & ...
                    synthDT.RF.sigMajor==3 & synthDT.RF.sigMinor==1);
    if length(pm)==1
        plot(pm.BOLD.tSamples, pm.BOLD.predictedWithNoise);
        grid on;xlabel('Time (sec)');ylabel('Relative amplitude');
    else
        error('Filter pm to be length one or assign to a specific row')
    end
%}

%% Convert to other formats
% To nifti

% To csv (for Python implementations?)

%% Run different PRF models
% In the first version run one model per row. 
% We will change the stimuli and most of the parameters, even the algorithm
% config, so it makes sense looping all the rows now. 

% analyzePRF
analyzePRF_estimates = calculateFit(synthDT, 'analyzePRF');

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























