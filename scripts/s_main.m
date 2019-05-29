%% MAIN SCRIPT
% Add to path
% winawerlab/analyzePRF

%% Create tables with different parameters and generate synthetic BOLD timeseries

% Generate a default seed table with default values
synthDT = forwardModelTableCreate();
% Visualize the default
%{
    mrvNewGraphWin('predicted with noise');
    plot(synthDT.pm.BOLD.tSamples, synthDT.pm.BOLD.predictedWithNoise);
    grid on;xlabel('Time (sec)');ylabel('Relative amplitude');
%}
              
%% Add rows with the combinations of parameters we want to check
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

%% Run different PRF models
% In the first version run one model per row. 
% We will change the stimuli and most of the parameters, even the algorithm
% config, so it makes sense looping all the rows now. 

% analyzePRF
analyzePRF_estimates = calculateFit(synthDT, 'analyzePRF');


%% Compare models
% Method 1



% Method 2