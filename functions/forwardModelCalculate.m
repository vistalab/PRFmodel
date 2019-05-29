function dt = forwardModelCalculate(dt)
% Creates a table with all parameters (defaults) required to perform a forward
% calculation
% 
%  Inputs: table with one or several rows of parameters to calculate bold series
%          in fwd model
% 
%  Outputs: returns the same table with the last column of pm-s calculated
% 
%  See also: forwardModelTableCreate
% 
%  GLU Vistalab 2019.05



for ii=1:height(dt)
    % Do it row to row and parameter to parameter first, for debuging
    
    %% Initialize the basic model with defaults
    pm = dt.pm(ii);

    %% Stimulus
    % If the stimulus have been calculated and exists in /data, load it. 
    % Otherwise, calculate it and save it in /data/stimulus
    % 
    % Read the parameters from the table and generate the file name. 
    pm.stimulus.ExpName         = dt.Stimulus.ExpName(ii);
    pm.stimulus.fieldofviewHorz = dt.Stimulus.fovHorz(ii);
    pm.stimulus.fieldofviewVert = dt.Stimulus.fovVert(ii);
    
    stimName = strcat('Exp-',pm.stimulus.ExpName, ...
                      '_binary-', choose(dt.Stimulus.Binary(ii),'true','false'), ...
                      '_size-', num2str(pm.stimulus.fieldofviewVert), 'x', ...
                               num2str(pm.stimulus.fieldofviewHorz), '.mat'); 
    stimNameWithPath = fullfile(pmRootPath,'data','stimulus',stimName);
    if exist(stimNameWithPath, 'file')
        % Load stimulus produced by s_pmStimulusInterface
        s = load(stimNameWithPath);
    else
        % TODO: convert s_pmStimulusInterface to pmStimulusInterface
        % Write here the routines that calls and saves stimuli based on the
        % parameters of the table.
        s.stim = pmStimulusGenerate('filename', stimNameWithPath);
    end
    % Add the image series of the stimuli, in numeric matrix form to the struct.
    % pm.stimulus.values  =  s.stim;
    % TODO: I think it would be more efficient if here we would just store the
    % link to the stimulus, a .mat file on file, that would be use on the
    % calculations and it would be reused. 
    pm.stimulus.values  =  stimNameWithPath;
    
    % Calculate the X and Y values as well. TODO: do it inline...
    pm = spatialSampleCompute(pm);
    
    
    %% Receptive field (RF)
    pm.RF.sigmaMajor = dt.RF.sigMajor;
    pm.RF.sigmaMinor = dt.RF.sigMinor;
    pm.RF.theta      = dt.RF.theta;
    pm.RF.center     = [dt.RF.x0, dt.RF.y0];
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
end

