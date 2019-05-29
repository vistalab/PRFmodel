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
    % Do it row to row and parameter to parameter first, for debugging
    
    %% Initialize the basic model with defaults
    pm = dt.pm(ii);
    
    %% TR
    pm.TR = dt.TR(ii);

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
    pm = pm.spatialSampleCompute;
    
    
    %% Receptive field (RF)
    % Read values from table
    pm.RF.sigmaMajor = dt.RF.sigMajor(ii);
    pm.RF.sigmaMinor = dt.RF.sigMinor(ii);
    pm.RF.theta      = dt.RF.theta(ii);
    pm.RF.center     = [dt.RF.x0(ii), dt.RF.y0(ii)];
    % Compute
    pm = pm.rfCompute;

    %% Obtain HRF
    % Read values from table
    pm.HRF.modelName = dt.HRF.Type(ii);
    pm.HRF.duration  = dt.HRF.duration(ii);
    pm.HRF.tSteps    = 0:(pm.TR):pm.HRF.duration; % Calculated
                   a = dt.HRF.Friston_a(ii,:);
                   b = dt.HRF.Friston_b(ii,:);
    pm.HRF.Friston_a = [a(1), a(2)];
    pm.HRF.Friston_b = [b(1), b(2)];
    pm.HRF.Friston_c = dt.HRF.Friston_c(ii);
    % Compute
    pm               = pm.HRFget;
    % TODO: Fix the whole parameter thing. Depending on the model we will have
    % different parameters. The output is creating params now but it shuold
    % update the parameters if they don't exist. 
    
    %% Synthetic time series
    % This function performs the Hadamard product between the RF and the stimuli,
    % and then the convolution with the hrf signal. 
    pm = pm.timeSeriesCompute;
    
    % DEBUGGING
    % Visualize the predicted time series
    %{
      tSteps = 0:size(pm.HRF.values,2)-1;
      mrvNewGraphWin; plot(tSteps, pm.HRF.values)
      grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');
      
      tSteps = pm.BOLD.tSamples;
      mrvNewGraphWin('predicted'); plot(tSteps, pm.BOLD.predicted)
      grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');
    %}    
    
    %% Apply different noise models
    % Read values from table
    pm.noise.Type    = dt.Noise.Type(ii);
    pm.noise.white_k = dt.Noise.white_k(ii);
    pm               = pm.noiseCompute;

    % DEBUGGING
    %{
      hold on; 
      plot(pm.BOLD.tSamples, pm.BOLD.predictedWithNoise)
      grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');
    %}
    
    %% Write back the updated volume
    dt.pm(ii) = pm;
    
end

