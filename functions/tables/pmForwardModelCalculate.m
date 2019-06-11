function dt = pmForwardModelCalculate(dt)
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
    pm.Stimulus.expName         = dt.Stimulus.ExpName(ii);
    pm.Stimulus.fieldofviewHorz = dt.Stimulus.fovHorz(ii);
    pm.Stimulus.fieldofviewVert = dt.Stimulus.fovVert(ii);
    
    stimName = strcat('Exp-',pm.Stimulus.expName, ...
                      '_binary-', choose(dt.Stimulus.Binary(ii),'true','false'), ...
                      '_size-', num2str(pm.Stimulus.fieldofviewVert), 'x', ...
                               num2str(pm.Stimulus.fieldofviewHorz), '.mat'); 
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
    % TODO: I think it would be more efficient if here we would just s  tore the
    % link to the stimulus, a .mat file on file, that would be use on the
    % calculations and it would be reused. 
    pm.Stimulus.values  =  stimNameWithPath;
    
    
    
    %% Receptive field (RF)
    % Read values from table
    pm.RF.sigmaMajor = dt.RF.sigMajor(ii);
    pm.RF.sigmaMinor = dt.RF.sigMinor(ii);
    pm.RF.Theta      = dt.RF.theta(ii);
    pm.RF.Center     = [dt.RF.x0(ii), dt.RF.y0(ii)];
    % This is computed on the fly once the values are set
    % TODO: do a profiling of the code, I think it is better to set all this
    %       values first, and then invoking pm.RF.compute

    %% Obtain HRF
    % Read values from table
    % TODO: check this assignment, removes values, add compute function
    % pm.HRF          = eval(strcat('pmHRF_', dt.HRF.Type(ii)));
    pm.HRF.Duration = dt.HRF.duration(ii);
    switch dt.HRF.Type(ii)
        case {'friston'}
            a = dt.HRF.params_a(ii,:);
            b = dt.HRF.params_b(ii,:);
            c = dt.HRF.params_c(ii);
            pm.HRF.params.a = [a(1), a(2)];
            pm.HRF.params.b = [b(1), b(2)];
            pm.HRF.params.c = c;
        case {'boynton'}
            disp('TODO')
        otherwise
            error('The %s HRF model not implemented yet',dt.HRF.Type(ii))
    end
    
    % Compute:  no need to compute HRF, calculated on the fly
    
    %% Synthetic time series (remove this)
    %  Compute
    % pm = pm.compute;
    
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
    
    %% Apply different noise models (edit this to add/edit all noise models)
    % Read values from table
    % pm.Noise.Type    = dt.Noise.Type(ii);
    % pm.noise.white_k = dt.Noise.white_k(ii);
    % pm               = pm.noiseCompute;

    % DEBUGGING
    %{
      hold on; 
      plot(pm.BOLD.tSamples, pm.BOLD.predictedWithNoise)
      grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');
    %}
    
    %% Compute the synthetic signal
    pm.compute;
    %% Write back the updated pm model
    dt.pm(ii) = pm;
    
end

