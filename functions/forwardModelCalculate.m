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
    if exist(fullfile(pmRootPath,'data','stimulus',stimName), 'file')
        % Load stimulus produced by s_pmStimulusInterface
        s = load(fullfile(pmRootPath,'data','stimulus',stimName));
    else
        % TODO: convert s_pmStimulusInterface to pmStimulusInterface
        % Write here the routines that calls and saves stimuli based on the
        % parameters of the table.
    end
                           
    
    
    
    
    
    

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
end

