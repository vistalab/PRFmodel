classdef prfModel < matlab.mixin.SetGet
    % Initialize a prf model object for the forward time series
    % calculation
    %
    % Syntax:
    %      pm = prfModel(varargin);
    %
    % Inputs:
    %
    % Outputs:
    %
    % Optional key/value pairs:
    %
    % Description
    %
    % See also
    %
    
    % Examples
    %{
       pm = prfModel;
       pm.rfCompute;
    %}
    
    properties (Access = private)
        % What do we need to create our synthetic bold series
        % The main model will have other subclasses with the required
        % components. Analogy: the main pm model is a car, and the subclasses are
        % components of the car (wheel, for example).
        % The color of the wheel is not independent of the color of the car, so
        % the pm model itself will be a property of the subclass as well.
        uniqueTR;          % This is part of the main model. It will be set in all subclasses
    end
    properties (GetAccess=public, SetAccess=private)
        Type            ;
        Stimulus        ;
        RF              ;
        HRF             ;
        Noise           ;
    end
    properties (Dependent)
        TR              ; % This one is derived and copied to all other classes
        timePointsN     ; % Defined by the stimulus.
        timePointsSeries; % Defined by the stimulus.
        defaultsTable   ;
        values          ; % Final value, composed of the BOLD + noise
    end
    
    methods
        % Constructor
        function pm = prfModel
            pm.TR       = 1;
            pm.Type     = 'basic';
            
            % Create the classes
            disp('Creating stimulus ...')
            pm.Stimulus = pmStimulus;     % Create an Stimulus object
            pm.Stimulus.PM = pm;          % Initialize the prfModel inside it
            disp('                  ... created')
            
            disp('Creating HRF ...')
            pm.HRF      = pmHRF_friston;  % Create an HRF object
            pm.HRF.PM   = pm;             % Initialize the prfModel inside it
            disp('                  ... created')
            
            disp('Creating RF ...')
            pm.RF       = pmRF;
            pm.RF.PM    = pm;             % Initialize the prfModel inside it
            disp('                  ... created')
            
            disp('Creating Noise ...')
            pm.Noise{1} = pmNoise_white;  % TODO: this should be able to take several noise models
            pm.Noise{1}.PM = pm;
            pm.Noise{2} = pmNoise_cardiac;
            pm.Noise{2}.PM = pm;
            pm.Noise{3} = pmNoise_respiratory;
            pm.Noise{3}.PM = pm;
            disp('                  ... created')
        end
        
        % Functions that apply the setting of main parameters to subclasses
        % TR: set and get
        function set.TR(pm, tr)
            pm.uniqueTR           = tr;
            % pm.HRF.setTR(tr);
            % pm.Stimulus.setTR(tr);
        end
        function v = get.TR(pm)
            v = pm.uniqueTR;
        end
        % timePointsN and timePointsSeries: Stimulus is the source for the
        % number of timePoints that will be used in the rest of places.
        % NOTE: in a real experiment I think that the number of dicoms should
        % lead this, because then we calculate the stimuli that was shown in
        % every dicom.
        %
        % When creating this synthetic data, we are creating the forward model,
        % where starting with the stimulus we calculate the synthetic BOLD.
        
        function v = get.timePointsN(pm)
            v = pm.Stimulus.timePointsN;
        end
        function v = get.timePointsSeries(pm)
            v  = pm.Stimulus.timePointsSeries;
        end
        
        function defaultsTable = get.defaultsTable(pm)
            % TODO
            % This function will obtain all the defaults from all the
            % subclasses, so that we can construct a parameter table
            defaultsTable = table();
            
            % defaultsTable.Stimulus = pm.Stimulus.defaultsTable;
            % defaultsTable.HRF      = pm.HRF.defaultsTable;
            % defaultsTable.RF       = pm.RF.defaultsTable;
            % defaultsTable.Noise    = pm.Noise.defaultsTable;
            
            
            % pmFindAttributes(pm,'Dependent',true)
            
        end
        
        
        
        function values = get.values(pm)
            disp('Creating synthetic BOLD with noise ...')
            sumOfNoise = zeros(size(pm.BOLD));
            for ii=1:length(pm.Noise)
                sumOfNoise = sumOfNoise + pm.Noise{ii}.values;
            end
            values = pm.BOLD + sumOfNoise;
            disp('                  ... created')
        end
        
        % Plot it
        function plot(pm, what)
            switch what
                case 'noNoise'
                    mrvNewGraphWin([pm.Type 'synthetic BOLD signal without noise']);
                    plot(pm.timePointsSeries, pm.BOLD);
                    grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');
                case 'withNoise'
                    mrvNewGraphWin([pm.Type 'synthetic BOLD signal with noise']);
                    plot(pm.timePointsSeries, pm.values);
                    grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');
                case 'both'
                    mrvNewGraphWin([pm.Type 'synthetic BOLD signal with and without noise']);
                    plot(pm.timePointsSeries, pm.values); hold on;
                    plot(pm.timePointsSeries, pm.BOLD);
                    grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');
                otherwise
                    error('Only noNoise, withNoise and both are acepted')
            end
            
        end
        
    end
    
end
