classdef prfModel < matlab.mixin.SetGet % handle % matlab.mixin.SetGet, matlab.mixin.Copyable is another possible
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
    properties (GetAccess=public, SetAccess=public) % Changed from SetAccess=private, check
        % Basic, CSS ... (char)
        Type            ; 
        % Components required to build the synthetic bold signal (classes)
        Stimulus        ;
        RF              ;
        HRF             ;
        Noise           ;
        % Other required options (double)
        BOLDMeanValue   ; % Mean value of the synthetic BOLD signal
        % The result: synthetic BOLD series (1 dim array of doubles)
        BOLDnoise       ; % Final value, composed of the BOLD + noise
    end
    properties (Dependent)
        TR              ; % This one is derived and copied to all other classes
        timePointsN     ; % Defined by the stimulus.
        timePointsSeries; % Defined by the stimulus.
        defaultsTable   ;
    end
    
    methods
        % Constructor
        function pm = prfModel
            pm.TR       = 1;
            
            pm.Type     = 'basic';
            % Required mean signal values. 
            pm.BOLDMeanValue = 1000;
            
            % Create the classes, and initialize a prfModel inside it
            pm.Stimulus = pmStimulus(pm);  
            pm.HRF      = pmHRF_friston(pm);  
            pm.RF       = pmRF(pm);       
                       
            % We should initialize the RNG here and return the seed
            
            % TODO: this should be able to take several noise models
            pm.Noise{1} = pmNoise_white(pm);  
            pm.Noise{2} = pmNoise_cardiac(pm);
            pm.Noise{3} = pmNoise_respiratory(pm);
            
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
        
        % 
        function compute(pm)
            % Computes the mean BOLD response and then adds noise.
            
            % Every sub-class has a computeBOLD function to compute
            % the mean response.
            pm.computeBOLD;
            
            % Here is how we add the noise terms.
            sumOfNoise = zeros(size(pm.BOLD));
            for ii=1:length(pm.Noise)
                pm.Noise{ii}.compute;
                sumOfNoise = sumOfNoise + pm.Noise{ii}.values;
            end
            
            pm.BOLDnoise = pm.BOLD + sumOfNoise;
        end
        
        % Plot it
        function plot(pm, what)
            what = mrvParamFormat(what);
            switch what
                case 'nonoise'
                    mrvNewGraphWin([pm.Type 'Synthetic BOLD signal (no noise)']);
                    plot(pm.timePointsSeries, pm.BOLD);
                    grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');
                case 'withnoise'
                    mrvNewGraphWin([pm.Type 'Synthetic BOLD signal (noise)']);
                    plot(pm.timePointsSeries, pm.BOLDnoise);
                    grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');
                case 'both'
                    mrvNewGraphWin([pm.Type 'Synthetic BOLD signals']);
                    plot(pm.timePointsSeries, pm.BOLD); hold on;
                    plot(pm.timePointsSeries, pm.BOLDnoise);
                    grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');
                    legend({'No Noise','With Noise'})
                otherwise
                    error("Only 'no noise', 'with noise' and 'both' are acepted")
            end
            
        end
        
    end
    
end
