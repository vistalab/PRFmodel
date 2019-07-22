classdef prfModel < matlab.mixin.SetGet & matlab.mixin.Copyable
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
        BOLDmeanValue   ; % Required mean value of the synthetic BOLD signal
        BOLDcontrast    ; % Required  contrast  of the synthetic BOLD signal
        timeSeries      ;
        % BOLD signal value (before noise)
        BOLD            ;
        % The result: synthetic BOLD series (1 dim array of doubles)
        BOLDnoise       ; % Final value, composed of the BOLD + noise
    end
    properties (Dependent)
        TR              ; % This one is derived and copied to all other classes
        timePointsN     ; % Defined by the stimulus.
        timePointsSeries; % Defined by the stimulus.
        defaultsTable   ;
    end
    
    methods (Static)
        function d = defaultsGet
            d.TR            = 1;
            d.Type          = 'basic';
            d.BOLDmeanValue = 10000;
            d.BOLDcontrast  = 8;
            % Convert to table and return
            d = struct2table(d,'AsArray',true);
        end
    end
    methods
        % Constructor
        function pm = prfModel(varargin)
            % Obtain defaults table. If a parameters is not passed, it will use
            % the default one defined in the static function
            d = pm.defaultsGet;
            % Make varargin lower case, remove white spaces...
            varargin = mrvParamFormat(varargin);
            % Parse the inputs/assign default values
            p = inputParser;
            p.addParameter('tr'           ,  d.TR           , @isnumeric);
            p.addParameter('type'         ,  d.Type{:}      , @ischar);
            p.addParameter('boldmeanValue',  d.BOLDmeanValue, @isnumeric);
            p.addParameter('boldcontrast' ,  d.BOLDcontrast , @isnumeric);
            
            
            p.parse(varargin{:});
            % Assign defaults/parameters to class/variables
            pm.TR            = p.Results.tr;
            pm.Type          = p.Results.type;
            pm.BOLDmeanValue = 10000;
            pm.BOLDcontrast  = 8;  % In percentage
            
            % Create the classes, and initialize a prfModel inside it
            pm.Stimulus      = pmStimulus(pm);
            pm.HRF           = pmHRF(pm); 
            pm.RF            = pmRF(pm);
            
            % We should initialize the RNG here and return the seed
            
            % pm.Noise takes one or several noise models that will be added
            % recursively and independently to the bold signal, i.e. each noise
            % is independent from each other.
            % Add parameters like this:,'params',struct('frequency',1.3)
            pm.Noise = {pmNoise(pm, 'Type','white'), ...
                        pmNoise(pm, 'Type','cardiac'), ...
                        pmNoise(pm, 'Type','respiratory')};
            % pm.Noise = {};
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
            % This function obtains all the defaults from all the
            % subclasses, so that we can construct a parameter table
            defaultsTable           = pm.defaultsGet;
            defaultsTable.HRF       = pm.HRF.defaultsGet;
            defaultsTable.RF        = pm.RF.defaultsGet;
            defaultsTable.Stimulus  = pm.Stimulus.defaultsGet;
            % Noise are several noise models by default, add them all
            % TODO: think about this. 
            %       right now, return only params
            % for nn = 1:length(pm.Noise)
            %     tmp                 = pm.Noise{nn}.defaultsGet;
            %     tmp.Type            = pm.Noise{nn}.Type;
            %     defaultsTable.(['Noise_' tmp.Type]) = tmp;
            % end
            defaultsTable.Noise  = pm.Noise{1}.defaultsGet;
        end
        % Compute synthetic BOLD without noise
        function computeBOLD(pm,varargin)
            % For this linear model we take the inner product of the
            % receptive field with the stimulus contrast. Then we
            % convolve that value with the HRF.
            
            varargin = mrvParamFormat(varargin);
            p = inputParser;
            p.addRequired('pm',@(x)(isa(x,'prfModel')));
            p.addParameter('randomseed',1000,@(x)(round(x) == x && x > 0));
            p.parse(pm,varargin{:});
            
            % First, compute the values in the required sub-classes
            % TODO: optimize this to not repeat operations
            pm.Stimulus.compute;
            pm.RF.compute;
            pm.HRF.compute;
                        
            % Load stimulus
            stimValues = pm.Stimulus.getStimValues;
            
            % Initialize timeSeries, it is the signal prior to convolution
            [r,c,t]   = size(stimValues);
            spaceStim = reshape(stimValues,r*c,t);
            
            switch pm.Type
                case 'basic'
                    % Calculate time series
                    pm.timeSeries  = spaceStim' * pm.RF.values(:);
                    
                    
                    % Initialize timeSeries
                    % timeSeries = zeros(1,pm.timePointsN);
                    % for tt = 1:pm.timePointsN
                    %     % This is called the hadamard product.  It is the pointwise
                    %     % multiplication of the RF with the stimulus.  The hadProduct is
                    %     % the same size as the stimulus
                    %     hadProduct = stimValues(:,:,tt) .* pm.RF.values;
                    %     
                    %     % Now, we add up all of the hadProduct values
                    %     timeSeries(tt) = sum(hadProduct(:));
                    % end
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    % Convolution between the timeSeries and the HRF
                    % The HRF is defined in its own time frame
                    % hrfValues   = pm.HRF.values;
                    % hrftSteps   = pm.HRF.tSteps;
                    % hrfDuration = pm.HRF.Duration;
                    % Resample the HRF to n timepoints corresponding to the TR and duration
                    % pm.TR has to be multiple of 0.1
                    % hrfResampled= 0:pm.TR:hrfDuration-pm.TR;
                    
                    % hrf = getcanonicalhrf(pm.TR,pm.TR);
                    % hrf = hrf();
                    % pm.BOLD = conv2(timeSeries',hrf,'same');
                    
                    pm.BOLD = conv(pm.timeSeries',pm.HRF.values);
                    pm.BOLD = pm.BOLD(1:(end+1-length(pm.HRF.values)));
                    % pm.BOLD = [0,pm.BOLD(1:(end-length(pm.HRF.values)))];
                case 'CSS'
                    % Import function from analyzePRF and use it to generate the predicted BOLD signal.
                    %                     aprff = ...
                    %                         @(pp,dd) ...  % Defines the inputs to the function, params (to be guessed) and data
                    %                         conv2run(...  % Defines the main function, a convolution. We will want to guess the params
                    %                         posrect(pp(4)) * ...  % FIRST part of the convolution. It has form A * B. This is A
                    %                         (...                  % B starts here. Separate it as well
                    %                         dd * ...               % Data matrix
                    %                         [vflatten( ...         % Vector. vflatten just returns a vertical vector. Same as kk(:).
                    %                         placematrix( ...    % Substitutes matrix 2 into matrix 1 depending on matrix 3
                    %                         zeros(res), ...    % Matrix 1
                    %                         makegaussian2d(resmx,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / ... % Matrix 2: element1
                    %                         (2*pi*abs(pp(3))^2) ...                                                 % Matrix 2: element2
                    %                         ) ...   % Close placematrix
                    %                         )...          % Close vflatten
                    %                         ;0]...               % Adds a 0 to the flattened matrix (vector) at the end
                    %                         ).^posrect(pp(5)),...% End of B part B of A*B part of the convolution. This is the parameter in the exponential
                    %                         options.hrf,...      % SECOND part of the convolution, the HRF signal
                    %                         dd(:,prod(res)+1) ...% THIRD part of the conv, according knk. Separates convs between runs. This third part is basically a categorization. See help conv2run
                    %                         );
                    % Remove most of his additionnal functions to make it simpler
                    % Obtain all the parameters required to run it.
                    % pp are the params and dd is the stimulus
                    % When we created signal with RF of
                    %       x0=0,y0=0,theta=0,sigmaMinor=1,sigmaMajor=1
                    % the parameters returned with a R2=57 where:
                    %       50.538, 50.278, 0.12756, 412.16, 0.015133
                    pp(1) = 50;  % number of row location in a grid
                    pp(2) = 50;  % number of col location in a grid
                    pp(3) = 1; % Sigma minor and sigma mayor. only circles for now
                    pp(4) = 400;
                    pp(5) = 0.15;
                    res   = 100;
                    [xx,yy] = calcunitcoordinates(res);
                    RF = makegaussian2d(res,pp(1),pp(2),abs(pp(3)),abs(pp(3)),xx,yy,0,0) / (2*pi*pp(3)^2);
                    dd = spaceStim;
                    timeSeries = pp(4) * (dd' * RF(:)).^posrect(pp(5));
                    
                    pm.BOLD = conv2(timeSeries, getcanonicalhrf(pm.TR,pm.TR), 'same')'; % This is the same now
                    
                otherwise
                    error('Model %s not implemented, select basic or CSS', pm.Type)
            end
            % Scale the signal so that it has the required mean and contrast
            pm.BOLD  = pm.scaleAmpMean(pm.BOLD, pm.BOLDcontrast, pm.BOLDmeanValue);
            
        end
        function v = scaleAmpMean(pm, signal, contrast, meanBOLD)
            % Normalize to 0-1, so that min(normBOLD) == 0
            normBOLD    = (signal - min(signal))/(max(signal) - min(signal));
            % Calculate the max value, so that the relation between the
            % meanValue and the amplitude is the one set in pm.BOLDcontrast
            maxValue    = (contrast/100) * meanBOLD;
            % No obtain the scaled value with the following formula
            minValue    = 0;
            scaledBOLD  = minValue + normBOLD .* (maxValue - minValue);
            % Now 100 * (max(scaledBOLD) - min(scaledBOLD)) / pm.BOLDmeanValue
            % is equal to pm.BOLDcontrast. Now we need to add the correct mean
            % to the signal so that pm.BOLD == pm.BOLDmeanValue
            v           = (scaledBOLD - mean(scaledBOLD)) + meanBOLD;
            
            
            
            
            
%             % Normalize to 0-1, so that min(normBOLD) == 0
%             normBOLD    = (pm.BOLD - min(pm.BOLD))/(max(pm.BOLD) - min(pm.BOLD));
%             % Calculate the max value, so that the relation between the
%             % meanValue and the amplitude is the one set in pm.BOLDcontrast
%             maxValue    = (pm.BOLDcontrast/100) * pm.BOLDmeanValue;
%             % No obtain the scaled value with the following formula
%             minValue    = 0;
%             scaledBOLD  = minValue + normBOLD .* (maxValue - minValue);
%             % Now 100 * (max(scaledBOLD) - min(scaledBOLD)) / pm.BOLDmeanValue
%             % is equal to pm.BOLDcontrast. Now we need to add the correct mean
%             % to the signal so that pm.BOLD == pm.BOLDmeanValue
%             pm.BOLD     = (scaledBOLD - mean(scaledBOLD)) + pm.BOLDmeanValue;
        end
        % Compute synthetic BOLD with noise in top of the BOLD signal
        function compute(pm)
            % Computes the mean BOLD response and then adds noise.
            
            % First, compute the values in the required sub-classes
            pm.Stimulus.compute;
            pm.RF.compute;
            pm.HRF.compute;
            
            % Every sub-class has a computeBOLD function to compute the mean response.
            pm.computeBOLD;
            
            % Compute and add noise
            sumOfNoise = zeros(size(pm.BOLD));
            for ii=1:length(pm.Noise)
                pm.Noise{ii}.compute;
                sumOfNoise = sumOfNoise + pm.Noise{ii}.values;
            end
            
            pm.BOLDnoise = pm.BOLD + sumOfNoise;
        end
        % Plot it
        function plot(pm, varargin)
            % Read the inputs
            varargin = mrvParamFormat(varargin);
            p = inputParser;
            p.addRequired ('pm'  ,  @(x)(isa(x,'prfModel')));
            p.addParameter('what', 'nonoise', @ischar);
            p.addParameter('window', true, @islogical);
            p.parse(pm,varargin{:});
            what = mrvParamFormat(p.Results.what);
            w  = p.Results.window;
            
            switch what
                case 'nonoise'
                    if w;mrvNewGraphWin([pm.Type 'Synthetic BOLD signal (no noise)']);end
                    plot(pm.timePointsSeries, pm.BOLD);
                case 'withnoise'
                    if w;mrvNewGraphWin([pm.Type 'Synthetic BOLD signal (noise)']);end
                    grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');
                case 'both'
                    if w;mrvNewGraphWin([pm.Type 'Synthetic BOLD signals']);end
                    plot(pm.timePointsSeries, pm.BOLD); hold on;
                    plot(pm.timePointsSeries, pm.BOLDnoise);
                    legend({'No Noise BOLD','With Noise BOLD'})
                case 'all'
                    if w;mrvNewGraphWin([pm.Type 'Synthetic BOLD signals']);end
                    plot(pm.timePointsSeries, pm.scaleAmpMean(pm.timeSeries, pm.BOLDcontrast, pm.BOLDmeanValue)); hold on;
                    plot(pm.timePointsSeries, pm.BOLD);
                    plot(pm.timePointsSeries, pm.BOLDnoise);
                    legend({'Time Series','No Noise BOLD','With Noise BOLD'})                    
                case 'withnoisetimeseries'
                    if w;mrvNewGraphWin([pm.Type 'Synthetic BOLD signals']);end
                    plot(pm.timePointsSeries, pm.scaleAmpMean(pm.timeSeries, pm.BOLDcontrast, pm.BOLDmeanValue)); hold on;
                    plot(pm.timePointsSeries, pm.BOLDnoise);
                    legend({'Time Series','With Noise BOLD'})
                case 'nonoisetimeseries'
                    if w;mrvNewGraphWin([pm.Type 'Synthetic BOLD signals']);end
                    plot(pm.timePointsSeries, pm.scaleAmpMean(pm.timeSeries, pm.BOLDcontrast, pm.BOLDmeanValue)); hold on;
                    plot(pm.timePointsSeries, pm.BOLD);
                    legend({'Time Series','No Noise BOLD'})
                otherwise
                    error("'no noise', 'with noise', 'both', 'all', 'with noise timeseries', 'no noise timeseries' are acepted")
            end
            grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');
            aa = gca;
            text(5,aa.YLim(1)*(1.0025),...
                sprintf('TR=%1.1fs, Center(x0,y0)=[%1.1f,%1.1f]deg, Theta=%1.1frad, \\sigmaMaj=%1.1fdeg, \\sigmaMin=%1.1fdeg, FoVh=%1.1fdeg, FoVv=%1.1fdeg' ,...
                pm.TR,        pm.RF.Centerx0, pm.RF.Centery0, pm.RF.Theta,    pm.RF.sigmaMajor,     pm.RF.sigmaMinor,     pm.Stimulus.fieldofviewHorz, pm.Stimulus.fieldofviewVert))
            
        end
        
    end
    
end
