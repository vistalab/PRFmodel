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
       pm.plot;
    %}
    %{ 
       % Check mean BOLD, range, contrast
       pm = prfModel;
       pm.BOLDcontrast = 8;
       pm.Noise.seed = 12345
       pm.HRF.Type = 'boynton'
       pm.compute;
       pm.plot;
       noiselessRange = (max(pm.BOLD)-min(pm.BOLD))/2;  % should be 0.16
       assert(pm.BOLDcontrast/100 == noiselessRange)
       pm.signalPercentage = bold; % Give real signal with a mean of BOLDmeanValue
       pm.compute;
       noiselessRange = (max(pm.BOLD)-min(pm.BOLD))/2;  % should be 0.16
       assert(pm.BOLDcontrast/100 == noiselessRange/pm.BOLDmeanValue)
       pm.plot
    %}
    %{ 
        % Test if the contrast scaling thing is working or not
        pm = prfModel;
        pm.Noise.seed = 'none';
        pm.BOLDcontrast  = 10;
        pm.RF.Centerx0   = 10;
        pm.RF.Centery0   = 10;
        pm.RF.sigmaMajor = 2;
        pm.RF.sigmaMinor = 2;
        pm.BOLDcontrast
        pm.plot
        pm.scaleContrast = true;
        pm.compute
        pm.BOLDcontrast
        pm.plot
        % If we add noise, we see how in the scaled case, there is no signal, just noise
        pm = prfModel;
        pm.Noise.seed = 12345;
        pm.BOLDcontrast  = 10;
        pm.RF.Centerx0   = 10;
        pm.RF.Centery0   = 10;
        pm.RF.sigmaMajor = 2;
        pm.RF.sigmaMinor = 2;
        pm.BOLDcontrast
        pm.plot
        pm.scaleContrast = true;
        pm.compute
        pm.BOLDcontrast
        pm.plot
    
        close all;
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
        Type             ;
        % Components required to build the synthetic bold signal (classes)
        Stimulus         ;
        RF               ;
        HRF              ;
        Noise            ;
        % Other required options (double)
        signalPercentage ; % Provide results in BOLD signal (bold), in signal percent change (spc) or unitless (none) (default bold)
        BOLDmeanValue    ; % Required mean value of the synthetic BOLD signal (default 10000)
        BOLDcontrast     ; % Contrast of the synthetic BOLD signal, in % (default 8%)
        scaleContrast    ; % If we want to scale contrast to the max possible value
        timeSeries       ; % No scaling
        timeSeriesConv   ; % Time series and convolution, no scaling
        computeSubclasses; % Logical
        % BOLD signal value (before noise)
        BOLD             ;
        % The result: synthetic BOLD series (1 dim array of doubles)
        BOLDnoise        ; % Final value, composed of the BOLD + noise
    end
    properties (Dependent)
        TR               ; % This one is derived and copied to all other classes
        timePointsN      ; % Defined by the stimulus.
        timePointsSeries ; 
        frequencySeriesHz; 
        defaultsTable    ;
    end
    
    methods (Static)
        function d = defaultsGet
            % This provides the defaults of this class, which is the only one at
            % the top level. 
            d.TR                = 1;
            d.Type              = 'basic';
            d.signalPercentage  = 'bold';
            d.BOLDcontrast      = 5;    % Percent. So this will be 0.08
            d.scaleContrast     = false;    % Logical
            d.BOLDmeanValue     = 10000; % Mean BOLD, set signalPercentage to bold to use this
            d.computeSubclasses = true;
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
            p.addParameter('tr'               , d.TR                  , @isnumeric);
            p.addParameter('type'             , d.Type{:}             , @ischar);
            p.addParameter('signalpercentage' , d.signalPercentage{:} , @islogical);
            p.addParameter('boldmeanvalue'    , d.BOLDmeanValue       , @isnumeric);
            p.addParameter('boldcontrast'     , d.BOLDcontrast        , @isnumeric);
            p.addParameter('scalecontrast'    , d.scaleContrast       , @islogical);
            p.addParameter('computesubclasses', d.computeSubclasses   , @islogical);
            
            p.parse(varargin{:});
            % Assign defaults/parameters to class/variables
            pm.TR                = p.Results.tr;
            pm.Type              = p.Results.type;
            pm.signalPercentage  = p.Results.signalpercentage;
            pm.BOLDmeanValue     = p.Results.boldmeanvalue;
            pm.BOLDcontrast      = p.Results.boldcontrast;  % In percentage
            pm.scaleContrast     = p.Results.scalecontrast;
            pm.computeSubclasses = p.Results.computesubclasses;
            
            % Create the classes, and initialize a prfModel inside it
            pm.Stimulus          = pmStimulus(pm);
            pm.HRF               = pmHRF(pm); 
            pm.RF                = pmRF(pm);
            pm.Noise             = pmNoise(pm);
        end
        % Functions that apply the setting of main parameters to subclasses
        function set.TR(pm, tr)
            pm.uniqueTR      = tr;
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
        function v = get.frequencySeriesHz(pm)
            v  = (1/pm.TR)*(0:(pm.timePointsN-1))/pm.timePointsN;
        end
        function defaultsTable = get.defaultsTable(pm)
            % This function obtains all the defaults from all the
            % subclasses, so that we can construct a parameter table
            defaultsTable           = pm.defaultsGet;
            defaultsTable.HRF       = pm.HRF.defaultsGet;
            defaultsTable.RF        = pm.RF.defaultsGet;
            defaultsTable.Stimulus  = pm.Stimulus.defaultsGet;
            defaultsTable.Noise     = pm.Noise.defaultsGet;
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
            p.addParameter('showconv',false,@islogical);
            p.parse(pm,varargin{:});
            showconv = p.Results.showconv;
            
            % First, compute the values in the required sub-classes
            if pm.computeSubclasses
                pm.Stimulus.compute;
                pm.RF.compute;
                pm.HRF.compute;
            end            
            % Load stimulus
            stimValues = pm.Stimulus.getStimValues;
            
            % Initialize timeSeries, it is the signal prior to convolution
            [r,c,t]    = size(stimValues);
            spaceStim  = reshape(stimValues,r*c,t);
            switch pm.Type
                case 'basic'
                    % Calculate time series
                    pm.timeSeries = spaceStim' * pm.RF.values(:);
                    convValues    = conv(pm.timeSeries',pm.HRF.values);
                    % Create the bold signal with the correct size
                    % TODO: make all vectors columns whenever possible. Time vertical
                    % The conv is longer (HRF size -1), we need to cut the end
                    % of the conv, or paste the results in a correct sized vect
                    pm.timeSeriesConv = zeros(size(pm.timeSeries))';
                    pm.timeSeriesConv = convValues(1:length(pm.timeSeriesConv));
                    pm.BOLD = pm.timeSeriesConv;
                    if showconv
                        pm.showConvolution
                        hold on
                    end
                case 'CSS'
                    error('Untested, use only basic')
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
                    %{
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
                    %}
                otherwise
                    error('Model %s not implemented, select basic or CSS', pm.Type)
            end
            % Scale the signal so that it has the required mean and contrast
            
            % We need to know if the user wants to scale it to the max possible
            % signal first. We will edit the contrast if it is the case,
            % therefore the value of the set contrast will change
            if pm.scaleContrast
               maxcenter            = pm.Stimulus.maxStimCenter;
               % Create an pmRF class but with the same values coming from the main pm
               pmMax = pmRF(pm);
               % Change the center and rfsize, to calculate the max, maintain the rest
               pmMax.Type           = pm.RF.Type;
               pmMax.Centerx0       = maxcenter(1);
               pmMax.Centery0       = maxcenter(2);
               pmMax.sigmaMajor     = pm.RF.sigmaMajor;
               pmMax.sigmaMinor     = pm.RF.sigmaMinor;
               pmMax.Theta          = pm.RF.Theta;
               pmMax.dog_sigmaMajor = pm.RF.dog_sigmaMajor;
               pmMax.dog_sigmaMinor = pm.RF.dog_sigmaMinor;
               pmMax.dog_Theta      = pm.RF.dog_Theta;
               pmMax.dog_Scale      = pm.RF.dog_Scale;
               
               pmMax.compute;
               
               maxTimeSeries    = max(spaceStim' * pmMax.values(:));
               myMaxTimeSeries  = max(pm.timeSeries);
               % Scale the contrast value
               pm.BOLDcontrast = pm.BOLDcontrast * myMaxTimeSeries / maxTimeSeries;
            end
            % Convert the output is requested in pm.signalPercentage 
            switch pm.signalPercentage
                 case {'spc'}
                    % The last true means that we want to center in 0, it is nicer for visualization
                    pm.timeSeries = pm.unitless2contrast(pm.timeSeries, pm.BOLDcontrast,true);
                    pm.BOLD       = pm.unitless2contrast(pm.BOLD, pm.BOLDcontrast,true);
                case {'bold'}
                    pm.timeSeries = pm.contrast2BOLD(pm.timeSeries, pm.BOLDcontrast,true);
                    pm.BOLD       = pm.contrast2BOLD(pm.BOLD, pm.BOLDcontrast, pm.BOLDmeanValue);
                case {'none'}
                    % Do nothing, pm.BOLD remains unaltered in this step
                otherwise
                   error('%s provided, valid values are spc, none and bold', pm.signalPercentage)
            end
        end
        function showConvolution(pm)
            a_fig = mrvNewGraphWin('show convolution')
            % %%%%%%%%%%%%%%%%%%%%% Function: acnv.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %  Title: Animation of Graphical Convolution
            %
            %  To run this script just type 'acnv' on the MatLab prompt: > acnv
            %
            %  Description:
            %   1. This is a simple MatLab demo to animate the process of convolution.
            %      It is meant to help student to visualize how convolution works.
            %
            %   2. When this script is run, two function f(t) and go(t) are convolved
            %      and the output figure will show animated graphical convolution.
            %
            %   3. The functions "f" and "go" and their range of interval can be changed
            %      by editing this script at line numbers around "48 to 64"
            %
            %   4. Note:  For a better scaled plots of the functions f(t) and go(t1),
            %             it is recommended to set the functions such that their
            %             maximum value remains comparable. e.g one can use appropriate
            %             scaling. Other functions are also given 'commented out'
            %
            %             Interger values are recommended for the intervals
            %
            %   5. The animation can be made faster or slower by changing the value of
            %      the pause function in the animation loop. (around line number 134)
            %
            %  Author:
            %      Laine Berhane Kahsay
            %      Uni-Ulm, Germany
            %
            %   email: kahsay_2004@yahoo.com
            %
            %     ver: 1.0, written in Matlab 6.5/7.0
            %
            %  To see this help - type on the Matlab Prompt: > help acnv
            %
            % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % help acnv;
            % color of axis constant
            axis_color= [0.5 0.5 0.5];
            % sampling interval constant
            s_int = pm.TR;
            % interval for function 'f(t)'
            % t = [ -10:s_int:10 ];
            t = pm.timePointsSeries;
            % definition of function 'f(t)'
            % f = 0.1*(t.^2);
            f = pm.timeSeries';
            HRFcontrast = (max(pm.HRF.values)-min(pm.HRF.values))/2;
            f = 100 * pm.unitless2contrast(f,HRFcontrast,true);
            if f(1) >= 0, f = f - f(1);end
            if f(1) <= 0, f = f + abs(f(1));end
            %  f = 5*ones(1, length(t));
            %  f = t;
            % interval for function 'go(t1)'
            % t1 = [-7:s_int:7];
            t1 = pm.HRF.tSteps;
            
            
            % definition of function 'go(t1)'
            % go = -0.1*(t1.^2);
            go = pm.HRF.values;
            
            % go = .1*(t1.^3);
            % go = 5*cos(2*pi*t1);
            % go = 5*ones(1, length(t1));
            % go = zeros(1, length(t1));go(1)=5;
            % convolve: note the multiplation by the sampling interval
            c = s_int * conv(f, go);
            % Animation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % flip 'go(t1)' for the graphical convolutions g = go(-t1)
            g = fliplr(go);
            tf = fliplr(-t1);
            % slide range of 'g' to discard non-ovelapping areas with 'f' in the convolution
            tf = tf + ( min(t)-max(tf) );
            % get the range of function 'c' which is the convolution of 'f(t)' and 'go(t1)'
            tc = [ tf t(2:end)];
            tc = tc+max(t1);
            % start graphical output with three subplots
            % a_fig = figure;
            set(a_fig, 'Name', 'Animated Convolution', 'unit', 'pixel', ...
                'Position', [300, 150, 600, 750]);
            % plot f(t) and go(t1)
            ax_1 = subplot(3,1,1);
            op = plot(t,f, 'b',  t1, go, 'r');
            if isempty(pm.BOLD); pm.computeBOLD; end
            hold on; grid on;
            set(ax_1, 'XColor', axis_color, 'YColor', axis_color, 'Color', 'w', 'Fontsize', 9);
            xlim( [ ( min(t)-abs(max(tf)-min(tf)) - 1 ) ( max(t)+abs(max(tf)-min(tf)) + 1 ) ] );
            title('Graph of f(t) and go(t)', 'Color', axis_color );
            legend({'f(t)' 'go(t)'});
            % initialize animation the plot of 'g' is slided over the plot of 'f'
            % plot f in the subplot number 2
            ax_2 = subplot(3,1,2);
            p = plot(t, f);
            hold on; grid on;
            title('Graphical Convolution: f(t) and g = go(-t1)', 'Color', axis_color );
            
            % plot g in the subplot number 2
            q = plot(tf, g, 'r');
            xlim( [ ( min(t)-abs(max(tf)-min(tf))-1 ) ( max(t)+abs(max(tf)-min(tf))+1 ) ] );
            u_ym = get(ax_2, 'ylim');
            % plot two vertical lines to show the range of ovelapped area
            s_l = line( [min(t) min(t)], [u_ym(1) u_ym(2)], 'color', 'g'  );
            e_l = line( [min(t) min(t)], [u_ym(1) u_ym(2)], 'color', 'g'  );
            hold on; grid on;
            set(ax_2, 'XColor', axis_color, 'YColor', axis_color, 'Color', 'w', 'Fontsize', 9);
            % put a yellow shade on ovelapped region
            sg = rectangle('Position', [min(t) u_ym(1) 0.0001 u_ym(2)-u_ym(1)], ...
                'EdgeColor', [0 0 1]*0.7); %, 'FaceColor', [0 0 1]*0.7); %, 'EraseMode', 'xor');
            drawnow
            
            % initialize the plot the convolution result 'c'
            ax_3 = subplot(3,1,3);
            r = plot(tc, c);
            grid on; hold on;
            set(ax_3, 'XColor', axis_color, 'YColor', axis_color, 'Fontsize', 9);
            % xlim( [ min(tc)-1 max(tc)+1 ] );
            xlim( [ ( min(t)-abs(max(tf)-min(tf)) - 1 ) ( max(t)+abs(max(tf)-min(tf)) + 1 ) ] );
            title('Convolutional Product c(t). Black is our BOLD', 'Color', axis_color );
            % animation block
            for i=1:length(tc)
                
                % control speed of animation minimum is 0, the lower the faster
                pause(0.005);
                drawnow;
                
                % update the position of sliding function 'g', its handle is 'q'
                tf=tf+s_int;
                drawnow % set(q,'EraseMode','xor');
                
                set(q,'XData',tf,'YData',g);
                % show overlapping regions
                
                % show a vertical line for a left boundary of overlapping region
                sx = min( max( tf(1), min(t) ), max(t) );
                sx_a = [sx sx];
                drawnow % set(s_l,'EraseMode','xor');
                set(s_l, 'XData', sx_a);
                % show a second vetical line for the right boundary of overlapping region
                ex = min( tf(end), max(t) );
                ex_a = [ex ex];
                drawnow % set(e_l,'EraseMode','xor');
                set(e_l, 'XData', ex_a);
                
                % update shading on ovelapped region
                rpos = [sx u_ym(1) max(0.0001, ex-sx) u_ym(2)-u_ym(1)];
                set(sg, 'Position', rpos);
                
                % update the plot of convolutional product 'c', its handle is r
                drawnow % set(r,'EraseMode','xor');
                set(r,'XData',tc(1:i),'YData',c(1:i) );
                
            end;
            tccontrast = (max(c)-min(c))/2;
            ourBOLD = 100 * pm.unitless2contrast(pm.BOLD,tccontrast,true);
            if ourBOLD(1) >= 0, ourBOLD = ourBOLD - ourBOLD(1);end
            if ourBOLD(1) <= 0, ourBOLD = ourBOLD + abs(ourBOLD(1));end
            plot(ourBOLD,'k-.')
            
            %
            % end of acnv %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        end
        function v = contrast2BOLD(pm, signal, contrast, meanBOLD)
            % Normalize to 0-1, so that min(normBOLD) == 0
            normBOLD    = (signal - min(signal))/(max(signal) - min(signal));
            % Calculate the max value, so that the relation between the
            % meanValue and the amplitude is the one set in pm.BOLDcontrast
            maxValue    = (2*contrast/100) * meanBOLD;
            % No obtain the scaled value with the following formula
            minValue    = 0;
            scaledBOLD  = minValue + normBOLD .* (maxValue - minValue);
            % Now 100 * (max(scaledBOLD) - min(scaledBOLD)) / pm.BOLDmeanValue
            % is equal to pm.BOLDcontrast. Now we need to add the correct mean
            % to the signal so that pm.BOLD == pm.BOLDmeanValue
            v           = (scaledBOLD - mean(scaledBOLD)) + meanBOLD;
        end
        function v = unitless2contrast(pm, signal, contrast, centerzero)
            
            % TODO: the contrast we set, is the max contrast we would get when
            % the stimuli hits the center of the RF.
            % If the stimuli doesn't hit it, then the unitless max values that we are
            % goint to obtain from the matrix multiplication, are going to be
            % smaller. Right now, in the normalization process, it creates the
            % same amplitude for all. First aproximation to solve this problem
            % is that we are going to use the scale that it is coming from the
            % unitless timeseries. 
            
            
            
            
            
            
            
            % Normalize to 0-1, so that min(normBOLD) == 0
            normBOLD    = (signal - min(signal))/(max(signal) - min(signal));
            % Calculate the max value, so that the relation between the
            % meanValue and the amplitude is the one set in pm.BOLDcontrast
            maxValue    = (2*contrast/100);
            % Now obtain the scaled value with the following formula
            minValue    = 0;
            scaledBOLD  = minValue + normBOLD .* (maxValue - minValue);
            v           = (scaledBOLD - mean(scaledBOLD));
            if centerzero
                if v(1) >= 0, v = v - v(1);end
                if v(1) <= 0, v = v + abs(v(1));end
            end
            
        end
        % Compute synthetic BOLD with noise in top of the BOLD signal
        function compute(pm)
            % Computes the mean BOLD response and then adds noise.
            if pm.computeSubclasses
                % First, compute the values in the required sub-classes
                pm.Stimulus.compute;
                pm.RF.compute;
                pm.HRF.compute;
                pm.Noise.compute;
                
            end
            % Compute BOLD signal
            pm.computeBOLD;
            % Add the noise component. We want them to be separated. 
            % Depending on how the pm.BOLD is going, we need to do a separate thing with the noise
            switch pm.signalPercentage
                case {'spc'} 
                    % Both are in the same contrast scale, with zero mean
                    % We can add noise directly
                    pm.BOLDnoise   = pm.BOLD + pm.Noise.values;
                case {'bold'}
                    % De-scale the BOLD to contrast, add noise, and re-scale it back
                    signal         = (pm.BOLD - mean(pm.BOLD)) ./ mean(pm.BOLD);
                    signalAndNoise = signal + pm.Noise.values;
                    pm.BOLDnoise   = pm.contrast2BOLD(signalAndNoise, ...
                                             pm.BOLDcontrast, pm.BOLDmeanValue); 
               case {'none'}
                   % Now this is the trickiest one. I need to scale the noise to the unitless arbitrary values coming from the convolution.
                   % Maybe not the most elegant, but I'll use the same approach as above
                   
                   % Calculate the "contrast" of the upcoming unitless signal,
                   % to revert it back. Min is zero. Contrast is range / 2.
                   
                   % The noise will be added in the contrast we set up in
                   % pm.BOLDcontrast, which is 5% by default, but we can change
                   % it
                   signal = pm.unitless2contrast(pm.BOLD, pm.BOLDcontrast, true); % it is fractional, not in %, but...
                   % Now, we can add noise as in the other cases, because it has
                   % the same scale. 
                   signalAndNoise = signal + pm.Noise.values; 
                   % And, now scale it back to whatever scale the original
                   % unitless signal had, but it will have noise. 
                   
                   
                   % Normalize to 0-1, so that min(normBOLD) == 0
                    normBOLD    = (signalAndNoise)/(max(signalAndNoise));
                    % Make the max value the same as the orig signal
                    maxValue    = max(pm.BOLD); % This will change in every case
                    % Now obtain the scaled value
                    v  = normBOLD .* maxValue;
                    % Make it start in zero
                    if v(1) >= 0, v = v - v(1);end
                    if v(1) <= 0, v = v + abs(v(1));end
                   
                    pm.BOLDnoise = v;
                   
               otherwise
                   error('%s provided, valid values are spc, none and bold', pm.signalPercentage)
            end
            
        end
        % Plot it
        function plot(pm, varargin)
            % Read the inputs
            varargin = mrvParamFormat(varargin);
            p = inputParser;
            p.addRequired ('pm'  ,  @(x)(isa(x,'prfModel')));
            p.addParameter('what', 'both', @ischar);
            p.addParameter('window', true, @islogical);
            p.addParameter('addtext', false, @islogical);
            p.addParameter('color', 'b');
            p.addParameter('centerzero', false);
            
            p.parse(pm,varargin{:});
            what = mrvParamFormat(p.Results.what);
            w  = p.Results.window;
            t  = p.Results.addtext;
            c  = p.Results.color;
            z  = p.Results.centerzero;
            
            switch what
                case 'timeseries'
                    pm.computeBOLD
                    if w;mrvNewGraphWin([pm.Type 'timeseries']);end
                    % plot(pm.timePointsSeries, pm.unitless2contrast(pm.timeSeries, ...
                    %      pm.BOLDcontrast,true),'-o','color',c,'LineWidth',2);
                    plot(pm.timePointsSeries, pm.timeSeries, ...
                         '-o','color','b','LineWidth',2);
                    legend({'Time Series'})
                case 'nonoisetimeseries'
                    pm.computeBOLD
                    if w;mrvNewGraphWin([pm.Type 'Synthetic BOLD signals']);end
                    % plot(pm.timePointsSeries, pm.unitless2contrast(pm.timeSeries, ...
                    %      pm.BOLDcontrast,true),'--','color','k','LineWidth',1); hold on;
                    plot(pm.timePointsSeries, pm.timeSeries,'--','color','k',...
                                                        'LineWidth',1); hold on;
                    plot(pm.timePointsSeries, pm.BOLD,'color','b');
                    legend({'Time Series','No Noise BOLD'})
                case {'nonoise','noiseless','noisefree'}
                    pm.computeBOLD
                    if w;mrvNewGraphWin([pm.Type 'Synthetic BOLD signal (no noise)']);end
                    plot(pm.timePointsSeries, pm.BOLD,'color',c);
                    grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');
                case 'withnoise'
                    pm.compute;
                    if w;mrvNewGraphWin([pm.Type 'Synthetic BOLD signal (noise)']);end
                    plot(pm.timePointsSeries, pm.BOLDnoise,'color',c);
                    grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');
                  case 'both'
                    pm.compute;
                    if w;mrvNewGraphWin([pm.Type 'Synthetic BOLD signals']);end
                    plot(pm.timePointsSeries, pm.BOLD,'--','color','k'); hold on;
                    plot(pm.timePointsSeries, pm.BOLDnoise,'color',c);
                    legend({'No Noise BOLD','With Noise BOLD'})
                case 'all'
                    pm.compute;
                    if w;mrvNewGraphWin([pm.Type 'Synthetic BOLD signals']);end
                    plot(pm.timePointsSeries, pm.unitless2contrast(pm.timeSeries, pm.BOLDcontrast,true),'color','k'); hold on;
                    plot(pm.timePointsSeries, pm.BOLD);
                    plot(pm.timePointsSeries, pm.BOLDnoise);
                    legend({'Time Series','No Noise BOLD','With Noise BOLD'})                    
                case 'withnoisetimeseries'
                    pm.compute;
                    if w;mrvNewGraphWin([pm.Type 'Synthetic BOLD signals']);end
                    plot(pm.timePointsSeries, pm.unitless2contrast(pm.timeSeries, pm.BOLDcontrast,true),'color','k'); hold on;
                    plot(pm.timePointsSeries, pm.BOLDnoise);
                    legend({'Time Series','With Noise BOLD'})
                otherwise
                    error('no noise, with noise, both, all, with noise timeseries, no noise timeseries are acepted')
            end
            grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');
            aa = gca;
            if t, text(5,aa.YLim(1)*(1.0025),...
                    sprintf('TR=%1.1fs, Center(x0,y0)=[%1.1f,%1.1f]deg, Theta=%1.1frad, \\sigmaMaj=%1.1fdeg, \\sigmaMin=%1.1fdeg, FoVh=%1.1fdeg, FoVv=%1.1fdeg' ,...
                    pm.TR,        pm.RF.Centerx0, pm.RF.Centery0, pm.RF.Theta,    pm.RF.sigmaMajor,     pm.RF.sigmaMinor,     pm.Stimulus.fieldofviewHorz, pm.Stimulus.fieldofviewVert))
            end
        end
    end
end
