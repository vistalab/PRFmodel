classdef pmNoise <  matlab.mixin.SetGet & matlab.mixin.Copyable
    % This is a superclass for Noise-s. 
    %
    % The noise returned has zero mean and an amplitude that makes sense
    % for the mean level of the simulated BOLD signal.  We are considering
    % returning a unit contrast for the noise and scaling its amplitude. At
    % this time, The mean BOLD is specified by a parameter (BOLDmeanValue).
    % The range of the BOLD signal is specified by a parameter
    % (BOLDcontrast).  These set the level and range of the BOLD signal.
    % The amplitude of the noise is related to these two values.
    %
    % To adjust the contrast of the noise, we set the parameter noiseMult
    % and the parameter noise2signal.
    %
    % Every particular instance of this class
    % will have different parameters, so it will be a children class. For
    % example:
    %   - White noise (white)
    %
    %
    %   - Motion related (translation and rotation)
    %
    %   - Cardiac Related
    %
    %   - Respiration related
    %
    %   - Low frequency physiological fluctuations
    %
    %   - Draining veins
    %
    %   - Low frequency drifts
    %     (slow head displacements, scanner related (e.g. heating...)
    %
    %   - Hardware related instabilities
    %
    % TO-DO:
    %    %   - Eye motion jitter
    %     eyeMotionJitter = 1;  % Deg
    % Syntax:
    %      noise = Noise;
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
    % At the end of this function we have the calculations in normal subjects to
    % have an idea of how a representative noise looks
    
    % Examples
    %{
      pm = prfModel;
      pm.Noise.compute;
      pm.Noise.plot
    
      pm.BOLDmeanValue=100;
      pm.Noise.compute;
      pm.Noise.plot
    
    
    %}
    %{
      pm = prfModel;
      pm.Noise.seed = 'none';
      pm.Noise.compute;
      pm.Noise.plot
    
    %}
    %{
      pm = prfModel;
    pm.Noise.seed=1;
      pm.Noise.compute;
      pm.Noise.plot
    
      pm.BOLDmeanValue=100;
    pm.Noise.seed=1;
      pm.Noise.compute;
      pm.Noise.plot
    %}
        %{
      pm = prfModel;
    pm.Noise.seed=1;
      pm.Noise.compute;
      pm.Noise.plot
    
      pm.BOLDmeanValue=100;
    pm.Noise.seed=2;
      pm.Noise.compute;
      pm.Noise.plot
    %}
   %{
      pm = prfModel;
    pm.Noise.lowfrequ_amplitude=0;
    
      pm.Noise.compute;
      pm.Noise.plot
    %}
     %{
      pm = prfModel;
      pm.Noise.compute;
      F = abs(fft(pm.Noise.values'));
      mrvNewGraphWin;plot(F);
      
    %}
    %{
      pm    = prfModel;
      pm.TR = 1.5;
      pm.Noise.plot;
    
      pm.Noise.white_amplitude=3;

      pm.Noise.plot;
    %}
    %{
      pm = prfModel;
      pm.Noise.white_amplitude  = 0;
      pm.Noise.respiratory_amplitude  = 0.5;
      pm.Noise.respiratory_frequency  = 0.2;
      pm.Noise.plot;
      pm.TR=2;
      pm.Noise.plot;
    %}
    %{
      pm = prfModel;
      pm.Noise.white_amplitude  = 0;
      pm.Noise.respiratory_amplitude  = 0.5;
      pm.Noise.respiratory_frequency  = 0.2;
      pm.Noise.cardiac_amplitude  = 0.3;
      pm.Noise.cardiac_frequency  = 1.17;
      pm.Noise.plot;
      pm.TR=2;
      pm.Noise.plot;
    %}

    %{
      pm = prfModel;
      pm.Noise.white_amplitude  = 0;
      pm.Noise.respiratory_amplitude  = 0;
      pm.Noise.respiratory_frequency  = 0.2;
      pm.Noise.cardiac_amplitude  = 0;
      pm.Noise.cardiac_frequency  = 1.17;
      pm.Noise.lowfrequ_amplitude = 0.25;
      pm.Noise.lowfrequ_frequ     = 120; % 1/120=0.0083
      pm.Noise.plot;
      pm.TR=2;
      pm.Noise.plot;
    %}    
    properties
        PM;                 % prfModel that has some of the variables we need, such as TR
        seed;               % 'none' for no noise, 'random' or numeric(seed) otherwise
        jitter;             % 0-1 numeric 2 vector, [freq_jitter, amplitude_jitter], [0 0] for no jitter
        white_amplitude;    % amplitude of white noise, the stdev of the gaussian, 0 for no white noise
        cardiac_amplitude;  % amplitude of cardiac noise, 0 for no cardiac noise
        cardiac_frequency;  % cardiac frequ in Hz
        respiratory_amplitude; % amplitude of resp. noise, 0 for no resp. noise
        respiratory_frequency; % resp. frequ. in Hz
        lowfrequ_frequ;     % Seconds
        lowfrequ_amplitude; % amplitude for lowfrequ noise, 0 for no lowfrequ noise
        values;
    end
    
    properties(Dependent= true, SetAccess = private, GetAccess = public)
        TR;            % Seconds, it will be read from the parent class pm
    end
    
    
    
    %%
    methods (Static)
        function d = defaultsGet(varargin)
            varargin = mrvParamFormat(varargin);
            p = inputParser;
            p.addParameter('voxel','mid', @ischar);
            p.parse(varargin{:});
            d.voxel = p.Results.voxel;
            
            d.seed   = 'random';
            d.jitter = [0,0];
            
            switch strrep(lower(d.voxel),' ','')
                case {'mid','midnoise'}
                    % White Noise
                    d.white_amplitude    = 0.032; % white noise constant to relate to random number generator
                    % Cardiac
                    d.cardiac_amplitude  = 0.01;
                    d.cardiac_frequency  = 1.05;
                    % Respiratory
                    d.respiratory_amplitude  = 0.01;
                    d.respiratory_frequency  = 0.3;
                    % Low frequency drift
                    d.lowfrequ_amplitude     = 0.01;
                    d.lowfrequ_frequ         = 120;  % Seconds
                case {'good','low','lownoise'}
                    % White Noise
                    d.white_amplitude    = 0.016; % white noise constant to relate to random number generator
                    % Cardiac
                    d.cardiac_amplitude  = 0.004;
                    d.cardiac_frequency  = 1.05;
                    % Respiratory
                    d.respiratory_amplitude  = 0.004;
                    d.respiratory_frequency  = 0.28;
                    % Low frequency drift
                    d.lowfrequ_amplitude     = 0.004;
                    d.lowfrequ_frequ         = 120;  % Seconds
                case {'bad','high','highnoise'}
                    % White Noise
                    d.white_amplitude    = 0.05; % white noise constant to relate to random number generator
                    % Cardiac
                    d.cardiac_amplitude  = 0.025;
                    d.cardiac_frequency  = 1.055;
                    % Respiratory
                    d.respiratory_amplitude  = 0.01;
                    d.respiratory_frequency  = 0.3;
                    % Low frequency drift
                    d.lowfrequ_amplitude     = 0.015;
                    d.lowfrequ_frequ         = 120;  % Seconds
                otherwise
                    error('Voxel type %s not recognized',d.voxel)
            end
            d = struct2table(d,'AsArray',true);
        end
        
    end
    methods
        % Constructor
        
        function noise = pmNoise(pm,varargin)
            % Obtain defaults table. If a parameters is not passed, it will use
            % the default one defined in the static function
            d = noise.defaultsGet;
            % Read the inputs
            varargin = mrvParamFormat(varargin);
            p = inputParser;
            p.addRequired('pm'                                , @(x)(isa(x,'prfModel')));
            p.addParameter('seed'                  ,d.seed                          );
            p.addParameter('jitter'                ,d.jitter               , @isnumeric);
            p.addParameter('white_amplitude'       ,d.white_amplitude      , @isnumeric);
            p.addParameter('cardiac_amplitude'     ,d.cardiac_amplitude    , @isnumeric);
            p.addParameter('cardiac_frequency'     ,d.cardiac_frequency    , @isnumeric);
            p.addParameter('respiratory_amplitude' ,d.respiratory_amplitude, @isnumeric);
            p.addParameter('respiratory_frequency' ,d.respiratory_frequency, @isnumeric);
            p.addParameter('lowfrequ_frequ'        ,d.lowfrequ_frequ       , @isnumeric);
            p.addParameter('lowfrequ_amplitude'    ,d.lowfrequ_amplitude   , @isnumeric);
            p.parse(pm,varargin{:});
            
            % Initialize the pm model and hrf model parameters
            noise.PM                     = pm;
            % params
            noise.seed                   = p.Results.seed;
            noise.jitter                 = p.Results.jitter;
            noise.white_amplitude        = p.Results.white_amplitude;
            noise.cardiac_amplitude      = p.Results.cardiac_amplitude;
            noise.cardiac_frequency      = p.Results.cardiac_frequency;
            noise.respiratory_amplitude  = p.Results.respiratory_amplitude;
            noise.respiratory_frequency  = p.Results.respiratory_frequency;
            noise.lowfrequ_frequ         = p.Results.lowfrequ_frequ;
            noise.lowfrequ_amplitude     = p.Results.lowfrequ_amplitude;
            
            
            % Check that jitter is between 0 and 1
            validateattributes(noise.jitter(1), {'numeric'},{'>=',0,'<=',1})
            validateattributes(noise.jitter(2), {'numeric'},{'>=',0})
        end
        
        function noise = setVoxelDefaults(noise,voxelType)
            % Obtain defaults table. If a parameters is not passed, it will use
            % the default one defined in the static function
            d = table2struct(noise.defaultsGet('voxel',voxelType));
            
            noise.white_amplitude       = d.white_amplitude;
            noise.cardiac_amplitude     = d.cardiac_amplitude;
            noise.cardiac_frequency     = d.cardiac_frequency;
            noise.respiratory_amplitude = d.respiratory_amplitude;
            noise.respiratory_frequency = d.respiratory_frequency;
            noise.lowfrequ_frequ        = d.lowfrequ_frequ;
            noise.lowfrequ_amplitude    = d.lowfrequ_amplitude;
        end
        
        % Methods available to this class and his childrens (friston, boynton... classes)
        function TR = get.TR(noise)
            TR = noise.PM.TR;
        end
        % COMPUTE
        function compute(noise)
            if iscell(noise.seed)
                noise.seed = noise.seed{:};
            end
            if length(noise.jitter)==1
                noise.jitter = [noise.jitter,noise.jitter]; 
                warning('Only one value for jitter provided, used it for amplitude and frequency')
            end
            if length(noise.jitter)~=2
                error('pm.Noise.jitter needs to be numeric of length 2')
            end
            % Set a seed so that the rands generates the same sequence
            switch noise.seed
                case {'none','nonoise'}
                    addnoise = false;
                case 'random'
                    % Initialize with a random seed
                    rng('shuffle','twister')
                    addnoise = true;
                otherwise
                    addnoise = true;
                    % Use this as the rng seed
                    if isnumeric(noise.seed)
                        % rng(noise.seed);
                        rng(noise.seed,'twister')
                    else
                        error('Unknown noise status %s\n',noise.seed);
                    end
            end
            
            tmpNoise     = {};
            % WHITE
            tmpNoise{1}  = {};
            if noise.white_amplitude > 0 && addnoise
                aNoise      = noise.white_amplitude;
                tmpNoise{1} = aNoise * randn(1,noise.PM.timePointsN);
            end
                  
            % CARDIAC
            tmpNoise{2}  = {};
            if noise.cardiac_amplitude > 0 && addnoise
                % Jitter represents the standard deviation of the Gaussian
                % noise added to the parameter, as a fraction of the parameter 
                % value.  A value of 1 means that
                % the sd is equal to the parameter value itself.  
                
                % Frequency [Hz]
                fNoise = noise.cardiac_frequency*(1 + noise.jitter(1)*randn(1,1));
                % Amplitude
                aNoise = noise.cardiac_amplitude*(1 + noise.jitter(2)*randn(1,1));
                % Time points
                t      = noise.PM.timePointsSeries;
                % Calculate the noise
                tmpNoise{2} = aNoise * sin(2 * pi .* t .* fNoise);                
            end
            
            % RESPIRATORY
            tmpNoise{3}  = {};
            if noise.respiratory_amplitude > 0 && addnoise
                % See cardiac for the jitter explanation
                
                % Frequency [Hz]
                fNoise = noise.respiratory_frequency*(1 + noise.jitter(1)*randn(1,1));
                % Amplitude
                aNoise = noise.respiratory_amplitude * (1 + noise.jitter(2)*randn(1,1));
                % Time points
                t      = noise.PM.timePointsSeries;
                % Calculate the noise
                tmpNoise{3} = aNoise * sin(2 * pi .* t .* fNoise);                
            end
            
            % LOW FREQU DRIFT
            tmpNoise{4}  = {};
            if noise.lowfrequ_amplitude && addnoise
                % This is the R function used in neuRosim
                % CALL: n.low <- lowfreqdrift(dim = 1, nscan = 100, TR = 2, freq = 120) 
                % RESULTS: R_result = [ 4.240198436, 4.220685725, 4.181794990, 4.123794530, 4.047084113, 3.952191782, 3.839769644, 3.710588669, 3.565532546, 3.405590654, 3.231850196, 3.045487573, 2.847759065, 2.639990899, 2.423568788, 2.199927027, 1.970537244, 1.736896889, 1.500517574, 1.262913341, 1.025588978, 0.790028456, 0.557683615, 0.329963153, 0.108222047,-0.106248532,-0.312230718,-0.508588788,-0.694277263,-0.868348215,-1.029957718,-1.178371387,-1.312968979,-1.433248000,-1.538826310,-1.629443706,-1.704962466,-1.765366865,-1.810761670,-1.841369635,-1.857528016,-1.859684152,-1.848390149,-1.824296721,-1.788146247,-1.740765108,-1.683055375,-1.615985929,-1.540583084,-1.457920814,-1.369110653,-1.275291381,-1.177618564,-1.077254066,-0.975355600,-0.873066441,-0.771505356,-0.671756872,-0.574861945,-0.481809115,-0.393526235,-0.310872826,-0.234633135,-0.165509952,-0.104119229,-0.050985562,-0.006538548, 0.028889935, 0.055067521, 0.071862196, 0.079242204, 0.077275002, 0.066125304, 0.046052212, 0.017405489,-0.019379004,-0.063784647,-0.115220131,-0.173026195,-0.236482966,-0.304817820,-0.377213707,-0.452817845,-0.530750696,-0.610115155,-0.690005853,-0.769518482,-0.847759065,-0.923853077,-0.996954326,-1.066253523,-1.130986450,-1.190441648,-1.243967567,-1.290979094,-1.330963404,-1.363485089,-1.388190496,-1.404811260,-1.413166969];
                % Matlab implementation
                % function spm_drift has been implemented as a function file
                % nscan = 100; TR=2; freq = 120;
                % n     = floor(2 * (nscan * TR) / freq + 1);
                % isclose(tmpNoise{4},R_result','tolerance',0.000000001)
                
                % Jitter: see cardiac and respiratory
                fNoise = noise.lowfrequ_frequ * (1 + noise.jitter(1)*randn(1,1));
                aNoise = noise.lowfrequ_amplitude * (1 + noise.jitter(2)*randn(1,1));
                
                n = floor(2 * (noise.PM.timePointsN * noise.PM.TR)/fNoise + 1);
                if (n < 3) 
                    error("Number of basis functions is too low. Lower frequency or longer scanning time should be provided.")
                end
                spmdrift    = spm_drift(noise.PM.timePointsN, n);
                driftbase   = spmdrift(:,2:end);  % Eliminate the DC term
                driftbase   = sum(driftbase,2)'; 
                % Now we need to scale this noise to the actual noise values
                tmpNoise{4} = aNoise * driftbase;
            end
            
            
            % SUM ALL NOISE MODEL
            vals = zeros([1,noise.PM.timePointsN]);
            for nn=1:length(tmpNoise)
                if ~isempty(tmpNoise{nn})
                    vals = vals + tmpNoise{nn};
                end
            end
            
            % RETURN VALUES
            noise.values = vals;

        end
        
        % Plot it
        function plot(noise,varargin)
            % Read the inputs
            varargin = mrvParamFormat(varargin);
            p = inputParser;
            p.addRequired ('noise'  ,  @(x)(isa(x,'pmNoise')));
            p.addParameter('window',true, @islogical);
            p.addParameter('windowfreq',true, @islogical);
            p.addParameter('color','b');
            p.addParameter('line','-');
            
            p.parse(noise,varargin{:});
            w = p.Results.window;
            wf = p.Results.windowfreq;
            c = p.Results.color;
            line = p.Results.line;
            
            if w; mrvNewGraphWin('Noise'); end
            noise.compute;
            % Plot it
            if wf; subplot(2,1,1); end
            plot(noise.PM.timePointsSeries, noise.values,'-','color',c,'LineStyle',line);hold on;
            grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');
            title(['Noise: time domain (not scaled to BOLD), TR=' num2str(noise.TR)])
            % If it is not even we cannot plot the half of it, so remove the last time point
            timePoints = noise.PM.timePointsN;
            timeValues = noise.values;
            if ~iseven(timePoints)
                timePoints = timePoints - 1;
                timeValues = timeValues(1:timePoints);
            end

            
            % Only plot frequ domain if in our own window
            if wf
                subplot(2,1,2)
                % Obtain the FFT of the signal
                P = abs(fft(timeValues)/timePoints);
                % Halve it
                % {
                    F = P(1:timePoints/2);
                    % The amplitude is divided in the two halfs, so we take one half and
                    % multiply the amplitude by two
                    F(2:end-1) = 2*F(2:end-1);
                    % Get the frequency vector (for the half points)
                    f = (1/noise.PM.TR)*(0:(timePoints/2)-1)./timePoints;
                %}

                % Do not halve it (but double it for the amplitude)
                %{
                    F = P;
                    F(2:end-1) = 2*F(2:end-1);
                    % Obtain the frequency vector
                    f = (1/noise.PM.TR)*(1:(noise.PM.timePointsN))/noise.PM.timePointsN;
                %}
                plot(f, F);hold on;
                grid on; xlabel('f[Hz]'); ylabel('Relative amplitude');
                title(['Noise: frequ domain (not scaled to BOLD), TR=' num2str(noise.TR)])
            end
        end
    end
    
end

    %{
      pm = prfModel;
      pm.TR = 1;
      pm.Noise.white_amplitude  = 0.0;
      pm.Noise.respiratory_amplitude  = 1;
      pm.Noise.respiratory_frequency  = 0.2;
      pm.Noise.cardiac_amplitude  = 1;
      pm.Noise.cardiac_frequency  = 1.2;
      pm.Noise.lowfrequ_amplitude = 0;
      pm.Noise.lowfrequ_frequ     = 120; % 1/120=0.0083
      pm.Noise.plot;
      
    %}
    %{ 
       % Check jitter
       pm = prfModel;
       pm.TR=1.5;
       pm.BOLDcontrast = 16;
       % Same seed, jitter should be the same every repetition
       pm.Noise.seed = 12346
       % Make all zero except respiratory
       pm.Noise.white_amplitude=0;
       pm.Noise.lowfrequ_amplitude=0;
       pm.Noise.cardiac_amplitude=0;
       % Set desired respiratory
       pm.Noise.respiratory_amplitude=1;
       pm.Noise.respiratory_frequency=0.2;
       % Plot
       pm.Noise.plot;title('no jitter')
       pm.Noise.jitter = 0.1;
       pm.Noise.plot;title('0.1 jitter')
       pm.Noise.jitter = 0.5;
       pm.Noise.plot;title('0.5 jitter')
       % close all
    %}
    %{ 
       % Check jitter
       pm = prfModel;
       pm.TR=1.5;
       pm.BOLDcontrast = 16;
       % Same seed, jitter should be the same every repetition
       pm.Noise.seed = 12345
       % Make all zero except lowfrequ
       pm.Noise.white_amplitude=0;
       pm.Noise.respiratory_amplitude=0;
       pm.Noise.cardiac_amplitude=0;
       % Plot
       pm.Noise.plot;title('no jitter')
       pm.Noise.jitter = 0.1;
       pm.Noise.plot;title('0.1 jitter')
       pm.Noise.jitter = 0.4;
       pm.Noise.plot;title('0.5 jitter')
       % close all
    %}

    


