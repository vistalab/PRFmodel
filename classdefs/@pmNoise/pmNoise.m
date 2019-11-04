classdef pmNoise <  matlab.mixin.SetGet & matlab.mixin.Copyable
    % This is a superclass for Noise-s. Every particular instance of this class
    % will have different parameters, so it will be a children class. For
    % example:
    %   - White noise (white)
    %
    %   - Eye motion jitter
    %     eyeMotionJitter = 1;  % Deg
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
       
    %}
    
    properties
        PM;           % prfModel that has some of the variables we need, such as TR
        seed;         % 'none' for no noise, 'random' or numeric(seed) otherwise
        jitter;       % 0-1 numeric, 0 for no jitter
        noisemult;    % noise amplitude multiplier
        white_noise2signal; % amplitude of white noise, 0 for no white noise
        cardiac_amplitude;  % amplitude of cardiac noise, 0 for no cardiac noise
        cardiac_frequency;  % cardiac frequ in Hz
        respiratory_amplitude; % amplitude of resp. noise, 0 for no resp. noise
        respiratory_frequency; % resp. frequ. in Hz
        lowfrequ_frequ;  % Seconds
        lowfrequ_amplitude;  % amplitude for lowfrequ noise, 0 for no lowfrequ noise
        values;
    end
    
    properties(Dependent= true, SetAccess = private, GetAccess = public)
        TR;            % Seconds, it will be read from the parent class pm
    end
    
    
    
    %%
    methods (Static)
        function d = defaultsGet
            d.seed               = 'random';
            d.jitter             = 0;
            d.noisemult          = 1;
            
            % White Noise
            d.white_noise2signal = 0.022; % white noise; constant to relate to signal
            
            % Cardiac
            d.cardiac_amplitude  = 0.024;
            d.cardiac_frequency  = 1.17; 
            
            % Respiratory
            d.respiratory_amplitude  = 0.024;
            d.respiratory_frequency  = 0.23; 
            
            % Low frequency drift
            d.lowfrequ_frequ         = 120;  % Seconds
            d.lowfrequ_amplitude     = 0.12;  % Seconds
            
            
            % Convert to table and return
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
            p.addParameter('noisemult'             ,d.noisemult            , @isnumeric);
            p.addParameter('white_noise2signal'    ,d.white_noise2signal   , @isnumeric);
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
            noise.noisemult              = p.Results.noisemult;
            noise.white_noise2signal     = p.Results.white_noise2signal;
            noise.cardiac_amplitude      = p.Results.cardiac_amplitude;
            noise.cardiac_frequency      = p.Results.cardiac_frequency;
            noise.respiratory_amplitude  = p.Results.respiratory_amplitude;
            noise.respiratory_frequency  = p.Results.respiratory_frequency;
            noise.lowfrequ_frequ         = p.Results.lowfrequ_frequ;
            noise.lowfrequ_amplitude     = p.Results.lowfrequ_amplitude;
            
            % Check that jitter is between 0 and 1
            validateattributes(noise.jitter, {'numeric'},{'>=',0,'<=',1})
        end
        
        % Methods available to this class and his childrens (friston, boynton... classes)
        function TR = get.TR(noise)
            TR = noise.PM.TR;
        end
        
        % COMPUTE
        function compute(noise)
            % Simplifying this code. Noise is a simple class, and the noise
            % models will be specified in the parameters. 
            if iscell(noise.seed)
                noise.seed = noise.seed{:};
            end
            % Set a seed so that the rands generates the same sequence
            switch noise.seed
                case {'none','nonoise'}
                    % No noise case
                    % noise.white_noise2signal     = 0;
                    % noise.cardiac_amplitude      = 0;
                    % noise.respiratory_amplitude  = 0;
                    % noise.lowfrequ_amplitude     = 0;
                    addnoise = false;
                case 'random'
                    % Initialize with a random seed, say time
                    rng(now)
                    addnoise = true;
                otherwise
                    addnoise = true;
                    % Use this as the rng seed
                    if isnumeric(noise.seed)
                        rng(noise.seed);
                    else
                        error('Unknown noise status %s\n',noise.seed);
                    end
            end
            
            
            % Obtain the BOLD signal without noise. 
            noise.PM.computeBOLD;
            signal       = noise.PM.BOLD; % Retrieved from the parent model
            tmpNoise     = {};
            
            % WHITE
            tmpNoise{1}  = {};
            if noise.white_noise2signal > 0 && addnoise
                aNoise      = noise.white_noise2signal * noise.noisemult;
                aNoise      = aNoise * (max(signal) - min(signal));
                tmpNoise{1} = aNoise * randn([1,noise.PM.timePointsN]);
            end
                  
            % CARDIAC
            tmpNoise{2}  = {};
            if noise.cardiac_amplitude > 0 && addnoise
                % Jitter represents the standard deviation of the Gaussian
                % noise added to the parameter, as a fraction of the parameter 
                % value.  A value of 1 means that
                % the sd is equal to the parameter value itself.  
                
                % Frequency [Hz]
                fNoise = noise.cardiac_frequency*(1 + noise.jitter*randn(1,1));
                % Amplitude
                aNoise = noise.cardiac_amplitude*(1 + noise.jitter*randn(1,1));
                aNoise = aNoise * noise.noisemult;
                aNoise = aNoise * (max(noise.PM.BOLD) - min(noise.PM.BOLD));
                
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
                fNoise = noise.respiratory_frequency*(1 + noise.jitter*randn(1,1));
                % Amplitude
                aNoise = noise.respiratory_amplitude * (1 + noise.jitter*randn(1,1));
                aNoise = aNoise * noise.noisemult;
                aNoise = aNoise * (max(noise.PM.BOLD) - min(noise.PM.BOLD));
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
                fNoise = noise.lowfrequ_frequ * (1 + noise.jitter*randn(1,1));
                aNoise = noise.lowfrequ_amplitude * (1 + noise.jitter*randn(1,1));
                aNoise = aNoise * noise.noisemult;
                
                n = floor(2 * (noise.PM.timePointsN * noise.PM.TR)/fNoise + 1);
                if (n < 3) 
                    error("Number of basis functions is too low. Lower frequency or longer scanning time should be provided.")
                end
                spmdrift    = spm_drift(noise.PM.timePointsN, n);
                driftbase   = spmdrift(:,2:end);  % Eliminate the DC term
                driftbase   = sum(driftbase,2)'; 
                % Now we need to scale this noise to the actual noise values
                aNoise      = aNoise * (max(signal) - min(signal));
                tmpNoise{4} = aNoise * driftbase;
            end
            
            
            % SUM ALL NOISE MODEL
            vals = zeros(size(signal));
            for nn=1:length(tmpNoise)
                if ~isempty(tmpNoise{nn})
                    vals = vals + tmpNoise{nn};
                end
            end
            
            % RETURN VALUES
            noise.values = vals;

        end
        
        
        % Plot it
        function plot(noise)
            % I think we don't want this to be stored in the object.
            % Calculate it and return every time we need it.
            mrvNewGraphWin('Noise to be added to signal');
            noise.compute;
            plot(noise.PM.timePointsSeries, noise.values);
            grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');
        end
    end
    
    end


    


