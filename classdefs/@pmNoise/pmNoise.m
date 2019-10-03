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
    %
    
    % Examples
    %{
       
    %}
    
    properties
        PM;           % prfModel that has some of the variables we need, such as TR
        white;
        white_noise2signal;
        cardiac;
        cardiac_amplitude;
        cardiac_frequency;
        respiratory;
        respiratory_amplitude;
        respiratory_frequency;
        values;
    end
    
    properties(Dependent= true, SetAccess = private, GetAccess = public)
        TR;            % Seconds, it will be read from the parent class pm
    end
    
    
    
    %%
    methods (Static)
        function d = defaultsGet
            % White Noise
            d.white              = false; 
            d.white_noise2signal = 0.1; % white noise; constant to relate to signal
            
            % Cardiac
            d.cardiac            = false; 
            d.cardiac_amplitude  = 0.1;
            d.cardiac_frequency  = 1.17; 
            
            % Respiratory
            d.respiratory            = false; 
            d.respiratory_amplitude  = 0.1;
            d.respiratory_frequency  = 0.2; 
            
            % Low frequency
            
            
            % Autocorrelation
            
            
            
            % Task related
            
            
            
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
            p.addRequired('pm'                                             , @(x)(isa(x,'prfModel')));
            p.addParameter('white'                 ,d.white                , @islogical);
            p.addParameter('white_noise2signal'    ,d.white_noise2signal   , @isnumeric);
            p.addParameter('cardiac'               ,d.cardiac              , @islogical);
            p.addParameter('cardiac_amplitude'     ,d.cardiac_amplitude    , @isnumeric);
            p.addParameter('cardiac_frequency'     ,d.cardiac_frequency    , @isnumeric);
            p.addParameter('respiratory'           ,d.respiratory          , @islogical);
            p.addParameter('respiratory_amplitude' ,d.respiratory_amplitude, @isnumeric);
            p.addParameter('respiratory_frequency' ,d.respiratory_frequency, @isnumeric);
            p.parse(pm,varargin{:});
            
            % Initialize the pm model and hrf model parameters
            noise.PM       = pm;
            % params
            noise.white    = p.Results.white;
            noise.white_noise2signal = p.Results.white_noise2signal;
            noise.cardiac            = p.Results.cardiac;
            noise.cardiac_amplitude  = p.Results.cardiac_amplitude;
            noise.cardiac_frequency  = p.Results.cardiac_frequency;
            noise.respiratory            = p.Results.respiratory;
            noise.respiratory_amplitude  = p.Results.respiratory_amplitude;
            noise.respiratory_frequency  = p.Results.respiratory_frequency;
        end
        
        % Methods available to this class and his childrens (friston, boynton... classes)
        function TR = get.TR(noise)
            TR = noise.PM.TR;
        end
        
        % COMPUTE
        function compute(noise)
            % Simplifying this code. Noise is a simple class, and the noise
            % models will be specified in the parameters. 
            
            % Obtain the BOLD signal without noise. 
            noise.PM.computeBOLD;
            signal       = noise.PM.BOLD; % Retrieved from the parent model
            tmpNoise     = {};
            
            % WHITE
            tmpNoise{1}  = {};
            if noise.white
                n = noise.white_noise2signal * (max(signal) - min(signal));
                tmpNoise{1} = n * randn([1,noise.PM.timePointsN]);
            end
            
            % CARDIAC
            tmpNoise{2}  = {};
            if noise.cardiac
                % Frequency [Hz]
                fNoise = noise.cardiac_frequency;
                % Amplitude
                aNoise = noise.cardiac_amplitude * (max(noise.PM.BOLD) - min(noise.PM.BOLD));
                t      = noise.PM.timePointsSeries;
                % Calculate the noise
                tmpNoise{2} = aNoise*sin(2*pi.*t.*fNoise);                
            end
            
            % RESPIRATORY
            tmpNoise{3}  = {};
            if noise.respiratory
                % Frequency [Hz]
                fNoise = noise.respiratory_frequency;
                % Amplitude
                aNoise = noise.respiratory_amplitude * (max(noise.PM.BOLD) - min(noise.PM.BOLD));
                t      = noise.PM.timePointsSeries;
                % Calculate the noise
                tmpNoise{2} = aNoise*sin(2*pi.*t.*fNoise);                
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
            
            
            
            % Old code, to be deleted
            %{
            % TODO: move the code to the .m functions of the same name
            switch noise.Type
                case 'white'  % Indstrumen
                case 'pink'  % fisiological
                case 'cardiac'
                    %                     % Make sure that when cardiac is selected, params comes with
                    %                     % the correct parameters, otherwise use defaults for cardiac
                    %                     if ~isfield(noise.params,'frequency')
                    %                         params.frequency = 1.25;% 1.25 Hz:75 beats/min
                    %                     end
                    %                     if ~isfield(noise.params,'amplitude')
                    %                         params.amplitude = 0.10;  % 0-1 proportion over mean BOLD signal
                    %                     end
                    %
                    %                     % Update the obj just in case
                    %                     noise.params = params;
                    
                    % Compute the noise
                    signal = noise.PM.BOLD;  % BOLD signal, add noise to this
                    
                case 'respiratory'
                    %                     % Make sure that when respiratory is selected, params comes with
                    %                     % the correct parameters, otherwise use defaults for respiratory
                    %                     if ~isfield(noise.params,'frequency')
                    %                         params.frequency = 0.3;  % 0.3 Hz : 18 breaths/min
                    %                     end
                    %                     if ~isfield(noise.params,'amplitude')
                    %                         params.amplitude = 0.1;  % 0-1 proportion over BOLD signal amplitude
                    %                     end
                    %
                    %                     % Update the obj just in case
                    %                     noise.params = params;
                    
                    % Compute the noise
                    signal = noise.PM.BOLD;  % BOLD signal, add noise to this
                    % Frequency [Hz]
                    fNoise = noise.params.frequency;
                    % Amplitude
                    aNoise = noise.params.amplitude * (max(noise.PM.BOLD) - min(noise.PM.BOLD));
                    t      = noise.PM.timePointsSeries;
                    % Calculate the noise
                    noise.values = aNoise*sin(2*pi.*t.*fNoise);
                case 'eyemovement'
                    warning('eyemovement noise model not implemented yet')
                otherwise
                    error('Noise model %s not implemented or valid.',noise.Type);
            end
            %}
        end
        
        
        % Plot it
        function plot(noise)
            % I think we don't want this to be stored in the object.
            % Calculate it and return every time we need it.
            mrvNewGraphWin('Noise to be added to signal');
            plot(noise.PM.timePointsSeries, noise.values);
            grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');
        end
    end
    
end


