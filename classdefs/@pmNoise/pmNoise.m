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
        Type;         % 'White', ...
        params;
        values;
    end
    
    properties(Dependent= true, SetAccess = private, GetAccess = public)
        TR;            % Seconds, it will be read from the parent class pm
    end
    
    
    
    %%
    methods (Static)
        function d = defaultsGet
            d.Type                = 'white';
            d.params.noise2signal = 0; % white noise; constant to relate to signal
            d.params.amplitude    = 0.1; % cardiac and respiratory
            d.params.frequency    = 1.25; % cardiac and respiratory
            % Convert to table and return
            d = struct2table(d.params,'AsArray',true);
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
            p.addRequired('pm'       ,             @(x)(isa(x,'prfModel')));
            % p.addParameter('type'    ,d.Type{:}  , @ischar);
            % p.addParameter('params'  ,d.params   , @isstruct);
            p.addParameter('type'    ,'white'  , @ischar);
            p.addParameter('params'  ,table2struct(d)   , @isstruct);
            p.parse(pm,varargin{:});
            
            % Initialize the pm model and hrf model parameters
            noise.PM       = pm;
            noise.Type     = p.Results.type;
            % Check if only some of the fields in params where set
            % noise.params   = pmParamsCompletenessCheck(p.Results.params,d.params);
            noise.params   = pmParamsCompletenessCheck(p.Results.params,table2struct(d));
            
        end
        
        % Methods available to this class and his childrens (friston, boynton... classes)
        function TR = get.TR(noise)
            TR = noise.PM.TR;
        end
        
        % COMPUTE
        function compute(noise)
            % TODO: move the code to the .m functions of the same name
            switch noise.Type
                case 'white'
                    % This is the amplitude of the white noise, the goal is to
                    % scale it so that it h
                    noise.PM.computeBOLD;
                    signal       = noise.PM.BOLD; % Retrieved from the parent model
                    n            = noise.params.noise2signal * (max(signal) - min(signal));
                    noise.values = n * randn([1,noise.PM.timePointsN]);
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
                    % Frequency [Hz]
                    fNoise = noise.params.frequency;
                    % Amplitude
                    aNoise = noise.params.amplitude * (max(noise.PM.BOLD) - min(noise.PM.BOLD));
                    t      = noise.PM.timePointsSeries;
                    % Calculate the noise
                    noise.values = aNoise*sin(2*pi.*t.*fNoise);
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


