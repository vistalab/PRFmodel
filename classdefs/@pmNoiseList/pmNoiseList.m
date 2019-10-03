classdef pmNoiseList <  matlab.mixin.SetGet & matlab.mixin.Copyable
    % This is a superclass for List of Noise-s. 
    % This class will invoke different noise models and add them. The result
    % will be added to the noiseless BOLD signal in the main model. 
    %
    % Syntax:
    %      ;
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
        noiseparams;
        noisevalues;
    end
    properties(Dependent = true, SetAccess = public, GetAccess = public)
        % List of all possible options
    end
    
    %%
    methods (Static)
        function D = defaultsGet
            d1.Type                = 'white';
            d1.params.noise2signal = 0; % white noise; constant to relate to signal
            d1.params.amplitude    = 0.1; % cardiac and respiratory
            d1.params.frequency    = 1.25; % cardiac and respiratory
            % Convert to table and return
            d1 = struct2table(d1,'AsArray',true);
            
            d2.Type                = 'cardiac';
            d2.params.noise2signal = 0; % white noise; constant to relate to signal
            d2.params.amplitude    = 0.1; % cardiac and respiratory
            d2.params.frequency    = 1.25; % cardiac and respiratory
            % Convert to table and return
            d2 = struct2table(d2,'AsArray',true);
            
            D = {d1,d2};
                                 
        end
    end
    methods
        % Constructor
        function noises = pmNoiseList(pm, varargin)
            % Obtain defaults table. If a parameters is not passed, it will use
            % the default one defined in the static function
            D = noises.defaultsGet;
            % Read the inputs
            varargin = mrvParamFormat(varargin);
            p = inputParser;
            p.addRequired('pm',             @(x)(isa(x,'prfModel')));
            p.addParameter('d', D , @iscell);
            p.parse(pm,varargin{:});
            
            % Initialize the pm model and hrf model parameters
            noises.PM             = pm;
            noises.noiselist      = p.Results.d;
        end
        
        function compute(noises)
            % This invokes the individual pmNoise models. 
            pm    = noises.PM;
            for nn=1:length(noises.noiselist)
                noises.noisevalues{nn} = pmNoise(pm, ...
                                   'Type', noises.noiseparams{nn}.Type{:}, ...
                                   'params',noises.noiseparams{nn}.params);
            end
            
        end
        % MAKE A NOISE PLOT, WITH OPTIONS TO PLOT THEM SEPARATELY OR COMBINED
        % % Plot it
        % function plot(noise)
        %   % I think we don't want this to be stored in the object.
        %   % Calculate it and return every time we need it.
        %   mrvNewGraphWin('Noise to be added to signal');
        %   plot(noise.PM.timePointsSeries, noise.values);
        %   grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');
        % end
    end
    
end


