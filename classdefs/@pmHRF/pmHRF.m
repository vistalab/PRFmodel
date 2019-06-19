classdef pmHRF <  matlab.mixin.SetGet & matlab.mixin.Copyable
    % This is a superclass for HRF-s. Every particular instance of this class
    % will have different parameters, so it will be a children class. For
    % example, "Friston" or "Boynton"
    %
    % Syntax:
    %      hrf = HRF;
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
        Type;          % Firston, Boyton, canonical...
        PM;            % prfModel that has some of the variables we need, such as TR
        Duration;      % Seconds
        params;        % Different values depending on the type of HRF
        values;
    end
    
    properties (Dependent = true, Access = public)
        tSteps;
    end
    properties(Dependent= true, SetAccess = private, GetAccess = public)
        TR;            % Seconds, it will be read from the parent class pm
    end
    
    %%
    methods (Static)
        function d = defaultsGet
            d.Type              = 'friston';
            d.Duration          = 20;
            d.params.a          = [6  ,  12];
            d.params.b          = [0.9, 0.9];
            d.params.c          = 0.35;
            % Convert to table and return
            d = struct2table(d,'AsArray',true);
        end
    end
    methods
        function hrf = pmHRF(pm,varargin)
            % Obtain defaults table. If a parameters is not passed, it will use
            % the default one defined in the static function
            d = hrf.defaultsGet;
            % Read the inputs
            varargin = mrvParamFormat(varargin);
            p = inputParser;
            p.addRequired('pm'       ,        @(x)(isa(x,'prfModel')));
            p.addParameter('type'    ,d.Type{:}  , @ischar);
            p.addParameter('params'  ,d.params   , @isstruct);
            p.addParameter('duration',d.Duration , @isnumeric);
            p.parse(pm,varargin{:});
            % Assign it
            params   = p.Results.params;
            Duration = p.Results.duration;
            
            % Initialize the pm model and hrf model parameters
            hrf.PM       = pm;
            hrf.Type     = p.Results.type;
            hrf.Duration = p.Results.duration;
            hrf.params   = p.Results.params;
        end
        
        function compute(hrf)
            switch hrf.Type
                case 'friston'
                    % TODO: use the function again, remove code from here. 
                    a = hrf.params.a;
                    b = hrf.params.b;
                    c = hrf.params.c;
                    
                    t = hrf.tSteps;
                    % Calculate d
                    for ii = 1:2
                        d(ii) = a(ii)*b(ii);
                    end
                    % Calculate actual values
                    hrf.values = (t/d(1)).^a(1)   .* exp(-(t - d(1))/b(1)) ...
                        - c*(t/d(2)).^a(2) .* exp(-(t-d(2))/b(2));
                case 'boynton'
                case 'canonical'
                otherwise
                    error('HRF method %s not implemented or valid.',hrf.Type);
            end
        end
        
        
        
        
        
        
        
        
        
        % Methods available to this class and his childrens (friston, boynton... classes)
        function v = get.TR(hrf)
            v = hrf.PM.TR;
        end
        function tSteps = get.tSteps(hrf)
            % I think we don't want this to be stored in the object.
            % Calculate it and return every time we need it.
            tSteps  = 0: hrf.TR: hrf.Duration;
        end
        
        
        
        
        % Plot it
        function plot(hrf)
            % Calculate it and return every time we need it.
            % Compute it just in case, to have the latest version
            hrf.compute;
            % Plot it
            mrvNewGraphWin([hrf.Type ' HRF']);
            plot(hrf.tSteps, hrf.values);
            grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');
        end

    end

    
end



