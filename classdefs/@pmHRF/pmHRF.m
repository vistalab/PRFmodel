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
        PM;            % prfModel that has some of the variables we need, such as TR
        Duration;      % Seconds
        params;        % Different values depending on the type of HRF
        values;
    end
    
    properties (Dependent = true, Access = public)
        tSteps;
    end
    properties(Dependent= true, SetAccess = private, GetAccess = public)
        Type;
        TR;            % Seconds, it will be read from the parent class pm
    end
    
    %%
    methods
        % Methods available to this class and his childrens (friston, boynton... classes)
        function v = get.TR(hrf)
            v = hrf.PM.TR;
        end
        function tSteps = get.tSteps(hrf)
            % I think we don't want this to be stored in the object.
            % Calculate it and return every time we need it.
            tSteps  = 0: hrf.TR: hrf.Duration;
        end
        
        
        
        
        
        % Compute it
        function hrf = pmHRF_friston(pm,varargin)
            % Create default parameter struct
            a  = [6  ,  12]; b  = [0.9, 0.9]; c  = 0.35;
            params.a = a; params.b = b; params.c = c;
            
            % Read the inputs
            varargin = mrvParamFormat(varargin);
            p = inputParser;
            p.addRequired('pm'       ,        @(x)(isa(x,'prfModel')));
            p.addParameter('params'  ,params, @isstruct);
            p.addParameter('duration',20    , @isnumeric);
            p.parse(pm,varargin{:});
            % Assign it
            params   = p.Results.params;
            Duration = p.Results.duration;
            
            % Initialize the pm model and hrf model parameters
            hrf.PM       = pm;
            hrf.Type     = 'friston';
            hrf.Duration = Duration;
            hrf.params   = params;
        end
        
        function compute(hrf)
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
        
        
        
        % Return the defaults that can be changed to generate different HRF
        function HRF = defaultsTable
            HRF = table();
        end
    end

    
end



