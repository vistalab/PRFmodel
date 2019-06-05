classdef pmHRF <  matlab.mixin.SetGet
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
        PM;           % prfModel that has some of the variables we need, such as TR
        Duration;     % Seconds
        Type;         % 'friston', 'boynton'
    end
    
    properties (Dependent = true, Access = public)
        tSteps;
    end
    properties(Dependent= true, SetAccess = private, GetAccess = public)
        TR;            % Seconds, it will be read from the parent class pm
    end
    
    %%
    methods
        % Constructor
        function hrf = pmHRF
            hrf.Duration = 20; % Seconds
            hrf.Type     = 'friston';
        end
        % Methods available to this class and his childrens (friston, boynton... classes)
        % Put here the values that are duplicated from prfModel but that we need
        % to work with them independently.
        function v = get.TR(hrf)
            v       = hrf.PM.TR;
        end
        function tSteps = get.tSteps(hrf)
            % I think we don't want this to be stored in the object.
            % Calculate it and return every time we need it.
            tSteps  = pmTimePointsSeries(hrf.TR, hrf.Duration);
        end
        % Plot it
        function plot(hrf)
            % I think we don't want this to be stored in the object.
            % Calculate it and return every time we need it.
            mrvNewGraphWin([hrf.Type ' HRF']);
            plot(hrf.tSteps, hrf.values);
            grid on; xlabel('Time (sec)'); ylabel('Relative amplitude');
        end
    end
    % Here the methods that are only allowed from the pm class
    %     methods (Access = {?prfModel,?prfModel_basic})
    %         function setTR(hrf, tr)
    %             hrf.TR       = tr;
    %         end
    %     end
    
end



