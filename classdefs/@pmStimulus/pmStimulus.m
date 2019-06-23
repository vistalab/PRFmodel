classdef pmStimulus <  matlab.mixin.SetGet & matlab.mixin.Copyable
    % This is a superclass for Stimulus-s. Every particular instance of this class
    % will have different parameters, so it will be a children class. For
    % example, "bars" or "wedges". 
    % TODO: Right now uses Kendrick's convention. Bars with words is 103. 
    %
    % Syntax:
    %      stim = Stimulus;
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
        PM             ;   % prfModel that has some of the variables we need, such as TR
        fieldofviewHorz;   % Degrees
        fieldofviewVert;   % Degrees
        expName        ;   
        Binary         ;   % Logical. If true only apertures are shown. 
        barWidth       ;   % Degrees
        values         ;   % char/string with path to the stimulus .mat
        videoFileName  ;
    end
    
    properties (Dependent = true, Access = public)
        spatialSampleHorz;
        spatialSampleVert;
        XY;
        timePointsN;
        timePointsSeries;
    end
    properties(Dependent= true, SetAccess = private, GetAccess = public)
        TR;            % Seconds, it will be read from the parent class pm
    end   
    
    
    %%
    
    methods (Static)
        function d = defaultsGet
            d.fieldofviewHorz = 20;    % Degrees
            d.fieldofviewVert = 20;    % Degrees
            d.expName         = "103"; % TODO: give it meaningful names
            d.Binary          = true;  % True: only apertures. False: images/words inside apertures
            d.barWidth        = 2;     % Degrees. TODO
            % Convert to table and return
            d = struct2table(d,'AsArray',true);
        end
    end
    % Constructor
    methods
        function stim = pmStimulus(pm, varargin)
            % Obtain defaults table. If a parameters is not passed, it will use
            % the default one defined in the static function
            d = stim.defaultsGet;
            % Read the inputs
            varargin = mrvParamFormat(varargin);
            p = inputParser;
            p.addRequired ('pm'             ,                   @(x)(isa(x,'prfModel')));
            p.addParameter('fieldofviewhorz',d.fieldofviewHorz, @isnumeric);
            p.addParameter('fieldofviewvert',d.fieldofviewVert, @isnumeric);
            p.addParameter('expname'        ,d.expName{:}     , @ischar);
            p.addParameter('binary'         ,d.Binary         , @islogical);
            p.addParameter('barwidth'       ,d.barWidth       , @isnumeric);
            p.parse(pm,varargin{:});
            
            % Initialize the PM model
            stim.PM              = pm;
            % Edit the parameters
            stim.fieldofviewHorz = p.Results.fieldofviewhorz;
            stim.fieldofviewVert = p.Results.fieldofviewvert;
            stim.expName         = p.Results.expname;
            stim.Binary          = p.Results.binary;
            stim.barWidth        = p.Results.barwidth;
            % If it does not exist, create the stim file.
            % Always store just the path and the name
            stimName = strcat('Exp-',stim.expName, ...
                '_binary-', choose(stim.Binary,'true','false'), ...
                '_size-', num2str(stim.fieldofviewVert), 'x', ...
                num2str(stim.fieldofviewHorz));
            stimNameWithPath = fullfile(pmRootPath,'data','stimulus',strcat(stimName,'.mat'));
            if ~exist(stimNameWithPath, 'file')
                pmStimulusGenerate('filename', stimNameWithPath);
            end
            stim.values        =  char(stimNameWithPath);
            % Default fileName if we want to write a video of the stimuli
            stim.videoFileName = char(fullfile(pmRootPath,'local',strcat(stimName,'.avi')));
        end
        
        % Methods
        function v = get.TR(hrf)
            v       = hrf.PM.TR;
        end
        function stimValues = getStimValues(stim)
            % Obtain the values if it is a path
            % TODO: does storing it as a property occupy more space? Right now
            % it is being used as a method that obtains the values on demand
            % reading from the file system
            stimValues      = pmStimulusRead(stim.values);
        end
        function spatialSampleHorz = get.spatialSampleHorz(stim)
            spatialSampleHorz = stim.fieldofviewHorz/size(stim.getStimValues,2);
        end
        function spatialSampleVert = get.spatialSampleVert(stim)
            % Obtain the values if it is a path
            spatialSampleVert = stim.fieldofviewHorz/size(stim.getStimValues,1);
        end
        function XY = get.XY(stim)
            x = (stim.spatialSampleVert:stim.spatialSampleVert:stim.fieldofviewVert);
            x = x - mean(x);
            y = (stim.spatialSampleHorz:stim.spatialSampleHorz:stim.fieldofviewHorz);
            y = y - mean(y);
            % Calculate the spatial sampling parameters
            [X,Y] = meshgrid(x,y);
            XY = [{X},{Y}];
        end
        function v = get.timePointsN(stim)
            v = size(stim.getStimValues,3);
        end
        function timePointsSeries = get.timePointsSeries(stim)
            timePointsSeries = pmTimePointsSeries(stim.TR, stim.timePointsN);
        end
        
        
        % COMPUTE
        function compute(stim)
            % If it does not exist, create the stim file.
            % Always store just the path and the name
            stimName = strcat('Exp-',stim.expName, ...
                '_binary-', choose(stim.Binary,'true','false'), ...
                '_size-', num2str(stim.fieldofviewVert), 'x', ...
                num2str(stim.fieldofviewHorz));
            stimNameWithPath = fullfile(pmRootPath,'data','stimulus',strcat(stimName,'.mat'));
            if ~exist(stimNameWithPath, 'file')
                fprintf('Computing and storing new stimulus file in %s',stimNameWithPath)
                pmStimulusGenerate('filename', string(stimNameWithPath));
            end
            % fprintf('Retrieving stimulus file in %s',stimNameWithPath)
            stim.values        =  char(stimNameWithPath);
            % Default fileName if we want to write a video of the stimuli
            stim.videoFileName = char(fullfile(pmRootPath,'local',strcat(stimName,'.avi')));
        end
        
        
        
        
        
        
        % VISUALIZATION
        % Plot it
        function plot(stim)
            mrvNewGraphWin('Stimulus file montage');
            img = makeMontage(pmStimulusRead(stim.values));
            imagesc(img); colormap(gray);
            grid off; axis equal off; 
        end
        % TODO: add here the code to make it a video out of it.
        function toVideo(stim)
            % Write .avi for visualization
            % TODO:
            %      - check/flip the images so that the words read well in the video
            % Prepare the new file.
            vidObj = VideoWriter(char(stim.videoFileName));
            open(vidObj);
            % Create an animation.
            axis equal tight off
            set(gca,'nextplot','replacechildren');
            stimValues = pmStimulusRead(stim.values);
            for k = 1:size(stimValues,3)
                imagesc(stimValues(:,:,k)); colormap gray;
                % Write each frame to the file.
                currFrame = getframe(gcf);
                writeVideo(vidObj,currFrame);
            end
            % Close the file.
            close(vidObj);
        end
    end
    
end

