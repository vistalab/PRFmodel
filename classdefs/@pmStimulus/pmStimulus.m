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
     % RFSIZE = 0.5;
     mrvNewGraphWin;
     window = false;
     nrow = 2; ncol = 4;
     pm               = prfModel;  
     pm.RF.sigmaMajor = 0.5;
     pm.RF.sigmaMinor = pm.RF.sigmaMajor;
     pm.Stimulus.frameduration = 4;
     
     subplot(nrow,ncol,1)
     pm.TR = 1;
     pm.Stimulus.durationSecs = 300;
     pm.Stimulus
     pm.Stimulus.compute
     % pm.Stimulus.plot
     pm.plot('what','nonoisetimeseries','window',window)
     
     subplot(nrow,ncol,2)
     pm.TR = 1;
     pm.Stimulus.durationSecs = 150;
     pm.Stimulus
     pm.Stimulus.compute
     % pm.Stimulus.plot
     pm.plot('what','nonoisetimeseries','window',window)
     
     subplot(nrow,ncol,3)
     pm.TR = 2;
     pm.Stimulus.durationSecs = 300;
     pm.Stimulus.compute
     % pm.Stimulus.plot
     pm.plot('what','nonoisetimeseries','window',window)
    
     subplot(nrow,ncol,4)
     pm.TR = 2;
     pm.Stimulus.durationSecs = 150;
     pm.Stimulus.compute
     % pm.Stimulus.plot
     pm.plot('what','nonoisetimeseries','window',window)
      
    % RFSIZE = 2;
     pm               = prfModel;  
     pm.RF.sigmaMajor = 2;
     pm.Stimulus.frameduration = 4;
     pm.RF.sigmaMinor = pm.RF.sigmaMajor;    
     % pm.RF.plot    
   
     subplot(nrow,ncol,5)
     pm.TR = 1;
     pm.Stimulus.durationSecs = 300;
     pm.Stimulus
     pm.Stimulus.compute
     % pm.Stimulus.plot
     pm.plot('what','nonoisetimeseries','window',window)
     
     subplot(nrow,ncol,6)
     pm.TR = 1;
     pm.Stimulus.durationSecs = 150;
     pm.Stimulus
     pm.Stimulus.compute
     % pm.Stimulus.plot
     pm.plot('what','nonoisetimeseries','window',window)
     
     subplot(nrow,ncol,7)
     pm.TR = 2;
     pm.Stimulus.durationSecs = 300;
     pm.Stimulus.compute
     % pm.Stimulus.plot
     pm.plot('what','nonoisetimeseries','window',window)
    
     subplot(nrow,ncol,8)
     pm.TR = 2;
     pm.Stimulus.durationSecs = 150;
     pm.Stimulus.compute
     % pm.Stimulus.plot
     pm.plot('what','nonoisetimeseries','window',window)
   %} 
   %{
     pm.TR = 4;
     pm.Stimulus.compute
     pm.Stimulus.plot
     pm.plot('what','nonoise')
     
     pm.TR = 1.5;
     pm.Stimulus.durationSecs = (6+128)*1.5;
     pm.Stimulus.compute
     pm.Stimulus.plot
     pm.plot('what','nonoise')

    
     pm.TR = 6;
     pm.Stimulus.durationSecs=300;
     pm.Stimulus.compute
     pm.Stimulus.plot
     pm.plot('what','nonoise')
   %}
    
   % Create a 15 second on/off stimulus, 5 off, 5 on, 5 off
   %{
     See stimtests.m
    
   %}
    
    properties
        PM               ;   % prfModel that has some of the variables we need, such as TR
        fieldofviewHorz  ;   % Degrees
        fieldofviewVert  ;   % Degrees
        expName          ;   
        Binary           ;   % Logical. If true only apertures are shown. 
        Resize           ;   % Logical. If true image resized according to ImageSideSize
        ResizedHorz      ;   % Numeric. Size in pixels of resized horizontal side
        ResizedVert      ;   % Numeric. Size in pixels of resized vertical side
        barWidth         ;   % Degrees
        durationSecs     ;   % Numeric, duration of the stimuli in secs, default 300secs
        frameduration    ;   % Numeric, how many frames we want a refresh to last
        Shuffle          ;   
        values           ;   % char/string with path to the stimulus .mat
        videoFileName    ;
        niftiFileName    ;
        DataPath         ;
        LocalPath        ;
        userVals         ;  % No calculation required, accept values directly. Empty by default. 
    end
    
    properties (Dependent = true, Access = public)
        spatialSampleHorz;
        spatialSampleVert;
        XY               ;
        timePointsN      ;
        timePointsSeries ;
        Name             ;
        maxStimCenter    ;
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
            d.Resize          = true;  % True by default. Resize to ResizedHorz x ResizedVert
            d.ResizedHorz     = 101;   % 101 is the default in mrVista
            d.ResizedVert     = 101;   % 101 is the default in mrVista
            d.barWidth        = 2;     % Degrees. TODO
            d.durationSecs    = 200;   % Seconds
            d.frameduration   = 4;   
            d.Shuffle         = false; % Shuffle bars or content
            
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
            p.addParameter('resize'         ,d.Resize         , @islogical);
            p.addParameter('resizedhorz'    ,d.ResizedHorz    , @isnumeric);
            p.addParameter('resizedvert'    ,d.ResizedVert    , @isnumeric);
            p.addParameter('barwidth'       ,d.barWidth       , @isnumeric);
            p.addParameter('durationsecs'   ,d.durationSecs   , @isnumeric);
            p.addParameter('frameduration'  ,d.frameduration  , @isnumeric);
            p.addParameter('shuffle'        ,d.Shuffle        , @islogical);
            p.addParameter('uservals'       ,[]               , @isnumeric);
            p.parse(pm,varargin{:});
            
            % Initialize the PM model
            stim.PM              = pm;
            % Edit the parameters
            stim.fieldofviewHorz = p.Results.fieldofviewhorz;
            stim.fieldofviewVert = p.Results.fieldofviewvert;
            stim.expName         = p.Results.expname;
            stim.Binary          = p.Results.binary;
            stim.Resize          = p.Results.resize;
            stim.ResizedHorz     = p.Results.resizedhorz;
            stim.ResizedVert     = p.Results.resizedvert;
            stim.barWidth        = p.Results.barwidth;
            stim.durationSecs    = p.Results.durationsecs;
            stim.frameduration   = p.Results.frameduration;
            stim.Shuffle         = p.Results.shuffle;
            stim.userVals        = p.Results.uservals;
            
            % If we pass uservales, override the calculations
            if isempty(stim.userVals)
                % If it does not exist, create the stim file.
                % Always store just the path and the name
                stim.LocalPath       = fullfile(pmRootPath,'local');
                stim.DataPath        = fullfile(pmRootPath,'data','stimulus');
                stimNameWithPath     = fullfile(stim.DataPath, [stim.Name '.mat']);
                if ~exist(stimNameWithPath, 'file')
                    % TODO: add all parameters, see .compute below
                    pmStimulusGenerate('filename', stimNameWithPath,...
                        'totalduration',stim.durationSecs, ...
                        'TR', stim.TR, ...
                        'frameduration',stim.frameduration);
                end
                stim.values        =  char(stimNameWithPath);
                % Default fileName if we want to write a video of the stimuli
                stim.videoFileName = fullfile(stim.LocalPath,[stim.Name '.avi']);
                % Default fileName if we want to write a nifti of the stimuli
                stim.niftiFileName = fullfile(stim.LocalPath,[stim.Name '.nii.gz']);
            end
        end
        function Name = get.Name(stim)
            Name = [...
                   'Exp-'          char(stim.expName) ...
                   '_bin-'         choose(stim.Binary,'true','false') ...
                   '_size-'        num2str(stim.fieldofviewVert) 'x' ...
                                   num2str(stim.fieldofviewHorz) ...
                   '_resize-'      choose(stim.Resize,'true', 'false') ...
                   '_Horz-'        num2str(stim.ResizedVert) 'x' ...
                                   num2str(stim.ResizedHorz) ...
                   '_barW-'        num2str(stim.barWidth) ...
                   '_dur-'         num2str(stim.durationSecs) ...
                   '_TR-'          num2str(stim.TR) ...
                   '_framedur-'    num2str(stim.frameduration) ...
                   '_Shuffle-'     choose(stim.Resize,'true', 'false') ...
                   ];
               assert(isa(Name, 'char'));
        end
        % Methods
        function v = get.TR(hrf)
            v       = hrf.PM.TR;
        end
        function stimValues = getStimValues(stim)
            % Obtain the values if it is a path
            if iscell(stim.userVals);
                sv = stim.userVals{:};
            else
                sv = stim.userVals;
            end
            if isempty(sv)
                stimValues = pmStimulusRead(stim.values);
            else
                stimValues = sv;
            end
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
        function imageCentroid = get.maxStimCenter(stim)
            % If we want to normalize the contrast value between voxels, we need
            % to know what it the possible max value for our rfsize
            % To maintain the voxel calculations independent, we can calculate
            % it on the fly every time (TODO: check performance issues...)
            % This function will return the center of the maximum stimuls time
            % point, so that we can center an RF and do a calculation there
            
            % Get stimuli
            vals = stim.getStimValues;
            % Get the image with the most nonzeros
            [~,maxindx] = max(squeeze(sum(sum(vals))));
            % Select that image and get it
            maximage = vals(:,:,maxindx);
            % Get the centroid of the image
            imageCentroid = regionprops(maximage,'centroid');
            % Round the values, these are cols and rows
            imageCentroid = round(imageCentroid.Centroid);
            % Change from rows and columns to degrees
            % First, center it
            if all(iseven(stim.ResizedHorz))
                imageCentroid(2) = imageCentroid(2) - stim.ResizedHorz/2;
            else
                imageCentroid(2) = imageCentroid(2) - ((stim.ResizedHorz-1)/2 + 1);
            end
            if iseven(stim.ResizedVert)
                imageCentroid(1) = imageCentroid(1) - stim.ResizedVert/2;
            else
                imageCentroid(1) = imageCentroid(1) - ((stim.ResizedVert-1)/2 + 1);
            end
            % Then convert it
            imageCentroid(2)     = (stim.spatialSampleHorz * imageCentroid(2));
            imageCentroid(1)     = (stim.spatialSampleVert * imageCentroid(1));
        end
       
        % COMPUTE
        function compute(stim)
            % If it does not exist, create the stim file.
            % Store just the path and the name
            
            % If the user passed its values, override this and maintain the default
            if iscell(stim.userVals)
                uv = stim.userVals{:};
            else
                uv = stim.userVals;
            end
            if isempty(uv)
                stimNameWithPath = fullfile(stim.DataPath, [stim.Name '.mat']);
                if ~exist(stimNameWithPath, 'file')
                    fprintf('Computing and storing new stimulus file in %s',stimNameWithPath)
                    % TODO: pass all the variables and make it more flexible
                    pmStimulusGenerate('filename', stimNameWithPath,...
                                        'totalduration',stim.durationSecs, ...
                                        'TR', stim.TR, ...
                                        'frameduration',stim.frameduration);
                end
                % fprintf('Retrieving stimulus file in %s',stimNameWithPath)
                stim.values        =  char(stimNameWithPath);
            end
        end
        
  
        
        % VISUALIZATION
        % Plot it
        function plot(stim, varargin)
            % Read the inputs
            varargin = mrvParamFormat(varargin);
            p = inputParser;
            p.addRequired ('stim'  ,  @(x)(isa(x,'pmStimulus')));
            p.addParameter('slicelist',[], @isnumeric);
            p.addParameter('window',true, @islogical);
            p.parse(stim,varargin{:});
            slicelist = p.Results.slicelist;
            w = p.Results.window;
            
            if w, mrvNewGraphWin('Stimulus file montage'); end
            image3D = pmStimulusRead(stim.values);
            [sz1,sz2,sz3] = size(image3D);
            img = pmMakeMontage(image3D,slicelist);
            imagesc(img); colormap(gray);
            grid off; axis equal off; 
            aa = gca;text(1,aa.YLim(2)*(1.05),sprintf('%3.0fx%3.0fx%3.0f',sz1,sz2,sz3));
        end
        function toVideo(stim,varargin)
             % Read the inputs
            varargin = mrvParamFormat(varargin);
            p = inputParser;
            p.addRequired ('stim'  ,  @(x)(isa(x,'pmStimulus')));
            p.addParameter('fname',stim.videoFileName, @ischar);
            p.parse(stim,varargin{:});
            fname = p.Results.fname;
            
            % Write .avi for visualization
            % TODO:
            %      - check/flip the images so that the words read well in the video
            % Prepare the new file.
            vidObj = VideoWriter(fname);
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
        function fname = toNifti(stim,varargin)
            % Read the inputs
            varargin = mrvParamFormat(varargin);
            p = inputParser;
            p.addRequired ('stim'  ,  @(x)(isa(x,'pmStimulus')));
            p.addParameter('fname',stim.niftiFileName, @ischar);
            p.parse(stim,varargin{:});
            fname = p.Results.fname;
            stimValues = pmStimulusRead(stim.values,'format','nifti');
            writeFileNifti(niftiCreate('data', stimValues, ...
                'fname',fname, 'tr',stim.TR));
            
        end
    end
    
end

