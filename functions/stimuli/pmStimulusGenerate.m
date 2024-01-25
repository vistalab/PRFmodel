function [stim, params] = pmStimulusGenerate(varargin)
% Wrapper function to generate stimuli
%{
expname         = "103";
onlymasks       = false;
checkimages     = true;
wantdownsample  = true;
wantresize      = false;
resizedvert     = 101;
resizedhorz     = 101;
normalize01     = true;
binarize        = false;
savestimmat     = true;
shuffle         = false;
shuffleseed     = 12345;
filename        = '/Users/glerma/soft/morphing/DATA/retWordsMagno/Words_test.mat';
barwidth        = 2;
totalduration   = 300;
tr              = 1.5;
frameduration   = 4;

stim = pmStimulusGenerate('expname', expname,...
                          'onlymasks', onlymasks, ...
                          'checkimages', checkimages, ...
                          'wantdownsample', wantdownsample, ...
                          'wantresize', wantresize, ...
                          'resizedvert', resizedvert, ...
                          'resizedhorz', resizedhorz, ...
                          'normalize01', normalize01, ...
                          'binarize', binarize, ...
                          'savestimmat', savestimmat, ...
                          'shuffle', shuffle, ...
                          'shuffleseed', shuffleseed, ...
                          'filename', filename, ...
                          'barwidth', barwidth, ...
                          'totalduration', totalduration, ...
                          'tr', tr, ...
                          'frameduration', frameduration);
%}

%% Read parameters
bgimages  = fullfile(pmRootPath,'data','words','spanish.mat');
masks = fullfile(pmRootPath,'data','images','maskimages.mat');
stimulusdir = fullfile(pmRootPath,'data');
params.tr   = 999;

varargin = mrvParamFormat(varargin);
            
p = inputParser;
p.addParameter('bgfile'        , bgimages, @isstring) ; % Path with background images that will be uncovered by bars
p.addParameter('masks'         , masks, @isstring) ; % Path with background images that will be uncovered by bars
p.addParameter('stimulusdir'   , stimulusdir, @isstring) ;
p.addParameter('expname'       , "103", @isstring) ; % See "doc pmShowmulticlass" for explanation of alternatives                
p.addParameter('onlymasks'     , true , @islogical); % Select if only mask or mask and images are required
p.addParameter('checkimages'   , false, @islogical); % To visualize a sample of the stimuli in each step, set to true
p.addParameter('wantdownsample', true , @islogical); % To downsample to the same number of volume images
p.addParameter('wantresize'    , true , @islogical); % To resize the images to 100x100 for example
p.addParameter('resizedhorz'   , 101  , @isnumeric); % Size of resized side. horz
p.addParameter('resizedvert'   , 101  , @isnumeric); % Size of resized side. vert
p.addParameter('normalize01'   , true , @islogical); % To normalize all values between 0 and 1
p.addParameter('binarize'      , true , @islogical); % To binarize (threshold 0.5)
p.addParameter('savestimmat'   , true , @islogical); % To save images in fileName
p.addParameter('shuffle'       , false, @islogical); % To save images in fileName
p.addParameter('shuffleseed'   , 12345);             % It can be 'shuffle' or any integer number
p.addParameter('barwidth'      , 2    , @isnumeric); % Width of bar in deg
p.addParameter('filename'      , './stimulus.mat' , @isstring); % filename
p.addParameter('saveparamstimfile'      , false , @islogical); 
p.addParameter('totalduration' , 300    , @isnumeric); % Duration in seconds
p.addParameter('tr'            , 2      , @isnumeric); % TR
p.addParameter('params'        , params , @isstruct); 
p.addParameter('frameduration' , 4      , @isnumeric); % 

p.parse(varargin{:});
            
bgfile          = p.Results.bgfile;
masks           = p.Results.masks;
expName         = p.Results.expname;
onlyMasks       = p.Results.onlymasks;
checkImages     = p.Results.checkimages;
wantDownsample  = p.Results.wantdownsample;
wantResize      = p.Results.wantresize;
ResizedVert     = p.Results.resizedvert;
ResizedHorz     = p.Results.resizedhorz;
normalize01     = p.Results.normalize01;
binarize        = p.Results.binarize;
saveStimMat     = p.Results.savestimmat;
Shuffle         = p.Results.shuffle;
shuffleSeed     = p.Results.shuffleseed;
fileName        = p.Results.filename;
barWidth        = p.Results.barwidth;
saveParamStimFile=p.Results.saveparamstimfile;
totalDuration   = p.Results.totalduration;
TR              = p.Results.tr;
params          = p.Results.params;
frameduration   = p.Results.frameduration;

%% Generate stimuli
% TODO: simplify pmShowmulticlass to take just the stimuli we are interested with,
%       for example "bars with words" (this is 103), or similar. 
expnum = str2num(expName);

% fileName       = 'Exp-103_binary-true_size-100x100';
matFileName    = [fileName '.mat'];

% To write a video file with the stimuli
createVideo    = false;
videoFileName  = [fileName '.avi'];

% directory where we have the spanish word matrices
load(bgfile,'images');

% Mask images and fixation are here
load(masks,'maskimages');

% loads the variable spatialoverlay
fixation = fullfile(pmRootPath,'data','images','fixationgrid.mat');
a1       = load(fixation);


% background color
grayval = uint8(127); 

% how many frames we want a refresh to last. default: 15. previously: 4. 
% the lower the faster
% frameduration = 4;

dres    = [];           % <dres> (optional) is
                        %   [A B] where this is the desired resolution to imresize the images to (using bicubic interpolation).
                        %     if supplied, the imresize takes place immediately after loading the images in, and this imresized 
                        %     version is what is cached in the output of this function.
                        %  -C where C is the <scfactor> input in ptviewmovie.m
                        %   default is [] which means don't do anything special.  
                        % seems to crash for values >= 1
                            
offset  = [0 0];        % [X Y] where X and Y are the horizontal and vertical
                        % offsets to apply.  for example, [5 -10] means shift 
                        % 5 pixels to right, shift 10 pixels up.
                           
%%  There is one file that contains many different mask images
ptres = [];  % display resolution. [] means to use current display resolution.

% fixation dot
fixationinfo = {uint8([255 0 0; 0 0 0; 0 255 0]) 0.5};  % dot colors and alpha value
fixationsize = 10;         % dot size in pixels
meanchange = 3;            % dot changes occur with this average interval (in seconds)
changeplusminus = 2;       % plus or minus this amount (in seconds)


movieflip = [0 0];         % [A B] where A==1 means to flip vertical dimension
                           % and B==1 means to flip horizontal dimension
% set rand state
rand('state',sum(100*clock));
randn('state',sum(100*clock));

% prepare inputs
trialparams = [];
ptonparams = {ptres,[],0};
iscolor = 1;
soafun = @() round(meanchange*(60/frameduration) + changeplusminus*(2*(rand-.5))*(60/frameduration));

if ~exist('images','var')
  images = [];
  maskimages = [];
end
filename = [];
triggerkey = [];
gridImage = [];
tfun = [];
[images,maskimages,frameorder] = ...
  pmShowmulticlass(filename,offset,movieflip,frameduration,fixationinfo,...
                   fixationsize,tfun, ...
                   ptonparams,soafun,0,images,expnum,[],grayval,iscolor,...
                   [],[],[],dres,triggerkey, ...
                   [],trialparams,[],maskimages,gridImage,stimulusdir);  
 images = images{1};            
if checkImages
    figure(1);
    image(images(:,:,:,8) .* maskimages(:,:,2000));colormap gray;axis equal tight off;
end

% This is what has been shown
ImgRefs  = frameorder(1,:);
maskRefs = frameorder(2,:);

% Create the new trAdjImages stimuli set
% Images
trAdjImages = zeros([size(images,1),size(images,2),size(ImgRefs,2)]);
trAdjImages(:,:,ImgRefs~=0) = images(:,:,ImgRefs(ImgRefs~=0));
% Masks
trAdjMaskImages = zeros([size(maskimages,1),size(maskimages,2),size(maskRefs,2)]);
trAdjMaskImages(:,:,maskRefs~=0) = maskimages(:,:,maskRefs(maskRefs~=0));


% We can just work with masks or create the final stimuli multplying the mask
% with the actual images


if onlyMasks
    stim = trAdjMaskImages;
else
    stim = trAdjImages .* trAdjMaskImages;
end
if checkImages
    figure(2); image(stim(:,:,1000)); colormap gray; axis equal tight off;
end

% Downsample
if wantDownsample
    % TODO: create a version depending on refresh rate, tr-s and number of images
    % Parameters
    numOfScans   = round(totalDuration/TR);    
    % Create index
    imgIndx      = round(linspace(1, size(frameorder,2), numOfScans));
    % Downsample
    stim         = stim(:,:,imgIndx);
end
if checkImages
    figure(3); image(stim(:,:,66)); colormap gray; axis equal tight off;
end


% If working only with masks maybe we want to resize and binarize
if wantResize
    temp = zeros(ResizedVert, ResizedHorz, size(stim,3));
    for p=1:size(stim, 3)
        temp(:,:,p) = imresize(stim(:,:,p),[ResizedVert ResizedHorz],'cubic');
    end
    % Make sure value go from 0 to 1
    %     normTemp = temp - min(temp(:));
    %     normTemp = normTemp ./ max(normTemp(:));
    %     stim = normTemp;
    stim = temp;
    
    % ensure that all values are between 0 and 1
    % stim(stim < 0) = 0;
    % stim(stim > 1) = 1;
end
if checkImages
    figure(4); image(stim(:,:,66)); colormap gray; axis equal tight off;
end

if normalize01
    % Normalize to 0>1 and binarize 
    nstim   = stim - min(stim(:));
    nstim   = nstim ./ max(nstim(:));
    stim    = nstim;
    % We cannot binarize if it has not been normalized. 
    % TODO: warning if normalize01 is false and binarize true?
    if binarize
        % Add it to the pm we just created. 
        stim    = imbinarize(stim,.5);
    end
end

if Shuffle
    % Detect which ones are empty at the beginning and the end
    indx        = ~(squeeze(sum(stim,[1,2])) == 0);
    beginAndEnd = (cumsum(indx) < 1) | flip(cumsum(flip(indx))< 1);
    
    % Obtain reduced file
    stimReduced = stim(:,:,~beginAndEnd);
    
    % Shuffle
    shuff = 1:size(stimReduced,3);
    rng(shuffleSeed,'twister')
    shuff = shuff(randperm(length(shuff)));
    
    % Reconstruct
    stim(:,:,~beginAndEnd) = stimReduced(:,:,shuff);
end

if saveStimMat
    [p,f,e] = fileparts(fileName);
    if ~exist(p,'dir')
        mkdir(p)
    end
    save(fileName, 'stim');
end

if saveParamStimFile
    % Now we recreate the whole thing to be read by the program
    [fp,pn,fe] = fileparts(fileName);
    params.loadMatrix = fullfile(fp,strcat(pn,"_ParamsStims",fe));
    stimulus = stimFromParamsImages(params,stim)
    save(params.loadMatrix, 'params','stimulus'); 
end

if createVideo
% Write .avi for visualization
% TODO: 
%      - flip the images so that the words read well in the video
    % Prepare the new file.
    vidObj = VideoWriter(fullfile(pmRootPath,'local',videoFileName));
    open(vidObj);
    % Create an animation.
    axis equal tight off
    set(gca,'nextplot','replacechildren');
    for k = 1:size(stim,3)
        image(stim(:,:,k)); colormap gray;
       % Write each frame to the file.
       currFrame = getframe(gcf);
       writeVideo(vidObj,currFrame);
    end
    % Close the file.
    close(vidObj);
end


end




%% Function to create stimulus struct from params and validate with images
function stimulus = stimFromParamsImages(params,justimages)
    images      = cell(1,1);
    images{1}   = justimages; 
    
    stimulus.images = images;
    stimulus.seq = [1:size(justimages,3)]';
    stimulus.seqtiming = (stimulus.seq * params.tr) - params.tr;
    stimulus.cmap      = params.display.gammaTable;
    N = params.numImages;
    Changes = 5;
    howManyVols = ceil(N/(Changes*2));
    stimulus.fixSeq    = repmat([ones(howManyVols,1);2*ones(howManyVols,1)],[Changes,1]);
    if size(stimulus.fixSeq,1) ~= N
        stimulus.fixSeq = stimulus.fixSeq(1:N);
    end
end
