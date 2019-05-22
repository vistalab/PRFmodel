%% wrapper script to call <call_retinotopyKnk>
clear all; close all; clc; 

%% modify here

% subject ID
subID = '5944'; 

% directory where we have the spanish word matrices
words  = fullfile(pmRootPath,'data','words','spanish.mat');
load(words,'images');

% Mask images and fixation are here
masks = fullfile(pmRootPath,'data','images','maskimages.mat');
load(masks,'maskimages');

% loads the variable spatialoverlay
fixation = fullfile(pmRootPath,'data','images','fixationgrid.mat');
a1 = load(fixation);

% provide the dir with the simuli, log is written here. 
stimulusdir = fullfile(pmRootPath,'data');

% background color
grayval = uint8(127); 

% whether we want the fixation grid. 0 = no, 1 = yes
fixationGrid = 1;

% stimulus starts when this key is detected
triggerkey = 's';          

% full path of image matrix
% pathImages = fullfile(prfModelPath,wordsDir,imagesName);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expnum      = 103;
runnum      = 3; 
% pathMasks = fullfile('./', 'maskimages.mat');

% how many frames we want a refresh to last. default: 15. previously: 4. 
% the lower the faster
frameduration = 4;

% directories
% path to directory that contains the stimulus .mat files
% stimulusdir = fullfile('./');

%tweaking (in the cni)
dres    = -.9;              % <dres> (optional) is
                            %   [A B] where this is the desired resolution to imresize the images to (using bicubic interpolation).
                            %     if supplied, the imresize takes place immediately after loading the images in, and this imresized 
                            %     version is what is cached in the output of this function.
                            %  -C where C is the <scfactor> input in ptviewmovie.m
                            %   default is [] which means don't do anything special.  
                            % seems to crash for values >= 1
                            

offset  = [0 -35];        % [X Y] where X and Y are the horizontal and vertical
                           % offsets to apply.  for example, [5 -10] means shift 
                           % 5 pixels to right, shift 10 pixels up.
     
                           
% dres    = [];
% offset  = [0 0];

%%  There is one file that contains many different mask images
% display. currently not sure when/whether we should specify or define it as empty
% ptres = [1024 768 60 32]; 
% OJO CON ESTO: YO ME ASEGURARIA DE PONER LO NUESTRO BIEN, COPIAR DE MI SCRIPT THE MINI
ptres = [];  % display resolution. [] means to use current display resolution.



% fixation dot
fixationinfo = {uint8([255 0 0; 0 0 0; 0 255 0]) 0.5};  % dot colors and alpha value
fixationsize = 10;         % dot size in pixels
meanchange = 3;            % dot changes occur with this average interval (in seconds)
changeplusminus = 2;       % plus or minus this amount (in seconds)

% trigger
tfun = @() fprintf('STIMULUS STARTED.\n');  % function to call once trigger is detected

% tweaking
% offset = [0 0];
movieflip = [0 0];         % [A B] where A==1 means to flip vertical dimension
                           % and B==1 means to flip horizontal dimension

%%%%%%%%%%%%%%%%%%%%%%%%%% DO NOT EDIT BELOW

% set rand state
rand('state',sum(100*clock));
randn('state',sum(100*clock));

% ask the user what to run
% if ~exist('subjnum','var') 
%   subjnum = input('What is the subj id? ','s')
% end
% expnum = input('What experiment (89=CCW, 90=CW, 91=expand, 92=contract, 93=multibar, 94=wedgeringmash)? ')
% runnum = input('What run number (for filename)? ')

% prepare inputs
trialparams = [];
ptonparams = {ptres,[],0};
iscolor = 1;
soafun = @() round(meanchange*(60/frameduration) + changeplusminus*(2*(rand-.5))*(60/frameduration));

% assign depending on whether or not we want fixation grid
if fixationGrid
    g = a1.specialoverlay;
    gscale = imresize(g, -dres);
    gridImage = gscale;
else
    gridImage = [];
end

% some prep
if ~exist('images','var')
  images = [];
  maskimages = [];
end
% filename = sprintf('%s_subj%d_run%02d_exp%02d.mat',gettimestring,subID,runnum,expnum);
filename = sprintf('%s_subj%s_run%i_exp%i.mat',gettimestring,subID,runnum,expnum);



[images,maskimages,frameorder] = ...
  pmShowmulticlass(filename,offset,movieflip,frameduration,fixationinfo,...
                 fixationsize,tfun, ...
                 ptonparams,soafun,0,images,expnum,[],grayval,iscolor,...
                 [],[],[],dres,triggerkey, ...
                 [],trialparams,[],maskimages,gridImage,stimulusdir);  
             
             
% Write .avi for visualization
% Read the mask images (the apertures, they saw words inside the apertures)

%     params{ns}.scanDuration = 300;
%     params{ns}.totalScans   = 165;
%     params{ns}.tr           = 1.82;
%     config{ns}.config.rescaleStimuliImages = true;
%     config{ns}.config.imageSideSize = 100;  
    
    

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
    
    
    % If working only with masks maybe we want to resize and binarize
    wantResize = 0
    if wantResize
        temp = zeros(config{ns}.config.imageSideSize, config{ns}.config.imageSideSize, size(stimFileBcbl,3));
        for p=1:size(trAdjImages{ns}, 3)
            temp(:,:,p) = imresize(trAdjImages{ns}(:,:,p),[config{ns}.config.imageSideSize config{ns}.config.imageSideSize],'cubic');
        end
        trAdjImages{ns} = temp;

        % ensure that all values are between 0 and 1
        trAdjImages{ns}(trAdjImages{ns} < 0) = 0;
        trAdjImages{ns}(trAdjImages{ns} > 1) = 1;
    end
    

% We can just work with masks or create the final stimuli multplying the mask
% with the actual images
stim = trAdjImages .* trAdjMaskImages;
% Check:
% figure(6);image(stim(:,:,1000));colormap gray;axis equal tight off;

% TODO: 
% flip the images so that the words read well

% Prepare the new file.
vidObj = VideoWriter(fullfile(pmRootPath,'local','stimuliFlip.avi'));
open(vidObj);
% Create an animation.
axis equal tight off
set(gca,'nextplot','replacechildren');
for k = 1:size(stim,3)
    image(flip(stim(:,:,k))); colormap gray;
   % Write each frame to the file.
   currFrame = getframe(gcf);
   writeVideo(vidObj,currFrame);
end
% Close the file.
close(vidObj);
