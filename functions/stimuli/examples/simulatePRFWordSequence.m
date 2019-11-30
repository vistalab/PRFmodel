%% MAKE A SIMULATED PRF STIMULUS SEQUENCE, WITH WORDS FLASHING ONE AT A TIME ALONG THE HORIZONTAL MERIDIAN 
% Alex L. White, November 2019
% Creates a matrix imageSeq, which can be used to simulate BOLD data from stimulated PRFs.
% (At least, Gari can use it to simulate PRFs.) 
% Specifically, this script simulates an experiment with words flashing one at a time at a range
% of locations along the horizontal meridian. 
% The "scan" is divided into "blocks" of 15 "trials." Each trial is 2s long,
% with a word flashing briefly at a random location at the beginning of the trial. 
% After each "block", there is a 12s blank period. 
%
% imageSeq is a [H x W x T] matrix, where H is the image height in pixels;
% W is the image width in pixels, and T is the total number of TRs; 
%
% This script generates two versions of imageSeq, one with the words
% appearing at 5 unique positions, and one with words appearing at 11
% unique positions (between -6 deg and 6 deg). 
% 
% I'm curious whether the PRF fits are better (less noisy) with 11
% positions, given the same total number of trials. 


clear; 

showMovie = false; %whether to play a 'movie' of the stimulus sequence

%% choices about stimuli 

imageSizePx  = 101; %number of pixels across the image (which is assumed to be square)
imageSizeDeg = 20;  %Number of degrees of visaul angle the image subtends, all the way across

pxPerDeg = imageSizePx/imageSizeDeg; %assumed pixels per degree: 

%Duration of the TR in s. Note that the time dimension in imageSeq is in units of TRs
TR = 1.4; %s

%stimulus duration in seconds, which is the rounded up to TRs
stimDurS = 0.120; % s

%trial duration in seconds
trialDurS = 2; %s

%blank duration in seconds
blankDurS = 12; %s

%number of trials in the scan
totalTrials = 180;
%and in each block 
trialsPerBlock = 15; 


%number of locations to try stimulating 
numLocations = [5 11];

 %max eccentricity to use. So for each version we will try numLocations between -maxEcc and +maxEcc
maxEcc = 6; %deg visual angle

% vertical position of the words, relative to screen center 
stimY = 0; %deg visual angle 

%dimensions of each word, in degrees visual angle
stimWidthDeg = 1.6; %deg
stimHeightDeg = 0.5; %deg


%% compute some things 

%number of simulations to make, each with a different number of unique
%positions
numNumLocs = length(numLocations);

%centerX and y of the screen, in pixels 
ctrX = round(imageSizePx/2);
ctrY = ctrX;

%stimulus dimensions in pixels
stimWidthPx = round(stimWidthDeg*pxPerDeg); % px
stimHeightPx = round(stimHeightDeg*pxPerDeg); %px

%durations in TRs
stimDur = ceil(stimDurS/TR); 
trialDur = ceil(trialDurS/TR);
blankDur = ceil(blankDurS/TR);

%scan segment durations 
numBlocks = floor(totalTrials/trialsPerBlock);
blockDur = (trialsPerBlock * trialDur) + blankDur; 
totalDur = blockDur*numBlocks; %in TRs

%% Create images 

%loop through the number of simulations, each with a different number of
%unique stimus locations 
for ni = 1:numNumLocs
    nLocs = numLocations(ni);

    %horizontal positions, in pixels, relative to left edge
    stimXs = round(pxPerDeg*linspace(-maxEcc, maxEcc, nLocs)) + ctrX;
    
    
    %Make an image of the word at each position 
    locImages = zeros(imageSizePx, imageSizePx, nLocs);
    for li=1:nLocs
        
        %word center
        wordCtrX = stimXs(li); 
        wordCtrY = round(ctrY + stimY*pxPerDeg);
       
        startX = (wordCtrX - floor(stimWidthPx/2));
        endX   = (wordCtrX + ceil(stimWidthPx/2) - 1);
        startY = (wordCtrY - floor(stimHeightPx/2));
        endY   = (wordCtrY + ceil(stimHeightPx/2) -1);
        
        locImages(startY:endY, startX:endX, li) = 1;
    end    
    
    %make full stimulus sequence 
    trialsPerLoc = ceil(totalTrials/nLocs);

    %first set trial order, which is just random 
    locOrder = repmat(1:nLocs, 1, trialsPerLoc);
    locOrder = locOrder(randperm(length(locOrder))); 
    locOrder = locOrder(1:totalTrials);
    
    %imageSeq: HxWxT
    % GLU: changing it to stim, for compatibility with prfModel
    stim = zeros(imageSizePx, imageSizePx, totalDur);
    
    %counter of which TR we're at
    trI = 1;
    
    %counter of which trial we're at
    trialNum = 0;
    
    %loop through blocks 
    for bi=1:numBlocks
        
        %loop through trials 
        for ti=1:trialsPerBlock
            trialNum = trialNum+1;
            %insert images of this trial's word at the right time
            stim(:,:,(trI:(trI+stimDur-1))) = locImages(:,:,locOrder(trialNum));
            
            %advance to next trial
            trI = trI+trialDur;
        end
        
        %advance to next block, inserting blank 
        trI = trI+blankDur;
    end
    
    %save 
    save(sprintf('imageSeq_%iLoctns.mat', nLocs), 'imageSeq');
    
    if showMovie
        for trI=1:size(stim,3)
            figure(1); clf;
            imshow(stim(:,:,trI));
            WaitSecs(TR/5); %speed it up by factor of 5
        end
    end
    
    %check by showing max across frames
    %figure(1+ni);
    %imshow(max(imageSeq,[],3));
end


    


