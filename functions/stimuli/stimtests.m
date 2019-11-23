[images4,maskimages4,frameorder4] = ...
  pmShowmulticlass(filename,offset,movieflip,4,fixationinfo,...
                 fixationsize,tfun, ...
                 ptonparams,soafun,0,images,expnum,[],grayval,iscolor,...
                 [],[],[],dres,triggerkey, ...
                 [],trialparams,[],maskimages,gridImage,stimulusdir); 
             

[images15,maskimages15,frameorder15] = ...
  pmShowmulticlass(filename,offset,movieflip,15,fixationinfo,...
                 fixationsize,tfun, ...
                 ptonparams,soafun,0,images,expnum,[],grayval,iscolor,...
                 [],[],[],dres,triggerkey, ...
                 [],trialparams,[],maskimages,gridImage,stimulusdir); 
             
             
             
imagesc(images(:,:,:,50))
imagesc(maskimages(:,:,1900))
imagesc(stim(:,:,1900))
imagesc(stim(:,:,1000))


isequal(images4,images15)
isequal(maskimages4,maskimages15)
isequal(frameorder4(1,:),frameorder15(1,:))  % FALSE
isequal(frameorder4(2,:),frameorder15(2,:))  % TRUE

cc = [frameorder4(1,:);frameorder15(1,:)];
cc(:,200:300)


mrvNewGraphWin;
doIt = false;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the middle bar one frame
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,1)
pm = prfModel;
pm.TR=1;
pm.Stimulus.durationSecs = 31;
pm.Stimulus.compute
% pm.Stimulus.plot('window',window)
if doIt
    % Edit it manually to make it a single on/off bar
    fName = fullfile(pmRootPath,'data','stimulus','Exp-103_bin-true_size-20x20_resize-true_Horz-101x101_barW-2_dur-31_TR-1_framedur-4.mat');
    kk = load(fName);
    % Repeat fourth bar, make rest zero, and save it with the same name, it will
    % not overwrite it.
    stim = kk.stim;
    stim(:,:,3) = repmat(kk.stim(:,:,4),[1,1,1]);
    stim(:,:,[1,2,4:end]) = zeros(size(kk.stim(:,:,[1,2,4:end])));
    % Check it
    montage(stim)
    % Save it
    save(fName,'stim');
    % Check again as part of class, it should just read it and plot what we want
    pm.Stimulus.plot
end
% Compute the time series
% pm.compute: we can't do this because the number of basis functions too low, for slow drift noise
pm.computeBOLD
pm.plot('what','nonoisetimeseries','window',window)
title('middle bar one frame')
ylim([-.1,.2]);xlim([0,35])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the bar in the middle on for 5 TRs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,2)
pm = prfModel;
pm.TR=1;
pm.Stimulus.durationSecs = 30;
pm.Stimulus.compute
% pm.Stimulus.plot('window',window)
if doIt
    % Edit it manually to make it a single on/off bar
    fName = fullfile(pmRootPath,'data','stimulus','Exp-103_bin-true_size-20x20_resize-true_Horz-101x101_barW-2_dur-30_TR-1_framedur-4.mat');
    kk = load(fName);
    % Repeat fourth bar, make rest zero, and save it with the same name, it will
    % not overwrite it.
    stim = kk.stim;
    stim(:,:,3:7) = repmat(kk.stim(:,:,4),[1,1,5]);
    stim(:,:,[1,2,8:30]) = zeros(size(kk.stim(:,:,[1,2,8:30])));
    % Check it
    montage(stim)
    % Save it
    save(fName,'stim');
    % Check again as part of class, it should just read it and plot what we want
    pm.Stimulus.plot
end
% Compute the time series
% pm.compute: we can't do this because the number of basis functions too low, for slow drift noise
pm.computeBOLD
pm.plot('what','nonoisetimeseries','window',window)
title('bar in the middle on for 5 TRs')
ylim([-.1,.2]);xlim([0,35])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the middle bar 10 frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,3)
pm = prfModel;
pm.TR=1;
pm.Stimulus.durationSecs = 32;
pm.Stimulus.compute
% pm.Stimulus.plot('window',window)
if doIt
    % Edit it manually to make it a single on/off bar
    fName = fullfile(pmRootPath,'data','stimulus','Exp-103_bin-true_size-20x20_resize-true_Horz-101x101_barW-2_dur-32_TR-1_framedur-4.mat');
    kk = load(fName);
    % Repeat fourth bar, make rest zero, and save it with the same name, it will
    % not overwrite it.
    stim = kk.stim;
    stim(:,:,3:12) = repmat(kk.stim(:,:,4),[1,1,10]);
    stim(:,:,[1,2,13:end]) = zeros(size(kk.stim(:,:,[1,2,13:end])));
    % Check it
    montage(stim)
    % Save it
    save(fName,'stim');
    % Check again as part of class, it should just read it and plot what we want
    pm.Stimulus.plot
end
% Compute the time series
% pm.compute: we can't do this because the number of basis functions too low, for slow drift noise
pm.computeBOLD
pm.plot('what','nonoisetimeseries','window',window)
title('middle bar 10 frames')
ylim([-.1,.2]);xlim([0,35])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Make the middle bar 20 frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,4)
pm = prfModel;
pm.TR=1;
pm.Stimulus.durationSecs = 40;
pm.Stimulus.compute
% pm.Stimulus.plot('window',window)
if doIt
    % Edit it manually to make it a single on/off bar
    fName = fullfile(pmRootPath,'data','stimulus','Exp-103_bin-true_size-20x20_resize-true_Horz-101x101_barW-2_dur-40_TR-1_framedur-4.mat');
    kk = load(fName);
    % Repeat fourth bar, make rest zero, and save it with the same name, it will
    % not overwrite it.
    stim = kk.stim;
    stim(:,:,3:22) = repmat(kk.stim(:,:,5),[1,1,20]);
    stim(:,:,[1,2,23:end]) = zeros(size(kk.stim(:,:,[1,2,23:end])));
    % Check it
    montage(stim)
    % Save it
    save(fName,'stim');
    % Check again as part of class, it should just read it and plot what we want
    pm.Stimulus.plot
end
% Compute the time series
% pm.compute: we can't do this because the number of basis functions too low, for slow drift noise
pm.computeBOLD
pm.plot('what','nonoisetimeseries','window',window)
title('middle bar 20 frames')
ylim([-.1,.2]);xlim([0,35])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select a side bar 5 frames
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2,3,5)
pm = prfModel;
pm.TR=1;
pm.Stimulus.durationSecs = 33;
pm.Stimulus.compute
pm.Stimulus.plot('window',window)
if doIt
    % Edit it manually to make it a single on/off bar
    fName = fullfile(pmRootPath,'data','stimulus','Exp-103_bin-true_size-20x20_resize-true_Horz-101x101_barW-2_dur-33_TR-1_framedur-4.mat');
    kk = load(fName);
    % Repeat fourth bar, make rest zero, and save it with the same name, it will
    % not overwrite it.
    stim = kk.stim;
    stim(:,:,3:7) = repmat(kk.stim(:,:,3),[1,1,5]);
    stim(:,:,[1,2,8:end]) = zeros(size(kk.stim(:,:,[1,2,8:end])));
    % Check it
    montage(stim)
    % Save it
    save(fName,'stim');
    % Check again as part of class, it should just read it and plot what we want
    pm.Stimulus.plot
end
% Compute the time series
% pm.compute: we can't do this because the number of basis functions too low, for slow drift noise
pm.computeBOLD
pm.plot('what','nonoisetimeseries','window',window)
title('side bar 5 frames')
ylim([-.1,.2]);xlim([0,35])


%% See if the time series and convolution is right
% Read the 5 bars in the middle
% close all
pmmid = prfModel;
pmmid.Stimulus.durationSecs = 30;
pmmid.Stimulus.compute
% pmmid.Stimulus.plot
pmmid.computeBOLD
pmmid.timeSeries
pmmid.plot('what','nonoisetimeseries','window',true);hold on;
ylim([-.06,.14]);xlim([0,25])


% Read the 5 bars on the side
pmside = prfModel;
pmside.Stimulus.durationSecs = 33;
pmside.Stimulus.compute
% pmside.Stimulus.plot
pmside.computeBOLD
pmside.timeSeries
pmside.plot('what','nonoisetimeseries','window',false)
ylim([-.06,.14]);xlim([0,25])



