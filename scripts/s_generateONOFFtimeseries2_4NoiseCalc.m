%% Read on-off stimuli to calculate real human noise values and select synth noise params
% Steps in this script:
% 1/ Read real time series, with on-off stimuli. If not in /local, it will
%    download from OSF
% 2/ Select three representative voxels based on coherence and noise levels
%    Coherence only is not enough because we can hace high coherence with low
%    and high noise. We will select a set of voxels with high coherence (so,
%    from the visual cortex), and then using noise distributions we will select
%    a high, low and mid noise voxel.
% 3/ Try to mimic the noise spectrum of these three voxels in our syntetic signal
% 

%% Read NYU niftis 
clear all; close all; clc;
tbUse prfmodel;
URLToOSFNoiseTests = 'https://files.osf.io/v1/resources/2dusf/providers/osfstorage/5dc5d6787f37e3000ecb477a?action=download&version=1&direct';
% rootpth = '/Volumes/server/Projects/7tUMinn/functional/';  % Old location in NYU server
rootpth = fullfile(pmRootPath, 'local','noisedata');
if ~exist(rootpth,'dir')
   websave(fullfile(pmRootPath,'local','noisedata.zip'),URLToOSFNoiseTests);
   unzip(fullfile(pmRootPath,'local','noisedata.zip'), fullfile(pmRootPath,'local'));
   delete(fullfile(pmRootPath,'local','noisedata.zip'));
end

which_session = 2;

switch which_session
    case 1
        projpth = 'KK20110623/vistasession/Inplane/';
        mnpath =  'Averages/TSeries/tSeriesScan2.nii.gz';
        funpth =  {'Original/TSeries/tSeriesScan5.nii.gz'};
        
        keepframes = 7:86;
        numcycles  = 5;
        tr         = 1.5;
    case 2
        projpth = 'JW20110714/vistasession/Inplane/';
        mnpath =  'Averages/TSeries/tSeriesScan1.nii.gz';
        funpth =  {'Original/TSeries/tSeriesScan3.nii.gz' ...
            'Original/TSeries/tSeriesScan6.nii.gz' ...
            'Original/TSeries/tSeriesScan9.nii.gz' ...
            'Original/TSeries/tSeriesScan12.nii.gz' };
        
        keepframes = 5:132;
        numcycles  = 8;
        tr         = 1.5;
    case 3
        projpth = 'JW20110526/vistasession/Inplane/';
        mnpath =  'Averages/TSeries/tSeriesScan1.nii.gz';
        funpth =  {'Original/TSeries/tSeriesScan1.nii.gz' ...
            'Original/TSeries/tSeriesScan2.nii.gz' ...
            'Original/TSeries/tSeriesScan3.nii.gz' ...
            'Original/TSeries/tSeriesScan4.nii.gz' };
        
        keepframes = 7:66;
        numcycles  = 5;
        tr         = 2;
end

%% Prepare data and calculate coherence

% Compute the coherence at the stimulus frequency
period   = length(keepframes) / numcycles;
stimidx  = numcycles + 1;
noiseidx = numcycles+(-1:3);

mn       = niftiread(fullfile(rootpth, projpth, mnpath));
mn       = mn(:,:,:,keepframes);
dim      = size(mn);
nVoxels  = dim(1)*dim(2)*dim(3);
nTimes   = dim(4);
% Make voxels that are not on the brain NaN. Obtain the mean BOLD and make 0 all
% the voxels that are below that (all the zeros will make the mean much smaller)
meanBOLD = mean(mn(:));
CSFindx  = mean(mn,4)<meanBOLD;
for ii = 1:dim(1),for jj=1:dim(2),for kk=1:dim(3)
    if mean(mn(ii,jj,kk,:)) < 0.75*meanBOLD
        mn(ii,jj,kk,:) = zeros(1,nTimes);
    end
end,end,end

% Calculate fft, the amplitude is divided in the two halfs, so we take one half
% and multiply the amplitude by two
MN      = abs(fft(mn, [], 4));
meanMN  = squeeze(mean(MN,1:3));
% Get the frequency vector
f = (1/tr)*(0:(nTimes/2)-1)./nTimes;
% Plot spectrum of all voxels and vertical line where the expected signal is
% Plot the freque range to calculate the coherence
mrvNewGraphWin;plot(f(2:end),meanMN(2:nTimes/2));hold on;
plot(f(stimidx)*[1,1],[0,max(meanMN(2:nTimes/2))],'r-.', 'LineWidth', 1.5)
plot(f(noiseidx(1))*[1,1],[0,max(meanMN(2:nTimes/2))],'r-', 'LineWidth', 0.75)
plot(f(noiseidx(end))*[1,1],[0,max(meanMN(2:nTimes/2))],'r-', 'LineWidth', 0.75)
xlabel(sprintf('f [Hz] (TR=%1.1fs)',tr));
ylabel('BOLD Amplitude')

% Plot the signal being compared for the coherence calculation
X = MN(:,:,:,stimidx);
Y = sum(MN(:,:,:,noiseidx),4);

coh = X ./ Y;

% check it
mrvNewGraphWin; montage(coh)

% Select voxels with coherence inside a higher band of values
mrvNewGraphWin;
[pdf, X_values] = ksdensity(coh(:));
a       = plot(X_values,pdf,'Color','k','LineStyle','-','LineWidth',3); hold on;
[peak, maxCohLoc] = max(pdf);
mincoh  = 0.80*max(coh(:));
maxcoh  = 0.99*max(coh(:));
[~,minCohii]  = min(abs(X_values - (mincoh  * ones(size(X_values)))));
[~,maxCohii]  = min(abs(X_values - (maxcoh  * ones(size(X_values)))));
hmin    = plot(mincoh*[1,1],[0 pdf(minCohii)],'Color','k','LineStyle','-.','LineWidth',2);hold on;
tmin    = text(mincoh,pdf(minCohii)+0.11,sprintf('Coh.:%0.3f',mincoh),'HorizontalAlignment','center','FontSize',14);
hmax    = plot(maxcoh*[1,1],[0 pdf(maxCohii)],'Color','k','LineStyle','-.','LineWidth',2);
tmax    = text(maxcoh,pdf(maxCohii)+0.15,sprintf('Coh.:%0.3f',maxcoh),'HorizontalAlignment','center','FontSize',14);
jbfill(X_values(minCohii:maxCohii), pdf(minCohii:maxCohii), 0*pdf(minCohii:maxCohii), 'k','k',1,0.4)
xlabel('Coherence Values'); ylabel('Probability density'); title('');

%% From all the coherent voxels, select 3 based on noise

% Load time series from individual scans, the mean value, and reshape it to 2D
ts = [];
% Load the data and store time series for selected voxels
for ii = 1:numel(funpth)
    % niftiread is a Matlab function, niftiRead is vistasoft
    data = niftiread(fullfile(rootpth, projpth, funpth{ii}));   
    data = data(:,:,:,keepframes);
    data = reshape(data, [], nTimes);
    ts(:,:,ii) = data;
end
% Load and reshape the mean as well
tsmn  = reshape(mn, [], nTimes);
tscoh = squeeze(reshape(coh, [], 1));

% Selects all the voxels with coherence between mincoh and maxcoh
selected = 1:length(tscoh);
selected = selected((tscoh >= mincoh) & (tscoh <= maxcoh));

% Filter only the  voxels we are interested in
ts   = ts(selected,:,:);
tsmn = tsmn(selected, :)';

% Convert them to contrast:
realNormts   = (ts-mean(ts,2)) ./ mean(ts,2);
realNormtsmn = (tsmn-mean(tsmn)) ./ mean(tsmn);

% Obtain the noise
% Fit a sinusoidal to each one of the repetitions and obtain
% the difference

% How many repetitions do we have in this session
nReps   = size(ts,3);
V1diff1 = zeros([length(selected),nTimes, nReps]);
for nv=1:length(selected)
    for thisRep = 1 :nReps
        oneRep                       = realNormts(nv,1:nTimes,thisRep);
        V1diff1(nv,1:nTimes,thisRep) = pm_fitSinusoidal(oneRep, period); 
        % To visualize the fit execute:
        % pm_fitSinusoidal(oneRep, period,'plotit',true);
    end
end 

% Calculate the noise probability distributions per all voxels
PDF   = zeros(length(selected), 100);
XVALS = zeros(length(selected), 100);
for nv=1:length(selected)
    allDiffValues = V1diff1(nv,:,:);
    [pdf,xvals] = ksdensity(allDiffValues(:));
    PDF(nv,:)   = pdf;
    XVALS(nv,:) = xvals;
end

% Now select the voxels with the higher peak (less noise), lower peak
% (highest noise) and middle peak (middle noise).
maxes  = max(PDF,[],2);  % All the max values, peaks of the distributions
maxmax = max(maxes);     % The absolute max value

% Based on the maxmax value, select low and high noise voxels
[~,maxNoiseIndex]  = min(abs(maxes - (0.15 * maxmax  * ones(size(maxes)))));
[~,midNoiseIndex]  = min(abs(maxes - (0.35 * maxmax  * ones(size(maxes)))));
[~,minNoiseIndex]  = min(abs(maxes - (0.85 * maxmax  * ones(size(maxes)))));

% Plot the noise distributions
mrvNewGraphWin('diff noise hist');
amin    = plot(XVALS(maxNoiseIndex,:),PDF(maxNoiseIndex,:),'Color','r','LineStyle','-','LineWidth',2); hold on;
amin    = plot(XVALS(midNoiseIndex,:),PDF(midNoiseIndex,:),'Color','b','LineStyle','-','LineWidth',2);
amin    = plot(XVALS(minNoiseIndex,:),PDF(minNoiseIndex,:),'Color','g','LineStyle','-','LineWidth',2);
xlabel('Relative noise values'); ylabel('Probability density'); title('');
legend({'High noise voxel','Mid noise voxel','Low noise voxel'})

% Plot time series and noise spectrum of the new selected 3 voxels
mrvNewGraphWin('Time series and noise spectrum of selected voxels');
subplot(2,1,1)
t = (1:size(data,2)) * tr;
plot(t, squeeze(ts(maxNoiseIndex,:,:)), 'r-', 'LineWidth', 1); hold on;
plot(t, squeeze(ts(midNoiseIndex,:,:)), 'b-', 'LineWidth', 1);
plot(t, squeeze(ts(minNoiseIndex,:,:)), 'g-', 'LineWidth', 1);
plot(t, tsmn(:,[maxNoiseIndex,midNoiseIndex,minNoiseIndex]), 'k-', 'LineWidth', 2);
set(gca, 'XTick', tr*keepframes(1:period:end), 'XGrid', 'on', 'FontSize', 16)
xlabel('Time (seconds)')
ylabel('BOLD response')
% legend({'High noise voxel','Mid noise voxel','Low noise voxel'},'location','best')

% To calculate the fft, use the first repetition time series always
realFallmin          = abs(fft(V1diff1(maxNoiseIndex,:,1))/nTimes);
realFallmin(2:end-1) = 2 * realFallmin(2:end-1);
realFallmin          = realFallmin(1:(nTimes/2))';

realFallmid          = abs(fft(V1diff1(midNoiseIndex,:,1))/nTimes);
realFallmid(2:end-1) = 2 * realFallmid(2:end-1);
realFallmid          = realFallmid(1:(nTimes/2))';

realFallmax          = abs(fft(V1diff1(minNoiseIndex,:,1))/nTimes);
realFallmax(2:end-1) = 2 * realFallmax(2:end-1);
realFallmax          = realFallmax(1:(nTimes/2))';

frequency            = ((1/tr)*(0:(nTimes/2)-1)/nTimes)';

subplot(2,1,2)
% mrvNewGraphWin('abs fft of all real diffs');
plot(frequency(2:end),realFallmin(2:end),'r');hold on;
plot(frequency(2:end),realFallmid(2:end),'b');
plot(frequency(2:end),realFallmax(2:end),'g');
legend({'High noise voxel','Mid noise voxel','Low noise voxel'})
xlabel('Hz'); ylabel('relative amplitude')


%% Fit the parameters in pmNoise until we have similar values
% It is already a signal percent change
% No need to create noise free signal and substract. Our noise is already fMRI
% BOLD signal percent change
% Create generic
pm                 = prfModel;
pm.TR              = tr;
pm.BOLDcontrast    = 8;
pm.Noise.seed      = 12345;
pm.Noise.jitter    = [0,0]; % [freq, ampl]

% PLOT IT
sfrequency      = ((1/pm.TR)*(0:(pm.timePointsN/2)-1)/pm.timePointsN)';

mrvNewGraphWin('synth vs real');

subplot(2,3,1)
% Change for min voxel
pm.Noise.setVoxelDefaults('high noise');
pm.Noise.compute;
Fsynth          = abs(fft(pm.Noise.values)/pm.timePointsN)';
Fsynth(2:end-1) = 2*Fsynth(2:end-1);
Fsynth          = Fsynth(1:(pm.timePointsN/2));
Fsynthmin       = Fsynth;
pm.plot('what','withnoise','window',false,'addtext',false,'color','r');hold on

subplot(2,3,2)
% Change for mid voxel
pm.Noise.setVoxelDefaults('mid noise');
pm.Noise.compute;
Fsynth          = abs(fft(pm.Noise.values)/pm.timePointsN)';
Fsynth(2:end-1) = 2*Fsynth(2:end-1);
Fsynth          = Fsynth(1:(pm.timePointsN/2));
Fsynthmid       = Fsynth;
pm.plot('what','withnoise','window',false,'addtext',false,'color','b')

subplot(2,3,3)
% Obtain the max voxel
pm.Noise.setVoxelDefaults('low noise');
pm.Noise.compute;
Fsynth          = abs(fft(pm.Noise.values)/pm.timePointsN)';
Fsynth(2:end-1) = 2*Fsynth(2:end-1);
Fsynth          = Fsynth(1:(pm.timePointsN/2));
Fsynthmax       = Fsynth;
pm.plot('what','withnoise','window',false,'addtext',false,'color','g')


subplot(2,3,4)
plot(sfrequency(2:end,:),Fsynthmin(2:end,:),'r-.');hold on;
plot(frequency(2:end,:),realFallmin(2:end,:),'r-');  
legend({'High noise (synth)','High noise (real)'},'location','best')
xlabel('Hz'); ylabel('Relative amplitude')

subplot(2,3,5)
plot(sfrequency(2:end,:),Fsynthmid(2:end,:),'b-.');hold on;
plot(frequency(2:end,:),realFallmid(2:end,:),'b-');  
legend({'Mid noise (synth)','Mid noise (real)'},'location','best')
xlabel('Hz'); ylabel('Relative amplitude')

subplot(2,3,6)
plot(sfrequency(2:end,:),Fsynthmax(2:end,:),'g-.');hold on;
plot(frequency(2:end,:),realFallmax(2:end,:),'g-');  
legend({'Low noise (synth)','Low noise (real)'},'location','best')
xlabel('Hz'); ylabel('Relative amplitude')
% title('Noise spectrum of the real and synth voxels')
