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
% Prepare and start always from scratch
clear all; close all; clc;
tbUse prfmodel;
% Save files to
saveTo = '/Users/glerma/gDrive/STANFORD/PROJECTS/2019_PRF_Validation_methods_(Gari)/__PUBLISH__/PAPER_SUBMISSION01/Figures/RAW';
% Three colors
rcolor = [1 0 0]*0.75;
gcolor = [0 1 0]*0.75;
bcolor = [0 0 1]*0.75;

%% Read NYU niftis 
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
CSFindx  = mean(mn,4) < meanBOLD;
for ii = 1:dim(1),for jj=1:dim(2),for kk=1:dim(3)
    if mean(mn(ii,jj,kk,:)) < 0.75*meanBOLD
        mn(ii,jj,kk,:) = zeros(1,nTimes);
    end
end,end,end

% Calculate fft, the amplitude is divided in the two halfs, so we take one half
% and multiply the amplitude by two
MN      = abs(fft(mn, [], 4));
meanMN  = squeeze(mean(MN,1:3));

% Plot the signal being compared for the coherence calculation
X = MN(:,:,:,stimidx);
Y = sum(MN(:,:,:,noiseidx),4);
coh = X ./ Y;

% Plot the max coherence voxels amplitude spectrum
% Get the frequency vector
f = (1/tr)*(0:(nTimes/2)-1)./nTimes;
% Plot spectrum of all voxels and vertical line where the expected signal is
% Plot the frequeNCY range to calculate the coherence
% THIS PLOTS THE MEAN AMPLITUDE FREQUENCY
%{
mrvNewGraphWin;plot(f(2:end),meanMN(2:nTimes/2));hold on;
plot(f(stimidx)*[1,1],[0,max(meanMN(2:nTimes/2))],'r-.', 'LineWidth', 1.5)
plot(f(noiseidx(1))*[1,1],[0,max(meanMN(2:nTimes/2))],'r-', 'LineWidth', 0.75)
plot(f(noiseidx(end))*[1,1],[0,max(meanMN(2:nTimes/2))],'r-', 'LineWidth', 0.75)
xlabel(sprintf('f [Hz] (TR=%1.1fs)',tr));set(gca,'FontSize',14)
ylabel('BOLD Amplitude')
%}
% THIS PLOTS JUST TWO VOXELS
% Obtain max and min coh, and indexes
[maxcohval,maxcohind] = max(coh(:)); [xma,yma,zma] = ind2sub(size(coh),maxcohind);
[mincohval,mincohind] = min(coh(:)); [xmi,ymi,zmi] = ind2sub(size(coh),mincohind);
mrvNewGraphWin;
v1 = plot(f(2:end),squeeze(MN(xma,yma,zma,2:nTimes/2)));hold on;
v2 = plot(f(2:end),squeeze(MN(xmi,ymi,zmi,2:nTimes/2)));
plot(f(stimidx)*[1,1],[0,3000],'K-.', 'LineWidth', 1.5)
plot(f(noiseidx(1))*[1,1],[0,3000],'k-', 'LineWidth', 0.75)
plot(f(noiseidx(end))*[1,1],[0,3000],'k-', 'LineWidth', 0.75)
xlabel(sprintf('f [Hz] (TR=%1.1fs)',tr));set(gca,'FontSize',14)
ylabel('BOLD Amplitude')
ylim([0,3000])
legend([v1,v2], {'Max local coherence voxel','Min local coherence voxel'},'location','best')
title('BOLD signal amplitude spectrum for two voxels')
fnameRoot = 'AmplitudeSpectrumOfTwoVoxels';
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.svg')),'svg');

% check the coherence with  montage
mrvNewGraphWin; montage(coh)
title('Coherence maps')
fnameRoot = 'CoherenceMapMontage';
colorbar('southoutside'); 
set(gca, 'FontSize',14)
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.svg')),'svg');


% Select voxels with coherence inside a higher band of values
mrvNewGraphWin('coherence distr','wide');
[pdf, X_values] = ksdensity(coh(:));
a       = plot(X_values,pdf,'Color','k','LineStyle','-','LineWidth',3); hold on;
[peak, maxCohLoc] = max(pdf);
mincoh  = 0.8*max(coh(:));
maxcoh  = max(coh(:));
[~,minCohii]  = min(abs(X_values - (mincoh  * ones(size(X_values)))));
[~,maxCohii]  = min(abs(X_values - (maxcoh  * ones(size(X_values)))));
hmin    = plot(mincoh*[1,1],[0 pdf(minCohii)],'Color','k','LineStyle','-.','LineWidth',2);hold on;
tmin    = text(mincoh,pdf(minCohii)+0.11,sprintf('Coh.:%0.3f',mincoh),'HorizontalAlignment','center','FontSize',14);
hmax    = plot(maxcoh*[1,1],[0 pdf(maxCohii)],'Color','k','LineStyle','-.','LineWidth',2);
tmax    = text(maxcoh,pdf(maxCohii)+0.15,sprintf('Coh.:%0.3f',maxcoh),'HorizontalAlignment','center','FontSize',14);
jbfill(X_values(minCohii:maxCohii), pdf(minCohii:maxCohii), 0*pdf(minCohii:maxCohii), 'k','k',1,0.4)
xlabel('Coherence Values'); ylabel('Probability density'); 
xlim([0,1])
title('Coherence density distribution');set(gca,'FontSize',16)
fnameRoot = 'Coherence_density_distribution';
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.svg')),'svg');

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

% Convert BOLD signal to contrast:
tsrealNorm   = (ts-mean(ts,2)) ./ mean(ts,2);
tsmnrealNorm = (tsmn-mean(tsmn,2)) ./ mean(tsmn,2);

% Calcualte the contrast value as well
tscontr      = ((max(tsmnrealNorm') -  min(tsmnrealNorm'))/2)'; 
% Select voxels with contrast inside a higher band of values
mrvNewGraphWin('ContrastDistr','wide');
[pdf, X_values] = ksdensity(tscontr(:));
a       = plot(X_values,pdf,'Color','k','LineStyle','-','LineWidth',3); hold on;
[peak, maxCohLoc] = max(pdf);
mincontr  = 0.08;
maxcontr  = 0.12;
[~,minContrii]  = min(abs(X_values - (mincontr  * ones(size(X_values)))));
[~,maxContrii]  = min(abs(X_values - (maxcontr  * ones(size(X_values)))));
hmin    = plot(mincontr*[1,1],[0 pdf(minContrii)],'Color','k','LineStyle','-.','LineWidth',2);hold on;
tmin    = text(mincontr,pdf(minContrii)+2,sprintf('Contr.:%0.3f',mincontr),'HorizontalAlignment','center','FontSize',14);
hmax    = plot(maxcontr*[1,1],[0 pdf(maxContrii)],'Color','k','LineStyle','-.','LineWidth',2);
tmax    = text(maxcontr,pdf(maxContrii)+2,sprintf('Contr.:%0.3f',maxcontr),'HorizontalAlignment','center','FontSize',14);
jbfill(X_values(minContrii:maxContrii), pdf(minContrii:maxContrii), 0*pdf(minContrii:maxContrii), 'k','k',1,0.4)
xlabel('Contrast Values'); ylabel('Probability density'); 
xlim([0,0.25])
title('Contrast density distribution');set(gca,'FontSize',16)
fnameRoot = 'Contrast_density_distribution';
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.svg')),'svg');





% figure; hist(tscontr,500);xlim([0.01,1])

% Selects all the voxels with coherence between mincoh and maxcoh
% and
% we want to have a high contrast in all voxels regardless of noise
selected = 1:length(tscoh);
selected = selected((tscoh >= mincoh) & (tscoh <= maxcoh) & ...
                    (tscontr >= mincontr) & (tscontr <= maxcontr));

% Filter only the  voxels we are interested in
ts           = ts(selected,:,:);
tsmn         = tsmn(selected, :)';
tsrealNorm   = tsrealNorm(selected,:,:);
tsmnrealNorm = tsmnrealNorm(selected,:)';
tscontr      = tscontr(selected,:)';
tscoh        = tscoh(selected,:)';

% Obtain the noise
% Fit a sinusoidal to each one of the repetitions and obtain the difference

% How many repetitions do we have in this session
nReps   = size(ts,3);
V1diff1 = zeros([length(selected),nTimes, nReps]);
for nv=1:length(selected)
    for thisRep = 1 :nReps
        oneRep                       = tsrealNorm(nv,1:nTimes,thisRep);
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
[~,maxNoiseIndex]  = min(abs(maxes - (0.1 * maxmax  * ones(size(maxes)))));
[~,midNoiseIndex]  = min(abs(maxes - (0.45 * maxmax  * ones(size(maxes)))));
[~,minNoiseIndex]  = min(abs(maxes - (0.95 * maxmax  * ones(size(maxes)))));

% Plot the noise distributions
mrvNewGraphWin('diff noise hist');
amax    = plot(XVALS(maxNoiseIndex,:),PDF(maxNoiseIndex,:),'Color',rcolor,'LineStyle','-','LineWidth',2); hold on;
amid    = plot(XVALS(midNoiseIndex,:),PDF(midNoiseIndex,:),'Color',bcolor,'LineStyle','-','LineWidth',2);
amin    = plot(XVALS(minNoiseIndex,:),PDF(minNoiseIndex,:),'Color',gcolor,'LineStyle','-','LineWidth',2);
xlabel('Relative noise values'); ylabel('Probability density'); title('');
legend({'High noise voxel','Mid noise voxel','Low noise voxel'})
set(gca,'FontSize',14)
fnameRoot = 'SelectedVoxels_NoiseDistribution_TOGETHER';
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.svg')),'svg');



% Make same plot but separated in three
mrvNewGraphWin('diff noise hist','wide');
subplot(1,3,1)
amax    = plot(XVALS(maxNoiseIndex,:),PDF(maxNoiseIndex,:),'Color',rcolor,'LineStyle','-','LineWidth',2); hold on;
xlabel('Relative noise values'); ylabel('Probability density'); 
title('Maximum noise voxel')
xlim([-0.25,0.25]); ylim([0,22]);set(gca,'FontSize',20)

subplot(1,3,2)
amid    = plot(XVALS(midNoiseIndex,:),PDF(midNoiseIndex,:),'Color',bcolor,'LineStyle','-','LineWidth',2);
xlabel('Relative noise values'); ylabel('Probability density'); 
title('Middle noise voxel')
xlim([-0.25,0.25]); ylim([0,22]);set(gca,'FontSize',20)

subplot(1,3,3)
amin    = plot(XVALS(minNoiseIndex,:),PDF(minNoiseIndex,:),'Color',gcolor,'LineStyle','-','LineWidth',2);
xlabel('Relative noise values'); ylabel('Probability density'); 
title('Lower noise voxel')
xlim([-0.25,0.25]); ylim([0,22]);set(gca,'FontSize',20)

fnameRoot = 'SelectedVoxels_NoiseDistribution';
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.svg')),'svg');


% Plot time series and noise spectrum of the new selected 3 voxels
mrvNewGraphWin('Time series and noise spectrum of selected voxels');
subplot(2,2,1)
t = (1:size(data,2)) * tr;
plot(t, squeeze(ts(maxNoiseIndex,:,:)), 'r-', 'LineWidth', 1); hold on;
plot(t, squeeze(ts(midNoiseIndex,:,:)), 'b-', 'LineWidth', 1);
plot(t, squeeze(ts(minNoiseIndex,:,:)), 'g-', 'LineWidth', 1);
plot(t, tsmn(:,[maxNoiseIndex,midNoiseIndex,minNoiseIndex]), 'k-', 'LineWidth', 2);
set(gca, 'XTick', tr*keepframes(1:period:end), 'XGrid', 'on', 'FontSize', 16)
xlabel('Time (seconds)')
ylabel('BOLD response');set(gca,'FontSize',14)
title('BOLD time series')
% legend({'High noise voxel','Mid noise voxel','Low noise voxel'},'location','best')

subplot(2,2,2)
t = (1:size(data,2)) * tr;
plot(t, tsmnrealNorm(:,maxNoiseIndex), '-','color',rcolor, 'LineWidth', 2); hold on;
[~,yhat] = pm_fitSinusoidal(tsmnrealNorm(:,maxNoiseIndex), period);
plot(t, yhat, '-.','color',rcolor, 'LineWidth', 1); 

plot(t, tsmnrealNorm(:,midNoiseIndex), '-','color',bcolor, 'LineWidth', 2);
[~,yhat] = pm_fitSinusoidal(tsmnrealNorm(:,midNoiseIndex), period);
plot(t, yhat, '-.','color',bcolor, 'LineWidth', 1);


plot(t, tsmnrealNorm(:,minNoiseIndex), '-','color',gcolor, 'LineWidth', 2);
[~,yhat] = pm_fitSinusoidal(tsmnrealNorm(:,minNoiseIndex), period);
plot(t, yhat, '-.','color',gcolor, 'LineWidth', 1);

set(gca, 'XTick', tr*keepframes(1:period:end), 'XGrid', 'on', 'FontSize', 16)
xlabel('Time (seconds)')
ylabel('BOLD response');set(gca,'FontSize',14)
title('contrast, signal percent change')
text(1,0.09,sprintf('Low (coh.:%0.2f, contr: %0.2f), Mid (coh.:%0.2f, contr: %0.2f), High (coh.:%0.2f, contr: %0.2f)',...
                         tscoh(minNoiseIndex),tscontr(minNoiseIndex),...
                         tscoh(midNoiseIndex),tscontr(midNoiseIndex),...
                         tscoh(maxNoiseIndex),tscontr(maxNoiseIndex)))

% To calculate the fft, use the first repetition time series always
subplot(2,2,3)
realAFmin          = abs(fft(squeeze(ts(maxNoiseIndex,:,1)))/nTimes);
realAFmin(2:end-1) = 2 * realAFmin(2:end-1);
realAFmin          = realAFmin(1:(nTimes/2))';

realAFmid          = abs(fft(squeeze(ts(midNoiseIndex,:,1)))/nTimes);
realAFmid(2:end-1) = 2 * realAFmid(2:end-1);
realAFmid          = realAFmid(1:(nTimes/2))';

realAFmax          = abs(fft(squeeze(ts(minNoiseIndex,:,1)))/nTimes);
realAFmax(2:end-1) = 2 * realAFmax(2:end-1);
realAFmax          = realAFmax(1:(nTimes/2))';

frequency            = ((1/tr)*(0:(nTimes/2)-1)/nTimes)';


plot(frequency(2:end),realAFmin(2:end),'r');hold on;
plot(frequency(2:end),realAFmid(2:end),'b');
plot(frequency(2:end),realAFmax(2:end),'g');
legend({'High noise voxel','Mid noise voxel','Low noise voxel'})
xlabel('Hz'); ylabel('relative amplitude');set(gca,'FontSize',14)
title('amplitude spectrum of the 3 selected voxels')

subplot(2,2,4)
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

plot(frequency(2:end),realFallmin(2:end),'r');hold on;
plot(frequency(2:end),realFallmid(2:end),'b');
plot(frequency(2:end),realFallmax(2:end),'g');
legend({'High noise voxel','Mid noise voxel','Low noise voxel'})
xlabel('Hz'); ylabel('relative amplitude');set(gca,'FontSize',14)
title('noise spectrum of the 3 selected voxels')

% Save just the second subplot but in three separated voxels
mrvNewGraphWin('Three voxels and their sinusoidal','wide');
t = (1:size(data,2)) * tr;
subplot(1,3,1)
plot(t, tsmnrealNorm(:,minNoiseIndex), '-','color',gcolor, 'LineWidth', 2);hold on;
[~,yhat] = pm_fitSinusoidal(tsmnrealNorm(:,minNoiseIndex), period);
plot(t, yhat, '-','color','k', 'LineWidth', 1);
ylim([-0.1,0.1])
ylabel('Signal percent change');
xlabel('Time [sec]');set(gca, 'XGrid', 'off', 'FontSize', 16)
title('Low noise voxel')

subplot(1,3,2)
plot(t, tsmnrealNorm(:,midNoiseIndex), '-','color',bcolor, 'LineWidth', 2);hold on;
[~,yhat] = pm_fitSinusoidal(tsmnrealNorm(:,midNoiseIndex), period);
plot(t, yhat, '-','color','k', 'LineWidth', 1);
ylim([-0.1,0.1])
xlabel('Time [sec]');set(gca, 'XGrid', 'off', 'FontSize', 16)
title('Middle noise voxel')

subplot(1,3,3)
plot(t, tsmnrealNorm(:,maxNoiseIndex), '-','color',rcolor, 'LineWidth', 2); hold on;
[~,yhat] = pm_fitSinusoidal(tsmnrealNorm(:,maxNoiseIndex), period);
plot(t, yhat, '-','color','k', 'LineWidth', 1); 
ylim([-0.1,0.1])
xlabel('Time [sec]');set(gca, 'XGrid', 'off', 'FontSize', 16)
title('High noise voxel')

fnameRoot = 'Three_voxels_and_their_sinusoidal';
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.svg')),'svg');


%% Fit the parameters in pmNoise until we have similar values
% It is already a signal percent change
% No need to create noise free signal and substract. Our noise is already fMRI
% BOLD signal percent change
% Create generic
pm                       = prfModel;
pm.TR                    = 1.5;
pm.BOLDcontrast          = 10;
pm.Noise.jitter          = [0,0]; % [freq, ampl]
pm.signalPercentage      = 'spc';
pm.Stimulus.durationSecs = 250;
pm.compute;

% PLOT IT
sfrequency      = ((1/pm.TR)*(0:(pm.timePointsN/2)-1)/pm.timePointsN)';

mrvNewGraphWin('synth vs real','wide');

subplot(2,3,1)
% Obtain the max voxel
pm.Noise.setVoxelDefaults('low noise');
pm.Noise.seed   = 12345;
pm.Noise.compute;
difflow = pm.Noise.values;
Fsynth          = abs(fft(pm.Noise.values)/pm.timePointsN)';
Fsynth(2:end-1) = 2*Fsynth(2:end-1);
Fsynth          = Fsynth(1:(pm.timePointsN/2));
Fsynthmax       = Fsynth;
pm.plot('what','both','window',false,'addtext',false,'color',gcolor)
xlabel('Time [sec]'),set(gca,'FontSize',16)
ylabel('BOLD contrast')
ylim([-0.2,0.3])
title('LOW NOISE VOXEL')

subplot(2,3,2)
% Change for mid voxel
pm.Noise.setVoxelDefaults('mid');
pm.Noise.seed   = 12347;
pm.Noise.compute;
diffmid = pm.Noise.values;
Fsynth          = abs(fft(pm.Noise.values)/pm.timePointsN)';
Fsynth(2:end-1) = 2*Fsynth(2:end-1);
Fsynth          = Fsynth(1:(pm.timePointsN/2));
Fsynthmid       = Fsynth;
pm.plot('what','both','window',false,'addtext',false,'color',bcolor)
set(gca,'FontSize',16);
ylim([-0.2,0.3])
xlabel('Time [sec]'), 
title('MID NOISE VOXEL')

subplot(2,3,3)
% Change for min voxel
pm.Noise.setVoxelDefaults('high noise');
pm.Noise.seed   = 12348;
pm.Noise.compute;
diffhigh = pm.Noise.values;
Fsynth          = abs(fft(pm.Noise.values)/pm.timePointsN)';
Fsynth(2:end-1) = 2*Fsynth(2:end-1);
Fsynth          = Fsynth(1:(pm.timePointsN/2));
Fsynthmin       = Fsynth;
pm.plot('what','both','window',false,'addtext',false,'color',rcolor);
ylim([-0.2,0.3])
xlabel('Time [sec]'), set(gca,'FontSize',16)
title('HIGH NOISE VOXEL')



subplot(2,3,4)
plot(sfrequency(2:end,:),Fsynthmax(2:end,:),'-','color',gcolor);hold on;
% plot(frequency(2:end,:),realFallmax(2:end,:),'-.','color',[0.8 1 0.8],'linewidth',2);  
% legend({'Low noise (synth)','Low noise (real)'},'location','best')
ylabel('Noise amplitude spectrum')
xlabel('f [Hz]'); set(gca,'FontSize',14)
ylim([0,0.03]); xlim([0,0.33])
% title('Noise spectrum of the real and synth voxels')

subplot(2,3,5)
plot(sfrequency(2:end,:),Fsynthmid(2:end,:),'-','color',bcolor);hold on;
% plot(frequency(2:end,:),realFallmid(2:end,:),'-.','color',[0.8 0.8 1],'linewidth',2);  
% legend({'Mid noise (synth)','Mid noise (real)'},'location','best')
xlabel('f [Hz]'); set(gca,'FontSize',14)
ylim([0,0.03]); xlim([0,0.33])

subplot(2,3,6)
plot(sfrequency(2:end,:),Fsynthmin(2:end,:),'-','color',rcolor);hold on;
% plot(frequency(2:end,:),realFallmin(2:end,:),'-.','color',[1 0.8 0.8],'linewidth',2);  
% legend({'High noise (synth)','High noise (real)'},'location','best')
xlabel('f [Hz]'); 
set(gca,'FontSize',14)
ylim([0,0.03]); xlim([0,0.33])
fnameRoot = 'SyntheticVoxelTimeSeriesAndSpectrum';
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.svg')),'svg');


% Plot the noise distribuion now
[pdflow,xvalslow] = ksdensity(difflow);
[pdfmid,xvalsmid] = ksdensity(diffmid);
[pdfhigh,xvalshigh] = ksdensity(diffhigh);

% Plot the noise distributions
mrvNewGraphWin('diff noise hist synth');
amax    = plot(xvalslow,pdflow,'Color',gcolor,'LineStyle','-','LineWidth',2); hold on;
amid    = plot(xvalsmid,pdfmid,'Color',bcolor,'LineStyle','-','LineWidth',2);
amin    = plot(xvalshigh,pdfhigh,'Color',rcolor,'LineStyle','-','LineWidth',2);
xlabel('Relative noise values'); ylabel('Probability density'); title('');
% legend({'High noise voxel','Mid noise voxel','Low noise voxel'})
set(gca,'FontSize',14)
fnameRoot = 'SyntheticVoxels_NoiseDistribution_TOGETHER';
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.svg')),'svg');

%% Create the fig for sup mat asked by jon
A = load('/Users/glerma/gDrive/STANFORD/PROJECTS/2019_PRF_Validation_methods_(Gari)/__PUBLISH__/PAPER_SUBMISSION01/Figures/RAW/sub-paper_ses-sess01-prf_acq-normal_run-01_bold.mat');
kk = mrvNewGraphWin('NoiselessCloudPoints3sizes','wide');
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  0.8  0.9]);

% Apply params to all
nslvl  = 'low';
addcihist = true;

% Plot individuals
subplot(3,4,1)
tools  = {'vista'};
useHRF = 'vista_twogammas';
pmCloudOfResults(A.compTable   , tools ,'onlyCenters',false ,'userfsize' , 0.25, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 2     ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
                ... 'xlims',[2.75, 3.25],'ylims',[2.75, 3.25], 'xtick',[2.5:.125:3.5],'ytick',[2.5:.125:3.5], ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')

subplot(3,4,2)
tools  = {'afni'};
useHRF = 'afni_spm';
% nslvl  = 'none';
pmCloudOfResults(A.compTable   , tools ,'onlyCenters',false ,'userfsize' , 0.25, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 2     ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
                 ... 'xlims',[2.75, 3.25],'ylims',[2.75, 3.25], 'xtick',[2.5:.125:3.5],'ytick',[2.5:.125:3.5], ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')

subplot(3,4,3)
tools  = {'popeye'};
useHRF = 'popeye_twogammas';
% nslvl  = 'none';
pmCloudOfResults(A.compTable   , tools ,'onlyCenters',false ,'userfsize' , 0.25, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 2     ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
                 ... 'xlims',[2.75, 3.25],'ylims',[2.75, 3.25], 'xtick',[2.5:.125:3.5],'ytick',[2.5:.125:3.5], ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')

subplot(3,4,4)
tools  = {'aprf'};
useHRF = 'canonical';
% nslvl  = 'none';
pmCloudOfResults(A.compTable   , tools ,'onlyCenters',false ,'userfsize' , 0.25, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 1.5   ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
                 ... 'xlims',[2.75, 3.25],'ylims',[2.75, 3.25], 'xtick',[2.5:.125:3.5],'ytick',[2.5:.125:3.5], ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')

subplot(3,4,5)
tools  = {'vista'};
useHRF = 'vista_twogammas';
% nslvl  = 'none';
pmCloudOfResults(A.compTable   , tools ,'onlyCenters',false ,'userfsize' , 1, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 2     ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
                ...  'xlims',[2, 4],'ylims',[2, 4], 'xtick',[2:.5:4],'ytick',[2:.5:4], ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')

subplot(3,4,6)
tools  = {'afni'};
useHRF = 'afni_spm';
% nslvl  = 'none';
pmCloudOfResults(A.compTable   , tools ,'onlyCenters',false ,'userfsize' , 1, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 2     ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
                ...  'xlims',[2, 4],'ylims',[2, 4], 'xtick',[2:.5:4],'ytick',[2:.5:4], ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')

subplot(3,4,7)
tools  = {'popeye'};
useHRF = 'popeye_twogammas';
% nslvl  = 'none';
pmCloudOfResults(A.compTable   , tools ,'onlyCenters',false ,'userfsize' , 1, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 2     ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
                ...  'xlims',[2, 4],'ylims',[2, 4], 'xtick',[2:.5:4],'ytick',[2:.5:4], ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')

subplot(3,4,8)
tools  = {'aprf'};
useHRF = 'canonical';
% nslvl  = 'none';
pmCloudOfResults(A.compTable   , tools ,'onlyCenters',false ,'userfsize' , 1, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 1.5   ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
                ...  'xlims',[2, 4],'ylims',[2, 4], 'xtick',[2:.5:4],'ytick',[2:.5:4], ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')
             
             
subplot(3,4,9)
tools  = {'vista'};
useHRF = 'vista_twogammas';
% nslvl  = 'none';
pmCloudOfResults(A.compTable   , tools ,'onlyCenters',false ,'userfsize' , 2, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 2     ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')

subplot(3,4,10)
tools  = {'afni'};
useHRF = 'afni_spm';
% nslvl  = 'none';
pmCloudOfResults(A.compTable   , tools ,'onlyCenters',false ,'userfsize' , 2, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 2     ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')

subplot(3,4,11)
tools  = {'popeye'};
useHRF = 'popeye_twogammas';
% nslvl  = 'none';
pmCloudOfResults(A.compTable   , tools ,'onlyCenters',false ,'userfsize' , 2, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 2     ,'noiselevel' ,nslvl ,'addcihist', addcihist, ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')

subplot(3,4,12)
tools  = {'aprf'};
useHRF = 'canonical';
% nslvl  = 'none';
pmCloudOfResults(A.compTable   , tools ,'onlyCenters',false ,'userfsize' , 2, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 1.5   ,'noiselevel' ,nslvl ,'addcihist', addcihist, ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')             
             
             
             
fnameRoot = 'Noisefree_accuracy_3sizes_lownoise';
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.svg')),'svg');