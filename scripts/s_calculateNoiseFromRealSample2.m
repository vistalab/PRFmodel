%% Calculate noise in a typical subject 
clear all;
if (0)
    % The data comes from the processed reading_prf/heb_pilot09/RetAndHebrewLoc/Gray subject
    % The data is now in /share/wandell/data but it will be in Flywheel too
    basedir = '/share/wandell/data/reading_prf/heb_pilot09/RetAndHebrewLoc/Gray';
    
    % N1 = niftiRead('Words_Hebrew1/TSeries/tSeriesScan1.nii.gz');
    % N2 = niftiRead('Words_Hebrew2/TSeries/tSeriesScan1.nii.gz');
    T1 = load(fullfile(basedir,'Words_Hebrew1/TSeries/Scan1/tSeries1.mat'));
    T2 = load(fullfile(basedir,'Words_Hebrew2/TSeries/Scan1/tSeries1.mat'));
    T3 = load(fullfile(basedir,'Words_English1/TSeries/Scan1/tSeries1.mat'));
    
    Tsum  = (T1.tSeries+T2.tSeries+T3.tSeries) / 3;
    diff1 = T1.tSeries - Tsum;
    signalmeanmaxmin = [mean(mean(Tsum)), mean(max(Tsum)), mean(min(Tsum)),...
                 mean(max(Tsum))-mean(min(Tsum)), mean(mean(diff1))];
    mrvNewGraphWin('Mean noise signal'); plot(mean(diff1,2))
    
    % Obtain the FFT
    F1 = fft(diff1);
    % The mean of all the voxels
    F1mean = mean(abs(F1),2);
    
    % Obtain a single voxel and try to imitate the same noise spectrum
    pm=prfModel;
    pm.BOLDmeanValue = signalmeanmaxmin(1);
    pm.BOLDcontrast  = (signalmeanmaxmin(4)/signalmeanmaxmin(1))*100;
    pm.Noise.white_noise2signal    = 0.002;
    pm.Noise.cardiac_amplitude     = 0.007;
    pm.Noise.respiratory_amplitude = 0.007;   
    pm.Noise.lowfrequ_amplitude    = 0.02;
    pm.Noise.compute;
    % pm.plot;
    % pm.Noise.plot;
    F = abs(fft(pm.Noise.values' * 5));
    plot(F(5:end-3,:),'b');hold on;
    plot(F1mean,'r');  set(gca,'yscale','log')
     % Save as data so that we can check noise levels elsewhere, we know where
    % this data is coming
    save(fullfile(pmRootPath,'data','noise','F1mean.mat'),'F1mean')
    save(fullfile(pmRootPath,'data','noise','diff1mean.mat'),'diff1mean')
    save(fullfile(pmRootPath,'data','noise','signalmeanmaxmin.mat'),'signalmeanmaxmin')
end     
%% Use another approach now. We want to generate unitless noise, and
% only in the last moment scale the noise and add it to the signal
% The data comes from the processed reading_prf/heb_pilot09/RetAndHebrewLoc/Gray subject
% The data is now in /share/wandell/data but it will be in Flywheel too
basedir = '/share/wandell/data/reading_prf/heb_pilot09/RetAndHebrewLoc/Gray';
    
T1 = load(fullfile(basedir,'Words_Hebrew1/TSeries/Scan1/tSeries1.mat'));
T2 = load(fullfile(basedir,'Words_Hebrew2/TSeries/Scan1/tSeries1.mat'));
T3 = load(fullfile(basedir,'Words_English1/TSeries/Scan1/tSeries1.mat'));

V1 = niftiRead(fullfile(pmRootPath,'data','examples','reallyGood_cube.nii.gz'));
V1 = V1.data;

normT1 = (T1.tSeries-mean(T1.tSeries)) ./ mean(T1.tSeries);
normT2 = (T2.tSeries-mean(T2.tSeries)) ./ mean(T2.tSeries);
normT3 = (T3.tSeries-mean(T3.tSeries)) ./ mean(T3.tSeries);
% Obtain mean from three signals: we asume this is the true value, and
% the noise is one of them - the true value
Tmean = (normT1 + normT2 + normT3) / 3;
% Obtain variables of interest
realBOLDmean     = mean(mean(Tmean));
realBOLDmax      = mean(max(Tmean));
realBOLDmin      = mean(min(Tmean));
realBOLDrange    = realBOLDmax - realBOLDmin;
if isclose(realBOLDmean,0,'tolerance',0.1), realBOLDcontrast = 100 * realBOLDrange;
else, error('Data not demeaned correctly, mean should be close to zero'), end
% Obtain difference from first time series
diff1 = normT1 - Tmean;
% Calculate SNR = meanSignal/stdDevNoise
% Neuroimage 2011 xxx polimeni wald physiiological noise and ...
% 
realSNR = mean(mean(T1.tSeries)) / mean(std(T1.tSeries - Tsum));
realSNR = mean(mean(T1.tSeries)) / mean(std(T1.tSeries));
mrvNewGraphWin('signal vs "true" value: mean of all voxels');plot(mean(normT1,2));hold on;plot(mean(Tmean,2))
mrvNewGraphWin('signal vs "true" value: just one random voxel');plot(normT1(:,5555));hold on;plot(Tmean(:,5555))
% Obtain the fft of the differences
F1 = fft(diff1);
% The mean of all the voxels abs (fft returns complex numbers)
F1mean = mean(abs(F1),2);


% Now do the same with our synthetic data.
pm              = prfModel;
pm.BOLDcontrast = 4;
pm.Noise.seed   = 'none';
pm.compute;
noisefree       = pm.BOLDnoise;

pm.Noise.seed      = 12345;
pm.Noise.noisemult = 1;
pm.Noise.jitter    = [0,0]; % [freq, ampl]
pm.Noise.white_noise2signal    = 0.02;
pm.Noise.cardiac_amplitude     = 0.05;
pm.Noise.respiratory_amplitude = 0.05;   
pm.Noise.lowfrequ_amplitude    = 1;

pm.compute;
withnoise       = pm.BOLDnoise;
pm.plot;hold on;
% Normalize both signals
normnf = (noisefree-mean(noisefree)) ./ mean(noisefree);
normwn = (withnoise-mean(withnoise)) ./ mean(withnoise);
% Obtain variables of interest
synthBOLDmean     = mean(mean(normwn));
synthBOLDmax      = mean(max(normwn));
synthBOLDmin      = mean(min(normwn));
synthBOLDrange    = synthBOLDmax - synthBOLDmin;
if isclose(synthBOLDmean,0,'tolerance',0.1), synthBOLDcontrast = 100 * synthBOLDrange;
else,error('Data not demeaned correctly, mean should be close to zero'),end
% Obtain difference between with noise and noise free
diffsynth = normwn - normnf;
synthSNR  = mean(withnoise) / std(withnoise-noisefree);
synthSNR  = mean(withnoise) / std(withnoise);
% Obtain the fft of the differences
Fsynth = fft(diffsynth);
% The mean of all the voxels abs (fft returns complex numbers)
F1meansynth = abs(Fsynth)';

mrvNewGraphWin('synth vs real');
plot(F1meansynth(5:end-3,:),'b');hold on;
plot(F1mean(2:end,:),'r');  
set(gca,'yscale','log')
    
    
    
%% Select high quality voxels with a common expected mean response from the retinotopy data
% We shuold find block data. Retinotopy data will be shifting in time
% because different voxels react differently

%  

V1       = niftiRead(fullfile(pmRootPath,'data','examples','reallyGood_cube.nii.gz'));
% Reshape it to make time x Nvoxels
nVoxels = V1.dim(1)*V1.dim(2)*V1.dim(3);
nTimePoints = V1.dim(4);
newShape = [nVoxels,1,1,nTimePoints];
V1       = squeeze(reshape(V1.data,newShape))';
% Plot to visualize how the time series look
mrvNewGraphWin('Real V1 data');plot(V1);
% Separate one of the time series (1st) and obtain the mean of the rest, plot it
% Select a random voxel
% thisVoxel = 9;
V1diff1 = zeros([nTimePoints, nVoxels]);
for thisVoxel=1:nVoxels
    z = ones(1,nVoxels);
    
    V1s1 = V1(:,thisVoxel);
    z(thisVoxel) = 0;
    V1s2_12  = V1(:,logical(z));
    V1mean   = mean(V1s2_12,2);
    
    % Obtain difference from the singled out time series
    V1diff1(:,thisVoxel) = V1s1 - V1mean;
    
    % mrvNewGraphWin('Real V1 data: 1 BOLD versus Mean');
    % plot(V1s1);hold on;plot(V1mean);legend({'Real V1','Mean V1 rest'});
end
% If the noise would be pointwise independent then this would be the
% noise distributio, it is about a sd (.3% of the BOLD contrast)
hist(V1diff1(:),100)

% Look for structure in the time series
realF1all = abs(fft(V1diff1));
TR        = 1;
lowestFreq = 1/(TR*(nTimePoints));
frequency = 2*(1:nTimePoints)*lowestFreq;
mrvNewGraphWin('abs fft of all real diffs');plot(frequency',realF1all);
xlabel('Hz');ylabel('relative amplitude')
% We find there is a little more low freq noise than high frequency


mrvNewGraphWin('abs fft of all real diffs');plot(frequency',mean(realF1all,2));
xlabel('Hz');ylabel('relative amplitude')
realF =mean(realF1all,2); 


%% 
mrvNewGraphWin('V1 noise');plot(V1diff1);

% Calculate SNR = meanSignal/stdDevNoise
% Neuroimage 2011 xxx polimeni wald physiological noise and ...
realSNRd = mean(V1s1) / std(V1diff1);
realSNR  = mean(V1s1) / std(V1s1);
% Obtain the fft of the differences
V1F1 = fft(V1diff1);
% The mean of all the voxels abs (fft returns complex numbers)
V1F1mean = abs(V1F1);


% Create a noise free time series
pm              = prfModel;
pm.BOLDcontrast = 8;
pm.BOLDmeanValue = 10000;
pm.Noise.seed   = 'none';
pm.compute;
noisefree       = pm.BOLDnoise;
% Add noise to the same time series
pm.Noise.seed      = 12345;
pm.Noise.noisemult = 1;
pm.Noise.jitter    = [0,0]; % [freq, ampl]
pm.Noise.white_noise2signal    = 1;
pm.Noise.cardiac_amplitude     = 0;
pm.Noise.respiratory_amplitude = 0;   
pm.Noise.lowfrequ_amplitude    = 0;
pm.compute;
withnoise       = pm.BOLDnoise;

% Normalize both signals
normnf = (noisefree-mean(noisefree)) ./ mean(noisefree);
normwn = (withnoise-mean(withnoise)) ./ mean(withnoise);
% Obtain variables of interest
synthBOLDmean     = mean(mean(normwn));
synthBOLDmax      = mean(max(normwn));
synthBOLDmin      = mean(min(normwn));
synthBOLDrange    = synthBOLDmax - synthBOLDmin;
if isclose(synthBOLDmean,0,'tolerance',0.1), synthBOLDcontrast = 100 * synthBOLDrange;
else,error('Data not demeaned correctly, mean should be close to zero'),end
% Obtain difference between with noise and noise free
diffsynth = normwn - normnf;
synthSNR  = mean(withnoise) / std(withnoise-noisefree);
synthSNR  = mean(withnoise) / std(withnoise);
% Obtain the fft of the differences
Fsynth = fft(diffsynth);
% The mean of all the voxels abs (fft returns complex numbers)
F1meansynth = abs(Fsynth)';

mrvNewGraphWin('synth vs real');
plot(F1meansynth(5:end-3,:),'b');hold on;
plot(realF(2:end,:),'r');  
set(gca,'yscale','log')
   
    

    


%% Generate synthetic data without noise
COMBINE_PARAMETERS               = [];HRF=[];Noise=[];
COMBINE_PARAMETERS.TR            = [1.5];
COMBINE_PARAMETERS.RF.Centerx0   = [-4,-3,-2,-1,0,1,2,3,4,5]; % [-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6];
COMBINE_PARAMETERS.RF.Centery0   = [-5,-4,-3,-2,-1,0,1,2,3,4];
COMBINE_PARAMETERS.RF.Theta      = [0]; %, deg2rad(45)];
COMBINE_PARAMETERS.RF.sigmaMajor = [1,4];
COMBINE_PARAMETERS.RF.sigmaMinor = 'same';
% HRF goes in groups, create 3
HRF(1).Type                      = 'canonical';
HRF(2).Type                      = 'popeye_twogammas';
HRF(3).Type                      = 'afni_spm';
COMBINE_PARAMETERS.HRF           = HRF;
% Noise goes in groups, create 1
Noise(1).noisemult               = 0;
COMBINE_PARAMETERS.Noise         = Noise;
% CREATE AND CALCULATE
synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS,'repeats',3);
no_noiseTest = pmForwardModelCalculate(synthDT);
% Save it in local
save(fullfile(pmRootPath,'local','no_noiseTest.mat'),'no_noiseTest')
    
%% Generate synthetic data with normal noise, with slow drift
    COMBINE_PARAMETERS               = [];HRF=[];Noise=[];
    COMBINE_PARAMETERS.TR            = [1.5];
    COMBINE_PARAMETERS.RF.Centerx0   = [-4,-3,-2,-1,0,1,2,3,4,5]; % [-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6];
    COMBINE_PARAMETERS.RF.Centery0   = [-5,-4,-3,-2,-1,0,1,2,3,4];
    COMBINE_PARAMETERS.RF.Theta      = [0]; %, deg2rad(45)];
    COMBINE_PARAMETERS.RF.sigmaMajor = [1,4];
    COMBINE_PARAMETERS.RF.sigmaMinor = 'same';
    % HRF goes in groups, create 3
    HRF(1).Type                    = 'canonical';
    HRF(2).Type                    = 'popeye_twogammas';
    HRF(3).Type                    = 'afni_spm';
    COMBINE_PARAMETERS.HRF         = HRF;
    % Noise goes in groups, create 1
    Noise(1).noisemult             = 1;
    Noise(1).white_noise2signal    = 0.12 ;
    Noise(1).cardiac_frequency     = 1.17;
    Noise(1).cardiac_amplitude     = 0.014 ;
    Noise(1).respiratory_frequency = 0.23;
    Noise(1).respiratory_amplitude = 0.014 ;
    Noise(1).lowfrequ_frequ        = 120;
    Noise(1).lowfrequ_amplitude    = 0.23 ;
    COMBINE_PARAMETERS.Noise       = Noise;
    % CREATE AND CALCULATE
    synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS,'repeats',3);
    same_noiseTest_jit0 = pmForwardModelCalculate(synthDT);
    % Save it in local
    save(fullfile(pmRootPath,'local','same_noiseTest_jit0.mat'),'same_noiseTest_jit0')    
      
%% Generate synthetic data with half the noise
%{
    COMBINE_PARAMETERS               = [];HRF=[];Noise=[];
    COMBINE_PARAMETERS.TR            = [1.5];
    COMBINE_PARAMETERS.RF.Centerx0   = [-4,-3,-2,-1,0,1,2,3,4,5]; % [-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6];
    COMBINE_PARAMETERS.RF.Centery0   = [-5,-4,-3,-2,-1,0,1,2,3,4];
    COMBINE_PARAMETERS.RF.Theta      = [0]; %, deg2rad(45)];
    COMBINE_PARAMETERS.RF.sigmaMajor = [1,4];
    COMBINE_PARAMETERS.RF.sigmaMinor = 'same';
    % HRF goes in groups, create 3
    HRF(1).Type                    = 'canonical';
    HRF(2).Type                    = 'popeye_twogammas';
    HRF(3).Type                    = 'afni_spm';
    COMBINE_PARAMETERS.HRF         = HRF;
    % Noise goes in groups, create 1
    Noise(1).noisemult             = 0.5;
    COMBINE_PARAMETERS.Noise       = Noise;
    % CREATE AND CALCULATE
    synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS,'repeats',3);
    half_noiseTest_jit0 = pmForwardModelCalculate(synthDT);
    % Save it in local
    save(fullfile(pmRootPath,'local','half_noiseTest_jit0.mat'),'half_noiseTest_jit0')
    %}
%% Generate synthetic data with double noise
%{
    COMBINE_PARAMETERS               = [];HRF=[];Noise=[];
    COMBINE_PARAMETERS.TR            = [1.5];
    COMBINE_PARAMETERS.RF.Centerx0   = [-4,-3,-2,-1,0,1,2,3,4,5]; % [-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6];
    COMBINE_PARAMETERS.RF.Centery0   = [-5,-4,-3,-2,-1,0,1,2,3,4];
    COMBINE_PARAMETERS.RF.Theta      = [0]; %, deg2rad(45)];
    COMBINE_PARAMETERS.RF.sigmaMajor = [1,4];
    COMBINE_PARAMETERS.RF.sigmaMinor = 'same';
    % HRF goes in groups, create 3
    HRF(1).Type                    = 'canonical';
    HRF(2).Type                    = 'popeye_twogammas';
    HRF(3).Type                    = 'afni_spm';
    COMBINE_PARAMETERS.HRF         = HRF;
    % Noise goes in groups, create 1
    Noise(1).noisemult             = 2;
    COMBINE_PARAMETERS.Noise       = Noise;
    % CREATE AND CALCULATE
    synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS,'repeats',3);
    double_noiseTest_jit0 = pmForwardModelCalculate(synthDT);
    % Save it in local
    save(fullfile(pmRootPath,'local','double_noiseTest_jit0.mat'),'double_noiseTest_jit0')
      %} 
    %% Compare synthetic to the real data
    % Load real
    load(fullfile(pmRootPath,'data','noise','F1mean.mat'),'F1mean')
    % load(fullfile(pmRootPath,'data','noise','diff1mean.mat'),'diff1mean')
    % load(fullfile(pmRootPath,'data','noise','signalmeanmaxmin.mat'),'signalmeanmaxmin')
    % 
    % % Load synthetic
    % % Noise free
    % load(fullfile(pmRootPath,'local','nonoise_Test1.mat'))
    % % With noise
    % load(fullfile(pmRootPath,'local','noiseTest4.mat'))

    % Generate matrices
    % Noise free
    PMnf       = no_noiseTest.pm;
    nfdata     = repmat(nan([1,PMnf(1).timePointsN]), [length(PMnf),1]);
    for ii=1:length(PMnf); nfdata(ii,:)=PMnf(ii).BOLDnoise;end
    
    % With noise: same
    PMs        = same_noiseTest_jit0.pm;
    sdata      = repmat(nan([1,PMs(1).timePointsN]), [length(PMs),1]);
    for ii=1:length(PMs); sdata(ii,:)=PMs(ii).BOLDnoise;end
    % Generate same values
    diff_s = sdata - nfdata;  
    Fs     = fft(diff_s');
    same_Fsmean = mean(abs(Fs(5:end-3,:)),2);

     % With noise: half
    PMs        = half_noiseTest_jit0.pm;
    sdata      = repmat(nan([1,PMs(1).timePointsN]), [length(PMs),1]);
    for ii=1:length(PMs); sdata(ii,:)=PMs(ii).BOLDnoise;end
    % Generate same values
    diff_s = sdata - nfdata;  
    Fs     = fft(diff_s');
    half_Fsmean = mean(abs(Fs(5:end-3,:)),2);
    
     % With noise: double
    PMs        = double_noiseTest_jit0.pm;
    sdata      = repmat(nan([1,PMs(1).timePointsN]), [length(PMs),1]);
    for ii=1:length(PMs); sdata(ii,:)=PMs(ii).BOLDnoise;end
    % Generate same values
    diff_s = sdata - nfdata;  
    Fs     = fft(diff_s');
    double_Fsmean = mean(abs(Fs(5:end-3,:)),2);   
    
    % Plot it
    mrvNewGraphWin('F1mean');plot(F1mean);hold on;
                             plot(same_Fsmean);
    legend({'human','synth'})
    mrvNewGraphWin('F1mean');plot(F1mean);hold on;
                             plot(same_Fsmean);
                             plot(half_Fsmean);
                             plot(double_Fsmean);
    legend({'human','synth','synth/2','synth*2'})
    title('Real noise spectrum vs synthetic noise (with jitter 0.0)')
