%% Calculate noise in a typical subject 
% The data comes from the processed reading_prf/heb_pilot09 subject
% The data is now in /share/wandell/data but it will be in Flywheel too


% N1 = niftiRead('Words_Hebrew1/TSeries/tSeriesScan1.nii.gz');
% N2 = niftiRead('Words_Hebrew2/TSeries/tSeriesScan1.nii.gz');
T1 = load('Words_Hebrew1/TSeries/Scan1/tSeries1.mat');
T2 = load('Words_Hebrew2/TSeries/Scan1/tSeries1.mat');
T3 = load('Words_English1/TSeries/Scan1/tSeries1.mat');

Tsum  = (T1.tSeries+T2.tSeries+T3.tSeries) / 3;
diff1 = T1.tSeries - Tsum;

signalmeanmaxmin = [mean(mean(Tsum)), mean(max(Tsum)), mean(min(Tsum))];

figure(99); plot(diff1(:,[300034,12332,234536]))
diff1mean = mean(diff1,2);
figure(999); plot(diff1mean)


F1 = fft(diff1);

figure(98);plot(abs(F1(2:end,[300034,12332,234536])))

F1mean = mean(abs(F1(2:end,:)),2);

figure(97);plot(F1mean)

size(F1mean)

% Save as data so that we can check noise levels elsewhere, we know where
% this data is coming
save(fullfile(pmRootPath,'data','noise','F1mean.mat'),'F1mean')
save(fullfile(pmRootPath,'data','noise','diff1mean.mat'),'diff1mean')
save(fullfile(pmRootPath,'data','noise','signalmeanmaxmin.mat'),'signalmeanmaxmin')


%% Check that we obtain the same results in the simulated data than in the real
% Generate aprox 1000 synthetic voxels combining several locations and
% parameters but using the same noise params
% WE will have 10*10*2*3*3
COMBINE_PARAMETERS               = [];

COMBINE_PARAMETERS.TR            = [1.5];
COMBINE_PARAMETERS.RF.Centerx0   = [-4,-3,-2,-1,0,1,2,3,4,5]; % [-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6];
COMBINE_PARAMETERS.RF.Centery0   = [-5,-4,-3,-2,-1,0,1,2,3,4];
COMBINE_PARAMETERS.RF.Theta      = [0]; %, deg2rad(45)];
COMBINE_PARAMETERS.RF.sigmaMajor = [1,4];
COMBINE_PARAMETERS.RF.sigmaMinor = 'same';
% HRF goes in groups, create 3
    HRF(1).Type                  = 'canonical';
    HRF(2).Type                  = 'popeye_twogammas';
    HRF(3).Type                  = 'afni_spm';
COMBINE_PARAMETERS.HRF           = HRF;
% Noise goes in groups, create 1
    Noise(1).white                 = true;
    Noise(1).white_noise2signal    = 0.1;
    Noise(1).cardiac               = true;
    Noise(1).cardiac_amplitude     = 0.1;
    Noise(1).respiratory           = true;
    Noise(1).respiratory_amplitude = 0.1;
    Noise(1).lowfrequ              = false;
COMBINE_PARAMETERS.Noise           = Noise;
% CREATE AND CALCULATE
synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS,'repeats',3);
synthDT = pmForwardModelCalculate(synthDT);
% Save the file locally just in case
save(fullfile(pmRootPath,'local','noiseTest1.mat'),'synthDT')


%% Compare the noise looks the same
% Load the data generated in the server with real values
load(fullfile(pmRootPath,'data','noise','F1mean.mat'),'F1mean')
load(fullfile(pmRootPath,'data','noise','diff1mean.mat'),'diff1mean')
load(fullfile(pmRootPath,'data','noise','signalmeanmaxmin.mat'),'signalmeanmaxmin')



