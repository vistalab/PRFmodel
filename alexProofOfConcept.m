pm = prfModel;
% Stimulus
pm.Stimulus.expName = 'alex5';
pm.TR = 1.4;
pm.Stimulus.durationSecs = 226.8;
pm.Stimulus.compute;
pm.Stimulus.plot
% Noise
pm.Noise.seed = 'none';
% RF
pm.RF.sigmaMajor = 0.25;
pm.RF.sigmaMinor = 0.25;
pm.compute
close all
pm.plot('what','nonoise','window',true,'color','k');hold on
% RF
pm.RF.sigmaMajor = .5;
pm.RF.sigmaMinor = .5;
pm.compute
pm.plot('what','nonoise','window',false,'color','r')
pm.RF.sigmaMajor = 1;
pm.RF.sigmaMinor = 1;
pm.compute
pm.plot('what','nonoise','window',false,'color','b')
pm.RF.sigmaMajor = 2;
pm.RF.sigmaMinor = 2;
pm.compute
pm.plot('what','nonoise','window',false,'color','g')



%% ADD A LITTLE BIT OF NOISE
% Noise
pm.Noise.seed = 'random';
pm.Noise.setVoxelDefaults('low')
what = 'withnoise';
% RF
pm.RF.sigmaMajor = 0.25;
pm.RF.sigmaMinor = 0.25;
pm.compute
snr025=pm.SNR;
pm.plot('what',what,'window',true,'color','k');hold on
% RF
pm.RF.sigmaMajor = .5;
pm.RF.sigmaMinor = .5;
pm.compute
snr05=pm.SNR;
pm.plot('what',what,'window',false,'color','r')
pm.RF.sigmaMajor = 1;
pm.RF.sigmaMinor = 1;
pm.compute
snr1=pm.SNR;
pm.plot('what',what,'window',false,'color','b')
pm.RF.sigmaMajor = 2;
pm.RF.sigmaMinor = 2;
pm.compute
snr2=pm.SNR;
pm.plot('what',what,'window',false,'color','g')
legend({['RF=0.25,SNR=' num2str(snr025)],['RF=0.5,SNR=' num2str(snr05)],...
        ['RF=1,SNR=' num2str(snr1)],['RF=2,SNR=' num2str(snr2)]})


