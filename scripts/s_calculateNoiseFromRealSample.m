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
