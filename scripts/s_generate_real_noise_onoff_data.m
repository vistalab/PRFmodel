%% Load the mean time series from the average of 4 scans to select some voxels
mn = niftiread('/Volumes/server/Projects/7tUMinn/functional/JW20110714/vistasession/Inplane/Averages/TSeries/tSeriesScan1.nii.gz');

% Compute the coherence at the stimulus frequency
keepframes = 5:132;
numcycles  = 8;
period     = length(keepframes) / numcycles;
stimidx    = numcycles + 1;
noiseidx   = numcycles+(-1:3);
MN = abs(fft(mn(:,:,:,keepframes), [], 4));

coh =  MN(:,:,:,stimidx) ./ sum(MN(:,:,:,noiseidx),4);

% check it
figure, montage(coh)

% Choose some voxels
[~, idx] = sort(coh(:), 'descend', 'MissingPlacement', 'last');
selected = idx((1:100)*50);


figure, 
subplot(2,1,1)
histogram(coh), xlim([0 1])
title('All voxels')
xlabel('Coherence'); ylabel('Counts')

subplot(2,1,2)
histogram(coh(selected)), xlim([0 1])
title('Selected voxels')
xlabel('Coherence'); ylabel('Counts')


%% Load time series from individual scans

% These are the 4 scans with checkerboards
ni{1} = '/Volumes/server/Projects/7tUMinn/functional/JW20110714/vistasession/Inplane/Original/TSeries/tSeriesScan3.nii.gz';
ni{2} = '/Volumes/server/Projects/7tUMinn/functional/JW20110714/vistasession/Inplane/Original/TSeries/tSeriesScan6.nii.gz';
ni{3} = '/Volumes/server/Projects/7tUMinn/functional/JW20110714/vistasession/Inplane/Original/TSeries/tSeriesScan9.nii.gz';
ni{4} = '/Volumes/server/Projects/7tUMinn/functional/JW20110714/vistasession/Inplane/Original/TSeries/tSeriesScan12.nii.gz';


% Load the data and store time series for selected voxels
for ii = 1:4
    data = niftiread(ni{ii});
    
    sz = size(data);
    
    data = reshape(data, [], sz(4));
    
    ts(:,:,ii) = data(selected, :);
    
end


% check it!
figure,

plot(squeeze(ts(1,:,:)), 'LineWidth', 2);     hold on;
plot(squeeze(ts(end,:,:)), 'LineWidth', 2);
set(gca, 'XTick', keepframes(1:period:end), 'XGrid', 'on')
legend({'Good voxel scan 1' 'Good voxel scan 2' 'Good voxel scan 3' 'Good voxel scan 4' ...
    'Bad voxel scan 1' 'Bad voxel scan 2' 'Bad voxel scan 3' 'Bad voxel scan 4'})
    
% save it    
save('~/Desktop/ts.mat', 'ts');

