which_session = 3;

rootpth = '/Volumes/server/Projects/7tUMinn/functional/';

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



%% Load the mean time series from the average of multiple scans to select some voxels

% Compute the coherence at the stimulus frequency

period     = length(keepframes) / numcycles;
stimidx    = numcycles + 1;
noiseidx   = numcycles+(-1:3);

mn  = niftiread(fullfile(rootpth, projpth, mnpath));
MN  = abs(fft(mn(:,:,:,keepframes), [], 4));
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


ts = [];

% Load the data and store time series for selected voxels
for ii = 1:numel(funpth)
    data = niftiread(fullfile(rootpth, projpth, funpth{ii}));
    
    sz = size(data);
    
    data = reshape(data, [], sz(4));
    
    ts(:,:,ii) = data(selected, :);
    
end


% check it!
figure,
t = (1:size(data,2)) * tr;
plot(t, squeeze(ts(1,:,:)), 'r-o', 'LineWidth', 2);     hold on;
plot(t, squeeze(ts(end,:,:)), 'b-o', 'LineWidth', 2);
set(gca, 'XTick', tr*keepframes(1:period:end), 'XGrid', 'on', 'FontSize', 16)
xlabel('Time (seconds)')
ylabel('BOLD response')
title('A good (red) and a bad (blue) voxel')
% save it    
save(sprintf('~/Desktop/ts_%d.mat', which_session), 'ts', 'keepframes', 'numcycles', 'tr');

