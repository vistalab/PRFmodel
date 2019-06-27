clear all;


%% Create demo data based on PRFmodel
COMBINE_PARAMETERS.RF.Centerx0   = [-6, 6];    % etc ...
COMBINE_PARAMETERS.RF.Centery0   = [-6,-3,3,6];
COMBINE_PARAMETERS.RF.sigmaMinor = [1];
COMBINE_PARAMETERS.RF.sigmaMajor = [1];
COMBINE_PARAMETERS.TR            = 2;
synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS);
synthDT = synthDT(61:76,:);
synthDT = pmForwardModelCalculate(synthDT);


% Write the sitimulus to Nifti
%{
pm1 = synthDT.pm(1);
stimNiftiFname = pm5.Stimulus.toNifti;
% Write the BOLDnoise to Nifti
niftiBOLDfile  = pmForwardModelToNifti(synthDT, 'fname',...
    fullfile(pm5.Stimulus.LocalPath,'synthDataExample_TR2.nii.gz'));
% Check that we can read the nifti files properly
NIdata     = niftiRead(niftiBOLDfile);
data       = reshape(NIdata.data, [NIdata.dim(1)*NIdata.dim(2),NIdata.dim(4)]);
NIstimulus = niftiRead(stimNiftiFname);
stimulus   = squeeze(NIstimulus.data);
%}

% Analyze it with analyzePRF
results_analyzePRF = pmModelFit(synthDT,'analyzePRF');

% Analyze it with AFNI
results_AFNI       = pmModelFit(synthDT,'AFNI');



















%% Obtain vistasofts data, right before the pRF model is applied
% One of the first try-s was to take the tSeries and try to go back to a 3D, but
% this cannot be done, since there are repeated positions
% TODO: I wonder how is this later visualized, since we are going to have
% several voxels with prfs calculated. 

% After talking to Rosemary and Brian, we are creating a 2x2x2x144 4D nifti (since 
% the pRF model does not care about the continuity of the voxels). We will try
% to run the two AFNI pRF (Dumoulin & Wandell, 2008 -DW- and their own) models on this data, 
% to see if:
% 1.- Applying D&W algorithm in Vistasoft and AFNI, we should obtain the same
%     1.1.- If not, understand what is different: algorithm? params?
% 2.- Applying AFNI-pRF to see if we obtain the same as DW or we start getting
%     the differences they get: circle or 2:1 ellpise for DW and 6:1 for Baker.
% 3.- [...]     






% Run this on black to get the data
% cd '/sni-storage/wandell/data/reading_prf/ab/tiledLoc_sizeRet/Gray/Words/TSeries/Scan1'

% VolumeMeanMap = load('/sni-storage/wandell/data/reading_prf/ab/tiledLoc_sizeRet/Volume/Original/meanMap.mat');
% InplaneMeanMap = load('/sni-storage/wandell/data/reading_prf/ab/tiledLoc_sizeRet/Inplane/Words/meanMap.mat');
% func_xformNIfti = readFileNifti('/sni-storage/wandell/data/reading_prf/ab/tiledLoc_sizeRet/13_1_fMRI_Ret_knk/func_xform.nii.gz');
% tSeriesinGrayST = load('/sni-storage/wandell/data/reading_prf/ab/tiledLoc_sizeRet/Gray/Words/TSeries/Scan1/tSeries1');
% Make the tSeries int16 before writing the nifti
% tSeries = int16(tSeriesinGrayST.tSeries);

% Now read the reallyGood data from Rosemary
load(fullfile(pmRootPath,'data','AFNI/DATA/LV3v_reallyGood.mat'));
% Now the big data
% load(fullfile(pmRootPath,'data','AFNI/DATA/WangAtlas_V1v_co0p2_tSeries.mat'));

% tSeries = int16(tSeries');
% There is no need on converting it to int16
tSeries = tSeries';



%% Get the grey matter coordinates
% This was a test if there was a unique solution to go back to 3D. 
% There is not one as some are repeated.
% But, we can maintain this information to select voxels of interest to 
% replicate in ANTs of whatever.

%{
cd '/sni-storage/wandell/data/reading_prf/ab/tiledLoc_sizeRet';
vw = initHiddenGray;  % same as mrVista but no GUI
greyCoords = vw.coords';
% Are they unique coordinates?
[C,IA,IC] = unique(greyCoords,'rows');
if ~isequal(IA,IC);
    size(IC,1) - size(IA,1) % No, they are not...
end

% Confirmed with Brian that we cannot go back.

Path2t1acpced = '/sni-storage/wandell/biac2/wandell2/data/anatomy/bugno';
t1 = readFileNifti(fullfile(Path2t1acpced, 't1.nii.gz'));

% We need to know the order of the coordinates
size(t1.data)
size(greyCoords)
size(tSeries)
[min(greyCoords(:,1)),min(greyCoords(:,2)),min(greyCoords(:,3)); ...
 max(greyCoords(:,1)),max(greyCoords(:,2)),max(greyCoords(:,3))]

% Create an empty volume of the same size of the t1 but 4D
newVolume = int16(zeros([size(t1.data), size(tSeriesinGray,1)]));
% Loop per every volume
for ii=1:size(newVolume,4)
    % Loop per every grey matter voxel and write it
    for jj=1:size(tSeries,2)
        coords = vw.coords(:,jj)';
        newVolume([coords,ii]) = tSeries(ii,jj);
    end
end



temp4D       = t1;
temp4D.data  = newVolume;
temp4D.dim   = size(newVolume);
temp4D.fname = '/data/localhome/glerma/PROJECTS/prf_afni_test.nii.gz';
writeFileNifti(temp4D);

prfTEST = readFileNifti(temp4D.fname);

% it is not equal because we are not writting all the values
isequal(max(prfTEST.data(:)), max(tSeries(:)))
%}

% We agreed that first try will be to create a 10*10*10 nifti and see if we
% can use that.

%% Create the 2*2*3*144 file  % Rosemay gave me 12 voxels so use them all
% find max value and center around it. 
% TODO: ask Rosemay for this subject which voxels I should use based on her
% results
% [rows,cols] = find(tSeries==max(tSeries(:));
% myValues    = tSeries(:, cols-6:cols+5)';
% myCube      = reshape(myValues, [2,2,3,144]);

% I understand this is the 12 voxel one
myCube      = reshape(tSeries, [2,2,3,144]);
% And this is the big data link from above
% myCube      = reshape(tSeries, [10,10,5,144]);

% writeFileNifti(niftiCreate('data', myCube, 'fname','/data/localhome/glerma/PROJECTS/prf_afni_cube.nii.gz' ));
writeFileNifti(niftiCreate('data', myCube, 'fname',fullfile(pmRootPath,'data','AFNI','DATA/reallyGood_cube.nii.gz' )));
% writeFileNifti(niftiCreate('data', myCube, 'fname','/home/glerma/PROJECTS/PRF/AFNI/DATA/Wang_cube.nii.gz' ));

% to convert to  BRIK format:
% delete('/home/glerma/PROJECTS/PRF/AFNI/DATA/reallyGood_cube+orig.BRIK')
% delete('/home/glerma/PROJECTS/PRF/AFNI/DATA/reallyGood_cube+orig.HEAD')
% system('3dcopy /home/glerma/PROJECTS/PRF/AFNI/DATA/reallyGood_cube.nii.gz /home/glerma/PROJECTS/PRF/AFNI/DATA/reallyGood_cube+orig')
system('3dcopy /home/glerma/PROJECTS/PRF/AFNI/DATA/Wang_cube.nii.gz /home/glerma/PROJECTS/PRF/AFNI/DATA/Wang_cube+orig')

% copyfile('/home/glerma/PROJECTS/PRF/AFNI/DATA/reallyGood_cube+orig.BRIK', ...
%          '/home/glerma/PROJECTS/PRF/AFNI/ANALYSIS/DW_01_default/reallyGood_cube+orig.BRIK')
% copyfile('/home/glerma/PROJECTS/PRF/AFNI/DATA/reallyGood_cube+orig.HEAD', ...
%          '/home/glerma/PROJECTS/PRF/AFNI/ANALYSIS/DW_01_default/reallyGood_cube+orig.HEAD')
copyfile('/home/glerma/PROJECTS/PRF/AFNI/DATA/Wang_cube+orig.BRIK', ...
         '/home/glerma/PROJECTS/PRF/AFNI/ANALYSIS/DW_01_default/Wang_cube+orig.BRIK')
copyfile('/home/glerma/PROJECTS/PRF/AFNI/DATA/Wang_cube+orig.HEAD', ...
         '/home/glerma/PROJECTS/PRF/AFNI/ANALYSIS/DW_01_default/Wang_cube+orig.HEAD')     
%% Get params for the pRF (Comments are from Rosemary's email)

% Load "images_knk_fliplr.mat"
% This loads a variable called “images”.
% “images” is a 3D matrix of size 768x768x2560.
% It is a matrix of images: 2560 images that are 768x768
load('images_knk_fliplr.mat')

% Load “params_knkfull_multibar_blank.mat"
% This loads a bunch of variables, one of which is called “params"
% params.seq is a vector of length 4500 (15 frames per second for 300 seconds).
% The values of params.seq range from 1 to 2560. It tells you what image was 
% shown at what time point. 
% **Note** that the time series that you are working with has 144 time points.
% Each TR is 2 seconds, so this amounts to 288 seconds. This is less than 300 
% seconds because the first 6 frames were clipped from the data (scanner data is 
% noisy at the start of each run), so 6*2secs = 12secs

load('params_knkfull_multibar_blank.mat')

% TODO: there is no seq or no structure which is 4500 long. It is called
% stimulus.seq


% Create the 2D+time dataset of stimulus images over TRs.:
imIndx = stimulus.seq;
imIndx = imIndx(15:15:4500);
imIndx = imIndx(2:2:300);
imIndx = imIndx(7:150);
for ii=1:length(imIndx)
    thisImage = images(:,:,imIndx(ii));
    % imwrite(thisImage, sprintf('images/IMG_%03d.jpg', ii-1))
    
    
    % Do we want to downsample?
    % I do not want to interpolate, so instead of imresize do it manually
    thisImageDS = thisImage(1:6:end, 1:6:end);
    imwrite(thisImageDS, sprintf('images/DSIMG_%03d.jpg', ii-1))
    
end



% From this file I can build the parameters required. 

% Test below and come back until here we only have the last working code: 


% From Cesar, do this in bash (modified slightly)
%{
3dcopy DSIMG_000.jpg DSIMAGE_MOVIE 
export NIMAGES=143 # NUMBER OF IMAGES

for NUM in $(count -digits 3 1 $NIMAGES)
do
  3dTcat -glueto DSIMAGE_MOVIE+orig. DSIMG_${NUM}.jpg 
done
%}


% I do not think this is right:
% - change TR to 2
% - make it 2D instead of 3D
cd(fullfile(pmRootPath,'data','AFNI','ANALYSIS/DW_01_default'))
system('3drefit -TR 2 ../DATA/images/DSIMAGE_MOVIE+orig')
% system('3drefit -TR 2 reallyGood_cube+orig')
system('3drefit -TR 2 Wang_cube+orig')
system('3drefit -space orig images/DSIMAGE_MOVIE+orig')
% system('3drefit -space orig reallyGood_cube+orig')
system('3drefit -space orig Wang_cube+orig')
% See if the changes are there
system('3dinfo images/DSIMAGE_MOVIE+orig')
system('3dinfo reallyGood_cube+orig')



%% AFNI command for DW:pRF
% MAke everything full paths...

% Now follow the steps on the 3dNLfim.help file
3dDeconvolve -nodata 10 2 -polort -1                \
                   -num_stimts 1 -stim_times 1 '1D:0' GAM \
                   -x1D conv2.ref.GAM.1D


export AFNI_CONVMODEL_REF=conv.ref.GAM.1D \
export AFNI_MODEL_PRF_STIM_DSET=IMAGE_MOVIE+orig \
export AFNI_MODEL_DEBUG=3 \
export AFNI_MODEL_PRF_ON_GRID=YES \
export AFNI_MODEL_PRF_RAM_STATS=Y
    
% This is the impulse created with the function above
load('/Users/glerma/soft/PRFmodel/data/AFNI/DATA/hrfresponse.mat')
plot([1:2:20 ],   [0
 0.089639335870743
 0.89834398031235
 0.75842666625977
 0.23252660036087
 0.040924623608589
 0
 0
 0
 0
]); hold on; plot([1:2:20 ],hrfresponse(1:10),'r-');

% Save the first time points of Rosemary's HRF
save('vistasoftHRF.1D', 'hrfresponse', '-ascii')



%% Run this thing in AFNI:
% This is a sample code that it 



% The cmd below is not working in GCP, I think that its because it runs out
% of memory (it asks for 32Gb of RAM.
% Changed this param and now it is working
export AFNI_MODEL_PRF_SIGMA_NSTEPS=10
% This is check is pre-downsampling images, pre-using rosemary's prf from
% vistasoft and pre-using Rosemarys realluGood.mat data. 
cmd = (['3dNLfim -input functional_cube+orig  ' ...
             ' -noise Zero ' ...
             ' -signal Conv_PRF ' ...
             ' -sconstr 0 -10.0 10.0 ' ...
             ' -sconstr 1 -1.0 1.0 ' ...
             ' -sconstr 2 -1.0 1.0 ' ...
             ' -sconstr 3 0.0 1.0 ' ...
             ' -BOTH ' ...
             ' -nrand 10000 ' ...
             ' -nbest 5 ' ...
             ' -bucket 0 buck.PRF' ...
             ' -snfit snfit.PRF' ...
             ' -TR 2 ' ...
             ' -jobs 8'])
% This command finished. Before checking it do the same using Rosemary's
% reallyGood.mat so that she can compare results. 

% Now we will use downsampled images to 128
% and we will set the parameter back to 100
export AFNI_MODEL_PRF_SIGMA_NSTEPS=100
% And now we are using reallyGood_cube with rosemarys data, and we will use
% the same HRF they use in vistasoft
export AFNI_CONVMODEL_REF=vistasoftHRF.1D
export AFNI_MODEL_PRF_STIM_DSET=DSIMAGE_MOVIE+orig


% Check printenv before running just in case
printenv | grep AFNI 
cmd = (['3dNLfim -input reallyGood_cube+orig  ' ...
             ' -noise Zero ' ...
             ' -signal Conv_PRF ' ...
             ' -sconstr 0 -10.0 10.0 ' ...
             ' -sconstr 1 -1.0 1.0 ' ...
             ' -sconstr 2 -1.0 1.0 ' ...
             ' -sconstr 3 0.0 1.0 ' ...
             ' -BOTH ' ...
             ' -nrand 10000 ' ...
             ' -nbest 5 ' ...
             ' -bucket 0 buck.PRF' ...
             ' -snfit snfit.PRF' ...
             ' -TR 2 ' ...
             ' -jobs 12'])
         
         
% Rosemary checked the data and saw that the results are off by a factor of 30
% Change the call to this:
cmd = (['3dNLfim -input reallyGood_cube+orig  ' ...
             ' -noise Zero ' ...
             ' -signal Conv_PRF ' ...
             ' -sconstr 0 -20.0 20.0 ' ...
             ' -sconstr 1 -16.0 16.0 ' ...
             ' -sconstr 2 -16.0 16.0 ' ...
             ' -sconstr 3 0.0 16.0 ' ...
             ' -BOTH ' ...
             ' -nrand 10000 ' ...
             ' -nbest 5 ' ...
             ' -bucket 0 buck30.PRF' ...
             ' -snfit snfit30.PRF' ...
             ' -TR 2 ' ...
             ' -jobs 4'])         
         
         % The fit is worse with this one
%try with the GAM HRF now
export AFNI_CONVMODEL_REF=conv.ref.GAM.1D
cmd = (['3dNLfim -input reallyGood_cube+orig  ' ...
             ' -noise Zero ' ...
             ' -signal Conv_PRF ' ...
             ' -sconstr 0 -10.0 10.0 ' ...
             ' -sconstr 1 -1.0 1.0 ' ...
             ' -sconstr 2 -1.0 1.0 ' ...
             ' -sconstr 3 0.0 1.0 ' ...
             ' -BOTH ' ...
             ' -nrand 10000 ' ...
             ' -nbest 5 ' ...
             ' -bucket 0 buck_GAM.PRF' ...
             ' -snfit snfit_GAM.PRF' ...
             ' -TR 2 ' ...
             ' -jobs 12'])

%try with the GAM HRF now
export AFNI_CONVMODEL_REF=sisar.conv.ref.SPMG1.1D
cmd = (['3dNLfim -input reallyGood_cube+orig  ' ...
             ' -noise Zero ' ...
             ' -signal Conv_PRF ' ...
             ' -sconstr 0 -10.0 10.0 ' ...
             ' -sconstr 1 -1.0 1.0 ' ...
             ' -sconstr 2 -1.0 1.0 ' ...
             ' -sconstr 3 0.0 1.0 ' ...
             ' -BOTH ' ...
             ' -nrand 10000 ' ...
             ' -nbest 5 ' ...
             ' -bucket 0 buck_SPM1.PRF' ...
             ' -snfit snfit_SPM1.PRF' ...
             ' -TR 2 ' ...
             ' -jobs 12'])


         
% I will launch Wang data now
export AFNI_CONVMODEL_REF=vistasoftHRF.1D
cmd = (['3dNLfim -input Wang_cube+orig  ' ...
             ' -noise Zero ' ...
             ' -signal Conv_PRF ' ...
             ' -sconstr 0 -10.0 10.0 ' ...
             ' -sconstr 1 -1.0 1.0 ' ...
             ' -sconstr 2 -1.0 1.0 ' ...
             ' -sconstr 3 0.0 1.0 ' ...
             ' -BOTH ' ...
             ' -nrand 10000 ' ...
             ' -nbest 5 ' ...
             ' -bucket 0 buck_Wang.PRF' ...
             ' -snfit snfit_Wang.PRF' ...
             ' -TR 2 ' ...
             ' -jobs 10'])
% Launch everything the same, but using using the 6 params models
printenv | grep AFNI 
cmd = (['3dNLfim -input Wang_cube+orig  ' ...
             ' -noise Zero ' ...
             ' -signal Conv_PRF_6 ' ...
             ' -sconstr 0 -10.0 10.0 ' ...
             ' -sconstr 1 -1.0 1.0 ' ...
             ' -sconstr 2 -1.0 1.0 ' ...
             ' -sconstr 3 0.0 1.0 ' ...
             ' -sconstr 4 1.0 4.0 ' ...
             ' -sconstr 5 -1.571 1.570 ' ...
             ' -BOTH ' ...
             ' -nrand 10000 ' ...
             ' -nbest 5 ' ...
             ' -bucket 0 buck_Wang_6.PRF' ...
             ' -snfit snfit_Wang_6.PRF' ...
             ' -TR 2 ' ...
             ' -jobs 10'])
         
         
% We have error with the sima ratios at zero, as rick recommended, making

printenv | grep AFNI 

cmd = (['3dNLfim -input Wang_cube+orig  ' ...
             ' -noise Zero ' ...
             ' -signal Conv_PRF_6 ' ...
             ' -sconstr 0 -10.0 10.0 ' ...
             ' -sconstr 1 -1.0 1.0 ' ...
             ' -sconstr 2 -1.0 1.0 ' ...
             ' -sconstr 3 0.0 1.0 ' ...
             ' -sconstr 4 1.0 4.0 ' ...
             ' -sconstr 5 -1.571 1.570 ' ...
             ' -BOTH ' ...
             ' -nrand 10000 ' ...
             ' -nbest 5 ' ...
             ' -bucket 0 buck_Wang_6nogrid.PRF' ...
             ' -snfit snfit_Wang_6nogrid.PRF' ...
             ' -TR 2 ' ...
             ' -jobs 10'])
         
         
         
         
         
         
% With 4 we have seen that there are most of the results. Now launch it twice, 
% once constrained to 1 and constrained to 6 and see what happens. 
export AFNI_MODEL_PRF_ON_GRID=NO \
export AFNI_MODEL_DEBUG=3 \
export AFNI_MODEL_PRF_RAM_STATS=Y \
export AFNI_MODEL_PRF_ON_GRID=NO  \
export AFNI_CONVMODEL_REF=vistasoftHRF.1D \
export AFNI_MODEL_PRF_SIGMA_NSTEPS=100 \
export AFNI_MODEL_PRF_STIM_DSET=DSIMAGE_MOVIE+orig
cmd = (['3dNLfim -input Wang_cube+orig  ' ...   
             ' -noise Zero ' ...
             ' -signal Conv_PRF_6 ' ...
             ' -sconstr 0 -10.0 10.0 ' ...
             ' -sconstr 1 -1.0 1.0 ' ...
             ' -sconstr 2 -1.0 1.0 ' ...
             ' -sconstr 3 0.0 1.0 ' ...
             ' -sconstr 4 1.0 1.0 ' ...
             ' -sconstr 5 -1.571 1.570 ' ...
             ' -BOTH ' ...
             ' -nrand 10000 ' ...
             ' -nbest 5 ' ...
             ' -bucket 0 buck_Wang_6nogrid_sigRat1.PRF' ...
             ' -snfit snfit_Wang_6nogrid_sigRat1.PRF' ...
             ' -TR 2 ' ...
             ' -jobs 10'])
         
         
cmd = (['3dNLfim -input Wang_cube+orig  ' ...
             ' -noise Zero ' ...
             ' -signal Conv_PRF_6 ' ...
             ' -sconstr 0 -10.0 10.0 ' ...
             ' -sconstr 1 -1.0 1.0 ' ...
             ' -sconstr 2 -1.0 1.0 ' ...
             ' -sconstr 3 0.0 1.0 ' ...
             ' -sconstr 4 1.0 6.0 ' ...
             ' -sconstr 5 -1.571 1.570 ' ...
             ' -BOTH ' ...
             ' -nrand 10000 ' ...
             ' -nbest 5 ' ...
             ' -bucket 0 buck_Wang_6nogrid_sigRat6.PRF' ...
             ' -snfit snfit_Wang_6nogrid_sigRat6.PRF' ...
             ' -TR 2 ' ...
             ' -jobs 10'])         
         


%% Norm bvalues for FW test 
%{
bvals = dlmread('~/Downloads/103818_dwi_1000.bval');
[uniqueValues, ~, uindex] = unique(bvals);
roundedBval  = 100 * round(uniqueValues/100);              
valnorm      = roundedBval(uindex);
dlmwrite('~/Downloads/103818_dwi_1000_norm.bval', valnorm);
%}



%% HREF CESAR
% Make this launch in command line with system
%{
TR=2

echo "Create GAM"
3dDeconvolve -nodata 50 $TR -polort -1                \
                   -num_stimts 1 -stim_times 1 '1D:0' GAM \
                   -x1D sisar.conv.ref.GAM.1D

# La GAM por defecto esta normalizada a pico=1

echo "Create SPMG1"
3dDeconvolve -nodata 50 $TR -polort -1                \
                   -num_stimts 1 -stim_times 1 '1D:0' 'SPMG1(0)' \
                   -x1D sisar.conv.ref.SPMG1.1D

# El valor entre parentesis de (0) es necesario para que esté normalizada a pico=1

1dplot -one sisar.conv.ref.SPMG1.1D sisar.conv.ref.GAM.1D
%}

%% Read the results in Matlab so that we can compare to vistasoft
% First we need to see that we are doing things properly
cd(fullfile(pmRootPath, 'data', 'AFNI/ANALYSIS/DW_01_default'))
addpath(genpath('~/soft/afni_matlab'));
[err, V, Info, ErrMessage] = BrikLoad('snfit_Wang_6nogrid.PRF+orig');
% Plot the fit and the signal of voxel 2,2,2
plot(squeeze(myCube(5,5,2,:)),'k-');hold on;plot(squeeze(V(5,5,2,:)),'r-');

[berr, bV, bInfo, bErrMessage] = BrikLoad('buck_Wang_6nogrid.PRF+orig');


% Write it back to matlab so that Rosemary can check it
fitSeries = reshape(V,  [500, 144]);
% results   = reshape(bV, [500, 12]);
results   = reshape(bV, [500, 14]); % Two more params
save('fitSeries_Wang_6nogrid.mat', 'fitSeries');
save('results_Wang_6nogrid.mat', 'results');



%% Now we want to compare the measured the time series with the predicted 
% time series in the voxels with the largest ratios. 
measured  = reshape(myCube, [500, 144]);
predicted = reshape(V , [500, 144]);
results   = reshape(bV, [500, 14]);
sigRat    = results(:,5);

% decimate for visualization
% decFact   = 10;
% measured  = measured(1:decFact:500,:);
% predicted = predicted(1:decFact:500,:);
% sigRat    = sigRat(1:decFact:500,:);

% Create bins for ploting
[Y,E] = discretize(sigRat,5);
measuredMinusPredicted = measured - predicted;

% Plot the differences in different colours
for e=1:(length(E)-1)
    eSet = measuredMinusPredicted(Y==e, :);
    % plot(repmat([1:144], [size(eSet,1),1])', eSet'); hold on;
    plot([1:144], mean(eSet)); hold on;
end
legend('1<>1.6','1.6<>2.2','2.2<>2.8','2.8<>3.4','3.4<>4')
title('Categories come from different sigRat values (sigRat goes from 1 to 4 per constraint in the model.)')
xlabel('Time Points (144, TR=2)')
ylabel('Mean(Measured - predicted) values (500 voxels)')





