%% Explain


%% Download the 7T data, just 3 subjects, 
% Install neurophyty and other software
%{
    conda activate base
    conda install s3fs -c conda-forge
    conda install h5py
    pip install --user cpython
    pip install --user neuropythy
    export SUBJECTS_DIR=$PWD/fs
    export HCP_SUBJECTS_DIR=$PWD/hcp
    export NPYTHY_DATA_CACHE_ROOT=~/tsemp/npythy_cache
    export HCP_CREDENTIALS=~/.hcp-passwd
    export HCP_AUTO_DOWNLOAD=true
    jupyter console
%}

% In python
%{ 
    import neuropythy as ny
    import nibabel as nib
    ny.hcp.download(536647)
    ny.hcp.download(115017)
    ny.hcp.download(164131)
%}

% Download time series from NYU. (1) vpn to NYU, and to local, (2) vpn to
% Stanford and upload to black

% From the OSF page, downnload the apertures file (https://osf.io/5sj3m/io/i)

%% Convert the nifti files before processing them
% Some code:



% Load the retinotopy bar time series
subjects = {'164131','115017','536647'};
ses      = '01';
run      = '01';
hcpdir   = '/data/localhome/glerma/toolboxes/PRFmodel/local/realdata/hcp';
bidsdir  = '/data/localhome/glerma/toolboxes/PRFmodel/local/realdata/BIDS';
basedir  = '/data/localhome/glerma/toolboxes/PRFmodel/local/realdata';
fsdir    = fullfile(bidsdir,'derivatives','freesurfer');

% load the resampling
load('/data/localhome/glerma/toolboxes/PRFmodel/scripts/readDataAnalysis/HCP7T/cifti_fsaverage_lookup.mat', 'ciftifsaverageix');
% nearest neighbor interpolation to fsaverage
n = length(ciftifsaverageix)/2;

% Read the fsaverage labels to be used as masks
setenv('SUBJECTS_DIR',fsdir)
L_V1 = read_label('fsaverage', 'lh.V1_exvivo');
R_V1 = read_label('fsaverage', 'rh.V1_exvivo');
L_V2 = read_label('fsaverage', 'lh.V2_exvivo');
R_V2 = read_label('fsaverage', 'rh.V2_exvivo');
% Concatenate indexes and add 1 (fs is 0 based)
L_Mask = [L_V1(:,1);L_V2(:,1)] + 1;
R_Mask = [R_V1(:,1);R_V2(:,1)] + 1;
    
for ii=1:length(subjects)
    sub      = subjects{ii};
    niisAP   = niftiRead(fullfile(hcpdir,sub,'ret','tfMRI_RETBAR1_7T_AP_Atlas_MSMAll_hp2000_clean.dtseries.nii'));
    niisPA   = niftiRead(fullfile(hcpdir,sub,'ret','tfMRI_RETBAR1_7T_AP_Atlas_MSMAll_hp2000_clean.dtseries.nii'));
    niisAP.data = squeeze(niisAP.data);
    niisPA.data = squeeze(niisPA.data);
    % Obtain the mean of both runs
    avgdata     = ((niisAP.data + niisPA.data)/2);
    % separate left and right
    fs_left     = avgdata(:, ciftifsaverageix(1:n))';
    fs_right    = avgdata(:, ciftifsaverageix(n+(1:n)))';
    % Mask the files with the V1 and V2 ROIs
    fs_left     = fs_left(L_Mask, :);
    fs_right    = fs_right(R_Mask, :);
    % Concatenate to get one file per subject
    ccat        = [fs_left; fs_right];
    % Reshape it to be the same shape as required nifti
    reshaped    = reshape(ccat,[size(ccat,1),1,1,size(ccat,2)]);
    % Write the nifti-s back
    writeFileNifti(niftiCreate('data', reshaped, ...
        'fname',fullfile(bidsdir,['sub-' sub],'ses-01','func',...
        ['sub-' sub '_ses-' ses '_task-prf_acq-normal_run-' run '_bold.nii.gz']), ...
        'tr',1));
end


%% Plot a freesurfer mesh
%{
% read in a freesurfer mesh
[vertices, faces] = read_surf(fullfile(hcpdir,'fsaverage/surf/lh.white'));
figure(1); clf
trimesh(faces+1, vertices(:,1), vertices(:,2), vertices(:,3), mean(fs_left{1},1)); 
axis equal

% read in a freesurfer mesh
[vertices, faces] = read_surf(fullfile(hcpdir,'fsaverage/surf/rh.white'));
figure(1); hold on
trimesh(faces+1, vertices(:,1), vertices(:,2), vertices(:,3), mean(fs_right{1},1)); 
axis equal
%}

%% Write the apertures, it is the same for all of them
apertures               = load(fullfile(stimdir,'RETBARsmall.mat'));
writeFileNifti(niftiCreate('data', apertures.stim, ...
                           'fname',fullfile(stimdir,'sub-536647_ses-01_task-prf_apertures.nii'), ...
                           'tr',1));
writeFileNifti(niftiCreate('data', apertures.stim, ...
                           'fname',fullfile(stimdir,'sub-164131_ses-01_task-prf_apertures.nii'), ...
                           'tr',1));
writeFileNifti(niftiCreate('data', apertures.stim, ...
                           'fname',fullfile(stimdir,'sub-115017_ses-01_task-prf_apertures.nii'), ...
                           'tr',1));

%% Write whole mgz for plotting 
%{
% Read a template file for a mgz file
L_tempmgh = MRIread(fullfile(hcpdir,'fsaverage','surf','lh.white.avg.area.mgh'));
R_tempmgh = MRIread(fullfile(hcpdir,'fsaverage','surf','rh.white.avg.area.mgh'));

MRIwrite(L_164131_mgh, fullfile(hcpdir, 'L_164131.mgz'));
MRIwrite(R_164131_mgh, fullfile(hcpdir, 'R_164131.mgz'));
MRIwrite(L_115017_mgh, fullfile(hcpdir, 'L_115017.mgz'));
MRIwrite(R_115017_mgh, fullfile(hcpdir, 'R_115017.mgz'));
MRIwrite(L_536647_mgh, fullfile(hcpdir, 'L_536647.mgz'));
MRIwrite(R_536647_mgh, fullfile(hcpdir, 'R_536647.mgz'));

% Read a valid nifti-2 and use it as a template, just in case
% niit = niftiRead('/black/localhome/glerma/TESTDATA/prfmodel/jon_box/BIDS/sub-01/ses-01/func/sub-01_ses-01_task-prf_acq-normal_run-01_bold.nii.gz');

% Write every file in the required position so that our Docker knows how to
% run it
% mkdir(fullfile(bidsdir,'sub-164131'))


% Now we only need to apply the label to reduce size, convert to nifti 2,
% create the support files and run the analysis
%}

%% Convert output to mgz
% Convert the output of analyze PRF, i.e. results struct, to an MGZ
%
%
% INPUTS
%   bidsfolder  : path to BIDS project
%   subject     : BIDS subject name
%   session     : BIDS session name
%   desc        : type of model [default = ''];
%
% OUTPUTS
%
% Example
%{
    bidsfolder =  bidsdir;
    subject    =  '536647'
    session    = '01';
    desc       = 'coarse';
    debug      = 1;
    run        = '01';
%}
% Path to analyze PRF results
resultsdir = '/data/localhome/glerma/toolboxes/PRFmodel/local/realdata/BIDS/derivatives/prfanalyze-vista6/sub-536647/ses-01';
vistaPRF2MAP(bidsfolder,resultsdir,subject, session, desc,debug,run)


%% Plot pngs for checking
runnums  = {'01','02'};
subjects = {'536647', '115017', '164131'};
sessions = '01';

runnum  = '01';
session = '01';
subject = '536647';

vistaMAP2PNG(bidsdir, subject, session, runnum);
% Visualize an MGZ ret data on a freesurfer surface in Matlab
%
% Maps2PNG(bidsfolder, subject, session, desc)
%
% INPUTS
%   bidsfolder     : path to BIDS project
%   subject        : BIDS subject name
%   session        : BIDS session name
%   desc           : name of subfolder reflecting model/design that was
%                    used for analysis;

% Example
% bidsfolder = '/Volumes/server/Projects/SampleData/BIDS/';
% subject    = 'wlsubj042';
% session    = '01';
% desc       = 'coarse';

