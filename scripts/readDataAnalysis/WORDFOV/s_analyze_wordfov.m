%% Script to analyze CNI GE 3T 
% Folder names
tr = 1.4;
projdir = '/black/localhome/glerma/toolboxes/PRFmodel/local/WORDSFOV';
bidsdir = '/black/localhome/glerma/toolboxes/PRFmodel/local/WORDSFOV/BIDS';
fsdir   = fullfile(bidsdir,'derivatives','freesurfer')
fmriprep= fullfile(bidsdir,'derivatives','fmriprep')
stimdir = fullfile(bidsdir,'stimuli');


%% Stimnulus: create nifti from fmriprep
subname = '003';
ses = '1';
aperturesdir = fullfile(fmriprep,['sub-' subname],['ses-' ses],'stimuli');
% mats = dir(fullfile(aperturesdir,'*Retinotopy*.mat'));
apertures  = load(fullfile(aperturesdir,'sub-003_ses-1_task-Retinotopy_run-1_stimMovie.mat'));
niistim    = niftiCreate('data', apertures.stimulus, ...
                           'fname',fullfile(stimdir,['sub-' subname '_ses-' ses '_task-prf_apertures.nii.gz']), ...
                           'tr',tr);
niistim.pixdim(end) = tr; 
niftiWrite(niistim);

% cdiA = niftiRead(fullfile(stimdir,['sub-' subname '_ses-' ses '_task-prf_apertures.nii.gz'])); 

%% Create the datafile  
run = '2';
% We manually edited the other required  files in time series directories
origdatadir = fullfile(fmriprep,['sub-' subname],['ses-' ses],'func');

% Convert from gifti to mgz
Lfname = ['sub-003_ses-1_task-Retinotopy_run-' run '_space-fsnative_hemi-L.func'];
cmd = ['mri_convert ' fullfile(origdatadir,[Lfname '.gii '])  fullfile(origdatadir,[Lfname '.mgz'])];
system(cmd)

Rfname = ['sub-003_ses-1_task-Retinotopy_run-' run '_space-fsnative_hemi-R.func'];
cmd = ['mri_convert ' fullfile(origdatadir,[Rfname '.gii '])  fullfile(origdatadir,[Rfname '.mgz'])];
system(cmd)

% Read the mgz files
L = MRIread(fullfile(origdatadir,[Lfname '.mgz']));
R = MRIread(fullfile(origdatadir,[Rfname '.mgz']));
% Reshape this
L.vol = permute(L.vol,[2,1,3,4]);
R.vol = permute(R.vol,[2,1,3,4]);

% Read the labels
setenv('SUBJECTS_DIR',fsdir)
L_V1 = read_label('sub-003', 'lh.V1_exvivo.thresh');
R_V1 = read_label('sub-003', 'rh.V1_exvivo.thresh');
L_V2 = read_label('sub-003', 'lh.V2_exvivo.thresh');
R_V2 = read_label('sub-003', 'rh.V2_exvivo.thresh');
% Concatenate labels to create mask
L_Mask = [L_V1(:,1);L_V2(:,1)] + 1;
R_Mask = [R_V1(:,1);R_V2(:,1)] + 1;
% Create masked volume
L.vol = L.vol(L_Mask,:,:,:);
R.vol = R.vol(R_Mask,:,:,:);
% Concatenate
LR = L;
LR.vol = [L.vol;R.vol];

% Write the niftis directly, they don't need to be nifti-2 now, as they are
% smaller
writeFileNifti(niftiCreate('data', LR.vol, ...
    'fname',fullfile(bidsdir,['sub-' subname],['ses-' ses],'func',...
    ['sub-' subname '_ses-' ses '_task-prf_acq-normal_run-' run '_bold.nii.gz']), ...
    'tr',tr));
                       
%% Create the config file to run the analyses
% 1./ create the json file 
% 2./ Run the run_prfanalyze.sh script to run the Docker containers.s

                       
                       
                       