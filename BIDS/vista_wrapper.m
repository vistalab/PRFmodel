% To run PRF analysis on one subject with BIDS formatting with vista solver (docker)

%% 1. open matlab and add paths:
tbUse PRFmodel;

session                 = 'nyu3t01';
subject                 = 'wlsubj064';

runnumber               = 99;
task                    = 'prf';


mainDir                 = sprintf('.'); % points to a folder were your BIDS formated folder is sitting
BidsDir                 = 'BIDS'; % name of the folder with derivatives
projectDir              = sprintf('%s/%s/',mainDir,BidsDir);
apertureFolder          = sprintf('%sderivatives/stim_apertures',projectDir);
dataFolder              = 'fmriprep';


filesDir                = sprintf('%sderivatives/%s/sub-%s/ses-%s/func',projectDir,dataFolder,subject,session);
averageFolName          = 'averageTCs';
averageFolDir           = sprintf('%sderivatives/%s',projectDir,averageFolName);
space                   = 'native';
estHRF                  = 0;

setenv('SUBJECTS_DIR',fullfile(projectDir, 'derivatives', 'freesurfer'))

%% path2configs


cfg.basename            = sprintf('sub-%s_ses-%s_task-%s_acq-normal_run-%i',subject,session,task,runnumber);
cfg.param               = sprintf('%s/sub-%s/ses-%s/func/%s_cfg.json',averageFolDir,subject,session,cfg.basename,runnumber);
cfg.parambold           = sprintf('%s/sub-%s/ses-%s/func/%s_bold.json',averageFolDir,subject,session,cfg.basename,runnumber);
cfg.events              = sprintf('%s/sub-%s/ses-%s/func/sub-%s_ses-%s_task-%s_events.json',averageFolDir,subject,session,subject,session,task);
cfg.events_tsv          = sprintf('%s/sub-%s/ses-%s/func/sub-%s_ses-%s_task-%s_events.tsv',averageFolDir,subject,session,subject,session,task);
cfg.load                = 0; % create default cfg file (NYU color retinotopy settings)


dockerscript            = 'prfanalyze_docker.sh';


debug.ifdebug           = 1; % fit pRFs only in rois specifed below
debug.roiname           = {'V1_exvivo';'V2_exvivo'}; % Roi or Rois from freesurfer label directory for the debug mode
% debug.ifdebug           = 2; % fit pRFs only in 10 voxels

makeReport              = 1;
%% convert to mgz using freescdurfer

% d = dir(sprintf('%s/*%s*.gii',filesDir,cfg.space));
d = dir(sprintf('%s/*%s*.mgz',filesDir,space));

for ii = 1:length(d)
    
    [~, fname] = fileparts(d(ii).name);
    str = sprintf('mri_convert %s/%s.gii %s/%s.mgz',filesDir, fname, filesDir,fname);
    
    if ~exist(sprintf('%s/%s.mgz',filesDir, fname'),'file') == 1
        system(str);
    end
    
end

% runnums                   =  1:length(d)/2; % / because there are 2 hemi
runnums                   =  1:6;
% runnums                   =  2; % / because there are 2 hemi

dataStr                   =  sprintf('%s*.mgz',space);


if makeReport
    
    vistaMAP2PNG(projectDir, subject,session,runnums);
    
else
    bidsVistaPRF(mainDir,projectDir,subject,session,task,runnums,...
        dataFolder,dataStr,apertureFolder,filesDir,debug,averageFolDir,cfg,...
        dockerscript,runnumber,estHRF,scriptDir);
end


