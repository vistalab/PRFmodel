function [] = vistaPRF2MAP(bidsfolder,subject, session, desc,debug,run)
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
% 2021-02: Heavily edited by Gari
% 
% Example:
% 
%{
close all; clear all;
 bidsfolder = '/black/localhome/glerma/TESTDATA/RETPIPE/BIDS'
 resultsdir = fullfile(bidsfolder,'results')  
 subject    = 'ex22163'
 session    = 'T01taskshuffledPRFspacefsnativehemiL';
 desc       = 'css';
 debug      = 0;
 run        = '01';
 
vistaPRF2MAP(bidsfolder, subject,session,desc,debug,run)



cd(fullfile(bidsfolder, 'derivatives','prfanalyze-vista4',['sub-' subject]))
sess = dir('ses*');
for ns = 1:length(sess)
    session = sess(ns).name(5:end);
    vistaPRF2MAP(bidsfolder, subject,session,desc,debug,run)
end

%}
% Path to analyze PRF results

resultDir = fullfile(bidsfolder, 'derivatives','prfanalyze-vista4',...
                                            ['sub-' subject],['ses-' session]);
if ~exist(resultDir,'dir');error("Can't find results directory %s",resultDir);end
cd(resultDir)

%% Read the results back to Matlab
% Read all the nifti-s for this type of tool 
niftis = dir('*.nii.gz');
if length(niftis) < 6
    error('Not enough nifti result files in output, at least there should be x0,y0,sigmamajor,sigmaminor,theta and modelpred')
end
filesRead = 0;

pmEstimates = table();
for nn=1:length(niftis)
    fname = niftis(nn).name;
    % Assume always it will be .nii.gz and that the . has not been used in the filename
    trunkname = split(fname,'.');
    if length(trunkname) ~= 3; error('%s should be a filename with no dots and .nii.gz',fname);end	
    trunkname = trunkname{1};
    resname   = split(trunkname, '_');
    resname   = resname{end};
    % read the nifti and squeeze the result matrix
    fprintf('Attempting to read %s\n',fullfile(pwd,fname))
    tmp       = niftiRead(fname);
    data      = squeeze(tmp.data);	
    % asign it to the table
    % Noah used other names in the filenames, leave this hardcoded here until we decide what to do
    % noahnames = {'modelpred', 'sigmamajor','sigmaminor', 'testdata', 'theta', 'x0' ,'y0'};
    % {'Centerx0','Centery0','Theta','sigmaMinor','sigmaMajor'}
    switch resname
        case 'sigmaminor', resname = 'sMin';
        case 'sigmamajor', resname = 'sMaj';
        case 'theta'     , resname = 'Theta';
        case 'centerx0'  , resname = 'x0';
        case 'centery0'  , resname = 'y0';
    end
    pmEstimates.(resname) = data;	        
end

H = session(end);
if H=='L';h='lh';elseif H=='R';h='rh';else error('hemi %s not recognized',H);end

% Freesurfer directory
fspth = fullfile(bidsfolder,'derivatives','freesurfer',['sub-' subject]);
curv  = read_curv(fullfile(fspth, 'surf', [h '.curv']));

assert(isequal(numel(curv), height(pmEstimates)), ...
    'The number of vertices in the aprf results and the l&r curv files do not match;');

mgz = MRIread(fullfile(fspth, 'mri', 'orig.mgz'));
idx  = 1:numel(curv);


% assign from the mat file
x0    = pmEstimates.x0;
y0    = pmEstimates.y0;
sigma = pmEstimates.sMaj;
% vexpl = 1 - (results.model{1}.rss ./ results.model{1}.rawrss);
vexpl = pmEstimates.r2;
[angle,eccen] = cart2pol(x0,y0);



mgz.vol = angle(idx);
MRIwrite(mgz, fullfile(resultDir, [h '.angle.mgz']));

mgz.vol = eccen(idx);
MRIwrite(mgz, fullfile(resultDir, [h '.eccen.mgz']));

mgz.vol = sigma(idx);
MRIwrite(mgz, fullfile(resultDir, [h '.sigma.mgz']));

% r2 (convert from percentage to fraction)
mgz.vol = vexpl(idx);
MRIwrite(mgz, fullfile(resultDir, [h '.vexpl.mgz']));

mgz.vol = x0(idx);
MRIwrite(mgz, fullfile(resultDir, [h '.x.mgz']));

mgz.vol = y0(idx);
MRIwrite(mgz, fullfile(resultDir, [h '.y.mgz']));

end


