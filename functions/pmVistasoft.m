function results = pmVistasoft(homedir, stimfile, datafile, stimradius, varargin)
% Converted to function from the script s_syntheticPRFvistasoft sent by Jon Winawer




%% User-specified variables > now they are input variables
%{
stimfile = '/Users/jonathanwinawer/Downloads/Exp-103_binary-true_size-20x20.nii.gz';
datafile = '/Users/jonathanwinawer/Downloads/synthDataExample3_1D_TR2.nii.gz';

stimradius = 10; % degrees
% name it
sessioncode = 'pRFsynthetic01';

% where to put it
homedir = fullfile(filesep, 'Volumes', 'server', 'Projects', 'PRF_Simulations', sessioncode);
%}

% First clear the workspace just in case
mrvCleanWorkspace


% Read the inputs
varargin = mrvParamFormat(varargin);
p = inputParser;
p.addRequired('homedir'                         , @ischar);
p.addRequired('stimfile'                        , @ischar);
p.addRequired('datafile'                        , @ischar);
p.addRequired('stimradius'                      , @isnumeric);
p.addParameter('sessioncode'  , 'pRFsynthetic01', @ischar);
p.addParameter('model'        , 'one gaussian'  , @ischar);
p.addParameter('grid'         , false           , @islogical);
p.addParameter('wsearch'      , 'coarse to fine', @ischar);
p.addParameter('detrend'      , 1               , @isnumeric);
p.addParameter('keepAllPoints', false           , @islogical);
p.addParameter('numberStimulusGridPoints', 50   , @isnumeric);
p.parse(homedir, stimfile, datafile, stimradius, varargin{:});
% Assign it
sessioncode   = p.Results.sessioncode;
model         = p.Results.model;
grid          = p.Results.grid;
wSearch       = p.Results.wsearch;
detrend       = p.Results.detrend;
keepAllPoints = p.Results.keepAllPoints;
numberStimulusGridPoints = p.Results.numberStimulusGridPoints;

% Disp the input files for debugging
fprintf('\n[pmVistasoft] This is homedir: %s\n',homedir)
fprintf('\n[pmVistasoft] This is stimfile: %s\n',stimfile)
fprintf('\n[pmVistasoft] This is datafile: %s\n',datafile)
fprintf('\n[pmVistasoft] This is stimradius: %i\n',stimradius)



% TODO: write all options for wSearch and model

%   wSearch     : 1 = grid search only ("coarse"),
%                 2 = minimization search only ("fine"),
%                 3 = grid followed by minimization search [default]
%                 4 = grid followed by two minimization searches, the first
%                     with temporal decimation, the second without.
%                 5 = grid followed by two minimization searches, followed
%                     by HRF search, followed by PRF search






%% Set up files and directories
% mkdir(homedir); 



cd(homedir);

if exist(fullfile(homedir,'Raw'),'dir');error('RAW DIR EXISTS');end
mkdir(fullfile(homedir, 'Raw'));
mkdir(fullfile(homedir, 'Stimuli'));

% move the data into the new directories
copyfile(datafile, fullfile(homedir, 'Raw'));
copyfile(stimfile, fullfile(homedir, 'Stimuli'));

%% convert stim file to .mat format that vistasoft can read
[~, f, e] = fileparts(stimfile);
ni = niftiRead(fullfile(homedir, 'Stimuli', sprintf('%s%s', f, e)));

images = squeeze(ni.data);
pixdim = niftiGet(ni, 'pixdim');
tr     = pixdim(end);
sprintf('/n/n USING TR:%2.2f/n/n',tr)
stimulus.seq = 1:size(images,3);
stimulus.seqtiming = (stimulus.seq-1) * tr;

stimfileMat = fullfile('.', 'Stimuli', 'images_and_params');
save(stimfileMat, 'images', 'stimulus');

%% create a pseudo inplane underlay, required by vistasoft, by averaging the
%   time series for each voxel
fmri        = niftiRead(datafile);
ippath      = fullfile('.', 'Raw', 'inplane.nii.gz');
ip          = fmri; 
ip.data     = mean(fmri.data, length(size(fmri.data)));
ip.dim(end) = 1
niftiWrite(ip, ippath);

A = niftiRead(ippath)


%% Set up the vistasoft session
params = mrInitDefaultParams;

params.sessionDir   = homedir;
params.vAnatomy     = [];
params.sessionCode  = sessioncode;
params.inplane      = ippath;

[p, f, e] = fileparts(datafile);
params.functionals  = fullfile('.','Raw', sprintf('%s%s', f,e));

% Run it:
ok = mrInit(params);
params
dir(fullfile('.', filesep,'Raw'))

params
dir(fullfile('.', filesep,'Raw'))



%% Check it
%{
vw = initHiddenInplane();
sz =  viewGet(vw,'anatsize');
[a, b, c] = ind2sub(sz, 1:prod(sz));
coords = [a; b; c];
vw = newROI(vw,'all',true,'w',coords, 'all voxels');

% plot raw signal
newGraphWin();
plotMeanTSeries(vw, viewGet(vw, 'current scan'), [], true);

% plot percent signal modulation
newGraphWin();
plotMeanTSeries(vw, viewGet(vw, 'current scan'), [], false);
%}

%% Set up prf model

vw = initHiddenInplane();
% edit GLU: dataTYPES is not found, but it was stablished as global in mrInit()
% Load mrSESSION in here to see if this solves it
load(fullfile(homedir,'mrSESSION.mat'))
disp(dataTYPES)
dataTYPES.scanParams
% Set default retinotopy stimulus model parameters
sParams = rmCreateStim(vw);
sParams
sParams.stimType   = 'StimFromScan'; % This means the stimulus images will
                                     % be read from a file.
sParams.stimSize   = stimradius;     % stimulus radius (deg visual angle)
sParams.nDCT       = detrend;        % detrending frequeny maximum (cycles
                                     % per scan): 1 means 3 detrending
                                     % terms, DC (0 cps), 0.5, and 1 cps
sParams.imFile     = stimfileMat;    % file containing stimulus images
sParams.paramsFile = stimfileMat;    % file containing stimulus parameters
% 'thresholdedBinary',  whenreading in images, treat any pixel value
%                       different from background as a 1, else 0
sParams.imFilter   = 'none';
% we switch from the default positive Boynton hRF to the biphasic SPM style
sParams.hrfType    = 'two gammas (SPM style)';
% pre-scan duration will be stored in frames for the rm, but was stored in
% seconds in the stimulus file
sParams.prescanDuration = 0;


dataTYPES = dtSet(dataTYPES, 'rm stim params', sParams);

saveSession();
% Check it
vw = rmLoadParameters(vw);

% edit GLU: this is opening a new window, hide it
% [~, M] = rmStimulusMatrix(viewGet(vw, 'rmparams'), [], [], 1, false);


%% Solve prf Models
vw = initHiddenInplane();

vw = rmMain(vw, [], wSearch, ...
            'model', {model}, ...
            'matFileName', 'tmpResults', ...
            'keepAllPoints', keepAllPoints, ...
            'numberStimulusGridPoints', numberStimulusGridPoints);

% Load the results        

if grid 
    d = dir(fullfile(dataDir(vw), sprintf('%s*', 'tmpResults-gFit')));
    results = load(fullfile(dataDir(vw), d.name));
else
    d = dir(fullfile(dataDir(vw), sprintf('%s*', 'tmpResults')));
    [~,newestIndex] = max([d.datenum]);
    results = load(fullfile(dataDir(vw), d(newestIndex).name));
end

% Examples below, delete at some point
%{
vw = rmMain(vw, [], 'coarse to fine and hrf', ...
    'model', {'one gaussian'}, 'matFileName', 'rm-linear-isotropic');

% note the hrf search does not currently work for css model: fix it!
vw = rmMain(vw, [], 'coarse to fine', ...
    'model', {'css'}, 'matFileName', 'rm-css-isotropic');
%}

%% Delete all global variables created by mrVista
mrvCleanWorkspace

% Delete the temp folder with all the tmp files, we only want the results
cd(homedir)
cd('../')
rmdir(homedir, 's')

%% Prepare results

% Right now return just results, we do the modifications in pmModelFit, decide
% later if we want to give back the tables with only the results we are going to
% need


end
