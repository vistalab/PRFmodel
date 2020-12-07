function [] = vistaPRF2MAP(bidsfolder,resultsdir,subject, session, desc,debug,run)
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
%  bidsfolder =  '/Volumes/server/Projects/SampleData/BIDS/'
%  subject    =  'wlsubj042'
%  session    = '01';
%  desc       = 'coarse';
%  aPRF2Maps(bidsfolder, subject, session, desc)

% Path to analyze PRF results



% JSON file with input arguments to analyzePRF
% opts = loadjson(fullfile(resultsdir, sprintf('sub-%s_ses-%s_%s_inputVar.json', subject, session, desc)));

% pix2deg = @(x) x * opts.stimwidthdeg / opts.stimwidthpix;

% AnalyzePRF results file
% load(fullfile(pth, sprintf('sub-%s_ses-%s_%s_results.mat', subject, session, desc)), 'results');

load(fullfile(resultsdir, sprintf('sub-%s_ses-%s_task-prf_acq-normal_run-%s_results.mat', subject, session,run)), 'results');

% Freesurfer directory
fspth = fullfile(bidsfolder, 'derivatives', 'freesurfer', 'fsaverage');

lcurv = read_curv(fullfile(fspth, 'surf', 'lh.curv'));
rcurv = read_curv(fullfile(fspth, 'surf', 'rh.curv'));

hemi{1} = zeros(length(lcurv),1);
hemi{2} = zeros(length(rcurv),1);
%%%%% debug
fields = {'x0';'y0';'sigma';'r2';cp};
results.model{1}.sigma = results.model{1}.sigma.major;


if  debug.ifdebug == 1
for f = 1 : length(fields)
    
    for h = 1 : length(hemi)
        
        tmp{h}.(fields{f}) = zeros(size(hemi{h}));
        
        if h  == 1
            
            tmp{h}.(fields{f})(opts{1}{h}.vert) = results.model{1}.(fields{f})(1:length(opts{1}{h}.vert));
        else
            
            tmp{h}.(fields{f})(opts{1}{h}.vert) = results.model{1}.(fields{f})(length(opts{1}{h-1}.vert)+1:end);
        end
    end
end

for t = 1 : length(tmp)
    
    for f = 1 : length(fields)
        
        results.model{1}.(fields{f}) = cat(1,tmp{1}.(fields{f}),tmp{2}.(fields{f}));
        
    end
    
end

end
assert(isequal(numel(lcurv) + numel(rcurv), numel(results.model{1}.x0)), ...
    'The number of vertices in the aprf results and the l&r curv files do not match;');

mgz = MRIread(fullfile(fspth, 'mri', 'orig.mgz'));


leftidx  = 1:numel(lcurv);
rightidx = (1:numel(rcurv))+numel(lcurv);

% assign from the mat file
x0    = results.model{1}.x0;
y0    = results.model{1}.y0;

angle = atan2(-results.model{1}.y0,results.model{1}.x0);
eccen = sqrt(results.model{1}.x0.^2+results.model{1}.y0.^2);

% sigma = results.model{1}.sigma.major;
sigma = results.model{1}.sigma;

vexpl = 1 - (results.model{1}.rss ./ results.model{1}.rawrss);

mgz.vol = angle(leftidx);
MRIwrite(mgz, fullfile(resultsdir, 'lh.angle.mgz'));
mgz.vol = angle(rightidx);
MRIwrite(mgz, fullfile(resultsdir, 'rh.angle.mgz'));

mgz.vol = eccen(leftidx);
MRIwrite(mgz, fullfile(resultsdir, 'lh.eccen.mgz'));
mgz.vol = eccen(rightidx);
MRIwrite(mgz, fullfile(resultsdir, 'rh.eccen.mgz'));

mgz.vol = sigma(leftidx);
MRIwrite(mgz, fullfile(resultsdir, 'lh.sigma.mgz'));
mgz.vol = sigma(rightidx);
MRIwrite(mgz, fullfile(resultsdir, 'rh.sigma.mgz'));

% r2 (convert from percentage to fraction)
mgz.vol = vexpl(leftidx);
MRIwrite(mgz, fullfile(resultsdir, 'lh.vexpl.mgz'));
mgz.vol = vexpl(rightidx);
MRIwrite(mgz, fullfile(resultsdir, 'rh.vexpl.mgz'));

mgz.vol = x0(leftidx);
MRIwrite(mgz, fullfile(resultsdir, 'lh.x.mgz'));
mgz.vol = x0(rightidx);
MRIwrite(mgz, fullfile(resultsdir, 'rh.x.mgz'));

mgz.vol = y0(leftidx);
MRIwrite(mgz, fullfile(resultsdir, 'lh.y.mgz'));
mgz.vol = y0(rightidx);
MRIwrite(mgz, fullfile(resultsdir, 'rh.y.mgz'));



% gain (in percent signal change)
% mgz.vol = results.gain(leftidx);
% MRIwrite(mgz, fullfile(pth, 'lh.gain.mgz'));
% mgz.vol = results.gain(rightidx);
% MRIwrite(mgz, fullfile(pth, 'rh.gain.mgz'));

% exponent
% mgz.vol = results.expt(leftidx);
% MRIwrite(mgz, fullfile(pth, 'lh.expon.mgz'));
% mgz.vol = results.expt(rightidx);
% MRIwrite(mgz, fullfile(pth, 'rh.expon.mgz'));

end


