function [compTable, tSeries] = pmNoiseFreeTests(prfImplementation, varargin)
% Try to create perfect solutions for evey tool using synthetic data.
% 
% Syntax:
%    prfImplementation = 'vista' %'afni' 'popeye' % 'vista' % 'aprf';
%    [compTable, tSeries] = pmNoiseFreeTests(prfImplementation);
%  
% Brief description:
%    Creates a example dataset and for the tool selected, it will recreate a
%    perfect solution
%
% Inputs:
%   prfImplementation - String defining the model
%
% Outputs: 
%   pmEstimates: Table format of the pRF model parameters in results
%   results:     The struct from analyzePRF
%
% Key/val parameters (Optional)
%   N/A
%
% GLU Vistalab 07.2019
%
% See also:
%     pmModelFit
%
% 

% Examples:
%{
  pmCompute...
%}

%% Read the inputs
% Make varargin lower case, remove white spaces...
prfimplementation = mrvParamFormat(prfImplementation);
varargin          = mrvParamFormat(varargin);
% Parse
p = inputParser;
p.addRequired('prfimplementation',@ischar);
p.addParameter('usenifti'   ,  false           , @islogical);
p.addParameter('plotit'     ,  false           , @islogical);
% Implementation specifics
% AnalyzePRF
options  = struct('seedmode', [0,1,2], 'display', 'off');
p.addParameter('options'    ,  options        , @isstruct);
% Vistasoft
p.addParameter('model'      , 'one gaussian'  , @ischar);
p.addParameter('grid'       , false           , @islogical);
p.addParameter('wsearch'    , 'coarse to fine', @ischar);
% AFNI
% p.addParameter('wsearch'  , 'coarse to fine', @ischar);

% Parse. Assign result inside each case
p.parse(prfimplementation,varargin{:});
useNifti = p.Results.usenifti;
plotit   = p.Results.plotit;


%% Create the test data
COMBINE_PARAMETERS.RF.Centerx0        = [0,3]; % [-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6];
COMBINE_PARAMETERS.RF.Centery0        = [0,3];
COMBINE_PARAMETERS.RF.Theta           = [0]; %, deg2rad(45)];
COMBINE_PARAMETERS.RF.sigmaMajor      = [1,2];
COMBINE_PARAMETERS.RF.sigmaMinor      = "same";

switch prfimplementation
    case {'aprf','analyzeprf','aprfcss'}
        COMBINE_PARAMETERS.TR                   = [1.5];
        HRF(1).Type = 'canonical';
    case {'afni_4','afni_6','afni'}
        COMBINE_PARAMETERS.TR                   = [1.5];
        % HRF(1).Type = 'afni_spm';
        HRF(1).Type = 'vista_twogammas';
    case {'vista','mrvista','vistasoft'}
        COMBINE_PARAMETERS.TR                   = [1.5];
        HRF(1).Type = 'vista_twogammas';
        % COMBINE_PARAMETERS.Stimulus.ResizedHorz = [100];
        % COMBINE_PARAMETERS.Stimulus.ResizedVert = [100];
    case {'popeye','pop','popnohrf','popeyenohrf'}
        COMBINE_PARAMETERS.TR                   = [1.5]; % before it had to be one because the hrf was hardcoded
        % amazing, TR:1 and 3, all ok, for TR:2, the last test fails and it is
        % not capable of predicting anything. 
        % HRF(1).Type = 'popeye_twogammas';
        HRF(1).Type = 'vista_twogammas';
        HRF(2).Type = 'canonical';
    otherwise
        error('%s not yet implemented',prfimplementation);
end

COMBINE_PARAMETERS.HRF           = HRF;
Noise(1).seed                    = 'none';
COMBINE_PARAMETERS.Noise         = Noise;
synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS, 'repeats', 2);
synthDT = pmForwardModelCalculate(synthDT);
% Visually check that all the combinations we specified are there
% [synthDT.RF(:,{'Centerx0','Centery0','Theta','sigmaMajor','sigmaMinor'}), ...
%  synthDT(:,'TR'), ...
%  synthDT.HRF(:,'Type')...
% ]

% Save the default niftis with different TR and HRF to be used as tests later on
% niftiBOLDfile = fullfile(pmRootPath,'local', ...
%     ['defaultSynth_TR' num2str(COMBINE_PARAMETERS.TR) '_HRF-' HRF(1).Type '.nii.gz']);
% if ~exist(niftiBOLDfile, 'file')
%     pmForwardModelToNifti(synthDT,'fname',niftiBOLDfile, 'demean',false);
% end

% jsonSynthFile = fullfile(pmRootPath,'local', ...
%     ['defaultSynth_TR' num2str(COMBINE_PARAMETERS.TR) '_HRF-' HRF(1).Type '.json']);
% if ~exist(jsonSynthFile, 'file')
%     % Encode json
%     jsonString = jsonencode(synthDT(:,1:(end-1)));
%     % Format a little bit
%     jsonString = strrep(jsonString, ',', sprintf(',\r'));
%     jsonString = strrep(jsonString, '[{', sprintf('[\r{\r'));
%     jsonString = strrep(jsonString, '}]', sprintf('\r}\r]'));
%     % Write it
%     fid = fopen(jsonSynthFile, 'w');if fid == -1,error('Cannot create JSON file');end
%     fwrite(fid, jsonString, 'char');fclose(fid);
%     % Read the json
%     %{
%     A = struct2table(jsonread(jsonSynthFile));
%     for na=1:width(A)
%         if isstruct(A{:,na})
%             A.(A.Properties.VariableNames{na}) = struct2table(A{:,na});
%         end
%     end
%     %}
% end

%{
stimNiftiFname = fullfile(pmRootPath,'local', ['defaultStim_TR' num2str(COMBINE_PARAMETERS.TR) '.nii.gz']);
if ~exist(stimNiftiFname, 'file')
    pm1            = prfModel; 
    pm1.TR         = COMBINE_PARAMETERS.TR;
    pm1.compute;
    stimNiftiFname = pm1.Stimulus.toNifti('fname',stimNiftiFname);
end
%}






if useNifti
    % This is for nifti in pmModelFit purposes
    input = {niftiBOLDfile, jsonSynthFile, stimNiftiFname};
else
    input = synthDT;
end


%% Launch the analysis
switch prfimplementation
    case {'aprf','analyzeprf'}
        options  = struct('seedmode',[0,1,2], 'display','off', 'maxpolydeg',0,'usecss',false);
        results = pmModelFit(input,'analyzePRF','options',options);
    case {'aprfcss'}
        options  = struct('seedmode',[0,1,2], 'display','off', 'maxpolydeg',0,'usecss',true);
        results = pmModelFit(input,'analyzePRF','options',options);
    case {'afni_4','afni_6','afni'}
        results    = pmModelFit(input,'afni_4','afni_hrf','SPM');
    case {'vista','mrvista','vistasoft'}
        results      = pmModelFit(input,'vistasoft', ...
            'model','one gaussian', ...
            'grid', false, ... % if true, returns gFit
            'wSearch', 'coarse to fine', ...
            'detrend', 0, ...
            'keepAllPoints', true); % , ...
            % 'numberStimulusGridPoints', 50);  %  We need to remove it otherwise it will find an average HRF for all of them
    case {'popeye','pop'}
        results  = pmModelFit(input,'popeye');
    case {'popnoherf','popeyenohrf'}
        results  = pmModelFit(input,'popeyenohrf');
    otherwise
        error('%s not yet implemented',prfimplementation);
end

%% Create and display the results
paramDefaults = {'Centerx0','Centery0','Theta','sigmaMinor','sigmaMajor'};
if ~strcmp(synthDT.Properties.VariableNames{end},'pm')
    synthDT = pmForwardModelCalculate(synthDT);
end
[compTable, tSeries] = pmResultsCompare(synthDT, ... % Defines the input params
                            {prfimplementation}, ... % Analysis names we want to see: 'aPRF','vista',
                            {results}, ...
                            'params', paramDefaults, ...
                            'shorten names',true, ...
                            'addIscloseCol', true); 
% Visualize with 2 digits after comma
format bank; disp(compTable); format

if plotit
    pmTseriesPlot(tSeries, synthDT(:,'TR'), ...
        'to compare', {'synth', prfimplementation}, ...
        'voxel',[1:height(synthDT)], ... % 'metric','RMSE', ...
        'newWin',true)
end
end




