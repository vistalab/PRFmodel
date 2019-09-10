function [compTable, tSeries] = pmNoiseFreeTests(prfImplementation, varargin)
% Try to create perfect solutions for evey tool using synthetic data.
% 
% Syntax:
%    prfImplementation = 'popeye' %'afni' 'popeye' % 'vista' % 'aprf';
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
% This options structs are defaults for analyzePRF
options  = struct('seedmode', [0 1], 'display' , 'off');
% Implementation specifics
% AnalyzePRF
p.addParameter('options'    ,  options        , @isstruct);
% Vistasoft
p.addParameter('model'      ,  'one gaussian' , @ischar);
p.addParameter('grid'       , false           , @islogical);
p.addParameter('wsearch'    , 'coarse to fine', @ischar);
% AFNI
% p.addParameter('wsearch'    , 'coarse to fine', @ischar);

% Parse. Assign result inside each case
p.parse(prfimplementation,varargin{:});


%% Create the test data
COMBINE_PARAMETERS.RF.Centerx0        = [0, -6, 6]; % [-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6];
COMBINE_PARAMETERS.RF.Centery0        = [0, -6, 6];
COMBINE_PARAMETERS.RF.Theta           = [0]; %, deg2rad(45)];
COMBINE_PARAMETERS.RF.sigmaMinor      = [4];
COMBINE_PARAMETERS.RF.sigmaMajor      = [4];
COMBINE_PARAMETERS.Noise.noise2signal = [0];  % By default only white noise added

switch prfimplementation
    case {'aprf','analyzeprf'}
        COMBINE_PARAMETERS.TR                   = [1];
        HRF(1).Type = 'canonical';
    case {'afni_4','afni_6','afni'}
        COMBINE_PARAMETERS.TR                   = [2];
        HRF(1).Type = 'afni_spm';
    case {'vista','mrvista','vistasoft'}
        COMBINE_PARAMETERS.TR                   = [2];
        HRF(1).Type = 'vista_twogammas';
        % COMBINE_PARAMETERS.Stimulus.ResizedHorz = [100];
        % COMBINE_PARAMETERS.Stimulus.ResizedVert = [100];
    case {'popeye','pop'}
        COMBINE_PARAMETERS.TR                   = [2]; % before it had to be one because the hrf was hardcoded
        % amazing, TR:1 and 3, all ok, for TR:2, the last test fails and it is
        % not capable of predicting anything. 
        HRF(1).Type = 'popeye_twogammas';
    otherwise
        error('%s not yet implemented',prfimplementation);
end

COMBINE_PARAMETERS.HRF           = HRF;
synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS);
synthDT = pmForwardModelCalculate(synthDT);
% Visually check that all the combinations we specified are there
% [synthDT.RF(:,{'Centerx0','Centery0','Theta','sigmaMajor','sigmaMinor'}), ...
%  synthDT(:,'TR'), ...
%  synthDT.HRF(:,'Type')...
% ]




%% Launch the analysis
switch prfimplementation
    case {'aprf','analyzeprf'}
        options  = struct('seedmode',[0,1], 'display','off', 'maxpolydeg',0);
        results = pmModelFit(synthDT,'analyzePRF','options',options);
    case {'afni_4','afni_6','afni'}
        results    = pmModelFit(synthDT,'afni_4','afni_hrf','SPM');
    case {'vista','mrvista','vistasoft'}
        results      = pmModelFit(synthDT,'vistasoft', ...
            'model','one gaussian', ...
            'grid', false, ... % if true, returns gFit
            'wSearch', 'coarse to fine', ...
            'detrend', false, ...
            'keepAllPoints', true, ...
            'numberStimulusGridPoints', 50);  %  We need to remove it otherwise it will find an average HRF for all of them
    case {'popeye','pop'}
        results  = pmModelFit(synthDT,'popeye_onegaussian');
    otherwise
        error('%s not yet implemented',prfimplementation);
end

%% Create and display the results
paramDefaults = {'Centerx0','Centery0','Theta','sigmaMinor','sigmaMajor'};
[compTable, tSeries] = pmResultsCompare(synthDT, ... % Defines the input params
                            {prfimplementation}, ... % Analysis names we want to see: 'aPRF','vista',
                            {results}, ...
                            'params', paramDefaults, ...
                            'shorten names',true, ...
                            'addIscloseCol', true); 
% Visualize with 2 digits after comma
format bank; disp(compTable); format


pmTseriesPlot(tSeries, synthDT(:,'TR'), ...
    'to compare', {'synth', prfimplementation}, ...
    'voxel',[1:height(synthDT)], ... % 'metric','RMSE', ...
    'newWin',true)
end




