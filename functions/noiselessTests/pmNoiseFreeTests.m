function [compTable, tSeries, results] = pmNoiseFreeTests(prfImplementation, varargin)
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

% Tests:
%{
pmNoiseFreeTests('afni')
%}

%{
pmNoiseFreeTests('aprf')
%}

%{
pmNoiseFreeTests('aprfcss')
%}

%{
pmNoiseFreeTests('vista')
%}

%{
pmNoiseFreeTests('afni6')
%}

%{
pmNoiseFreeTests('vista6')
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
p.addParameter('ellipse'    ,  false           , @islogical);
p.addParameter('stimshuffle',  false           , @islogical);
% Implementation specifics
    options       = struct();
    options.aprf  = struct('seedmode'     , [0 1 2], ...
                           'display'      , 'off'  , ...
                           'usecss'       , true  );
    options.vista = struct('model'        ,'one gaussian'   , ...
                           'grid'         , false           , ...
                           'wsearch'      , 'coarse to fine', ...
                           'detrend'      , 1               , ...
                           'keepAllPoints', false           , ...
                           'numberStimulusGridPoints',   50);
    options.afni  = struct('model'        , 'afni4', ...
                           'hrf'          , 'SPM');
    options.mlr  = struct('quickFit'      , 0, ...
                           'doParallel'   , 0, ...
                           'rfType'       , 'gaussian');
p.addParameter('options'    ,  options    , @isstruct);

% Parse. Assign result inside each case
p.parse(prfimplementation,varargin{:});
useNifti    = p.Results.usenifti;
plotit      = p.Results.plotit;
ellipse     = p.Results.ellipse;
stimshuffle = p.Results.stimshuffle;

allOptions  = p.Results.options;
% We need to be sure that if only some of the params are passed, the rest will
% be taken from the defaults 
allOptions  = pmParamsCompletenessCheck(allOptions, options);

%% Create the test data
COMBINE_PARAMETERS                        = struct();
if ellipse
    COMBINE_PARAMETERS.TR                   = [1.5];
    COMBINE_PARAMETERS.RF                 = struct();
    COMBINE_PARAMETERS.RF.Centerx0        = 3;%[3]; 
    COMBINE_PARAMETERS.RF.Centery0        = 3;%[3];
    COMBINE_PARAMETERS.RF.Theta           = deg2rad(45);%[deg2rad(135)]; 
    COMBINE_PARAMETERS.RF.sigmaMajor      = [1,2,3];
    COMBINE_PARAMETERS.RF.sigmaMinor      = [1,2]; 
    % COMBINE_PARAMETERS.RF.sigmaMajor      = [1,2,4]/2;
    % COMBINE_PARAMETERS.RF.sigmaMinor      = [1,2]/2; 
else
    COMBINE_PARAMETERS.TR                 = [1.5];
    COMBINE_PARAMETERS.RF                 = struct();
    COMBINE_PARAMETERS.RF.Centerx0        = [0,3]; 
    COMBINE_PARAMETERS.RF.Centery0        = [0,3];
    COMBINE_PARAMETERS.RF.Theta           = [0]; %, deg2rad(45)];
    COMBINE_PARAMETERS.RF.sigmaMajor      = [1,2];% [1,2];
    COMBINE_PARAMETERS.RF.sigmaMinor      = "same";
end
if stimshuffle
    COMBINE_PARAMETERS.Stimulus.Shuffle   = true;
end
switch prfimplementation
    case {'aprf','analyzeprf','aprfcss'}
        HRF(1).Type = 'canonical';
    case {'afni_4','afni_6','afni','afni4','afni6'}
        HRF(1).Type = 'afni_spm';
    case {'vista','mrvista','vistasoft','vistaoval','vista4','vista6'}
        HRF(1).Type = 'vista_twogammas';
    case {'popeye','pop','popnohrf','popeyenohrf'}
        HRF(1).Type = 'popeye_twogammas';
     case {'mrtools','mlrtools','mlr'}
        HRF(1).Type = 'vista_twogammas';
    otherwise
        error('%s not yet implemented',prfimplementation);
end

COMBINE_PARAMETERS.HRF                   = HRF;
COMBINE_PARAMETERS.Stimulus.durationSecs = 300;
Noise(1).seed                    = 'none';
% Noise(1).seed                    = 'random';
COMBINE_PARAMETERS.Noise         = Noise;
synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS, 'repeats', 1);
synthDT = pmForwardModelCalculate(synthDT);

if useNifti
    % This is for nifti in pmModelFit purposes
    input = {niftiBOLDfile, jsonSynthFile, stimNiftiFname};
else
    input = synthDT;
end


%% Launch the analysis
switch prfimplementation
    case {'aprf','analyzeprf'}
        options.aprf            = allOptions.aprf;
        options.aprf.maxpolydeg = 0;
        options.aprf.usecss     = false;
        results                 = pmModelFit(input,'analyzePRF','options',options);
    case {'aprfcss'}
        options.aprf            = allOptions.aprf;
        options.aprf.maxpolydeg = 0;
        options.aprf.usecss     = true;
        results                 = pmModelFit(input,'analyzePRF','options',options);
    case {'afni_4','afni4','afni'}
        options.afni            = allOptions.afni;
        results                 = pmModelFit(input,'afni','options',options);
    case {'afni_6','afni6'}
        options.afni            = allOptions.afni;
        options.afni.model      = 'afni6';
        results                 = pmModelFit(input,'afni','options',options);
    case {'vista','mrvista','vistasoft','vista4'}
        options.vista            = allOptions.vista;
        options.vista.model      = 'one gaussian';
        options.vista.grid       = false;  % if true, returns gFit
        options.vista.wSearch    = 'coarse to fine'; 
        options.vista.detrend    = 0;
        options.vista.keepAllPoints            = true; 
        options.vista.numberStimulusGridPoints =  50;  
        results                  = pmModelFit(input,'vistasoft','options',options);    
    case {'vistaoval','vista6'}
        options.vista            = allOptions.vista;
        options.vista.model      = 'one oval gaussian';
        options.vista.grid       = false;  % if true, returns gFit
        options.vista.wSearch    = 'coarse to fine'; 
        options.vista.detrend    = 0;
        options.vista.keepAllPoints            = true; 
        options.vista.numberStimulusGridPoints =  50;  
        results                  = pmModelFit(input,'vistasoft','options',options);
    case {'popeye','pop'}
        results  = pmModelFit(input,'popeye');
    case {'popnoherf','popeyenohrf'}
        results  = pmModelFit(input,'popeyenohrf');
    case {'mrtools','mlrtools','mlr'}
        options.mlr            = allOptions.mlr;
        results  = pmModelFit(input,'mlr','options',options);        
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




