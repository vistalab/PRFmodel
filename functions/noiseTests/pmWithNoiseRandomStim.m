function [compTable_noshuffle, tSeries_noshuffle, compTable_withshuffle, tSeries_withshuffle] = pmWithNoiseRandomStim(prfImplementation, varargin)
% Try to create perfect solutions for evey tool using synthetic data.
% 
% Syntax:
%    prfImplementation = 'vista' %'afni' 'popeye' % 'vista' % 'aprf';
%    [compTable, tSeries] = pmWithNoiseTests(prfImplementation);
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
pmWithNoiseTests('afni')
%}

%{
pmWithNoiseTests('aprf')
%}

%{
pmWithNoiseTests('aprfcss')
%}

%{
pmWithNoiseTests('vista')
%}

%{
pmWithNoiseTests('afni6')
%}

%{
pmWithNoiseTests('vista6')
%}






%% Read the inputs
% Make varargin lower case, remove white spaces...
prfimplementation = mrvParamFormat(prfImplementation);
varargin          = mrvParamFormat(varargin);
% Parse
p = inputParser;
p.addRequired('prfimplementation',@ischar);
p.addParameter('plotit'     ,  false            , @islogical);
p.addParameter('plotts'     ,  true           , @islogical);
p.addParameter('seed'       ,  'random');
p.addParameter('voxel'      ,  'low'           , @ischar);
p.addParameter('jitter'     ,  [0.1, 0.1]      , @isnumeric);
p.addParameter('repeats'    ,  10              , @isnumeric);
p.addParameter('signalperc' ,  'none'          , @ischar);
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
                          'doParallel'    , 0, ...
                          'rfType'        , 'gaussian');
p.addParameter('options'    ,  options    , @isstruct);

% Parse. Assign result inside each case
p.parse(prfimplementation,varargin{:});
plotts      = p.Results.plotts;
plotit      = p.Results.plotit;
allOptions  = p.Results.options;
seed        = p.Results.seed;
voxel       = p.Results.voxel;
jitter      = p.Results.jitter;
repeats     = p.Results.repeats;
signalperc  = string(p.Results.signalperc);
% We need to be sure that if only some of the params are passed, the rest will
% be taken from the defaults 
allOptions  = pmParamsCompletenessCheck(allOptions, options);

%% Create the test data
COMBINE_PARAMETERS                        = struct();

    COMBINE_PARAMETERS                       = struct();
    COMBINE_PARAMETERS.signalPercentage      = signalperc;
    COMBINE_PARAMETERS.RF.Centerx0           = [3];
    COMBINE_PARAMETERS.RF.Centery0           = [3];  
    COMBINE_PARAMETERS.RF.sigmaMajor         = [2];  
    COMBINE_PARAMETERS.RF.sigmaMinor         = 'same';
    COMBINE_PARAMETERS.TR                    = 1;
    COMBINE_PARAMETERS.Stimulus.durationSecs = 200;
    COMBINE_PARAMETERS.Stimulus.Shuffle      = false;
    COMBINE_PARAMETERS.Stimulus.shuffleSeed  = 12345;
    Noise                                    = struct();
    Noise(1).seed                            = seed; 
    Noise(1).voxel                           = voxel;
    Noise(1).jitter                          = jitter;
    COMBINE_PARAMETERS.Noise                 = Noise;

    HRF                                      = struct();
    HRF(1).Type                              = 'boynton';  
    HRF(1).normalize                         = 'height'; 
    HRF(1).params.n = 3;
    HRF(1).params.tau = 1.08;
    HRF(1).params.delay = 2.05;
        
    HRF(2).Type                              = 'boynton';
    HRF(2).normalize                         = 'height'; 
    HRF(2).params.n = 3;
    HRF(2).params.tau = 1.38;
    HRF(2).params.delay = 2;
    
    HRF(3).Type                              = 'boynton';
    HRF(3).normalize                         = 'height'; 
    HRF(3).params.n = 3;
    HRF(3).params.tau = 1.68;
    HRF(3).params.delay = 1.75;
    
    HRF(4).Type                              = 'boynton';
    HRF(4).normalize                         = 'height'; 
    HRF(4).params.n = 3;
    HRF(4).params.tau = 1.935;
    HRF(4).params.delay = 1.65;    
    COMBINE_PARAMETERS.HRF                   = HRF;
    % This is the same one as before, but now we want to do the slow stimuli version
    % by Jon's suggestion
    synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS, 'repeats',repeats);
    synthDT = pmForwardModelCalculate(synthDT);
    sDT_noshuffle = synthDT;

    
    
    
    
    
    
    % CREATE ANOTHER STIM RANDOMIZING THE BARS
    COMBINE_PARAMETERS                       = struct();
    COMBINE_PARAMETERS.signalPercentage      = signalperc;
    COMBINE_PARAMETERS.RF.Centerx0           = [3];
    COMBINE_PARAMETERS.RF.Centery0           = [3];  
    COMBINE_PARAMETERS.RF.sigmaMajor         = [2];  
    COMBINE_PARAMETERS.RF.sigmaMinor         = 'same';
    COMBINE_PARAMETERS.TR                    = 1;
    COMBINE_PARAMETERS.Stimulus.durationSecs = 200;
    COMBINE_PARAMETERS.Stimulus.Shuffle      = true;
    COMBINE_PARAMETERS.Stimulus.shuffleSeed  = 12345;
    Noise                                    = struct();
    Noise(1).seed                            = seed; 
    Noise(1).voxel                           = voxel;
    Noise(1).jitter                          = jitter;
    COMBINE_PARAMETERS.Noise                 = Noise;

    HRF                                      = struct();
    HRF(1).Type                              = 'boynton';  
    HRF(1).normalize                         = 'height'; 
    HRF(1).params.n = 3;
    HRF(1).params.tau = 1.08;
    HRF(1).params.delay = 2.05;
        
    HRF(2).Type                              = 'boynton';
    HRF(2).normalize                         = 'height'; 
    HRF(2).params.n = 3;
    HRF(2).params.tau = 1.38;
    HRF(2).params.delay = 2;
    
    HRF(3).Type                              = 'boynton';
    HRF(3).normalize                         = 'height'; 
    HRF(3).params.n = 3;
    HRF(3).params.tau = 1.68;
    HRF(3).params.delay = 1.75;
    
    HRF(4).Type                              = 'boynton';
    HRF(4).normalize                         = 'height'; 
    HRF(4).params.n = 3;
    HRF(4).params.tau = 1.935;
    HRF(4).params.delay = 1.65;    
    COMBINE_PARAMETERS.HRF                   = HRF;
    % This is the same one as before, but now we want to do the slow stimuli version
    % by Jon's suggestion
    synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS, 'repeats',repeats);
    synthDT = pmForwardModelCalculate(synthDT);
    sDT_withshuffle = synthDT;
   

%% Launch the analysis
switch prfimplementation
    case {'aprf','analyzeprf'}
        % Options
        options.aprf            = allOptions.aprf;
        % options.aprf.maxpolydeg = 0;
        options.aprf.usecss     = false;
        
        % Launch
        results_noshuffle       = pmModelFit(sDT_noshuffle,'analyzePRF','options',options);
        results_withshuffle     = pmModelFit(sDT_withshuffle,'analyzePRF','options',options);
    case {'aprfcss'}
        options.aprf            = allOptions.aprf;
        % options.aprf.maxpolydeg = 0;
        options.aprf.usecss     = true;
        
        % Launch
        results_noshuffle                 = pmModelFit(sDT_noshuffle,'analyzePRF','options',options);
        results_withshuffle                 = pmModelFit(sDT_withshuffle,'analyzePRF','options',options);
    case {'afni_4','afni4','afni'}
        options.afni            = allOptions.afni;
        
        % Launch
        results_noshuffle                 = pmModelFit(sDT_noshuffle,'afni','options',options);
        results_withshuffle                 = pmModelFit(sDT_withshuffle,'afni','options',options);
    case {'afni_6','afni6'}
        options.afni            = allOptions.afni;
        options.afni.model      = 'afni6';
        
        % Launch
        results_noshuffle                 = pmModelFit(sDT_noshuffle,'afni','options',options);
        results_withshuffle                 = pmModelFit(sDT_withshuffle,'afni','options',options);
    case {'vista','mrvista','vistasoft','vista4'}
        options.vista            = allOptions.vista;
        options.vista.model      = 'one gaussian';
        options.vista.grid       = false;  % if true, returns gFit
        options.vista.wSearch    = 'coarse to fine'; 
        % options.vista.detrend    = 0;
        options.vista.keepAllPoints            = true; 
        options.vista.numberStimulusGridPoints =  50;  
        
        % Launch
        results_noshuffle                  = pmModelFit(sDT_noshuffle,'vistasoft','options',options);    
        results_withshuffle                  = pmModelFit(sDT_withshuffle,'vistasoft','options',options);    
    case {'vistaoval','vista6'}
        options.vista            = allOptions.vista;
        options.vista.model      = 'one oval gaussian';
        options.vista.grid       = false;  % if true, returns gFit
        options.vista.wSearch    = 'coarse to fine'; 
        % options.vista.detrend    = 0;
        options.vista.keepAllPoints            = true; 
        options.vista.numberStimulusGridPoints =  50;  
        
        % Launch
        results_noshuffle                  = pmModelFit(sDT_noshuffle,'vistasoft','options',options);
        results_withshuffle                  = pmModelFit(sDT_withshuffle,'vistasoft','options',options);
    case {'popeye','pop'}
        
        % Launch
        results_noshuffle  = pmModelFit(sDT_noshuffle,'popeye');
        results_withshuffle  = pmModelFit(sDT_withshuffle,'popeye');
    case {'popnoherf','popeyenohrf'}
        
        % Launch
        results_noshuffle  = pmModelFit(sDT_noshuffle,'popeyenohrf');
        results_withshuffle  = pmModelFit(sDT_withshuffle,'popeyenohrf');
    case {'mrtools','mlrtools','mlr'}
        options.mlr            = allOptions.mlr;
        options.mlr.quickFit   = 0;
        options.mlr.doParallel = 1;
        
        % Launch
        results_noshuffle  = pmModelFit(sDT_noshuffle,'mlr','options',options);        
        results_withshuffle  = pmModelFit(sDT_withshuffle,'mlr','options',options);        
    otherwise
        error('%s not yet implemented',prfimplementation);
end

%% Create and display the results
paramDefaults = {'Centerx0','Centery0','Theta','sigmaMinor','sigmaMajor'};

[compTable_noshuffle, tSeries_noshuffle] = pmResultsCompare(sDT_noshuffle, ... % Defines the input params
                            {prfimplementation}, ... % Analysis names we want to see: 'aPRF','vista',
                            {results_noshuffle}, ...
                            'params', paramDefaults, ...
                            'shorten names',true, ...
                            'addIscloseCol', true); 

[compTable_withshuffle, tSeries_withshuffle] = pmResultsCompare(sDT_withshuffle, ... % Defines the input params
                            {prfimplementation}, ... % Analysis names we want to see: 'aPRF','vista',
                            {results_withshuffle}, ...
                            'params', paramDefaults, ...
                            'shorten names',true, ...
                            'addIscloseCol', true); 
                        
                        
                        

if plotts
    hh = mrvNewGraphWin('HRF comparison');
    set(hh,'Position',[0.007 0.62  .8  0.3]);
    nrows = 1; ncols = 2;
    
    subplot(nrows,ncols,1)
    pmTseriesPlot(tSeries_noshuffle, sDT_noshuffle(:,'TR'), ...
       'to compare', {'synth'}, ... %, prfimplementation}, ...
       'voxel',1, ... %[1:height(synthDT)], ... % 'metric','RMSE', ...
       'newWin',false)
   
    subplot(nrows,ncols,2)
    pmTseriesPlot(tSeries_withshuffle, sDT_withshuffle(:,'TR'), ...
       'to compare', {'synth'}, ... %, prfimplementation}, ...
       'voxel',1, ... %[1:height(synthDT)], ... % 'metric','RMSE', ...
       'newWin',false)
   
end

if plotit
    %%  Plot it
    hh = mrvNewGraphWin('HRF comparison');
    set(hh,'Position',[0.007 0.62  0.6  0.6]);
    nrows = 2; ncols = 4;
    if strcmp(seed,'none')
        nslvl = 'none';
    else
        nslvl  = voxel;
    end
    
    Cs  = 0.65 * distinguishable_colors(6,'w');
    
    
    % Create the fit plots with the ground truth
    tools  = {prfimplementation}; 
    HRFs   = {'boynton','boynton','boynton','boynton'}; % ,'canonical'};
    for ii=1:4 % height(noshufflecompTable)
        subplot(nrows,ncols,ii)
        useHRF = HRFs{ii};
        % ttable = compTable_noshuffle(ii,:);  
        ttable = compTable_noshuffle((ii:4:(4*repeats-4)+ii)',:);
        pmCloudOfResults(ttable  , tools ,'onlyCenters',false ,'userfsize' , 2, ...
            'centerPerc', 90     , 'useHRF'     , useHRF,'lineStyle' , '-','color',Cs(ii+1,:), ...
            'lineWidth' , .7      , 'noiselevel' , nslvl , ...
            'newWin'    , false  , 'saveTo'     , '','saveToType','svg')
        axis equal
    end
    
    % Create the fit plots with the ground truth
    tools  = {prfimplementation}; 
    HRFs   = {'boynton','boynton','boynton','boynton'}; % ,'canonical'};
    for ii=1:4  % height(withshufflecompTable)
        subplot(nrows,ncols,ncols + ii)
        useHRF = HRFs{ii};
        % ttable = compTable_withshuffle(ii,:);  
        ttable = compTable_withshuffle((ii:4:(4*repeats-4)+ii)',:);
        pmCloudOfResults(ttable   , tools ,'onlyCenters',false ,'userfsize' , 2, ...
            'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-','color',Cs(ii+1,:), ...
            'lineWidth' , .7     ,'noiselevel' ,nslvl , ...
            'newWin'    , false ,'saveTo'     ,'','saveToType','svg')
        axis equal
    end
end



end




