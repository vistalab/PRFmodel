function [compTable_noshuffle, tSeries_noshuffle, sDT_noshuffle, ...
          compTable_withshuffle, tSeries_withshuffle, sDT_withshuffle] = ...
                              pmWithNoiseRandomStim(prfImplementation, varargin)
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
% pmWithNoiseRandomStim('aprf','plotit',true,'plotts',true,'voxel','mid','signalperc','bold')
% 
%{
    % TEST IT
    [compTable_noshuffle, tSeries_noshuffle, sDT_noshuffle, ...
     compTable_withshuffle, tSeries_withshuffle, sDT_withshuffle] = ...
     pmWithNoiseRandomStim('aprfcss','plotit',true,'plotts',false,...
                           'signalperc','bold','repeats',1,'shuffleseed',12345,...
                           'radius',2, 'seed','none','voxel','mid',...
                           'hrfnorm','height','hrftype','boynton', ...
                           'boldcontrast',4)




[compTable_noshuffle, tSeries_noshuffle, sDT_noshuffle, ...
     compTable_withshuffle, tSeries_withshuffle, sDT_withshuffle] = ...
     pmWithNoiseRandomStim('vista','plotit',true,'plotts',false,...
                           'signalperc','bold','repeats',10,'shuffleseed',12345,...
                           'radius',10, 'seed','random','voxel','high',...
                           'hrfnorm','norm','hrftype','vista_twogammas', ...
                           'boldcontrast',4,'onlybarshuffleorboth','noshuffle')
    

    % LAUNCH WITH THREE DIFFERENT SHUFFLLINGS
    [compTable_noshuffle, tSeries_noshuffle, sDT_noshuffle, ...
     compTable_withshuffle, tSeries_withshuffle, sDT_withshuffle] = ...
     pmWithNoiseRandomStim('aprf','plotit',true,'plotts',false,'voxel','mid',...
                           'signalperc','bold','repeats',10,'shuffleseed',12345)

    [compTable_noshuffle, tSeries_noshuffle, sDT_noshuffle, ...
     compTable_withshuffle, tSeries_withshuffle, sDT_withshuffle] = ...
     pmWithNoiseRandomStim('aprf','plotit',true,'plotts',false,'voxel','mid',...
                           'signalperc','bold','repeats',10,'shuffleseed',54321)

    [compTable_noshuffle, tSeries_noshuffle, sDT_noshuffle, ...
     compTable_withshuffle, tSeries_withshuffle, sDT_withshuffle] = ...
     pmWithNoiseRandomStim('aprf','plotit',true,'plotts',false,'voxel','mid',...
                           'signalperc','bold','repeats',10,'shuffleseed',12453)






%}
% 
% Default tests:
%{
    pmWithNoiseRandomStim('afni')
%}
%{
    pmWithNoiseRandomStim('aprf')
%}
%{
    pmWithNoiseRandomStim('aprfcss')
%}
%{
    pmWithNoiseRandomStim('vista')
%}
%{
    pmWithNoiseRandomStim('afni6')
%}
%{
    pmWithNoiseRandomStim('vista6')
%}






%% Read the inputs
% Make varargin lower case, remove white spaces...
prfimplementation = mrvParamFormat(prfImplementation);
varargin          = mrvParamFormat(varargin);
% Parse
p = inputParser;
p.addRequired('prfimplementation',@ischar);
p.addParameter('plotit'      ,  true            , @islogical);
p.addParameter('plotts'      ,  false           , @islogical);
p.addParameter('seed'        ,  'random');
p.addParameter('voxel'       ,  'low'           , @ischar);
p.addParameter('jitter'      ,  [0.1, 0.1]      , @isnumeric);
p.addParameter('repeats'     ,  10              , @isnumeric);
p.addParameter('signalperc'  ,  'bold'          , @ischar);
p.addParameter('shuffleseed' ,  12345           , @isnumeric);
p.addParameter('boldcontrast',  5               , @isnumeric);
p.addParameter('radius'      ,  1               , @isnumeric);
p.addParameter('hrftype'     , 'boynton'        , @ischar);
p.addParameter('hrfnorm'     , 'norm'           , @ischar);
p.addParameter('scanduration', 200              , @isnumeric);
p.addParameter('tr'          , 1                , @isnumeric);
p.addParameter('window'      , false            , @islogical);
p.addParameter('addtext'     , true             , @islogical);
p.addParameter('onlybarshuffleorboth', 'both'   , @ischar);
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
shuffleseed = p.Results.shuffleseed;
boldcontrast= p.Results.boldcontrast;
size        = p.Results.radius;
hrftype     = p.Results.hrftype;
hrfnorm     = p.Results.hrfnorm;
scanduration= p.Results.scanduration;
tr          = p.Results.tr;
window      = p.Results.window;
addtext     = p.Results.addtext;
onlybarshuffleorboth = p.Results.onlybarshuffleorboth;
signalperc  = string(p.Results.signalperc);
% We need to be sure that if only some of the params are passed, the rest will
% be taken from the defaults 
allOptions  = pmParamsCompletenessCheck(allOptions, options);

switch onlybarshuffleorboth
    case {'noshuffle'}
        doNS = true;
        doS  = false;
    case {'shuffle'}
        doNS = false;
        doS  = true;
    case {'both'}
        doNS = true;
        doS  = true;
    otherwise
        error('%s not recognized, options: shuffle, noshuffle, both',onlybarshuffleorboth)
end


%% Create the test data
if doNS
    COMBINE_PARAMETERS                       = struct();
    COMBINE_PARAMETERS.signalPercentage      = signalperc;
    COMBINE_PARAMETERS.BOLDcontrast          = boldcontrast;
    COMBINE_PARAMETERS.RF.Centerx0           = [3];
    COMBINE_PARAMETERS.RF.Centery0           = [3];  
    COMBINE_PARAMETERS.RF.sigmaMajor         = [size];  
    COMBINE_PARAMETERS.RF.sigmaMinor         = 'same';
    COMBINE_PARAMETERS.TR                    = tr;
    COMBINE_PARAMETERS.Stimulus.durationSecs = scanduration;
    COMBINE_PARAMETERS.Stimulus.Shuffle      = false;
    Noise                                    = struct();
    Noise(1).seed                            = seed; 
    Noise(1).voxel                           = voxel;
    Noise(1).jitter                          = jitter;
    COMBINE_PARAMETERS.Noise                 = Noise;

    HRF                                      = struct();
    HRF(1).Type                              = hrftype;  
    HRF(1).normalize                         = hrfnorm; 
    HRF(1).params.n = 3;
    HRF(1).params.tau = 1.08;
    HRF(1).params.delay = 2.05;
        
    HRF(2).Type                              = hrftype;
    HRF(2).normalize                         = hrfnorm; 
    HRF(2).params.n = 3;
    HRF(2).params.tau = 1.38;
    HRF(2).params.delay = 2;
    
    HRF(3).Type                              = hrftype;
    HRF(3).normalize                         = hrfnorm; 
    HRF(3).params.n = 3;
    HRF(3).params.tau = 1.68;
    HRF(3).params.delay = 1.75;
    
    HRF(4).Type                              = hrftype;
    HRF(4).normalize                         = hrfnorm; 
    HRF(4).params.n = 3;
    HRF(4).params.tau = 1.935;
    HRF(4).params.delay = 1.65;    
    COMBINE_PARAMETERS.HRF                   = HRF;
    % This is the same one as before, but now we want to do the slow stimuli version
    % by Jon's suggestion
    synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS, 'repeats',repeats);
    synthDT = pmForwardModelCalculate(synthDT,'useparallel', false);
    sDT_noshuffle = synthDT;
end
    
if doS
    sDT_withshuffle = table();
        % CREATE ANOTHER STIM RANDOMIZING THE BARS
        COMBINE_PARAMETERS                       = struct();
        COMBINE_PARAMETERS.signalPercentage      = signalperc;
        COMBINE_PARAMETERS.BOLDcontrast          = boldcontrast;
        COMBINE_PARAMETERS.RF.Centerx0           = [3];
        COMBINE_PARAMETERS.RF.Centery0           = [3];  
        COMBINE_PARAMETERS.RF.sigmaMajor         = [size];  
        COMBINE_PARAMETERS.RF.sigmaMinor         = 'same';
        COMBINE_PARAMETERS.TR                    = 1;
        COMBINE_PARAMETERS.Stimulus.durationSecs = 200;
        COMBINE_PARAMETERS.Stimulus.Shuffle      = true;
        COMBINE_PARAMETERS.Stimulus.shuffleSeed  = shuffleseed;
        Noise                                    = struct();
        Noise(1).seed                            = seed; 
        Noise(1).voxel                           = voxel;
        Noise(1).jitter                          = jitter;
        COMBINE_PARAMETERS.Noise                 = Noise;

        HRF                                      = struct();
        HRF(1).Type                              = hrftype;  
        HRF(1).normalize                         = hrfnorm; 
        HRF(1).params.n = 3;
        HRF(1).params.tau = 1.08;
        HRF(1).params.delay = 2.05;

        HRF(2).Type                              = hrftype;
        HRF(2).normalize                         = hrfnorm; 
        HRF(2).params.n = 3;
        HRF(2).params.tau = 1.38;
        HRF(2).params.delay = 2;

        HRF(3).Type                              = hrftype;
        HRF(3).normalize                         = hrfnorm; 
        HRF(3).params.n = 3;
        HRF(3).params.tau = 1.68;
        HRF(3).params.delay = 1.75;

        HRF(4).Type                              = hrftype;
        HRF(4).normalize                         = hrfnorm; 
        HRF(4).params.n = 3;
        HRF(4).params.tau = 1.935;
        HRF(4).params.delay = 1.65;    
        COMBINE_PARAMETERS.HRF                   = HRF;
        % This is the same one as before, but now we want to do the slow stimuli version
        % by Jon's suggestion
        synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS, 'repeats',repeats);
        synthDT = pmForwardModelCalculate(synthDT,'useparallel', false);
        sDT_withshuffle = synthDT;
end

%% Launch the analysis
% Separate the shuffled and non-shuffled. 
% The shuffled needs to be analyzed independently, because the stimulus will be
% different. 

% NON-SHUFFLED
if doNS
    switch prfimplementation
        case {'aprf','analyzeprf'}
            % Options
            options.aprf            = allOptions.aprf;
            if strcmp(seed,'none')
                options.aprf.maxpolydeg = 0;
            end
            options.aprf.usecss     = false;
            % Launch
            results_noshuffle       = pmModelFit(sDT_noshuffle,'analyzePRF','options',options);
        case {'aprfcss'}
            options.aprf            = allOptions.aprf;
            if strcmp(seed,'none')
                options.aprf.maxpolydeg = 0;
            end
            options.aprf.usecss     = true;
            % Launch
            results_noshuffle       = pmModelFit(sDT_noshuffle,'analyzePRF','options',options);
        case {'afni_4','afni4','afni'}
            options.afni            = allOptions.afni;
            % Launch
            results_noshuffle       = pmModelFit(sDT_noshuffle,'afni','options',options);
        case {'afni_6','afni6'}
            options.afni            = allOptions.afni;
            options.afni.model      = 'afni6';
            % Launch
            results_noshuffle                 = pmModelFit(sDT_noshuffle,'afni','options',options);
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
        case {'popeye','pop'}
            % Launch
            results_noshuffle  = pmModelFit(sDT_noshuffle,'popeye');
        case {'popnoherf','popeyenohrf'}
            % Launch
            results_noshuffle  = pmModelFit(sDT_noshuffle,'popeyenohrf');
        case {'mrtools','mlrtools','mlr'}
            options.mlr            = allOptions.mlr;
            options.mlr.quickFit   = 0;
            options.mlr.doParallel = 1;
            % Launch
            results_noshuffle  = pmModelFit(sDT_noshuffle,'mlr','options',options);
        otherwise
            error('%s not yet implemented',prfimplementation);
    end
end

% SHUFFLED
if doS
    tmpShuffleResults = {};
    for ns=1:length(shuffleseed)
        sseed  = shuffleseed(ns);
        tmpsDT = sDT_withshuffle(sDT_withshuffle.Stimulus.shuffleSeed==sseed,:);
        
        switch prfimplementation
        case {'aprf','analyzeprf'}
            % Options
            options.aprf            = allOptions.aprf;
            % options.aprf.maxpolydeg = 0;
            % options.aprf.usecss     = false;
            
            % Launch
            tmpShuffleResults{ns}     = pmModelFit(tmpsDT,'analyzePRF','options',options);
            
        case {'aprfcss'}
            options.aprf            = allOptions.aprf;
            % options.aprf.maxpolydeg = 0;
            options.aprf.usecss     = true;
            
            % Launch
            tmpShuffleResults{ns}     = pmModelFit(tmpsDT,'analyzePRF','options',options);
        case {'afni_4','afni4','afni'}
            options.afni            = allOptions.afni;
            
            % Launch
            tmpShuffleResults{ns}     = pmModelFit(tmpsDT,'afni','options',options);
        case {'afni_6','afni6'}
            options.afni            = allOptions.afni;
            options.afni.model      = 'afni6';
            
            % Launch
            tmpShuffleResults{ns}                 = pmModelFit(tmpsDT,'afni','options',options);
        case {'vista','mrvista','vistasoft','vista4'}
            options.vista            = allOptions.vista;
            options.vista.model      = 'one gaussian';
            options.vista.grid       = false;  % if true, returns gFit
            options.vista.wSearch    = 'coarse to fine';
            % options.vista.detrend    = 0;
            options.vista.keepAllPoints            = true;
            options.vista.numberStimulusGridPoints =  50;
            
            % Launch
            tmpShuffleResults{ns}                  = pmModelFit(tmpsDT,'vistasoft','options',options);  
        case {'vistaoval','vista6'}
            options.vista            = allOptions.vista;
            options.vista.model      = 'one oval gaussian';
            options.vista.grid       = false;  % if true, returns gFit
            options.vista.wSearch    = 'coarse to fine';
            % options.vista.detrend    = 0;
            options.vista.keepAllPoints            = true;
            options.vista.numberStimulusGridPoints =  50;
            
            % Launch
            tmpShuffleResults{ns}                  = pmModelFit(tmpsDT,'vistasoft','options',options);
        case {'popeye','pop'}
            
            % Launch
            tmpShuffleResults{ns}  = pmModelFit(tmpsDT,'popeye');
        case {'popnoherf','popeyenohrf'}
            
            % Launch
            tmpShuffleResults{ns}  = pmModelFit(tmpsDT,'popeyenohrf');
        case {'mrtools','mlrtools','mlr'}
            options.mlr            = allOptions.mlr;
            options.mlr.quickFit   = 0;
            options.mlr.doParallel = 1;
            
            % Launch
            tmpShuffleResults{ns}  = pmModelFit(tmpsDT,'mlr','options',options); 
        otherwise
            error('%s not yet implemented',prfimplementation);
        end
    end
    
    
    % Build it back 
    results_withshuffle = table();
    for ns=1:length(shuffleseed)
        results_withshuffle = [results_withshuffle;tmpShuffleResults{ns}];
    end
    
    
    
end


    %% Create and display the results
    
    if doNS;
        paramDefaults = {'Centerx0','Centery0','Theta','sigmaMinor','sigmaMajor','tau'};
        shortnames    = {'x0','y0','Th','sMin','sMaj','tau'};
        [compTable_noshuffle, tSeries_noshuffle] = pmResultsCompare(sDT_noshuffle, ... % Defines the input params
            {prfimplementation}, ... % Analysis names we want to see: 'aPRF','vista',
            {results_noshuffle}, ...
            'params', paramDefaults, ...
            'shorten names',true, ...
            'shortnames',shortnames,...
            'addIscloseCol', false, ...
            'addsnrcol',true);
    end
    
    if doS;
        paramDefaults = {'Centerx0','Centery0','Theta','sigmaMinor','sigmaMajor','tau','shuffleseed'};
        shortnames    = {'x0','y0','Th','sMin','sMaj','tau','sseed'};
        [compTable_withshuffle, tSeries_withshuffle] = pmResultsCompare(sDT_withshuffle, ... % Defines the input params
            {prfimplementation}, ... % Analysis names we want to see: 'aPRF','vista',
            {results_withshuffle}, ...
            'params', paramDefaults, ...
            'shorten names',true, ...
            'shortnames',shortnames,...
            'addIscloseCol', false, ...
            'addsnrcol',true);
    end
    
    
    
    
    if plotts
        if window;
            hh = mrvNewGraphWin('HRF comparison');
            set(hh,'Position',[0.007 0.62  .8  0.3]);
        end
        nrows = 1; ncols = 2;
        
        subplot(nrows,ncols,1)
        pms = sDT_noshuffle.pm;
        noisevals = zeros(length(pms), size(tSeries_noshuffle.synth.BOLDnoise(1,:),2));
        for ii =1:length(pms)
            noisevals(ii,:) = pms(ii).Noise.values;
        end
        pmTseriesPlot(tSeries_noshuffle, sDT_noshuffle(:,'TR'), ...
            'to compare', {'synth'}, ... %, prfimplementation}, ...
            'voxel',1, ... %[1:height(synthDT)], ... %
            'metric','SNR', ...
            'noisevals',noisevals,...
            'newWin',false)
        
        subplot(nrows,ncols,2)
        pms = sDT_withshuffle.pm;
        noisevals = zeros(length(pms), size(tSeries_withshuffle.synth.BOLDnoise(1,:),2));
        for ii =1:length(pms)
            noisevals(ii,:) = pms(ii).Noise.values;
        end
        pmTseriesPlot(tSeries_withshuffle, sDT_withshuffle(:,'TR'), ...
            'to compare', {'synth'}, ... %, prfimplementation}, ...
            'voxel',1, ... %[1:height(synthDT)], ... %
            'metric','SNR', ...
            'noisevals',noisevals,...
            'newWin',false)
        
    end
    
    if plotit
        %%  Plot it
        if window
            hh = mrvNewGraphWin('HRF comparison');
            if strcmp(onlybarshuffleorboth,'both');
                set(hh,'Position',[0.007 0.62  0.6  0.6]);
            else
                set(hh,'Position',[0.007 0.62  0.6  0.3]);
            end
        end
        if strcmp(onlybarshuffleorboth,'both');nrows = 2; ncols = 4;
        else
            nrows = 1; ncols = 4;
        end
        if strcmp(seed,'none')
            nslvl = 'none';
        else
            nslvl  = voxel;
        end
        
        Cs  = 0.65 * distinguishable_colors(6,'w');
        
        y0 = COMBINE_PARAMETERS.RF.Centery0;
        x0 = COMBINE_PARAMETERS.RF.Centerx0;
        
        % Create the fit plots with the ground truth
        tools  = {prfimplementation};
        HRFs   = {hrftype,hrftype,hrftype,hrftype};
        if doNS
            taus   = unique(compTable_noshuffle.synth.tau);
            if ncols ~= length(taus)
                error('Number of cols should be the same as the taus for Boynton HRF')
            end
            for ii=1:length(taus)
                tau = taus(ii);
                subplot(nrows,ncols,ii)
                useHRF = HRFs{ii};
                % ttable = compTable_noshuffle(ii,:);
                % ttable = compTable_noshuffle((ii:4:(4*repeats-4)+ii)',:);
                ttable = compTable_noshuffle(compTable_noshuffle.synth.tau==tau,:);
                pmCloudOfResults(ttable  , tools ,'onlyCenters',false ,...
                    'userfsize' , size, ...
                    'centerPerc', 90     , 'useHRF'     , useHRF,...
                    'addsnr',true,'adddice',true,...
                    'lineStyle' , '-','color',Cs(ii+1,:), ...
                    'lineWidth' , 2      , 'noiselevel' , nslvl , ...
                    'xlims', [x0-2*size+1   , x0+2*size-1],...
                    'ylims', [y0-2*size+1   , y0+2*size-1],...
                    'xtick', [x0-2*size+2 : x0+2*size-2],...
                    'ytick', [y0-2*size+2 : y0+2*size-2],...
                    'addtext', addtext, ...
                    'newWin'    , false  , 'saveTo'     , '','saveToType','svg')
            end
        end
        if doS
            % Plot the shuffled one
            taus   = unique(compTable_withshuffle.synth.tau);
            if ncols ~= length(taus)
                error('Number of cols should be the same as the taus for Boynton HRF')
            end
            for ii=1:length(taus)
                tau = taus(ii);
                addcolnum = 0;
                if strcmp(onlybarshuffleorboth,'both');addcolnum = ncols;end
                subplot(nrows,ncols,addcolnum + ii)
                useHRF = HRFs{ii};
                % ttable = compTable_withshuffle(ii,:);
                % ttable = compTable_withshuffle((ii:4:(4*repeats-4)+ii)',:);
                ttable = compTable_withshuffle(compTable_withshuffle.synth.tau==tau,:);
                pmCloudOfResults(ttable   , tools ,'onlyCenters',false ,...
                    'userfsize' , size, ...
                    'centerPerc', 90    ,'useHRF'     ,useHRF,...
                    'addsnr',true,'adddice',true,...
                    'lineStyle' , '-','color',Cs(ii+1,:), ...
                    'lineWidth' , 2     ,'noiselevel' ,nslvl , ...
                    'xlims', [x0-2*size+1  , x0+2*size-1],...
                    'ylims', [y0-2*size+1  , y0+2*size-1],...
                    'xtick', [x0-2*size+2:x0+2*size-2],...
                    'ytick', [y0-2*size+2:y0+2*size-2],...
                    'addtext', addtext, ...
                    'newWin'    , false ,'saveTo'     ,'','saveToType','svg')
            end
        end
    end
    
    
    
end




