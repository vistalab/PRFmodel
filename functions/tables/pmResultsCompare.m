function [newDT, tSeries] = pmResultsCompare(synthDT, resNames, resDT, varargin)
% Takes the original table that generated results and compares to the
%          analysisTool
% 
%  Inputs: table with one or several rows of parameters to calculate bold series
%          in fwd model
% 
%  Outputs: returns a table with the comparisons
% 
%  See also: 
% 
%  GLU Vistalab 2019.07


%% Read the inputs
varargin = mrvParamFormat(varargin);
p = inputParser;
p.addRequired ('synthDT'        , @istable);
p.addRequired ('resNames'       , @iscell);
p.addRequired ('resDT'          , @iscell);
paramDefaults = {'Centerx0','Centery0','Theta','sigmaMinor','sigmaMajor'};
shortDefaults = {'x0','y0','Th','sMin','sMaj'};
p.addParameter('params'         , paramDefaults , @iscell);
p.addParameter('shortennames'   , false         , @islogical);
p.addParameter('addisclosecol'  , false         , @islogical);
p.addParameter('tolerance'      , 0.0001        , @isnumeric);
p.addParameter('dotseries'      , true         , @islogical);
p.parse(synthDT, resNames, resDT, varargin{:});

params       = p.Results.params;
shortenNames = p.Results.shortennames;
addIscloseCol= p.Results.addisclosecol;
tolerance    = p.Results.tolerance;
dotSeries    = p.Results.dotseries;


%% Calculate
% Keep adding columns that will be the measures we will be using

% Decide what is the best structure to visualize this calculations
% Right now use it just to concatenate
% Later decide if this function returns another table with just the results

newDT   = table();
newDT.synth = synthDT.RF(:,params);
if shortenNames && isequal(params,paramDefaults)
    newDT.synth.Properties.VariableNames = shortDefaults;
end
% Add the HRF type
newDT.HRFtype   = synthDT.HRF.Type;
% Add the TR type
newDT.TR        = synthDT.TR;
% Add the White Noise param
newDT.noise2sig = synthDT.Noise.noise2signal;



for ii=1:length(resNames)
    newDT.(resNames{ii}) = resDT{ii}(:,params);
    newDT.(resNames{ii}).Properties.VariableNames = shortDefaults;
end

% Add the isclose col for the perfect prediction case if required.
if addIscloseCol
    newDT.([resNames{ii} '_test']) = newDT.(resNames{ii});
    varNames = newDT.([resNames{ii} '_test']).Properties.VariableNames;
    for nn=1:length(varNames)
        X = newDT.synth.(varNames{nn});
        Y = newDT.(resNames{ii}).(varNames{nn});
        newDT.([resNames{ii} '_test']).(varNames{nn}) = isclose(X,Y, ...
                                                            'tolerance',0.001, ...
                                                            'returnvector',true);
    end
end




%% Create the tSeries table
tSeries           = table();
if dotSeries
    % Extract first pm for information
    pm1               = synthDT.pm(1);
    % Create synth table and add ones
    tSeries.synth     = table();
    tSeries.synth.BOLDnoise = repmat(ones([1,pm1.timePointsN]), [height(synthDT),1]);
    % Populate it with the real content in the pm-s
    PMs               = synthDT.pm;
    for ii=1:height(tSeries); tSeries.synth{ii,'BOLDnoise'}=PMs(ii).BOLDnoise;end
    
    
    
    % Now add the analysis results that we have
    for ii=1:length(resNames)
        tSeries.(resNames{ii}) = resDT{ii}(:,{'testdata','modelpred'});
    end
end

