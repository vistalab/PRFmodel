function newDT = pmResultsCompare(paramsDT, resNames, resDT, varargin)
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
p.addRequired ('paramsDT'       , @istable);
p.addRequired ('resNames'       , @iscell);
p.addRequired ('resDT'          , @iscell);
paramDefaults = {'Centerx0','Centery0','Theta','sigmaMinor','sigmaMajor'};
shortDefaults = {'x0','y0','Th','sMin','sMaj'};
p.addParameter('params', paramDefaults, @iscell);
p.addParameter('shortennames', false, @islogical);
p.parse(paramsDT, resNames, resDT, varargin{:});

params       = p.Results.params;
shortenNames = p.Results.shortennames;


%% Calculate
% Keep adding columns that will be the measures we will be using

% Decide what is the best structure to visualize this calculations
% Right now use it just to concatenate
% Later decide if this function returns another table with just the results

newDT = table();
newDT.synth = paramsDT.RF(:,params);
if shortenNames && isequal(params,paramDefaults)
    newDT.synth.Properties.VariableNames = shortDefaults;
end
for ii=1:length(resNames)
    newDT.(resNames{ii}) = resDT{ii}(:,params);
    newDT.(resNames{ii}).Properties.VariableNames = shortDefaults;
end


                        
