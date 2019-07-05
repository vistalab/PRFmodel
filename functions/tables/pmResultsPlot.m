function newDT = pmResultsPlot(DT, varargin)
% Takes the aggregated table of results and plots them, showing optionaly a metric
% 
%  Inputs: table 
% 
%  Outputs: returns a table with the comparisons
% 
%  Example:
%{ 
    pmResultsPlot(compTable, 'to compare', {'synth','aPRF'}, ...
                             'result','Centerx0', ...
                             'metric','RMSE', ...
                             'newWin',true)
%}
% 
%  See also: 
% 
%  GLU Vistalab 2019.07


%% Read the inputs
varargin = mrvParamFormat(varargin);
p = inputParser;
p.addRequired ('DT'           , @istable);
p.addParameter('tocompare', {}, @iscell);  % Default is do them all
p.addParameter('result'   , 'Centerx0', @ischar);  
p.addParameter('metric'   , 'RMSE'    , @ischar);  
p.addParameter('newwin'   , true      , @islogical); % For subplot, set this to false
% Parse
p.parse(DT, varargin{:});
% Assign
toCompare = p.Results.tocompare;
result    = p.Results.result;
metric    = p.Results.metric;
w         = p.Results.newwin;
% Check if required
if isempty(toCompare)
    toCompare = DT.Properties.VariableNames;    
end
if length(toCompare) <= 2
   error('You need at least two analysis results to be able to compare them') 
end
% Check if all names are ok
if sum(ismember(toCompare,DT.Properties.VariableNames)) ~= length(toCompare)
    error('Some of the analysis names in <to compare> are not correct')
end

% Check if the result is present in all tables, otherwise, send error
for ii=1:length(toCompare)
    if ~ismember(result, DT.(toCompare{ii}).Properties.VariableNames)
        error('%s  is not part of the results, at least in %s table. Check if the names were shortened in the previous step.', result, toCompare{ii})
    end
end

%% Calculate
% Keep adding columns that will be the measures we will be using

if w;mrvNewGraphWin('decide plot name');end
scatter(DT.(toCompare{1}).(result), DT.(toCompare{2}).(result),40,'r'); 
if length(toCompare) >= 2
    hold on;
    for ii=3:length(toCompare)
        scatter(DT.(toCompare{1}).(result), DT.(toCompare{ii}).(result),40,'b');
    end
end

                        
