function newDT = pmTseriesPlot(DT, varargin)
% Takes the time series (predicted and original) and plots them, 
% It has option to show some metrics in the plot itself
% 
%  Inputs: table 
% 
%  Outputs: returns a table with the comparisons
% 
%  Example:
%{ 
    pmTseriesPlot(tSeries, 'to compare', {'synth','aPRF','vista'}, ...
                           'voxel',[1:5], ...
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
p.addRequired ('DT'                   , @istable);
p.addParameter('tocompare', {}        , @iscell);  % Default is do them all
p.addParameter('metric'   , 'RMSE'    , @ischar);  
p.addParameter('voxel'    , [1]       , @isnumeric);  
p.addParameter('newwin'   , true      , @islogical); % For subplot, set this to false
% Parse
p.parse(DT, varargin{:});
% Assign
toCompare = p.Results.tocompare;
metric    = p.Results.metric;
voxel     = p.Results.voxel;
w         = p.Results.newwin;

% If a metric was not introduced, not show it in the plot, otherwise, do
addMetric = true;
if contains('metric',p.UsingDefaults) 
    addMetric = false;
end

% Check if required
if isempty(toCompare)
    toCompare = DT.Properties.VariableNames;    
end
if length(toCompare) <= 1
   error('You need at least two analysis results to be able to compare them') 
end
% Check if all names are ok
if sum(ismember(toCompare,DT.Properties.VariableNames)) ~= length(toCompare)
    error('Some of the analysis names in <to compare> are not correct')
end


%% Calculate
if w;mrvNewGraphWin('decide plot name');end
if length(voxel) == 1
    plot(DT.synth.BOLDnoise(voxel,:), 'r'); 
    if length(toCompare) >= 2
        hold on;
        for ii=2:length(toCompare)
            plot(DT.(toCompare{ii}).modelpred(voxel,:));
            if addMetric
                aa = gca; text(5, (aa.YLim(1)*(1.0025)) + ii, ...
                    sprintf('RMSE(%s):%2.2f',toCompare{ii},99.99));
            end
        end
    end
    legend(toCompare);
    xlabel('time'); ylabel('relative amplitude')
    title(sprintf('Fit for voxel %i', voxel))
else
    numAcross = ceil(sqrt(length(voxel)));
    numDown   = ceil(length(voxel)/numAcross);
    for nv=1:length(voxel)
        subplot(numAcross,numDown,nv) % as they are time series, better to see them in vertical
        plot(DT.synth.BOLDnoise(voxel(nv),:), 'r');
        if length(toCompare) >= 2
            hold on;
            for ii=2:length(toCompare)
                plot(DT.(toCompare{ii}).modelpred(voxel(nv),:));
                if addMetric
                    aa = gca; text(5,aa.YLim(1)*ii*(1.0025), ...
                        sprintf('RMSE(%s):%2.2f',toCompare{ii},99.99));
                end
            end
        end
        legend(toCompare);
        xlabel('time'); ylabel('relative amplitude')
        title(sprintf('Fit for voxel %i', voxel(nv)))
    end
end
                        
