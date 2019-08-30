function [newDT, outlierReport] = pmResultsCleanOutliers(DT, varargin)
% outlier cleaning tool, that will output:
%   - Cleaned result table. 
%   - A table with outlier counts, that will evaluate quality of analysis 
% 
%  Inputs: compTable generated with pmResultsCompare()
% 
%  Outputs: returns the new table and the report on the outliers
% 
%  See also: 
% 
%  GLU Vistalab 2019.08


%% Read the inputs
varargin = mrvParamFormat(varargin);
p = inputParser;
p.addRequired ('DT'          , @istable);

p.addParameter('tools'       , {'all'}, @iscell);
p.addParameter('rfsizes'     , 9999   , @isnumeric);
p.addParameter('noiselevels' , 9999   , @isnumeric);
p.addParameter('outlierlevel', 20     , @isnumeric);
% Parse
p.parse(DT, varargin{:});
% Read the params
tools        = p.Results.tools;
RFsizes      = p.Results.rfsizes;
NoiseLevels  = p.Results.noiselevels;
OutlierLevel = p.Results.outlierlevel;


% Validate the inputs
if strcmp(tools, {'all'})
    tools = {DT.Properties.VariableNames{5:end}};
else
    if nnz(ismember(tools,{DT.Properties.VariableNames{5:end}})) < length(tools)
        error('Some of the tools are not in the database, values can be %s',...
            cell2str({DT.Properties.VariableNames{5:end}}))
    end
end
if RFsizes == 9999
    RFsizes     = unique(DT.synth.sMaj);
else
    if nnz(ismember(RFsizes,unique(DT.synth.sMaj))) < length(RFsizes)
        error('Some of the RF sizes are not in the database, values can be %s',cell2str({unique(DT.synth.sMaj)}))
    end
end

if NoiseLevels == 9999
    NoiseLevels = unique(DT.noise2sig); 
else
    if nnz(ismember(NoiseLevels,unique(DT.noise2sig))) < length(NoiseLevels)
        error('Some of the noise levels are not in the database, values can be %s',cell2str({unique(DT.noise2sig)}))
    end
end

%% Calculate

% Create a copy of the DT that will be updated:
newDT = DT;
% Create an empty datatable for the report
[p,q]  = meshgrid(RFsizes, NoiseLevels);
alto   = length(p(:));
ancho  = 2 + length(tools);
outlierReport = array2table(nan(alto, ancho),'VariableNames',...
                     {'RFsize','n2sig',tools{:}});
outlierReport.RFsize = p(:);    
outlierReport.n2sig  = q(:);


for checkTool = tools
    for rf = RFsizes'
        for nlvl = NoiseLevels'
            % Select the values of interest
            values = DT.(checkTool{:}).sMaj(DT.synth.sMaj==rf & DT.noise2sig==nlvl);
            % There are 3 values way off, at exactly 214.8481. Check what must have happened here.
            % Identify the outliers (randomly selected 20)
            indGreaterThan  = DT.(checkTool{:}).sMaj(DT.synth.sMaj==rf & DT.noise2sig==nlvl) > OutlierLevel;
            % Update the report
            outlierReport.(checkTool{:})(outlierReport.RFsize==rf & outlierReport.n2sig==nlvl) = sum(indGreaterThan);
            % If there at least one outlier, update the data
            if sum(indGreaterThan) > 0
                % Obtain the mean removing the outliers
                meanNoOutl      = mean(values(~indGreaterThan));
                % Create new values vector
                values(indGreaterThan) = repelem(meanNoOutl,sum(indGreaterThan));
                % Assign values
                newDT.(checkTool{:}).sMaj(newDT.synth.sMaj==rf & newDT.noise2sig==nlvl) = values;
            end
        end
    end
end

end






