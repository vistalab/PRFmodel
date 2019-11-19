function [result,sliceVals,sds] = pmNoise2AccPrec(DT, accorprec, varargin)
%pmAccPrecDistributions Return a distribution plot for several noise values
%   Return a distribution plot for several noise values

%% Read the inputs
% Make varargin lower case, remove white spaces...
varargin = mrvParamFormat(varargin);
% Parse
p = inputParser;
p.addRequired('DT'                               , @istable);
p.addRequired('accorprec'                        , @ischar);
p.addParameter('fnameroot'  , 'defaultFigNameRoot' , @ischar);
p.addParameter('plotit'     , true                 , @islogical);
p.addParameter('tool'       , 'first'              , @ischar);
p.addParameter('separatehrf', false                , @islogical);
p.addParameter('medianci'   , false                , @islogical);
p.addParameter('sorthrf'  , {'same'}               , @iscell);
p.addParameter('usemetric', 'rfsize'               , @ischar);
p.addParameter('cirange'  , 50                     , @isnumeric);



% Parse. Assign result inside each case
p.parse(DT, accorprec, varargin{:});
% Read here only the generic ones
fnameRoot   = p.Results.fnameroot;
plotIt      = p.Results.plotit;
tool        = p.Results.tool;
separatehrf = p.Results.separatehrf;
medianCI    = p.Results.medianci;
sortHRF     = p.Results.sorthrf;
usemetric   = p.Results.usemetric;
CIrange     = p.Results.cirange;


%% Do the thing
% Select the values/categories we are going to filter to generate the different distributions
noiseVals = unique(DT.noiseLevel)';

if length(sortHRF)==1 && strcmp(sortHRF{1},'same')
    HRFs      = unique(DT.HRFtype)';
else
    % TODO: compare that we are passing exactly the same HRF names, none missing
    HRFs      = sortHRF;
end

if separatehrf
    sliceVals = HRFs;
else
    sliceVals = noiseVals;
end

% Select what variable we are using, the rest will come removed
switch usemetric
    case {'rfsize'}
        usemetriccol = 'sMaj';
    case {'polarangle'}
        usemetriccol = 'angle';
    case {'eccentricity'}
        usemetriccol = 'eccentricity';
    otherwise
        error('Metric %s not defined', usemetric)
end

 
% What is the value it should be? Create a horiz line
objectiveValue = unique(DT.synth.(usemetriccol));

% There should be only one
if length(objectiveValue) > 1
    error('There should only one value, check the filtering of the table')
end

% Select the tool
if strcmp(tool,'first')
    DT.Properties.VariableNames{5} = 'tool';
elseif ismember(tool, DT.Properties.VariableNames)
    loc = find(ismember(DT.Properties.VariableNames, tool), 1);
    DT.Properties.VariableNames{loc} = 'tool';
else
    error('%s tool not found',tool)
end

% Loop over the noise values
result = [];
sds    = [];
for sc=1:length(sliceVals)
    % This is what it should be:
    objectiveValue = unique(DT.synth.(usemetriccol));
    
    % Create the x value and the distributions
    if separatehrf
        x = DT.tool.(usemetriccol)(strcmp(DT.HRFtype,sliceVals{sc}));
    else
        x = DT.tool.(usemetriccol)(DT.noiseLevel == sliceVals(sc));
    end
    
    % Calculate the summary stats
    [x_values, mu, sigma, mn, mx, med, ci] = dr_distPlottingVals(x,CIrange);
    
    % Obtain precision or accuracy
    accorprec = mrvParamFormat(accorprec);
    switch accorprec
        case {'both'}
            % Create the results    
            if medianCI
                % Create the vector of results
                result    = [result, med];
                % Create the vector of SDs
                sds       = [sds   ; ci];
            else
                % Create the vector of results
                result    = [result, mu];
                % Create the vector of SDs
                sds       = [sds   , sigma];
            end
            
            % Create strings for the plot
            xTitle    = 'Noise levels';
            yTitle    = 'Abs difference from the mean (deg)';
            plotTitle = ['ACCURACY (rfSize:' num2str(objectiveValue) ')'];

        case {'accuracy','acc'}
            % Obtain the difference from the mean
            meanDiff  = abs(objectiveValue - mu);
            % Create the vector of results
            result    = [result, meanDiff];
            % Create strings for the plot
            xTitle    = 'Noise levels';
            yTitle    = 'Abs difference from the mean (deg)';
            plotTitle = ['ACCURACY (rfSize:' num2str(objectiveValue) ')'];
        case {'precision','prec'}
            % Create the vector of results
            result = [result, sigma];
            % Create strings for the plot
            xTitle    = 'Noise levels';
            yTitle    = 'SD of the results (deg)';
            plotTitle = ['PRECISION (rfSize:' num2str(objectiveValue) ')'];
        otherwise
            error('%s not recognized',accorprec)
    end 
end


if plotIt
    % Create the plot
    plot(sliceVals, result, 'k');
    xlabel(xTitle); ylabel(yTitle); title(plotTitle);
    % xlim(myXlim); ylim(myYlim);
end


%{

% Plot the vertical line with the correct RFsize
h1 = plot(objectiveValue*[1,1],[0 1],'Color','k','LineStyle','-.','LineWidth',3);
hold on;

% Plot the noiseless case, it will need to be another vertical line
x = DT.tool.(usemetriccol)(DT.(varToCompare)==cats(1));
[x_values, mu, sigma, mn, mx] = dr_distPlottingVals(x);
h1 = plot(objectiveValue*[1,1],[0 1],'Color','b','LineStyle','-.','LineWidth',2);

% Loop over the other noise values
a = [];
for sc=2:length(cats)
    % Create the x value and the distributions
    x = DT.tool.(usemetriccol)(DT.(varToCompare)==cats(sc));
    [x_values, mu, sigma, mn, mx] = dr_distPlottingVals(x);

    % Plot the distributions, normal or ksdensity
    if strcmp(normalKsdensity, 'normal')
        pd = fitdist(x,'normal');
        indiv_pdf = pdf(pd, x_values);
    else
        [indiv_pdf, x_values] = ksdensity(x);
    end
    a = [a;plot(x_values, indiv_pdf, 'LineWidth',2)]; % , 'color', catcolors{sc,:})];
    % hArrw = (length(cats)-sc+1)*(0.9*(ylims(2)/length(cats)));
    % fixed = 0.9*(ylims(2)/length(cats));
    % h1=plot(mu*[1,1],[0 hArrw],'LineStyle','-.','LineWidth',1);  % 'Color',catcolors{sc,:}

end

%}

%{
    xlabel('\DeltaFA', 'FontWeight','bold');
    % ylim([0.1, 0.75]); yticks([0.2,.4,.6])
    set(gca,'FontSize',18);
    title(sprintf('%s',tn));ylim(ylims);xlim(xlims);
    cats      = categories(unsMeans.SliceCats);
    catcolors = unique(unsMeans.SliceCatsRGB);
    for sc=1:length(cats)
        


        % Now calculate the effect size with bootstrapping and mark it in plot
        stats = mes( X, x, ...
            'hedgesg','isDep',0,'nBoot',nBoot, ...
            'doPlot',0,'missVal','listwise','confLevel',.95,'ROCtBoot',false);
        if stats.t.p > 0.05; texto = 'n.s.';
        else;                texto = sprintf('d=%6.3f',stats.hedgesg); end
        drawArrow([mu MU],[hArrw,hArrw],{'Color',catcolors{sc,:},'LineWidth',.7, ...
            'string',texto, 'FontSize', 14});
    end
    % legend(a,[{'All Samples'};strrep(cats,'_','\_')],'Location','NorthEast'); hold off;
    

    
    % {
    % Same plot for TEST-RETEST
    HCPunsMeans=unsMeans(unsMeans.Proj=='HCP' & unsMeans.TRT=='TEST',:);
    HCPunsMeans.SliceCats = removecats(HCPunsMeans.SliceCats);
    HCPunsMeansrt=unsMeans(unsMeans.Proj=='HCP' & unsMeans.TRT=='RETEST',:);
    HCPunsMeansrt.SliceCats = removecats(HCPunsMeansrt.SliceCats);
    
    cats      = categories(HCPunsMeans.SliceCats);
    catcolors = unique(HCPunsMeansrt.SliceCatsRGB);
    catsrt      = categories(HCPunsMeansrt.SliceCats);
    catcolorsrt = catcolors;
    
    
    
    nBoot     = 5000;
    normalKsdensity = 'normal';
    figure('Name',fnameRoot, ...
        'NumberTitle','off', ...
        'visible',   'on', ...
        'color','w', ...
        'Units','pixel', ...
        'Position',[0 0 1900 1100]);
    nrow=nrowcol(1); ncol=nrowcol(2);
    for nt = 1: length(tracts)
        subplot(nrow,ncol,nt)
        tn = char(tracts());
        a = [];
        l = {};
        for sc=1:length(cats)
            x   = HCPunsMeans.(tn)(HCPunsMeans.SliceCats==cats{sc});
            xrt = HCPunsMeansrt.(tn)(HCPunsMeansrt.SliceCats==catsrt{sc});
            [x_values, mu, sigma, mn, mx] = dr_distPlottingVals(x);
            [x_valuesrt, murt, sigmart, mnrt, mxrt] = dr_distPlottingVals(xrt);
            if strcmp(normalKsdensity, 'normal')
                pd = fitdist(x,'normal');
                indiv_pdf = pdf(pd, x_values);
                pdrt = fitdist(xrt,'normal');
                indiv_pdfrt = pdf(pdrt, x_valuesrt);
                mySupTitle = 'Normal distribution of mean FA-s per project/tract.';
            else
                [indiv_pdf, x_values] = ksdensity(x);
                [indiv_pdfrt, x_valuesrt] = ksdensity(xrt);
                mySupTitle = 'Density distribution of mean FA-s per project/tract.';
            end
            at  = plot(x_values, indiv_pdf, 'LineWidth',2, 'color',catcolors{sc,:});hold on;
            art = plot(x_valuesrt, indiv_pdfrt, 'LineWidth',2, 'color', catcolors{sc,:}, 'LineStyle',':');
            xlim([0, .8]);ylim([0, 20]);
            xticks([0:0.2:.8]);yticks([0:5:20]);
            title(tn)
            a = [a;at;art];  % Concatenate plots to be used with legend
            l = [l; strcat(strrep(cats{sc},'_','\_'),'TEST'); ...
                strcat(strrep(cats{sc},'_','\_'),'RETEST')];
            % Now calculate the effect size with bootstrapping and mark it in plot
            %         stats = mes( x, xrt, ...
            %                      'hedgesg','isDep',0,'nBoot',nBoot, ...
            %                      'doPlot',0,'missVal','listwise','confLevel',.95,'ROCtBoot',false);
            %         if stats.t.p > 0.05; sig = 'n.s.';
            %         else;                sig = sprintf('d=%6.2f',stats.t.p); end
            %         texto = sprintf('d=%6.3f (%s)',stats.hedgesg,sig);
            %         hArrw = max(indiv_pdf);
            %         hArrwrt = max(indiv_pdfrt);
            %         h1=plot(mu*[1,1],[0 hArrw],'Color',catcolors{sc,:},'LineStyle','-','LineWidth',1);
            %            plot(murt*[1,1],[0 hArrwrt],'Color',catcolors{sc,:},'LineStyle',':','LineWidth',1);
            %         drawArrow([mu murt],[hArrw,hArrw],{'Color',catcolors{sc,:},'LineWidth',1,'string',texto});
        end
        xlabel('FA', 'FontWeight','bold');
        set(gca,'FontSize',18)
        title(sprintf('%s',tn))
        
        % legend(a, l,'Location','SouthWest'); hold off;
    end


if HCPTRT
    suptitle({strrep(mySupTitle,'_','\_'), 'HCP TEST-RETEST'})
else
    suptitle({strrep(mySupTitle,'_','\_')})
end



%% Save it
if ~exist(saveItHere,'dir')
    mkdir(saveItHere)
end
if saveSvg
    saveas(gcf,fullfile(saveItHere, strcat(fnameRoot,'.svg')),'svg');
end
if savePng
    saveas(gcf,fullfile(saveItHere, strcat(fnameRoot,'.png')),'png');
end

%}

end

