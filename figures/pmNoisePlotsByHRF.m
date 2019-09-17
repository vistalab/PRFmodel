function pmNoisePlotsByHRF(compTable, tools,varargin)
%pmNoiseSizePlots Return a distribution plot for several noise values
%   Return a distribution plot for several noise values

%% Read the inputs
% Make varargin lower case, remove white spaces...
varargin = mrvParamFormat(varargin);
% Parse
p = inputParser;
p.addRequired('compTable');
p.addRequired('tools'                            , @iscell);
p.addParameter('fnameroot'  , 'defaultFigNameRoot' , @ischar);
p.addParameter('randomhrf'  , false                , @islogical);
p.addParameter('x0y0'       , [0,0]                , @isnumeric);
p.addParameter('sorthrf'    , {'same'}             , @iscell);
p.addParameter('usemetric'  , 'rfsize'             , @ischar);
p.addParameter('userfsize'  , [999]                    , @isnumeric);
p.addParameter('noisevalues', [999]                , @isnumeric);

% Parse. Assign result inside each case
p.parse(compTable, tools, varargin{:});
% Read here only the generic ones
fnameRoot   = p.Results.fnameroot;
randomhrf   = p.Results.randomhrf;
x0y0        = p.Results.x0y0;
sortHRF     = p.Results.sorthrf;
usemetric   = p.Results.usemetric;
userfsize   = p.Results.userfsize;
noisevalues = p.Results.noisevalues;

%% Data selection
if istable(compTable)
    % Calculate the polar angle and eccentricity
    [TH,R] = cart2pol(compTable.synth.x0, compTable.synth.y0);
    compTable.synth.angle = rad2deg(TH);
    compTable.synth.eccentricity = R;
    
    if noisevalues==999
        noises   = unique(compTable.noise2sig);
    else
        noises   = noisevalues;
    end
    switch usemetric
        case {'rfsize'}
            if userfsize == 999
                metrics = unique(compTable.synth.sMaj);
                % Reduce the table removing the locations
                compTable = compTable(compTable.synth.x0 == x0y0(1) & ...
                            compTable.synth.y0 == x0y0(2), :);
            else
                metrics = userfsize;
                % Reduce the table removing the locations
                compTable = compTable(compTable.synth.x0 == x0y0(1) & ...
                                      compTable.synth.y0 == x0y0(2) & ...
                                      ismember(compTable.synth.sMaj,metrics), :);
            end
        case {'polarangle'}
            metrics = unique(compTable.synth.angle);
            % Reduce the table removing the rfSize
            compTable = compTable(compTable.synth.sMaj == userfsize , :);
        case {'eccentricity'}
            metrics = unique(compTable.synth.eccentricity);
            % Reduce the table removing the rfSize
            compTable = compTable(compTable.synth.sMaj == userfsize, :);
        otherwise
            error('Metric %s not defined', usemetric)
    end
else
    error('Input compTable is not right, it needs to be a table')
end



%% PLOTS
colors   = distinguishable_colors(length(tools),'w');
mrvNewGraphWin([usemetric ' Accuracy and precision']);
plotNum = 0;
for np=1:length(metrics)
    for nn=1:length(noises)
        plotNum = plotNum +1;
        subplot(length(metrics), length(noises), plotNum);
        metric     = metrics(np);
        noiselvl    = noises(nn);
        a = [];b = [];
        for nt=1:length(tools)
            tool = tools{nt};
            % Calculate the polar angle and eccentricity by tool
            [TH,R] = cart2pol(compTable.(tool).x0, compTable.(tool).y0);
            compTable.(tool).angle = rad2deg(TH);
            compTable.(tool).eccentricity = R;
            
            % Reduce the table for only the tools and metric we want
            switch usemetric
                case {'rfsize','size'}
                    DT = compTable(compTable.synth.sMaj==metric & ...
                        compTable.noise2sig==noiselvl & ...
                        compTable.synth.x0==x0y0(1) & ...
                        compTable.synth.y0==x0y0(2) ,:);
                    [result,HRFs,sds] = pmNoise2AccPrec(DT, 'both', ...
                        'plotIt',false,'tool',tool, ...
                        'separatehrf',true, ...
                        'medianCI', true, ...
                        'sortHRF',sortHRF, ...
                        'usemetric', 'rfsize');
                case {'angle','polarangle'}
                    DT = compTable(compTable.synth.angle==metric & ...
                        compTable.noise2sig==noiselvl & ...
                        compTable.synth.sMaj==userfsize,:);
                    [result,HRFs,sds] = pmNoise2AccPrec(DT, 'both', ...
                        'plotIt',false,'tool',tool, ...
                        'separatehrf',true, ...
                        'medianCI', true, ...
                        'sortHRF',sortHRF, ...
                        'usemetric', 'polarangle');
                case {'eccen','eccentricity'}
                    DT = compTable(compTable.synth.eccentricity==metric & ...
                        compTable.noise2sig==noiselvl & ...
                        compTable.synth.sMaj==userfsize,:);
                    [result,HRFs,sds] = pmNoise2AccPrec(DT, 'both', ...
                        'plotIt',false,'tool',tool, ...
                        'separatehrf',true, ...
                        'medianCI', true, ...
                        'sortHRF',sortHRF, ...
                        'usemetric', 'eccentricity');
                otherwise
                    error('Metric %s not defined', usemetric)
            end
            % ADD SOME JITTER
            xvaluesCenter = 1:length(HRFs);
            xvalues = xvaluesCenter + 0.07*(nt-1);
            a = [a;scatter(xvalues, result, 30, colors(nt,:), 'filled')];hold on;
            % Add the confidence intervals in every point now
            % xxx = [1:length(HRFs);1:length(HRFs)]';
            xxx = [xvalues;xvalues]';
            for ns=1:length(sds)
                plot(xxx(ns,:),sds(ns,:),'Color',colors(nt,:),'LineStyle','-','LineWidth',1);
            end
        end

        h1 = plot([1,length(HRFs)],metric*[1,1],'Color','k','LineStyle','-.','LineWidth',1);
        
        grid
        % xticks([])
        set(gca,'box','off','color','none')
        % ylabel('Acc. & prec. (deg)','FontSize',14,'FontWeight','bold','Color','k');
        ylabel(usemetric,'FontSize',14,'FontWeight','bold','Color','k');
        if (metric==2 && noiselvl==0)
            title(sprintf('(HRF SORT) %s: %1.1f | Noise: %0.1f',usemetric,metric,noiselvl));
        else
            title(sprintf('%s: %1.1f | Noise: %0.1f',usemetric,metric,noiselvl));
        end
        if plotNum == 1
            legend([h1;a],['synth',tools],'Location','westoutside');
        end
        
        % In the bottom row, add the HRF names
        if np==length(metrics)
            xticks(xvaluesCenter)
            xticklabels(strrep(HRFs,'_','\_'));xtickangle(45)
            % set(gca,'TickLength',[0 0])
            ax = gca;
            % ax.XGrid = 'off';
            xlabel('HRFs','FontSize',18,'FontWeight','bold','Color','k');
        end
        
    end
    


end

% Add title to the whole thing
Theta = unique(compTable.synth.Th);
TR    = unique(compTable.TR);
if strcmp(usemetric, 'rfsize')
    dr_suptitle(sprintf('%s HRF Comparison | Location=[%i,%i] | Theta=%i | TR=%1.2f | CI = 50%', usemetric, x0y0(1), x0y0(2), Theta, TR));
else
    dr_suptitle(sprintf('%s HRF Comparison | Theta=%i | TR=%1.2f | CI = 50%', usemetric, Theta, TR));
end











%% Review and delete
%{

% Plot the vertical line with the correct RFsize
h1 = plot(objectiveValue*[1,1],[0 1],'Color','k','LineStyle','-.','LineWidth',3);
hold on;

% Plot the noiseless case, it will need to be another vertical line
x = DT.tool.sMaj(DT.(varToCompare)==cats(1));
[x_values, mu, sigma, mn, mx] = dr_distPlottingVals(x);
h1 = plot(objectiveValue*[1,1],[0 1],'Color','b','LineStyle','-.','LineWidth',2);

% Loop over the other noise values
a = [];
for sc=2:length(cats)
    % Create the x value and the distributions
    x = DT.tool.sMaj(DT.(varToCompare)==cats(sc));
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






