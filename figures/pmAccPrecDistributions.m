function pmAccPrecDistributions(DT, varToCompare, varargin)
%pmAccPrecDistributions Return a distribution plot for several noise values
%   Return a distribution plot for several noise values

%% Read the inputs
% Make varargin lower case, remove white spaces...
varargin = mrvParamFormat(varargin);
% Parse
p = inputParser;
p.addRequired('DT'                                      , @istable);
p.addRequired('varToCompare'                            , @ischar);
p.addParameter('fnameroot'       , 'defaultFigNameRoot' , @ischar);
p.addParameter('normalksdensity' , 'normal'             , @ischar);

% Parse. Assign result inside each case
p.parse(DT, varToCompare, varargin{:});
% Read here only the generic ones
fnameRoot = p.Results.fnameroot;
normalKsdensity = p.Results.normalksdensity;



%% Do the thing
% Select the values/categories we are going to filter to generate the different distributions
cats = unique(DT.(varToCompare));
c    = distinguishable_colors(length(cats));

% What is the value it should be? Create a vertical line
objectiveValue = unique(DT.synth.sMaj);
% There should be only one
if length(objectiveValue) > 1
    error('There should only one value, check the filtering of the table')
end

% Plot the vertical line with the correct RFsize
h1 = plot(objectiveValue*[1,1],[0 1],'Color','k','LineStyle','-.','LineWidth',3);
hold on;

% Plot the noiseless case, it will need to be another vertical line
x = DT.tool.sMaj(DT.(varToCompare)==cats(1));
[x_values, mu, sigma, mn, mx] = dr_distPlottingVals(x);
h2 = plot(mu*[1,1],[0 1],'Color',c(1,:),'LineStyle','-.','LineWidth',2);

% Loop over the other noise values
a = [];
for sc=2:2:length(cats)
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
    a = [a;plot(x_values, indiv_pdf, 'LineWidth',2, 'color', c(sc,:))];
    % hArrw = (length(cats)-sc+1)*(0.9*(ylims(2)/length(cats)));
    % fixed = 0.9*(ylims(2)/length(cats));
    h=plot(mu*[1,1],[0 max(indiv_pdf)],'LineStyle','-.','LineWidth',1,'Color',c(sc,:));

end
legend([h1;h2;a],['Real';num2str(cats)])

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

