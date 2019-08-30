function pmNoiseSizePlots(compTable, tools,varargin)
%pmNoiseSizePlots Return a distribution plot for several noise values
%   Return a distribution plot for several noise values

%% Read the inputs
% Make varargin lower case, remove white spaces...
varargin = mrvParamFormat(varargin);
% Parse
p = inputParser;
p.addRequired('compTable');
p.addRequired('tools'                            , @iscell);
p.addParameter('fnameroot', 'defaultFigNameRoot' , @ischar);
p.addParameter('randomhrf', false                , @islogical);


% Parse. Assign result inside each case
p.parse(compTable, tools, varargin{:});
% Read here only the generic ones
fnameRoot = p.Results.fnameroot;
randomhrf = p.Results.randomhrf;


%% Do the thing
if istable(compTable)
    doDiff   = false;
    prfsizes = unique(compTable.synth.sMaj);
elseif iscell(compTable) && length(compTable)==2
    doDiff   = true;
    prfsizes = unique(compTable{1}.synth.sMaj);
else
    error('Input compTable is not right, it needs to be a table or a cell with two tables')
end


colors   = distinguishable_colors(length(tools),'w');

mrvNewGraphWin(['Accuracy and precision']);

for np=1:length(prfsizes)
    subplot(length(prfsizes), 1, np);
    prfsize      = prfsizes(np);
   a = [];b = [];
   for nt=1:length(tools)
       tool = tools{nt};
       % Reduce the table for only the RF size-s we want, and the tool we want
       if doDiff
           DT1           = compTable{1}(compTable{1}.synth.sMaj==prfsize,:);
           DT2           = compTable{2}(compTable{2}.synth.sMaj==prfsize,:);
           [result1,noiseVals1,sds1] = pmNoise2AccPrec(DT1, 'both','plotIt',false,'tool',tool);
           [result2,noiseVals2,sds2] = pmNoise2AccPrec(DT2, 'both','plotIt',false,'tool',tool);
           result        = result1 - result2;
           sds           = abs(sds1    - sds2);
           if isequal(noiseVals1, noiseVals2)
               noiseVals     = noiseVals1;
           else
               error('Difference amount of noise values in each table')
           end
           
       else
           DT           = compTable(compTable.synth.sMaj==prfsize,:);
           [result,noiseVals,sds] = pmNoise2AccPrec(DT, 'both','plotIt',false,'tool',tool);
       end
       
       
       a = [a;plot(noiseVals,result+sds,'Color',colors(nt,:),'LineStyle','-','LineWidth',2)];hold on;
       b = [b;plot(noiseVals,result-sds,'Color',colors(nt,:),'LineStyle','-','LineWidth',2)];
       
       jbfill(noiseVals,result+sds,result-sds,colors(nt,:),colors(nt,:),1,.05); hold on;
   end
    if doDiff
        h1 = plot([noiseVals(1),noiseVals(end)],[0,0],'Color','k','LineStyle','-.','LineWidth',1);
    else
        h1 = plot([noiseVals(1),noiseVals(end)],prfsize*[1,1],'Color','k','LineStyle','-.','LineWidth',1);
    end
    
    grid
    xticks([])
    set(gca,'box','off','color','none')
    ylabel('Acc. & prec. (deg)','FontSize',14,'FontWeight','bold','Color','k');
    % title(['ACCURACY AND PRECISION (vs noise and RF size)']);
    legend([h1;a],['synth',tools],'Location','southwest');
    
end
xticks(noiseVals)
% set(gca,'TickLength',[0 0])
ax = gca;
ax.XGrid = 'off';
xlabel('Noise levels','FontSize',18,'FontWeight','bold','Color','k');
if doDiff;compTable = compTable{1};
x0    = unique(compTable.synth.x0);
y0    = unique(compTable.synth.x0);
Theta = unique(compTable.synth.Th);
TR    = unique(compTable.TR);
if doDiff
        dr_suptitle(sprintf('Ideal - Random HRF | Location=[%i,%i] | Theta=%i | TR=%1.2f | Band= 1 SD', x0, y0, Theta, TR));
else
    if randomhrf
        dr_suptitle(sprintf('Random HRF | Location=[%i,%i] | Theta=%i | TR=%1.2f | Band= 1 SD', x0, y0, Theta, TR));
    else
        dr_suptitle(sprintf('Ideal HRF | Location=[%i,%i] | Theta=%i | TR=%1.2f | Band= 1 SD', x0, y0, Theta, TR));
    end
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






