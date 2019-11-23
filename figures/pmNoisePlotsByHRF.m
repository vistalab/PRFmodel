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
p.addParameter('userfsize'  , [999]                , @isnumeric);
p.addParameter('noisevalues', {'all'}              , @iscell);
p.addParameter('ylims'      , [0,0]                , @isnumeric);
p.addParameter('cirange'    , 50                   , @isnumeric);
p.addParameter('saveto'     , ''                   , @ischar);
p.addParameter('savetotype' , 'png'                , @ischar);
p.addParameter('fontsize'   , 14                   , @isnumeric);

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
ylims       = p.Results.ylims;
CIrange     = p.Results.cirange; 
saveTo      = p.Results.saveto;
saveToType  = p.Results.savetotype; 
fontSize    = p.Results.fontsize; 

%% Data selection
if istable(compTable)
    % Calculate the polar angle and eccentricity
    [TH,R] = cart2pol(compTable.synth.x0, compTable.synth.y0);
    compTable.synth.angle = rad2deg(TH);
    compTable.synth.eccentricity = R;
    if length(noisevalues)==1
        if strcmp(noisevalues{:}, 'all')
            noises   = unique(compTable.noiseLevel);
        else
            noises   = noisevalues;
        end
    else
        noises   = noisevalues;
    end
    if length(noises)==4,
        noises={"none","low","mid","high"};
        noisetxt = strcat(noises{:});
    end
    if length(noises)==3,
        noises={"low","mid","high"};
        noisetxt = strcat(noises{:});
    end
    if length(noises)==1
        noisetxt = noises{:};
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
            if isempty(compTable)
                warning('No values, set up the values for the rfSize')
            end
        case {'eccentricity'}
            metrics = unique(compTable.synth.eccentricity);
            % Reduce the table removing the rfSize
            compTable = compTable(compTable.synth.sMaj == userfsize, :);
            if isempty(compTable)
                warning('No values, set up the values for the rfSize')
            end                       
        otherwise
            error('Metric %s not defined', usemetric)
    end
else
    error('Input compTable is not right, it needs to be a table')
end



%% PLOTS
colors   = distinguishable_colors(length(tools),'w');
colors   = colors .* 0.75;
if length(noisevalues)==1
    bigfighandler = mrvNewGraphWin([usemetric ' Accuracy and precision'],'tall');
else
    bigfighandler = mrvNewGraphWin([usemetric ' Accuracy and precision'],'wide');
    set(bigfighandler,'Position',[0.007 0.62  0.7  0.75]);
end
plotNum = 0;
for np=1:length(metrics)
    for nn=1:length(noises)
        plotNum = plotNum +1;
        subplot(length(metrics), length(noises), plotNum);
        metric     = metrics(np);
        noiselvl    = string(noises(nn));
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
                        compTable.noiseLevel==noiselvl & ...
                        compTable.synth.x0==x0y0(1) & ...
                        compTable.synth.y0==x0y0(2) ,:);
                    [result,HRFs,sds] = pmNoise2AccPrec(DT, 'both', ...
                        'plotIt',false,'tool',tool, ...
                        'separatehrf',true, ...
                        'medianCI', true, ...
                        'sortHRF',sortHRF, ...
                        'usemetric', 'rfsize', ...
                        'CIrange',CIrange);
                case {'angle','polarangle'}
                    DT = compTable(compTable.synth.angle==metric & ...
                        compTable.noiseLevel==noiselvl & ...
                        compTable.synth.sMaj==userfsize,:);
                    [result,HRFs,sds] = pmNoise2AccPrec(DT, 'both', ...
                        'plotIt',false,'tool',tool, ...
                        'separatehrf',true, ...
                        'medianCI', true, ...
                        'sortHRF',sortHRF, ...
                        'usemetric', 'polarangle', ...
                        'CIrange',CIrange);
                case {'eccen','eccentricity'}
                    DT = compTable(compTable.synth.eccentricity==metric & ...
                        compTable.noiseLevel==noiselvl & ...
                        compTable.synth.sMaj==userfsize,:);
                    [result,HRFs,sds] = pmNoise2AccPrec(DT, 'both', ...
                        'plotIt',false,'tool',tool, ...
                        'separatehrf',true, ...
                        'medianCI', true, ...
                        'sortHRF',sortHRF, ...
                        'usemetric', 'eccentricity', ...
                        'CIrange',CIrange);
                otherwise
                    error('Metric %s not defined', usemetric)
            end
            % ADD SOME JITTER
            xvaluesCenter = 1:length(HRFs);
            xvalues = xvaluesCenter + 0.07*(nt-1);
            a = [a;scatter(xvalues, result, 50, colors(nt,:),'filled')];hold on;
            % Add the confidence intervals in every point now
            % xxx = [1:length(HRFs);1:length(HRFs)]';
            xxx = [xvalues;xvalues]';
            for ns=1:length(sds)
                plot(xxx(ns,:),sds(ns,:),'Color',colors(nt,:),'LineStyle','-','LineWidth',1);
            end
        end

        h1 = plot([1,length(HRFs)+1],metric*[1,1],'Color','k','LineStyle','-.','LineWidth',1);
        
        grid
        % xticks([])
        set(gca,'box','off','color','none','FontSize',fontSize)
        % ylabel('Acc. & prec. (deg)','FontSize',14,'FontWeight','bold','Color','k');
        ylabel(usemetric,'FontSize',fontSize,'FontWeight','bold','Color','k');
        if (metric==2 && noiselvl=="none")
            title(sprintf('%s: %1.1f | Noise: %s',usemetric,metric,noiselvl));
        else
            title(sprintf('%s: %1.1f | Noise: %s',usemetric,metric,noiselvl));
        end
        if plotNum == 1
            legend([h1;a],[{'synth'};tools],'Location','best') % 'westoutside');
        end
        
        % In the bottom row, add the HRF names
        if np==length(metrics)
            xticks(xvaluesCenter)
            HRFlabels = strrep(HRFs,'popeye_twogammas','pop_twogammas');
            xticklabels(strrep(HRFlabels,'_','\_'));xtickangle(25)
            % set(gca,'TickLength',[0 0])
            ax = gca;
            % ax.XGrid = 'off';
            xlabel('HRFs','FontSize',18,'FontWeight','bold','Color','k');
        end
        % Modify ylims
        if ylims(1) ~= ylims(2)
            ylim(ylims);
        end
    end
end

% Add title to the whole thing
Theta = unique(compTable.synth.Th);
TR    = unique(compTable.TR);
set(gca,'FontSize',fontSize)
if strcmp(usemetric, 'rfsize')
    dr_suptitle(sprintf('%s HRF Comparison | Location=[%i,%i] | Theta=%i | TR=%1.2f | CI = %i%', usemetric, x0y0(1), x0y0(2), Theta, TR, CIrange));
    fnameRoot = sprintf('Fig_%s_Noise-%s_HRF_Comparison_Location%i%i_Theta%i_TR=%1.2f_CI%i', usemetric, noisetxt, x0y0(1), x0y0(2), Theta, TR, CIrange);
else
    dr_suptitle(sprintf('%s HRF Comparison | Theta=%i | TR=%1.2f | CI = %i%', usemetric, Theta, TR, CIrange));
    fnameRoot = sprintf('Fig_%s_Noise-%s_HRF_Comparison_Theta%i_TR%1.2f_CI%i', usemetric, noisetxt, Theta, TR, CIrange);
end

%% Save it
if ~isempty(saveTo)
    if ~exist(saveTo,'dir')
        mkdir(saveTo)
    end
    switch saveToType
        case 'svg'
            saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.svg')),'svg');
        case 'png'
            saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.png')),'png');
        otherwise
            error('File type %s does not exist, use png or svg',saveToType)
    end

end

end






