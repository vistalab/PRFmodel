function pmCloudOfResults(compTable, tools,varargin)
%pmCloudOfResults Return a distribution plot for several noise values
%   It return only one plot, it need to be composited with another script

%% Read the inputs
% Make varargin lower case, remove white spaces...
varargin = mrvParamFormat(varargin);
% Parse
p = inputParser;
p.addRequired('compTable');
p.addRequired('tools'                    , @iscell);
p.addParameter('centerdistr'  , true    , @islogical);
p.addParameter('onlycenters'  , false    , @islogical);
p.addParameter('location'     , 'all');
p.addParameter('userfsize'    , 1        , @isnumeric);
p.addParameter('userfsizemin' , 1        , @isnumeric);
p.addParameter('centerperc'   , 50       , @isnumeric);
p.addParameter('usehrf'       , 'mix'    , @ischar);
p.addParameter('linestyle'    , '-'      , @ischar);
p.addParameter('linewidth'    , .7       , @isnumeric);
p.addParameter('fontsize'     , 14       , @isnumeric);
p.addParameter('newwin'       , true     , @islogical);
p.addParameter('noiselevel'   , 'none'   , @ischar);
p.addParameter('saveto'       , ''       , @ischar);
p.addParameter('savetotype'   , 'png'    , @ischar);
p.addParameter('color'        , 'old'             );
p.addParameter('addcibar'     , false    , @islogical);
p.addParameter('addcihist'    , false    , @islogical);
p.addParameter('useellipse'   , false    , @islogical);
p.addParameter('addtext'      , true     , @islogical);
p.addParameter('adddice'      , false     , @islogical);
p.addParameter('addsnr'       , false     , @islogical);
p.addParameter('xlims'        , [.8,5]    , @isnumeric);
p.addParameter('ylims'        , [.8,5]    , @isnumeric);
p.addParameter('xtick'        , [2:4]    , @isnumeric);
p.addParameter('ytick'        , [2:4]    , @isnumeric);

% Parse. Assign result inside each case
p.parse(compTable, tools, varargin{:});
% Read here only the generic ones
onlyCenters   = p.Results.onlycenters;
location      = p.Results.location;
centerDistr   = p.Results.centerdistr;
userfsize     = p.Results.userfsize;
if contains('userfsizemin',p.UsingDefaults)
    userfsizemin  = userfsize;
else
   userfsizemin  = p.Results.userfsizemin; 
end
centerPerc    = p.Results.centerperc;
useHRF        = p.Results.usehrf;
lineStyle     = p.Results.linestyle;
lineWidth     = p.Results.linewidth;
fontsize      = p.Results.fontsize;
newWin        = p.Results.newwin;
noiseLevel    = p.Results.noiselevel; noiseLevel  = string(noiseLevel);  
saveTo        = p.Results.saveto;
saveToType    = p.Results.savetotype; 
color         = p.Results.color; 
addcibar      = p.Results.addcibar; 
addcihist     = p.Results.addcihist; 
useEllipse    = p.Results.useellipse; 
addtext       = p.Results.addtext; 
adddice       = p.Results.adddice; 
addsnr        = p.Results.addsnr; 
xlims         = p.Results.xlims; 
ylims         = p.Results.ylims; 
xtick         = p.Results.xtick; 
ytick         = p.Results.ytick; 

%% Do the thing
if color=='old'
    Cs               = distinguishable_colors(1+length(tools),'w');
    Cs(2,:) = Cs(2,:) * 0.65;
else
    Cs               = distinguishable_colors(1+length(tools),'w');
    Cs(2,:)          = color;
end

if newWin
    mrvNewGraphWin(sprintf('Clouds, %s',useHRF),'wide');
end
% Inside function operations
% Check that only one size value has been passed
if length(userfsize)>1
    error('Only one size can be plotted at a time')
end
% Check percentage is 100 based
if centerPerc < 1; centerPerc = centerPerc*100; end
% Define the required confidence intervals as two percentiles
twoTailedRange = (100 - centerPerc) / 2;

% Filter the results table
dt = compTable(compTable.synth.sMaj == userfsize , :);
dt = dt(dt.synth.sMin == userfsizemin , :);

% Filter only one noise level
dt = dt(dt.noiseLevel == noiseLevel, :);

% Filter the selected location
if strcmp(location, 'all')
    % Think what to do and how to split
else
    dt = dt(dt.synth.x0 == location(1), :);
    dt = dt(dt.synth.y0 == location(2), :);
end

% Filter the selected HRF, if it is 'mix', use all HRFs together
if ~strcmp(useHRF,'mix')
    dt = dt(dt.HRFtype == string(useHRF), :);
end
toolLegend= {'synth'};
if length(tools) > 1
    error('Use only one tool, more than one tool are not clear to see')
end
for nt=1:length(tools)
    tool = tools{nt};
    % Remove the extremes based on the size, the location is not that bad...
    % We could obtain the outliers from all of them and then do and AND as well
    X0      = dt.(tool).x0;
    Y0      = dt.(tool).y0;
    Sizes   = dt.(tool).sMaj;
    Sizemin = dt.(tool).sMin;
    Thetas  = dt.(tool).Th;
    
    B       = prctile(Sizes, [twoTailedRange, 100 - twoTailedRange]);
    inRange = Sizes>=B(1) & Sizes<=B(2);
    % Apply
    X0      = X0(inRange);
    Y0      = Y0(inRange);
    Sizes   = Sizes(inRange);
    Sizemin = Sizemin(inRange);
    Thetas  = Thetas(inRange);
    if addsnr
        SNR     = dt.snr;
        SNR     = SNR(inRange);
        % Calculate mean and std of SNR
        meanSNR = mean(SNR);
        stdSNR  = std(SNR);
    end
    
    if onlyCenters
        % This will plot just the centers
        scatter(X0, Y0, 20, Cs(nt+1,:));
    else
        if adddice
            % Plot circunferences or ellipses
            Esynth.center = [unique(dt.synth.x0), unique(dt.synth.y0)]; 
            Esynth.sigma  = [unique(dt.synth.sMaj), unique(dt.synth.sMin)]; 
            Esynth.theta  = unique(dt.synth.Th);
            diceVals = [];
            for ne = 1:length(Sizes)   
                E2.center = [X0(ne),Y0(ne)]; 
                E2.sigma = [Sizes(ne),Sizemin(ne)]; 
                E2.theta = Thetas(ne);
                diceVals = [diceVals, ellipseDice(Esynth, E2, 'show', false)];
            end
            % Calculate mean and std of DICE
            meanDiceVal = mean(diceVals);
            stdDiceVal  = std(diceVals);
        end
        
        if useEllipse
            for ne = 1:length(Sizes)                
                h = drawellipse(X0(ne),Y0(ne),Thetas(ne),Sizes(ne),Sizemin(ne));
                set(h,'LineWidth',lineWidth,'LineStyle',lineStyle,'Color',Cs(nt+1,:));
                hold on
            end
        else
            centers = [X0,Y0];
            radii   = Sizes; % Viscircles needs radius and sigma-s are diameters
            viscircles(centers,radii,'LineWidth',lineWidth,'LineStyle',lineStyle,'Color',Cs(nt+1,:));
        end
    end
    if adddice && addsnr && addtext
        text(1.1*xlims(1),1.1*ylims(1),sprintf('DICE:%.2f(±%.2f) | SNR:%.2f(±%.2f)',...
                          meanDiceVal, stdDiceVal,meanSNR, stdSNR), ...
             'FontWeight','bold','FontSize',12)
    end
    hold on; grid on;axis equal
    % Edit the legend to understand how many values went in
    if strcmp(noiseLevel,"none")
        toolLegend = [toolLegend,{sprintf('%s',tool)}];
    else
        toolLegend = [toolLegend,{sprintf('%s (%i/%i)',tool,sum(inRange),length(inRange))}];
    end
end

% Plot Centers in top of the circles
a = [scatter(unique(dt.synth.x0), unique(dt.synth.y0),40,Cs(1,:),'filled')]; hold on
if onlyCenters
    for nt=1:length(tools)
        tool = tools{nt};
        a = [a; scatter(median(dt.(tool).x0), median(dt.(tool).y0),60,Cs(nt+1,:),'filled')]; hold on
    end
else
    if useEllipse
        if (length(unique(dt.synth.x0)) == 1 && length(unique(dt.synth.y0)) == 1)
            h = drawellipse(unique(dt.synth.x0),unique(dt.synth.y0),unique(dt.synth.Th),...
                        userfsize/2, userfsizemin/2);
            set(h,'LineWidth',2.5,'LineStyle','--','Color',Cs(1,:));
        else
            % xys = unique(dt.synth(:,{'x0','y0'}));
            % for nx=1:height(xys)
            %     viscircles([xys{nx,1}, xys{nx,2}],userfsize/2,...
            %     'LineWidth',2.5,'Color',Cs(1,:),'LineStyle','--');
            % end
        end
    else
        if (length(unique(dt.synth.x0)) == 1 && length(unique(dt.synth.y0)) == 1)
            viscircles([unique(dt.synth.x0),unique(dt.synth.y0)],userfsize,...
                'LineWidth',2.5,'Color',Cs(1,:),'LineStyle','--');
        else
            xys = unique(dt.synth(:,{'x0','y0'}));
            for nx=1:height(xys)
                viscircles([xys{nx,1}, xys{nx,2}],userfsize,...
                'LineWidth',2.5,'Color',Cs(1,:),'LineStyle','--');
            end
        end
    end
end
if centerDistr
    % Instead of the centers (see above), plot the ellipses showing the distribution of the centers
    for nt=1:length(tools)
        tool = tools{nt};
        
        % We need this for the legend, but only when the legend is shown
        a = [a; scatter(median(dt.(tool).x0), median(dt.(tool).y0),20,Cs(nt+1,:),'filled')]; hold on
        
        if length(X0) >1
            fitEllipse(X0,Y0,[0,0,0],centerPerc);
        end
    end
end

% If only plotting one circle, we want it to be on top, redraw it
if strcmp(noiseLevel, "none")
    if useEllipse
        h = drawellipse(unique(dt.(tool).x0),unique(dt.(tool).y0),unique(dt.(tool).Th),...
                        unique(dt.(tool).sMaj), unique(dt.(tool).sMin));
        set(h,'LineWidth',lineWidth,'LineStyle',lineStyle,'Color',Cs(nt+1,:));
    else
        viscircles(centers,radii,'LineWidth',lineWidth,'LineStyle',lineStyle,'Color',Cs(nt+1,:));
    end
end

% If we are plotting noise and addcibar is true, plot a bar showing the size
% confidence interval in the same plot
if addcibar && ~strcmp(noiseLevel, "none")
    xgt  = unique(dt.synth.x0);
    ygt  = unique(dt.synth.y0);
    smgt = unique(dt.synth.sMaj);
    x0   = median(dt.(tool).x0);
    y0   = median(dt.(tool).y0);
    sm   = median(dt.(tool).y0);
    % Select the min and max of the inRange (inside the selected CI) sizes
    smmin = min(Sizes);
    smmax = max(Sizes);
    % METHOD 1
    %{
        linestarts = x0 + smmin/2;
        lineends   = x0 + smmax/2;
        plot([linestarts, lineends],[y0,y0],'k-','LineWidth',4)
    %}
    
    % METHOD 2: 3 lines in X
    %{
        % The ground truth size
        linestarts = xgt - smgt/2;
        lineends   = xgt + smgt/2;
        plot([linestarts, lineends],[1,1],'-','color',Cs(1,:),'LineWidth',4);hold on
        % The longest CI
        linestarts = x0 - smmax/2;
        lineends   = x0 + smmax/2;
        plot([linestarts, lineends],[1.15,1.15],'-','color',color*0.75,'LineWidth',4);hold on
        % The shortest CI
        linestarts = x0 - smmin/2;
        lineends   = x0 + smmin/2;
        plot([linestarts, lineends],[1.3,1.3],'-','color',color*0.1,'LineWidth',4);hold on
     %}
% METHOD 3: 3 lines in Y
    % {
        % The ground truth size
        linestarts = ygt - smgt/2;
        lineends   = ygt + smgt/2;
        plot([xlims(1)+0.1,xlims(1)+0.1],[linestarts, lineends],'-','color',Cs(1,:),'LineWidth',4);hold on
        % The longest CI
        linestarts = y0 - smmax/2;
        lineends   = y0 + smmax/2;
        plot([xlims(1)+0.25,xlims(1)+0.25],[linestarts, lineends],'-','color',color*0.75,'LineWidth',4);hold on
        % The shortest CI
        linestarts = y0 - smmin/2;
        lineends   = y0 + smmin/2;
        plot([xlims(1)+0.4,xlims(1)+0.4],[linestarts, lineends],'-','color',color*0.1,'LineWidth',4);hold on

    %}
end


if addcihist && ~strcmp(noiseLevel, "none")
    % Create new axes
    bigax = xlims(2)-xlims(1);
    smax  = 0.3 * bigax;
    startx = xlims(1)+0.12*smax;
    starty = ylims(1)+0.2*smax;
    endx   = xlims(1)+smax;
    endy   = ylims(1)+smax;
    % Find equivalencies
    realxbeg = 0; % degrees
    realxend = 6; % degrees
    gt       = unique(dt.synth.sMaj);
    rgt      = rescale([realxbeg, gt, realxend],startx,endx);
    % Where to write text
    textx  = mean([startx,endx]);
    texty  = mean([starty,ylims(1)]);
    hax = plot([startx, endx],starty*[1,1],'Color','k','LineStyle','-','LineWidth',1);hold on;
    vax = plot(startx*[1,1],[starty, endy],'Color','k','LineStyle','-','LineWidth',1);
    % Calculate the size density function
    [pdfsizes,xsizes] = ksdensity(Sizes,'Support','positive','BoundaryCorrection','reflection');
    % Rescale the values to be inside the new axes
    tmp = rescale([realxbeg, xsizes(end), realxend],startx,endx);
    rxsizes   = rescale(xsizes,startx,tmp(2));
    rpdfsizes = rescale(pdfsizes,starty,endy);
    % Plot it
    plot(rxsizes,rpdfsizes,'Color','k','LineStyle','-','LineWidth',1.5); hold on;
    % Add line where the ground truth is
    hmin = plot(rgt(2)*[1,1],[starty endy],'Color','b','LineStyle','--','LineWidth',1.5);hold on;
    % Add a text with the ground truth value
    %  (use the percentile min and max)
    text(startx, texty, sprintf('%i',realxbeg), 'FontSize',10);
    text(endx, texty, sprintf('%i deg',realxend), 'FontSize',10);
    text(rgt(2), texty, sprintf('%i',gt), 'FontSize',10);
end

if addtext
    % legend(a,toolLegend,'location','northeast')
    xlabel('degrees')
    ylabel('degrees')
else
    % set(gca, 'xtick', [],'ytick', [])
end

xticks(xtick)
yticks(ytick)
xlim(xlims); 
ylim(ylims);  

set(gca, 'FontSize', fontsize)
if addtext, title({['Tool: ' tools{:} ', HRF: ' strrep(useHRF,'_','\_')],...
                   ['Noise: ' char(noiseLevel) ', Ratio: ' num2str(userfsize/userfsizemin) ]}),end
fnameRoot = [tools{:} '_' useHRF '_noise_' char(noiseLevel) ];

%% Save it
set(0, 'DefaultFigureRenderer', 'painters');
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

% title({sprintf('HRF: %s. Filled circle is the median value', useHRF), ...
%        sprintf('Only the center %i%% percentile of values used', centerPerc)})


% Check what we are plotting
% a = compTable;
% a = a(a.synth.sMaj==2 & a.noiseLevel==0.2 & a.HRFtype=="friston",{'synth','aprfcss','vista'})


end






