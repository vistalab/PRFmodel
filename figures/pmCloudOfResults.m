function pmCloudOfResults(compTable, tools,varargin)
%pmCloudOfResults Return a distribution plot for several noise values
%   It return only one plot, it need to be composited with another script

%% Read the inputs
% Make varargin lower case, remove white spaces...
varargin = mrvParamFormat(varargin);
% Parse
p = inputParser;
p.addRequired('compTable');
p.addRequired('tools'                  , @iscell);
p.addParameter('onlycenters', false    , @islogical);
p.addParameter('userfsize'  , 2        , @isnumeric);
p.addParameter('centerperc' , 50       , @isnumeric);
p.addParameter('usehrf'     , 'friston', @ischar);
p.addParameter('linestyle'  , '-'      , @ischar);
p.addParameter('linewidth'  , .7       , @isnumeric);
p.addParameter('newwin'     , true     , @islogical);
p.addParameter('noiselevel' , 'none'   , @ischar);
p.addParameter('saveto'     , ''       , @ischar);
p.addParameter('savetotype' , 'png'    , @ischar);
% Parse. Assign result inside each case
p.parse(compTable, tools, varargin{:});
% Read here only the generic ones
onlyCenters = p.Results.onlycenters;
userfsize   = p.Results.userfsize;
centerPerc  = p.Results.centerperc;
useHRF      = p.Results.usehrf;
lineStyle   = p.Results.linestyle;
lineWidth   = p.Results.linewidth;
newWin      = p.Results.newwin;
noiseLevel  = p.Results.noiselevel; noiseLevel  = string(noiseLevel);  
saveTo      = p.Results.saveto;
saveToType  = p.Results.savetotype; 

%% Do the thing

Cs               = distinguishable_colors(1+length(tools),'w');
Cs(2,:) = Cs(2,:) * 0.65;
if newWin
    mrvNewGraphWin(sprintf('Clouds, %s',useHRF));
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
dt        = compTable(compTable.synth.sMaj == userfsize , :);
% Only one noise level
dt        = dt(dt.noiseLevel == noiseLevel, :);
% Use only friston
dt        = dt(dt.HRFtype == string(useHRF), :);
toolLegend= {'synth'};
for nt=1:length(tools)
    tool = tools{nt};
    % Remove the extremes based on the size, the location is not that bad...
    % We could obtain the outliers from all of them and then do and AND as well
    X0      = dt.(tool).x0;
    Y0      = dt.(tool).y0;
    Sizes   = dt.(tool).sMaj;
    B       = prctile(Sizes,[twoTailedRange, 100-twoTailedRange]);
    inRange = Sizes>=B(1) & Sizes<=B(2);
    % Apply
    X0      = X0(inRange);
    Y0      = Y0(inRange);
    Sizes   = Sizes(inRange);

    if onlyCenters
        % This will plot just the centers
        scatter(X0, Y0, 20, Cs(nt+1,:));
    else
        % Plot circunferences
        centers = [X0,Y0];
        radii   = Sizes/2; % Viscircles needs radius and sigma-s are diameters
        viscircles(centers,radii,'LineWidth',lineWidth,'LineStyle',lineStyle,'Color',Cs(nt+1,:));
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
a = [scatter(unique(dt.synth.x0), unique(dt.synth.y0),200,Cs(1,:),'filled')]; hold on
if onlyCenters
    for nt=1:length(tools)
        tool = tools{nt};
        a = [a; scatter(median(dt.(tool).x0), median(dt.(tool).y0),150,Cs(nt+1,:),'filled')]; hold on
    end
else
    viscircles([unique(dt.synth.x0),unique(dt.synth.y0)],userfsize/2,...
        'LineWidth',4,'Color',Cs(1,:),'LineStyle','--');
    % Instead of the centers (see above), plot the ellipses showing the distribution of the
    % centers
    for nt=1:length(tools)
        tool = tools{nt};
        % a = [a; scatter(median(dt.(tool).x0), median(dt.(tool).y0),80,Cs(nt+1,:),'filled')]; hold on
        a = [a; scatter(median(dt.(tool).x0), median(dt.(tool).y0),20,Cs(nt+1,:),'filled')]; hold on
        fitEllipse(X0,Y0,[0,0,0]);
    end
end


% If only plotting one circle, we want it to be on top, redraw it
if strcmp(noiseLevel, "none")
    viscircles(centers,radii,'LineWidth',lineWidth,'LineStyle',lineStyle,'Color',Cs(nt+1,:));
end


legend(a,toolLegend,'location','northeast')
xlabel('degrees'); xlim([1.5,4.5]); xticks([2:.5:4])
ylabel('degrees'); ylim([1.5,4.5]); yticks([2:.5:4])

set(gca, 'FontSize', 16)
title({['Tool: ' tools{:} ', HRF: ' strrep(useHRF,'_','\_')],['Noise: ' char(noiseLevel) ]})
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






