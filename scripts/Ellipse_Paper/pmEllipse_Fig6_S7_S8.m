function pmEllipse_Fig6_S7_S8
% Make Figures 6A-B-C, S7 and S8
% 
% TODO: separate it in sub-scripts that use the same dataset that we load here
% 
% See also
%  s00_MainFiguresScript
% 
%% Plotting parameters
clear all
ext  = 'svg'; % Can be 'png', 'svg'
saveTo = fullfile(pmRootPath,'local','figures');  % Folder path
if ~exist(saveTo,'dir'), mkdir(saveTo); end


%% READ: Real Data 7T

proj   = 'realdata';
tools  = {'vista6','vista4'};
subs   = {'115017','164131','536647'}; 
ses    = '01';
run    = '01';

[compTable,bylabelsums] = pmEllipse_loadExpData(proj,tools,subs,ses,run);
nonfilteredbylabelsums = bylabelsums;

%% READ THE SYNTH DATA
% Read the synthetic data as well, this is the eccenv2 dataset, with mid and low
% noise levels, with TR=1 and 2, duration 400, and the ground truth aspect ratio
% limited to 1

% Generated TR=1, Dur=300 data to plot alongside with the real data
fprintf('\n\nLoading synthetic TR=1 300sec data')

sub = 'ellipse'; ses = 'tr1dur300v3';
p = fullfile(pmRootPath,'local',sub,'BIDS','derivatives','prfreport',['sub-' sub],['ses-' ses]);
f = ['sub-' sub '_ses-' ses '-prf_acq-normal_run-01_bold.mat'];
% tools = {'synth','vista4','vista6'};
tools = {'synth','vista6'};
C = load(fullfile(p,f));
dt = C.compTable;
for nt=1:length(tools)
    dt.(tools{nt}).aspect = dt.(tools{nt}).sMaj ./ dt.(tools{nt}).sMin;
    [TH,R] = cart2pol(dt.(tools{nt}).x0, dt.(tools{nt}).y0);
    dt.(tools{nt}).angle = rad2deg(TH);
    dt.(tools{nt}).eccen = R;
    dt.(tools{nt}).area  = pmEllipseArea(dt.(tools{nt}).sMaj,dt.(tools{nt}).sMin);
end
% GT aspect ratio is always one
dt   = dt(dt.synth.aspect==1,:);
A1A2 = dt;
disp ('... done with load')

%% DO THE FILTERING
% Obtain the same eccentricities as in the simulations
eccenvalues = linspace(2,7,6);


% v2 only goes to size 1.5, v3 goes all the way to 4, based on the calculations
% of VISTA4 over the HCP 7T data, see pmEllipse_CalculateVista4Vals
% The results from that script these
%    Eccen at 25 = 2.5 and 75 = 6.5 percentiles. 
%    Radius at 25 = 1.4 and 75 = 3.1 percentiles. 
%    Area min = 6.5, max = 30



% Some plot options
doSave     = true;
centerPerc = 95;
eccenInGT  = true;
xlims      = [0,10];
ylims      = [0,10];
tools      = {'vista6'}; 
useLabels  = {    'V1d', 'V2d', 'V3d','V1v', 'V2v', 'V3v'};
Cs         = .65*[1 0 0; 0 1 0; 0 0 1;1 0 0; 0 1 0; 0 0 1];
marks      =     [  '*',   '*',   '*',  'o',   'o',   'o',];
lstyle     =     { '-.',  '-.',  '-.',  '-',   '-',   '-'};

% useLabels  = {'V1','V2','V3'};
duration   = 300;
tr         = 1;
% Filter results (ALL FILTERING SHOULD HAPPEN HERE)
sMajMIN    = 1.4; % Not used, uses area
sMinMIN    = 1.4; % Not used, uses area
sMajMAX    = 3.1; % Not used, uses area
eccenMIN   = 2.5; % Comes from pmEllipse_CalculateVista4Vals
eccenMAX   = 6.5; % Comes from pmEllipse_CalculateVista4Vals
areaMIN    = 6.5; % Comes from pmEllipse_CalculateVista4Vals
areaMAX    = 30;  % Comes from pmEllipse_CalculateVista4Vals
minR2      = 0.25;
% Percentiles for conf intervals
lowerprct  = (100-centerPerc)/2;
upperprct  = 100 - lowerprct;
% How many bins
NeccenBins = 7;
NareaBins  = NeccenBins;
% Close all



% Filter the synthetic data
tools  = {'vista6'};
for nt=1:length(tools)
    tool = tools{nt};
    A1A2 = A1A2(...
            A1A2.(tool).area  >= areaMIN  & ...
            A1A2.(tool).area  <= areaMAX & ...
            A1A2.(tool).eccen >= eccenMIN & ...
            A1A2.(tool).eccen <= eccenMAX & ...
            A1A2.HRFtype=="vista_twogammas", :);    
end


% Filter the HCP 7T data
tools  = {'vista4','vista6'};
subs   = {'115017','164131','536647'};
bylabelsums = nonfilteredbylabelsums;
% Apply the restrictions
for nt=1:length(tools)
    tool = tools{nt};
    for nl = 1:length(useLabels)
        lab = useLabels{nl};
        [TH,R]      = cart2pol(bylabelsums.(tool).(lab).x0, bylabelsums.(tool).(lab).y0);
        bylabelsums.(tool).(lab).angle = rad2deg(TH);
        bylabelsums.(tool).(lab).eccen = R;
        bylabelsums.(tool).(lab).area  = pmEllipseArea(bylabelsums.(tool).(lab).sMaj, ...
                                                       bylabelsums.(tool).(lab).sMin);
        bylabelsums.(tool).(lab) = bylabelsums.(tool).(lab)(...
                                   bylabelsums.(tool).(lab).area  >= areaMIN & ...
                                   bylabelsums.(tool).(lab).area  <= areaMAX & ...
                                   bylabelsums.(tool).(lab).eccen >= eccenMIN & ...
                                   bylabelsums.(tool).(lab).eccen <= eccenMAX & ...
                                   bylabelsums.(tool).(lab).r2    >= minR2,:);
        
        fprintf('\n%s %.2g',lab,min(bylabelsums.(tool).(lab).sMin))
        % Theta can only be [-90,+90]
        % Vista and Afni treat it differently it seems
        % I added 90 deg to AFNI, but I still don't know if I need it or not. Remove it
        theta            = rad2deg(bylabelsums.(tool).(lab).Th - deg2rad(90));
        theta(theta>180) = theta(theta>180) -180;
        theta(theta>90)  = theta(theta>90) -180;
        bylabelsums.(tool).(lab).Th = theta;
        % We can express the theta in the same way, because we only care about the
        % radiality, not the exact angle
        angle            = bylabelsums.(tool).(lab).angle;
        angle(angle>180) = angle(angle>180) - 180;
        angle(angle>90)  = angle(angle>90) - 180;
        bylabelsums.(tool).(lab).angle = angle;
        % Calculate aspect ratio and area
        bylabelsums.(tool).(lab).aspect = bylabelsums.(tool).(lab).sMaj  ./ bylabelsums.(tool).(lab).sMin;
        
        
    end
end


        


%% PLOT 6 (A,B,C):Real Data 7T, ONLY VISTA

tool = 'vista6';
% Discretize by label, to bin the eccentricities
% SYNTHETIC DATA
A1A2.(tool).Y = zeros(size(A1A2.(tool).aspect));
A1A2.(tool).Y(A1A2.noiseLevel=="mid") = discretize(A1A2.(tool).eccen(A1A2.noiseLevel=="mid"),eccenvalues); 
A1A2.(tool).Y(A1A2.noiseLevel=="low") = discretize(A1A2.(tool).eccen(A1A2.noiseLevel=="low"),eccenvalues); 
A1A2low = A1A2(A1A2.noiseLevel=="low", :);
A1A2mid = A1A2(A1A2.noiseLevel=="mid", :);

aspectlow = A1A2low.(tool).aspect;
aspectmid = A1A2mid.(tool).aspect;

% Apply percentiles and plot individually
% Create the vectors and then plot all together
vistaMedLowEcc = zeros(1,length(eccenvalues)-1);
vistaMedMidEcc = zeros(1,length(eccenvalues)-1);
vista25LowEcc = zeros(1,length(eccenvalues)-1);
vista25MidEcc = zeros(1,length(eccenvalues)-1);
vista975LowEcc = zeros(1,length(eccenvalues)-1);
vista975MidEcc = zeros(1,length(eccenvalues)-1);
for ne=1:(length(eccenvalues)-1)
    aspecclow = A1A2low.(tool).aspect(A1A2low.(tool).Y==ne);
    aspeccmid = A1A2mid.(tool).aspect(A1A2mid.(tool).Y==ne);
    % Median
    vistaMedLowEcc(ne) = median(aspecclow);
    vistaMedMidEcc(ne) = median(aspeccmid);
    
    vista25LowEcc(ne) = prctile(aspecclow,lowerprct);
    vista25MidEcc(ne) = prctile(aspeccmid,lowerprct);
    vista975LowEcc(ne) = prctile(aspecclow,upperprct);
    vista975MidEcc(ne) = prctile(aspeccmid,upperprct);
end

% HCP 7T DATA
for nl  = 1:length(useLabels)
    lab = useLabels{nl};
    % Discretize by label, to bin the eccentricities
    bylabelsums.(tool).(lab).Y = discretize(bylabelsums.(tool).(lab).eccen,eccenvalues); 
    % Apply percentiles and plot individually
    % Create the vectors and then plot all together
    aspectmedecc.(lab) = zeros(1,length(eccenvalues)-1);
    aspectN.(lab) = zeros(1,length(eccenvalues)-1);
    aspectminecc.(lab) = zeros(1,length(eccenvalues)-1);
    aspectmaxecc.(lab) = zeros(1,length(eccenvalues)-1);
    
    for ne=1:(length(eccenvalues)-1)
        % ECC - ASPECT
        aspecc = bylabelsums.(tool).(lab).aspect(bylabelsums.(tool).(lab).Y==ne);
        % Median and std
        if isempty(aspecc)
            aspectmedecc.(lab)(ne) = 0;
            aspectN.(lab)(ne) = 0;
            aspectminecc.(lab)(ne) = 0;
            aspectmaxecc.(lab)(ne) = 0;
        else
            aspectmedecc.(lab)(ne) = median(aspecc);
            aspectN.(lab)(ne) = length(aspecc);
            aspectminecc.(lab)(ne) = prctile(aspecc,lowerprct);  % They use SEM, check
            aspectmaxecc.(lab)(ne) = prctile(aspecc,upperprct);
        end
    end
end


%{
% STAT ANOVA COMPARE EXP RESULTS
    kktab  = table();
    tool = 'vista6';
    
    for nt=1:length(useLabels)
        lab = useLabels{nt}
        kk      = [bylabelsums.(tool).(lab)(:,'aspect'), bylabelsums.(tool).(lab)(:,'Y')];
        kk.lab  = categorical(repmat(string(lab),[height(bylabelsums.(tool).(lab)(:,'aspect')),1]));
        kktab   = [kktab;kk];
    end
    

    anovan(kktab{:,'aspect'}, {kktab{:,'lab'}, string(kktab{:,'Y'}) })
    anova1(kktab{:,'aspect'}, string(kktab{:,'Y'}) )
    anova1(kktab{:,'aspect'}, kktab{:,'lab'})
%}










% PLOT 6A
fnameBegin = 'Fig6-A_RealData_Ecc&Size';
% Create main plot with the ground truth lines
fnameEnd = sprintf('TR-%i_Dur-%is_C.I.-%i',tr,duration,centerPerc);
fnameRoot = strcat(fnameBegin,'-', fnameEnd);
% disp(fnameRoot) 
kk = mrvNewGraphWin(fnameRoot);
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  .5 0.4]);
% ECCEN vs ASPECT
% Plot it 
E = eccenvalues;
Emidpoints = mean([E(2:end);E(1:end-1)]);
as = [];
for nl  = 1:length(useLabels)
    lab = useLabels{nl};
    as = [as;plot(Emidpoints,aspectmedecc.(lab),'Color',Cs(nl,:), ...
              'LineStyle',lstyle{nl},'LineWidth',3)];hold on
%     plot([Emidpoints+(nl-3)/40;Emidpoints+(nl-3)/40] ,...
%                [aspectminecc.(lab)  ; aspectmaxecc.(lab)], ...
%                'Color',Cs(nl,:),'LineStyle',lstyle{nl},'LineWidth',2);  % 0.75*[0 1 0]
end
% Plot noise values as bands
% lowplot = plot(Emidpoints, vistaMedLowEcc,'k--','LineWidth',3);
% midplot = plot(Emidpoints, vistaMedMidEcc,'k:','LineWidth',3);
% xpoints,upper,lower,color,edge,add,transparency
% jbfill(Emidpoints,vista75LowEcc,vista25LowEcc,[.25,.25,.25],[.5,.5,.5],1,0.5);hold on
sh = jbfill(Emidpoints,vista975MidEcc,vista25MidEcc,[.45,.45,.45],[.3,.3,.3],1,0.25);

legend([as;sh],[useLabels,'Synth Mid Noise'], 'Location','eastoutside')
title(strrep(sprintf('%s_TR-%i_Dur-%is_C.I.-%i',...
    tool,tr,duration,centerPerc),'_','\_'))
grid on
xlabel('Eccentricity (deg)')
ylabel('pRF aspect ratio')
xlim([Emidpoints(1)-.2,Emidpoints(end)+.2]);
ylim([1,5]);
xticks(Emidpoints);
set(gca, 'FontSize', 16);

fname = fullfile(saveTo, strcat(fnameRoot,['.' ext]));
saveas(gcf,fname,ext);
fprintf('\nSaved %s\n', fname)

% PLOT S7
aspects = [];
fnameBegin = 'FigS7_RealData_AspectHistogram_Separated';
% Create main plot with the ground truth lines
fnameEnd = sprintf('TR-%i_Dur-%is_C.I.-%i',tr,duration,centerPerc);
fnameRoot = strcat(fnameBegin,'-', fnameEnd);
% disp(fnameRoot) 
kk = mrvNewGraphWin(fnameRoot);
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  .5  .5]);
% SET UP DATA
% Here is the aspect we want to plot
% bylabelsums.(tool).(lab).aspect
tool = 'vista6';
for nl  = 1:length(useLabels)
    subplot(2,3,nl)
    lab = useLabels{nl};
    % Obtain aspect
    aspectvista = bylabelsums.(tool).(lab).aspect;
    % aspectvista = aspectvista(aspectvista < 5);
    aspects     = [aspects;aspectvista];
    
    % Plot it
    h = histogram(aspectvista,25,'Normalization','probability');
    set(h,'LineWidth',2,'EdgeColor',[.5 .5 .5],'EdgeAlpha',0,'FaceAlpha',1,'FaceColor',[.5 .5 .5]);hold on
    plot(median(aspectvista)*[1,1],[0,max(h.Values)],'r-')    
    
    % Add the low noise and mid noise lines now

    
    xlim([1,5])
    ylim([0,0.35])
    % legend({lab,'median','Synth Low Noise','Synth Mid Noise'});
    title(lab)
    xlabel('Aspect Ratio')
end
fname = fullfile(saveTo, strcat(fnameRoot,['.' ext]));
saveas(gcf,fname,ext);
fprintf('\nSaved %s\n', fname)

% PLOT 6B
fnameBegin = 'Fig6-B_RealData_AspectHistogram_Combined';
% Create main plot with the ground truth lines
fnameEnd = sprintf('TR-%i_Dur-%is_C.I.-%i',tr,duration,centerPerc);
fnameRoot = strcat(fnameBegin,'-', fnameEnd);
% disp(fnameRoot) 
kk = mrvNewGraphWin(fnameRoot);
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  .5  .5]);

binWidth = 0.05;

% hmid = histogram(aspectmid,'DisplayStyle','stairs','BinWidth',h.BinWidth,'Normalization','probability');
hmid = histogram(aspectmid,'BinWidth',binWidth,'Normalization','probability');
% set(hmid,'LineWidth',2,'EdgeColor',[.5 .5 .5 ],'LineStyle','-','EdgeAlpha',.5,'FaceAlpha',.5,'FaceColor',[.5 .5 .5 ]);
set(hmid,'LineWidth',2,'EdgeColor','k','FaceAlpha',1,'FaceColor','k');hold on
blow = plot(median(aspectmid)*[1,1],[0,.1],'LineWidth',2,'Color','k','LineStyle','--'); 


h = histogram(aspects,35,'Normalization','probability','BinWidth',binWidth);hold on
% set(h,'LineWidth',2,'EdgeColor','k','FaceAlpha',1,'FaceColor','k');
set(h,'LineWidth',2,'EdgeColor',[.5 .5 .5 ],'LineStyle','-','EdgeAlpha',0,'FaceAlpha',.75,'FaceColor',[.5 .5 .5 ]);
a = plot(median(aspects)*[1,1],[0,.1],'Color',[.5 .5 .5 ],'LineStyle','--');


% Add the low noise and mid noise lines now
% hlow = histogram(aspectlow,'DisplayStyle','stairs','BinWidth',h.BinWidth,'Normalization','probability');
% hlow = histogram(aspectlow,'BinWidth',h.BinWidth,'Normalization','probability');
% set(hlow,'LineWidth',2,'EdgeColor',[1 .5 .5 ],'LineStyle','-','EdgeAlpha',0,'FaceAlpha',.5,'FaceColor',[1 .5 .5 ]);
% alow = plot(median(aspectlow)*[1,1],[0,.1],'LineWidth',2,'Color',[1 .5 .5 ],'LineStyle','-'); 

xlim([1,5]);

% legend([h;a;hlow;hmid],{'Experimental Data','Median of Exp. Data','Synth Low Noise','Synth Mid Noise'});

legend([h;a;hmid;blow],{'Experimental Data','Median of Exp. Data','Synth Mid Noise','Median of Synth. Data'});
fname = fullfile(saveTo, strcat(fnameRoot,['.' ext]));
saveas(gcf,fname,ext);
fprintf('\nSaved %s\n', fname)

%% PLOT S8  % THETA vs ANGLE
% prepare data
tool = 'vista6';
thetas = [];
angles = [];
for nl  = 1:length(useLabels)
    lab = useLabels{nl};
    thetas = [thetas;bylabelsums.(tool).(lab).Th];
    angles = [angles;bylabelsums.(tool).(lab).angle];
end

fnameBegin = 'FigS8_RealData_AnglevsTheta';
% Create main plot with the ground truth lines
fnameEnd = sprintf('TR-%i_Dur-%is_C.I.-%i',tr,duration,centerPerc);
fnameRoot = strcat(fnameBegin,'-', fnameEnd);
% disp(fnameRoot) 
kk = mrvNewGraphWin(fnameRoot);
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  1  0.5]);
subplot(1,2,1)
% PLOT 2b
h = histogram(theta-angle,25,'Normalization','probability');
set(h,'LineWidth',2,'EdgeColor',[.5 .5 .5],'EdgeAlpha',0,'FaceAlpha',1,'FaceColor',[.5 .5 .5]);hold on
xlabel('Theta - Angle')

subplot(1,2,2)
plot(thetas, angles,'ko');xlabel('\Theta (deg)');ylabel('Angle (deg)');hold on
identityLine(gca);
xlim([-90,90]);
ylim([-90,90]);
xticks(-90:15:90);
yticks(-90:15:90);
fname = fullfile(saveTo, strcat(fnameRoot,['.' ext]));
saveas(gcf,fname,ext);
fprintf('\nSaved %s\n', fname)

%% PLOT 6C: Compare r2 values vista4/vista6
useLabels = {'V1','V2','V3'};   
fnameRoot = 'Fig6-C_RealData_R2diff_histogram_filteredAsTheRest';

% Create intermediate variables
% R2 in perc
v6 = 100*compTable.vista6.r2;
v4 = 100*compTable.vista4.r2;
% Eccentricity values for filtering
[~,R6] = cart2pol(compTable.vista6.x0, compTable.vista6.y0);
[~,R4] = cart2pol(compTable.vista4.x0, compTable.vista4.y0);
compTable.vista6.eccen = R6;
compTable.vista4.eccen = R4;


% Filter by variance explained
v6ind = (v6 > 100*minR2) & ...
        (compTable.vista6.sMaj  > sMajMIN)  & (compTable.vista6.sMaj < sMajMAX) & ...
        (compTable.vista6.eccen > eccenMIN) & (compTable.vista6.eccen < eccenMAX);
v4ind = (v4 > 100*minR2) & ...
        (compTable.vista4.sMaj  > sMajMIN)  & (compTable.vista4.sMaj < sMajMAX) & ...
        (compTable.vista4.eccen > eccenMIN) & (compTable.vista4.eccen < eccenMAX);
vind  = v6ind & v4ind;
% Create filtered version
v6f   = v6(vind);
v4f   = v4(vind);

v6m   = median(v6);
v4m   = median(v4);
v6fm  = median(v6f);
v4fm  = median(v4f);

v64   = 100 * (v6m - v4m)/v4m;
v64f  = 100 * (v6fm - v4fm)/v4fm;


kk = mrvNewGraphWin('R2 vista4/vista6');
set(kk,'Position',[0.007 0.62  0.5 0.5]);
h = histogram(v6f - v4f,'Normalization','probability');  % 'DisplayStyle','stairs'
set(h,'LineWidth',2,'EdgeColor',[.5 .5 .5],'EdgeAlpha',0,'FaceAlpha',1,'FaceColor',[.5 .5 .5]);hold on
a = plot(median(v6f - v4f)*[1,1],[0,max(h.Values)],'r-');
xlabel('Delta R2 (Elliptical - Circular; in %)')
legend(a,'Median of the difference')
xlim([-2,5])
set(gca,'FontSize',20)

% Print the variance explained for Insub
%{
kk = mrvNewGraphWin('R2 vista6 and vista4');
set(kk,'Position',[0.007 0.62  0.5 0.5]);
h6 = histogram(v6f,'Normalization','probability');  % 'DisplayStyle','stairs'
set(h6,'LineWidth',2,'EdgeColor',[.5 .5 .5],'EdgeAlpha',0,'FaceAlpha',1,'FaceColor','k');hold on
h4 = histogram(v4f,'Normalization','probability');  % 'DisplayStyle','stairs'
set(h4,'LineWidth',2,'EdgeColor',[.5 .5 .5],'EdgeAlpha',0,'FaceAlpha',.65,'FaceColor',[.5 .5 .5]);
xlabel('R2 (%)')
legend('Elliptical Fit','Circular fit')
xlim([20,100])
set(gca,'FontSize',20)
%}


% title(sprintf('Elliptical median variance explained is %1.2g%% larger than Circular',v64f))
fname = fullfile(saveTo, strcat(fnameRoot,['.' ext]));
saveas(gcf,fname,ext);
fprintf('\nSaved %s\n', fname)

end

