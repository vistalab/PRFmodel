function pmEllipse_FigAfniHCP
% Make Figures AFNI HCP for reviewers
% 
% TODO: separate it in sub-scripts that use the same dataset that we load here
% 
% See also
%  s00_MainFiguresScript

%% Plotting parameters
ext  = 'png'; % Could be svg
saveTo = fullfile(pmRootPath,'local','figuresTEST');  % Folder path
if ~exist(saveTo,'dir'), mkdir(saveTo); end


%% Create the r2 values for afni
%{
proj   = 'realdata';
tool   = 'afni6';
subs   = {'115017','164131','536647'}; 
ses    = '01';
run    = '01';

for ns=1:length(subs)
    sub   = subs{ns};
    path2 = fullfile(pmRootPath,'local',proj,'BIDS','derivatives',...
                     ['prfanalyze-' tool],['sub-' sub],['ses-' ses]);
    fname     = fullfile(path2,['sub-' sub '_ses-' ses '_task-prf_acq-normal_run-' run '_modelpred.nii.gz']);
    tmp       = niftiRead(fname);
    modelpred = squeeze(tmp.data);

    fname     = fullfile(path2,['sub-' sub '_ses-' ses '_task-prf_acq-normal_run-' run '_testdata.nii.gz']);
    tmp       = niftiRead(fname);
    testdata  = squeeze(tmp.data);
    
    r2        = calccod(modelpred, testdata, 2, 1, 1) / 100;
    
    % Write it as a file
    fname     = fullfile(path2,['sub-' sub '_ses-' ses '_task-prf_acq-normal_run-' run '_a.nii.gz']);
    tmp       = niftiRead(fname);
    tmp.data  = r2;
    tmp.fname = fullfile(path2,['sub-' sub '_ses-' ses '_task-prf_acq-normal_run-' run '_r2.nii.gz']);
    niftiWrite(tmp);
end
%}


%% READ: Real Data 7T

proj   = 'realdata';
tools  = {'afni6','vista6'};
subs   = {'115017','164131','536647'}; 
ses    = '01';
run    = '01';

[compTable,bylabelsums] = pmEllipse_loadExpData(proj,tools,subs,ses,run);
nonfilteredbuylabelsums = bylabelsums;

%% Plot some tests
fnameRoot = 'AFNI6 vs mrVISTA6 R2 plots';
kk = mrvNewGraphWin(fnameRoot);
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  .5 0.8]);

% EXPERIMENTAL HCP 7T
af6v1v2v3    = [bylabelsums.afni6.V1;  bylabelsums.afni6.V2;  bylabelsums.afni6.V3];
vi6v1v2v3    = [bylabelsums.vista6.V1; bylabelsums.vista6.V2; bylabelsums.vista6.V3];
% Calculate eccentricity
[~,R]       = cart2pol(af6v1v2v3.x0, af6v1v2v3.y0);
af6v1v2v3.eccen = R;
[~,R]       = cart2pol(vi6v1v2v3.x0, vi6v1v2v3.y0);
vi6v1v2v3.eccen = R;


subplot(2,2,1)
% Histograms of R2 values (differently calculated...)
% NOTE: 
% VISTA: Obtain the variance explained: r2 = 1-results.model{1}.rss./results.model{1}.rawrss;
% AFNI:      r2        = calccod(modelpred, testdata, 2, 1, 1) / 100;

% FILTER 1
% {
    afind = af6v1v2v3.r2 > 0.25;
    viind = vi6v1v2v3.r2 > 0.25;
%}

% FILTER 2
%{
    % Filter results
    sMajMIN    = 1 ; 
    sMinMIN    = 1 ; 
    sMajMAX    = 3 ; 
    eccenMIN   = 2 ; 
    eccenMAX   = 6 ; 
    minR2      = .1;
    afind = af6v1v2v3.sMaj  > sMajMIN  & af6v1v2v3.sMin  > sMinMIN & ...
            af6v1v2v3.sMaj  < sMajMAX  & af6v1v2v3.eccen > eccenMIN & ...
            af6v1v2v3.eccen < eccenMAX & af6v1v2v3.r2    > minR2;
    viind = vi6v1v2v3.sMaj  > sMajMIN  & vi6v1v2v3.sMin  > sMinMIN & ...
            vi6v1v2v3.sMaj  < sMajMAX  & vi6v1v2v3.eccen > eccenMIN & ...
            vi6v1v2v3.eccen < eccenMAX & vi6v1v2v3.r2    > minR2;
%}



h = histogram(100*af6v1v2v3.r2(afind),'FaceColor','b','FaceAlpha',0.5,'NumBins',50);hold on;
plot(median(100*af6v1v2v3.r2(afind))*[1,1],[0,.9*max(get(h,'Values'))],'Color','b','LineStyle','--');

histogram(100*vi6v1v2v3.r2(viind),'FaceColor',[.5 .5 .5],'NumBins',50);
plot(median(100*vi6v1v2v3.r2(viind))*[1,1],[0,.9*max(get(h,'Values'))],'Color',[.5 .5 .5 ],'LineStyle','--');
legend({'HCP 7T Afni6','median','HCP 7T Vista6','median'})
title('HCP R2 for Afni6 and vista6 (V1-V2-V3, R2>25%)')

subplot(2,2,3)

afni6aspect  = af6v1v2v3.sMaj(afind) ./ af6v1v2v3.sMin(afind);
afni6aspect  = afni6aspect(afni6aspect < 5);
vista6aspect = vi6v1v2v3.sMaj(viind) ./ vi6v1v2v3.sMin(viind);
vista6aspect = vista6aspect(vista6aspect < 5);
h = histogram(afni6aspect,'FaceColor','b','FaceAlpha',0.5,'NumBins',50);hold on;
xlim([1,5])
plot(median(afni6aspect)*[1,1],[0,.9*max(get(h,'Values'))],'Color','b','LineStyle','--');
h = histogram(vista6aspect,'FaceColor',[.5 .5 .5],'NumBins',50);
plot(median(vista6aspect)*[1,1],[0,.9*max(get(h,'Values'))],'Color',[.5 .5 .5 ],'LineStyle','--');
legend({'HCP 7T Afni6','median','HCP 7T Vista6','median'})
title('HCP Aspect Ratio for Afni6 and vista6 (V1-V2-V3, R2>25%)')
xlabel('Aspect Ratio')

% SYNTHETIC
bdir ='/Users/glerma/toolboxes/PRFmodel/local/ellipse/BIDS/derivatives/prfreport/sub-ellipse/ses-tr1dur300v2';
KK   = load(fullfile(bdir,'sub-ellipse_ses-tr1dur300v2-prf_acq-normal_run-01_bold.mat'));

subplot(2,2,2)
% I did the same in the server with the vista6's r2 results that are calculated
% by default but they are not in the prfreports
% Then using the same method as above, I calculated the r2 for AFNI
% Here I read them and plot the histograms
bdir ='/Users/glerma/toolboxes/PRFmodel/local/ellipse/BIDS/derivatives/prfreport/sub-ellipse/ses-tr1dur300v2';
afni6r2  = niftiRead(fullfile(bdir,'sub-ellipse_ses-tr1dur300v2_task-prf_acq-normal_run-01_r2calccod_AFNI6.nii.gz'));
vista6r2 = niftiRead(fullfile(bdir,'sub-ellipse_ses-tr1dur300v2_task-prf_acq-normal_run-01_r2calccod_VISTA6.nii.gz'));
afni6r2  = afni6r2.data;
vista6r2 = vista6r2.data;

afni6r2  = afni6r2(KK.compTable.HRFtype=="afni_spm");
vista6r2 = vista6r2(KK.compTable.HRFtype=="vista_twogammas");

% FILTER 1
% afind    = afni6r2  > .25;
% viind    = vista6r2 > .25;
% afni6r2  = afni6r2(afind);
% vista6r2 = vista6r2(viind);

histogram(100*afni6r2,'FaceColor','b','FaceAlpha',0.5,'NumBins',50);hold on;
plot(median(100*afni6r2)*[1,1],[0,600],'Color','b','LineStyle','--');
histogram(100*vista6r2,'FaceColor',[.5 .5 .5],'NumBins',75);
plot(median(100*vista6r2)*[1,1],[0,600],'Color',[.5 .5 .5 ],'LineStyle','--');
legend({'Synth Afni6','median','Synth Vista6','median'})
title('Synth R2 for Afni6 and vista6')
xlabel('R2 in %')



% SYNTHETIC ASPECT
subplot(2,2,4)


afni6aspect  = KK.compTable.afni6.sMaj ./ KK.compTable.afni6.sMin;
vista6aspect = KK.compTable.vista6.sMaj ./ KK.compTable.vista6.sMin;
afni6aspect  = afni6aspect(afni6aspect < 5);
vista6aspect = vista6aspect(vista6aspect < 5);
histogram(afni6aspect,'FaceColor','b','FaceAlpha',0.5,'NumBins',50);hold on;
xlim([1,5])
plot(median(afni6aspect)*[1,1],[0,600],'Color','b','LineStyle','--');
histogram(vista6aspect,'FaceColor',[.5 .5 .5],'NumBins',50);
plot(median(vista6aspect)*[1,1],[0,600],'Color',[.5 .5 .5 ],'LineStyle','--');
legend({'Synth Afni6','median','Synth Vista6','median'})
title('Synth Aspect Ratio for Afni6 and vista6 (V1-V2-V3)')
xlabel('Aspect Ratio (Ground Truth is 1)')


%% PLOT 6 (A,B,C):Real Data 7T, ONLY VISTA
% restart every time
bylabelsums = nonfilteredbuylabelsums;

doSave     = true;
centerPerc = 90;
eccenInGT  = true;
xlims      = [0,10];
ylims      = [0,10];
tools      = {'afni6','vista6'}; 
subs       = {'115017','164131','536647'};
useLabels  = {    'V1d', 'V2d', 'V3d','V1v', 'V2v', 'V3v'};
Cs         = .65*[1 0 0; 0 1 0; 0 0 1;1 0 0; 0 1 0; 0 0 1];
marks      =     [  '*',   '*',   '*',  'o',   'o',   'o',];
lstyle     =     { '-.',  '-.',  '-.',  '-',   '-',   '-'};

% Obtain the same eccentricities as in the simulations
eccenvalues = linspace(1.5,6.5,6);

% useLabels  = {'V1','V2','V3'};
duration   = 300;
tr         = 1;
% Filter results
sMajMIN    = -inf %1 ; % .5; % 1;
sMinMIN    = -inf %1 ; % .5; % .75;
sMajMAX    = inf %5 %3 ; % 1.5 % 3; % 4;
eccenMIN   = -inf %2 ; % 1 % 2;
eccenMAX   = inf %6 ; % 6;
minR2      = 0.1;
% How many bins
NeccenBins = 6;
NareaBins  = NeccenBins;



bylabelsums = nonfilteredbuylabelsums;
% Apply the restrictions
for nt=1:length(tools)
    tool = tools{nt};
    for nl = 1:length(useLabels)
        lab = useLabels{nl};
        [TH,R]      = cart2pol(bylabelsums.(tool).(lab).x0, bylabelsums.(tool).(lab).y0);
        bylabelsums.(tool).(lab).angle = rad2deg(TH);
        bylabelsums.(tool).(lab).eccen = R;
        bylabelsums.(tool).(lab).area  = pmEllipseArea(2*bylabelsums.(tool).(lab).sMaj, 2*bylabelsums.(tool).(lab).sMin);
        bylabelsums.(tool).(lab) = bylabelsums.(tool).(lab)(...
                                        bylabelsums.(tool).(lab).sMaj  > sMajMIN & ...
                                        bylabelsums.(tool).(lab).sMin  > sMinMIN & ...
                                        bylabelsums.(tool).(lab).sMaj  < sMajMAX & ...
                                        bylabelsums.(tool).(lab).eccen > eccenMIN & ...
                                        bylabelsums.(tool).(lab).eccen < eccenMAX & ...
                                        bylabelsums.(tool).(lab).r2    > minR2,:);
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
        
        bylabelsums.(tool).(lab).aspect = bylabelsums.(tool).(lab).sMaj  ./ bylabelsums.(tool).(lab).sMin;
        
        
    end
end
% Read the synthetic data as well, this is the eccenv2 dataset, with mid and low
% noise levels, with TR=1 and 2, duration 400, and the ground truth aspect ratio
% limited to 1



% Generated TR=1, Dur=300 data to plot alongside with the real data
fprintf('\n\nLoading synthetic TR=1 300sec data')

sub = 'ellipse'; ses = 'tr1dur300v2';
p = fullfile(pmRootPath,'local',sub,'BIDS','derivatives','prfreport',['sub-' sub],['ses-' ses]);
f = ['sub-' sub '_ses-' ses '-prf_acq-normal_run-01_bold.mat'];
tools = {'synth','afni6'};

C = load(fullfile(p,f));
dt = C.compTable;
for nt=1:length(tools)
    dt.(tools{nt}).aspect = dt.(tools{nt}).sMaj ./ dt.(tools{nt}).sMin;
    [TH,R] = cart2pol(dt.(tools{nt}).x0, dt.(tools{nt}).y0);
    dt.(tools{nt}).angle = rad2deg(TH);
    dt.(tools{nt}).eccen = R;
    dt.(tools{nt}).area  = pmEllipseArea(2*dt.(tools{nt}).sMaj, 2*dt.(tools{nt}).sMin);
end
% GT aspect ratio is always one
dt   = dt(dt.synth.aspect==1,:);
A1A2 = dt;
tool = 'afni6';
A1A2 = A1A2(A1A2.(tool).sMaj >= sMajMIN  & ...
             A1A2.(tool).sMaj <= sMajMAX & ...
             A1A2.(tool).eccen >= eccenMIN & ...
             A1A2.(tool).eccen <= eccenMAX , :); % & ...
         % A1A2.HRFtype=="vista_twogammas", :);
%}
%{
A1A2 = A1A2(A1A2.HRFtype=="vista_twogammas", :);
%}
% {
A1A2 = A1A2(A1A2.HRFtype=="afni_spm", :);
%}

%{
unique(A1A2.synth.sMaj)
unique(A1A2.synth.sMin)
unique(A1A2.synth.eccen)
unique(A1A2.HRFtype)
unique(A1A2.noiseLevel)
%}

disp ('... done with load')

aspect1  = A1A2.(tool).aspect(A1A2.noiseLevel=="low");
B1=prctile(aspect1, [5, 95]);inRange1 = aspect1 >= B1(1) & aspect1 <= B1(2);
aspect1  = aspect1(inRange1);
% sprintf('Low noise: Min aspect ratio for vista 6 is %g and max is %g', min(aspect1),max(aspect1))

aspect2  = A1A2.(tool).aspect(A1A2.noiseLevel=="mid");
B2=prctile(aspect2, [5, 95]);inRange2 = aspect2 >= B2(1) & aspect2 <= B2(2);
aspect2  = aspect2(inRange2);
% sprintf('Mid noise: Min aspect ratio for vista 6 is %g and max is %g', min(aspect2),max(aspect2))
mediansyntheticdatalow = median(aspect1);
mediansyntheticdatamid = median(aspect2);

% Discretize by label, to bin the eccentricities
A1A2.(tool).Y = zeros(size(A1A2.(tool).aspect));
A1A2.(tool).Y(A1A2.noiseLevel=="mid") = discretize(A1A2.(tool).eccen(A1A2.noiseLevel=="mid"),eccenvalues); 
A1A2.(tool).Y(A1A2.noiseLevel=="low") = discretize(A1A2.(tool).eccen(A1A2.noiseLevel=="low"),eccenvalues); 
A1A2low = A1A2(A1A2.noiseLevel=="low", :);
A1A2mid = A1A2(A1A2.noiseLevel=="mid", :);
% Apply percentiles and plot individually
% Create the vectors and then plot all together
vistaMedLowEcc = zeros(1,length(eccenvalues)-1);
vistaMedMidEcc = zeros(1,length(eccenvalues)-1);
for ne=1:(length(eccenvalues)-1)
    aspecclow = A1A2low.(tool).aspect(A1A2low.(tool).Y==ne);
    aspeccmid = A1A2mid.(tool).aspect(A1A2mid.(tool).Y==ne);
    % Median
    vistaMedLowEcc(ne) = median(aspecclow);
    vistaMedMidEcc(ne) = median(aspeccmid);
end

% prepare data
tool = 'afni6';
for nl  = 1:length(useLabels)
    lab = useLabels{nl};
    % Discretize by label, to bin the eccentricities
    bylabelsums.(tool).(lab).Y = discretize(bylabelsums.(tool).(lab).eccen,eccenvalues); 
    % Apply percentiles and plot individually
    % Create the vectors and then plot all together
    aspectmedecc.(lab) = zeros(1,length(eccenvalues)-1);
    aspectmeanecc.(lab) = zeros(1,length(eccenvalues)-1);
    aspectminecc.(lab) = zeros(1,length(eccenvalues)-1);
    aspectmaxecc.(lab) = zeros(1,length(eccenvalues)-1);
    
    for ne=1:(length(eccenvalues)-1)
        % ECC - ASPECT
        aspecc = bylabelsums.(tool).(lab).aspect(bylabelsums.(tool).(lab).Y==ne);
        % Median and std
        if isempty(aspecc)
            aspectmedecc.(lab)(ne) = 0;
            aspectmeanecc.(lab)(ne) = 0;
            aspectminecc.(lab)(ne) = 0;
            aspectmaxecc.(lab)(ne) = 0;
        else
            aspectmedecc.(lab)(ne) = median(aspecc);
            aspectmeanecc.(lab)(ne) = mean(aspecc);
            aspectminecc.(lab)(ne) = min(aspecc);  % They use SEM, check
            aspectmaxecc.(lab)(ne) = max(aspecc);
        end
    end
end


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
    % as = [as;plot(Emidpoints,aspectmeanecc.(lab),'Color',Cs(nl,:), ... 
    as = [as;plot(Emidpoints,aspectmedecc.(lab),'Color',Cs(nl,:), ...
              ... % marks(nl),'MarkerSize',12, ...
              'LineStyle',lstyle{nl},'LineWidth',2)];hold on
    % a  = plot([Emidpoints;Emidpoints] ,...
    %           [aspectminecc  ; aspectmaxecc], ...
    %           'Color','k','LineStyle','-','LineWidth',3);  % 0.75*[0 1 0]
end

lowplot = plot(Emidpoints, vistaMedLowEcc,'k--','LineWidth',3);
midplot = plot(Emidpoints, vistaMedMidEcc,'k:','LineWidth',3);

legend([useLabels,'Synth Low Noise','Synth Mid Noise'], 'Location','eastoutside')
title(strrep(sprintf('%s_TR-%i_Dur-%is_C.I.-%i',...
    tool,tr,duration,centerPerc),'_','\_'))
grid on
xlabel('Eccentricity (deg)')
ylabel('pRF aspect ratio')
xlim([Emidpoints(1)-.2,Emidpoints(end)+.2]);
ylim([0,4.5]);
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
tool = 'afni6';
for nl  = 1:length(useLabels)
    subplot(2,3,nl)
    lab = useLabels{nl};
    % Obtain aspect
    aspectvista = bylabelsums.(tool).(lab).aspect;
    aspects = [aspects;aspectvista];
    
    % Plot it
    h = histogram(aspectvista,20,'Normalization','probability');
    set(h,'LineWidth',2,'EdgeColor',[.5 .5 .5],'EdgeAlpha',0,'FaceAlpha',1,'FaceColor',[.5 .5 .5]);hold on
    plot(median(aspectvista)*[1,1],[0,max(h.Values)],'r-')    
    
    % Add the low noise and mid noise lines now

    
    xlim([1,5])
    ylim([0,0.28])
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

aspect2  = aspect2(aspect2 < 5);
% hmid = histogram(aspect2,'DisplayStyle','stairs','BinWidth',h.BinWidth,'Normalization','probability');
hmid = histogram(aspect2,'BinWidth',binWidth,'Normalization','probability');
% set(hmid,'LineWidth',2,'EdgeColor',[.5 .5 .5 ],'LineStyle','-','EdgeAlpha',.5,'FaceAlpha',.5,'FaceColor',[.5 .5 .5 ]);
set(hmid,'LineWidth',2,'EdgeColor','k','FaceAlpha',1,'FaceColor','k');hold on
blow = plot(median(aspect2)*[1,1],[0,.1],'LineWidth',2,'Color','k','LineStyle','--'); 


h = histogram(aspects,35,'Normalization','probability','BinWidth',binWidth);hold on
% set(h,'LineWidth',2,'EdgeColor','k','FaceAlpha',1,'FaceColor','k');
set(h,'LineWidth',2,'EdgeColor',[.5 .5 .5 ],'LineStyle','-','EdgeAlpha',0,'FaceAlpha',.75,'FaceColor',[.5 .5 .5 ]);
a = plot(median(aspects)*[1,1],[0,.1],'Color',[.5 .5 .5 ],'LineStyle','--');

tool = 'afni6';
% Add the low noise and mid noise lines now
%  aspect1  = aspect1(aspect1 < 4);
% hlow = histogram(aspect1,'DisplayStyle','stairs','BinWidth',h.BinWidth,'Normalization','probability');
% hlow = histogram(aspect1,'BinWidth',h.BinWidth,'Normalization','probability');
% set(hlow,'LineWidth',2,'EdgeColor',[1 .5 .5 ],'LineStyle','-','EdgeAlpha',0,'FaceAlpha',.5,'FaceColor',[1 .5 .5 ]);
% alow = plot(median(aspect1)*[1,1],[0,.1],'LineWidth',2,'Color',[1 .5 .5 ],'LineStyle','-'); 

xlim([1,4]);

% legend([h;a;hlow;hmid],{'Experimental Data','Median of Exp. Data','Synth Low Noise','Synth Mid Noise'});

legend([h;a;hmid;blow],{'Experimental Data','Median of Exp. Data','Synth Mid Noise','Median of Synth. Data'});
fname = fullfile(saveTo, strcat(fnameRoot,['.' ext]));
saveas(gcf,fname,ext);
fprintf('\nSaved %s\n', fname)



%% PLOT S8  % THETA vs ANGLE
% prepare data
tool = 'afni6';
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
v4 = 100*compTable.afni6.r2;
% Eccentricity values for filtering
[~,R6] = cart2pol(compTable.vista6.x0, compTable.vista6.y0);
[~,R4] = cart2pol(compTable.afni6.x0, compTable.afni6.y0);
compTable.vista6.eccen = R6;
compTable.afni6.eccen = R4;


% Filter by variance explained
v6ind = (v6 > 100*minR2) & ...
        (compTable.vista6.sMaj  > sMajMIN)  & (compTable.vista6.sMaj < sMajMAX) & ...
        (compTable.vista6.eccen > eccenMIN) & (compTable.vista6.eccen < eccenMAX);
v4ind = (v4 > 100*minR2) & ...
        (compTable.afni6.sMaj  > sMajMIN)  & (compTable.afni6.sMaj < sMajMAX) & ...
        (compTable.afni6.eccen > eccenMIN) & (compTable.afni6.eccen < eccenMAX);
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


kk = mrvNewGraphWin('R2 afni6/vista6');
set(kk,'Position',[0.007 0.62  0.5 0.5]);
h = histogram(v6f - v4f,'Normalization','probability');  % 'DisplayStyle','stairs'
set(h,'LineWidth',2,'EdgeColor',[.5 .5 .5],'EdgeAlpha',0,'FaceAlpha',1,'FaceColor',[.5 .5 .5]);hold on
a = plot(median(v6f - v4f)*[1,1],[0,max(h.Values)],'r-');
xlabel('Delta R2 (vista6 - afni6; in %)')
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

