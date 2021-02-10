function pmEllipse_Fig3
% Make Figure 3
%
%
% See also
%  s00_MainFiguresScript

%% EDIT IF REQUIRED
ext  = 'png'; % Could be svg
saveTo = fullfile(pmRootPath,'local','figures');  % Folder path
if ~exist(saveTo,'dir'), mkdir(saveTo); end

%% LOAD THE DATA
disp('Loading Figure 3 data')

sub = 'ellipse'; ses = 'eccsv2';
p = fullfile(pmRootPath,'local',sub,'BIDS','derivatives','prfreport',['sub-' sub],['ses-' ses]);
f = ['sub-' sub '_ses-' ses '-prf_acq-normal_run-01_bold.mat'];
A = load(fullfile(p,f));

% It seems that AFNI-s theta was not correctly corrected. This has been fixed now.
% It only affects to results in sub-ellipse/ses-*v2
% if strcmp(tool,'afni6');A.compTable.afni6.Th = A.compTable.afni6.Th + deg2rad(90);end
% Add the SNR values (this will come from prfreport in the future)
sub = 'ellipse'; ses = 'eccsv2SNR';
p = fullfile(pmRootPath,'local',sub,'BIDS','derivatives','prfsynth',['sub-' sub],['ses-' ses]);
f = ['sub-' sub '_ses-' ses '_task-prf_acq-normal_run-01_bold.json'];
B = struct2table(jsondecode(fileread(fullfile(p,f))));
A.compTable.SNR = B.SNR;

tools   = {'afni6','vista6'};
set(0,'defaultAxesFontName', 'Arial')
set(0,'defaultTextFontName', 'Arial')

disp('Done with load.')

%% PLOT

% Generic values coming from the config.json
onlyCenters = false;
userfsize   = 2;
location    = [3.1315,3.1315];  % [3,3]; %
useHRF      = {};
centerPerc  = 90;
lineStyle   = '-';
lineWidth   = 0.7;
fontsize    = 14;
noiselevel  = {'low'};
addtext     = true;
useellipse  = true;
color       = [0.5,0.5,0.5];
xlims       = [0,5.5];
ylims       = [0,5.5];
xtick       = [1,2,3,4,5];
ytick       = [1,2,3,4,5];
addcihist   = false;
addcibar    = false;
newWin      = false;

numanalysis = length(tools);

useHRFs = cell(1,numanalysis);
for nj=1:numanalysis
    tool = tools{nj};
    switch tool
        case {'vista','mrvista','vistasoft','vista4','vista6'}
            useHRF = 'vista_twogammas';
        case {'pop','popeye'}
            useHRF = 'popeye_twogammas';
        case {'afni','afni4','afni6','afnidog'}
            useHRF = 'afni_spm';
        case {'aprf','analyzeprf'}
            useHRF = 'canonical';
        otherwise
            warning('%s not recorded, using vista_twogammas as default',tool)
    end
    useHRFs{nj} = useHRF;
end

for nslvl = noiselevel
    fnameRoot = sprintf('Fig3_R1_CloudPlots4x4_Noise-%s_Size-%i', nslvl{:}, userfsize);
    mm        = mrvNewGraphWin(fnameRoot,[]);  % Turn off Visible to run it in the server or Docker container
    
    set(mm,'Units','centimeters','Position',[0 0 10*numanalysis 10*numanalysis]);
    np      = 0;
    
    % AFNI
    subplot(numanalysis,numanalysis,1)
    tool = {'afni6'}; useHRF = 'afni_spm';
    pmCloudOfResults(A.compTable   , tool ,'onlyCenters',onlyCenters ,...
        'userfsize' , userfsize, ...
        'centerPerc', centerPerc    ,'useHRF'     ,useHRF,...
        'lineStyle' , lineStyle, ...
        'lineWidth' , lineWidth     ,'noiselevel' ,nslvl{:} , ...
        'useellipse', useellipse, 'addsnr',true,...
        'location',location,...
        'addtext',addtext, 'adddice',false,'addsnr',true,...
        'color', color, 'xlims',xlims,'ylims',ylims,'fontsize', fontsize, ...
        'xtick',xtick,'ytick',ytick, 'addcibar', addcibar,'addcihist', addcihist,  ...
        'newWin'    , newWin ,'saveTo','','saveToType',ext)
    
    
    % Select the ground-truth
    afnitt = A.compTable.afni6(A.compTable.noiseLevel==string(nslvl{:}) & ...
        A.compTable.HRFtype==string(useHRF) & ...
        A.compTable.synth.sMaj==2 & ...
        A.compTable.synth.sMin==2 & ...
        A.compTable.synth.x0==3.1315 & ...
        A.compTable.synth.y0==3.1315,:);
    
    % Calculate the aspect ratios and plot
    subplot(numanalysis,numanalysis,2)
    afnitt.aspect    = afnitt.sMaj  ./ afnitt.sMin;
    h = histogram(afnitt.aspect,100);hold on
    
    % Plots a red line to show the median value of the aspect ratios
    medaspect  = median(afnitt.aspect);
    plot(medaspect*[1,1],[0,18],'r-','LineWidth',1)
    set(h,'LineWidth',2,'EdgeColor','k','FaceAlpha',1,'FaceColor','k');hold on
    
    title('AFNI Elliptical')
    xlabel('Aspect Ratio')
    set(gca,'FontName', 'Arial','FontSize',16)
    set(gca,'xlim',[1 5]); grid on;
    
    % VISTA
    subplot(numanalysis,numanalysis,3)
    tool = {'vista6'}; useHRF = 'vista_twogammas';
    pmCloudOfResults(A.compTable   , tool ,'onlyCenters',onlyCenters ,...
        'userfsize' , userfsize, ...
        'centerPerc', centerPerc    ,'useHRF'     ,useHRF,...
        'lineStyle' , lineStyle, ...
        'lineWidth' , lineWidth     ,'noiselevel' ,nslvl{:} , ...
        'useellipse', useellipse, 'addsnr',true,...
        'location',location,...
        'addtext',addtext, 'adddice',false,'addsnr',true,...
        'color', color, 'xlims',xlims,'ylims',ylims,'fontsize', fontsize, ...
        'xtick',xtick,'ytick',ytick, 'addcibar', addcibar,'addcihist', addcihist,  ...
        'newWin'    , newWin ,'saveTo'     ,'','saveToType',ext)
        
    % Select the aspect ratio data
    vistatt = A.compTable.vista6(A.compTable.noiseLevel==string(nslvl{:}) & ...
        A.compTable.HRFtype==string(useHRF) & ...
        A.compTable.synth.sMaj==2 & ...
        A.compTable.synth.sMin==2 & ...
        A.compTable.synth.x0==3.1315 & ...
        A.compTable.synth.y0==3.1315,:);
    
    % Calculate the aspect ratio and plot the histogram
    subplot(numanalysis,numanalysis,4)
    vistatt.aspect    = vistatt.sMaj  ./ vistatt.sMin;
    h = histogram(vistatt.aspect,100); hold on
    set(h,'LineWidth',2,'EdgeColor','k','FaceAlpha',1,'FaceColor','k');hold on
    
    % Plot the red line showing the median
    medaspect = median(vistatt.aspect);
    plot(medaspect*[1,1],[0,18],'r-','LineWidth',1)
    title('mrVista Elliptical')
    xlabel('Aspect Ratio')
    set(gca,'FontName', 'Arial','FontSize',16)
    set(gca,'xlim',[1 5]); grid on;
    fname = fullfile(saveTo, strcat(fnameRoot,['.' ext]));
    saveas(gcf,fname,ext);
    fprintf('\nSaved %s\n', fname)
end



%% STATS SECTION

synth        = afnitt;
synth.x0     = 3.1315*ones(height(afnitt),1);
synth.y0     = 3.1315*ones(height(afnitt),1);
synth.sMaj   = 2*ones(height(afnitt),1);
synth.sMin   = 2*ones(height(afnitt),1);
% See pmEllipse_FigS6.m to understand how this error is calculated. Basically,
% there is error, and aspect ratio cant go below 1
synth.aspect = 1.2*ones(height(afnitt),1);

% Differences in aspect ratios
afniAspect  = abs(synth.aspect - afnitt.aspect);
vistaAspect = abs(synth.aspect - vistatt.aspect);

length(afnitt.aspect)
sum(afnitt.aspect<1.05)
sum(afnitt.aspect<1.1)
std(afnitt.aspect)
mean(vistatt.aspect)
std(vistatt.aspect)

% Calculate Euclidian distances for the centers
% afniCenterDist  = sqrt((synth.x0 - afnitt.x0).^2 + (synth.y0 - afnitt.y0).^2);
% vistaCenterDist = sqrt((synth.x0 - vistatt.x0).^2 + (synth.y0 - vistatt.y0).^2);


  

%{
fprintf('\nASPECT RATIO GROUND TRUTH = 1')
fprintf('\npRF Center')
fprintf('\n[Figure 3: mean(std) vs G.T.] AFNI6: %.2g(%.2g), VISTA6: %.2g(%.2g)',...)
    mean(afniCenterDist), std(afniCenterDist), ...
    mean(vistaCenterDist), std(vistaCenterDist))
[H,P,CI,STATS] = ttest(afniCenterDist,vistaCenterDist);
fprintf('\n[Figure 3: paired t-test AFNIvsVISTA] t: %.2g, p: %.2g\n',STATS.tstat,P)
%}






fprintf('\npRF Aspect Ratio')
fprintf('\n[Figure 3: mean(std) vs G.T.] AFNI6: %.2g(%.2g), VISTA6: %.2g(%.2g)\n',...)
    mean(afniAspect), std(afniAspect), ...
    mean(vistaAspect), std(vistaAspect))
[H,P,CI,STATS] = ttest(afniAspect,vistaAspect);
fprintf('[Figure 3: paired t-test AFNIvsVISTA] t: %.2g, p: %.2g\n',STATS.tstat,P)














end

