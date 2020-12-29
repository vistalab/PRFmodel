function pmEllipse_Fig2
% Make Figure 2
%
%
% See also
%  s00_MainFiguresScript

%%
ext  = 'png'; % Could be svg
saveTo = fullfile(pmRootPath,'local','figures');  % Folder path
if ~exist(saveTo,'dir'), mkdir(saveTo); end

%% Noiseless: eccentricity plots, but noiseless
sub = 'ellipse'; ses = 'noiselessv2';
p = fullfile(pmRootPath,'local',sub,'BIDS','derivatives','prfreport',['sub-' sub],['ses-' ses]);
f = ['sub-' sub '_ses-' ses '-prf_acq-normal_run-01_bold.mat'];

theFitFile = fullfile(p,f);
if isfile(theFitFile)
    load(theFitFile,'compTable')
else
    disp('Calculate noise free tests')
    afnicompTable  = pmNoiseFreeTests('afni6' ,'eccen',true);
    vistacompTable = pmNoiseFreeTests('vista6','eccen',true);
    % Save it so that we don't need to generate every time
    compTable        = afnicompTable;
    compTable.vista6 = vistacompTable.vista6;
    
    disp('Saving noise free tests')
    if ~exist(p,'dir'); mkdir(p); end
    save(theFitFile, 'compTable')
end

%% SUMMARY PLOT: RATIO 1 and 2
% As a summary, plot the eccen vs aspect ratio
fnameBegin = 'Fig2_NoiselessEccSimRatio1and2';
% ext        = 'svg';
nlvl       = "none";
eccenInGT  = true;
checksizes = [0.5    ,     1,       2,    3];
ellipsizes = {[1,0.5], [2,1], [3,1.5], [4,2]};
tools      = {'afni6'          , 'vista6'};  % 'vista6' 'afni6' 'vista4' 'afni4'
% for vista is vista_twogammas, but only one value in synth, see B table
useHRFs    = {'afni_spm'       , 'afni_spm' };
duration   = 400;
tr         = 2;
nrow = 2; ncol = 4;
% Create main plot with the ground truth lines
fnameEnd   = sprintf('TR-%i_Dur-%is_Noise-%s',tr,duration,nlvl);
fnameRoot  = strcat(fnameBegin,'-', fnameEnd);

disp(fnameRoot)
kk = mrvNewGraphWin(fnameRoot);
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  1  0.6]);
np=0;

% ct = array2table([zeros(32,1),zeros(32,1),zeros(32,1)]);
% ct.Properties.VariableNames = {'aspect','x0','y0'};
CT1 = table();
CT2 = table();
for nt=1:length(tools)
    tool   = tools{nt};
    useHRF = useHRFs{nt};
    % Initialize tables for storing plotted values    
    CT1tool   = table();
    CT2tool   = table();
        
    for ns=1:length(checksizes)
        np=np+1; subplot(nrow, ncol, np)
        checksize  = checksizes(ns);
        dt         = compTable;
        % MAKE THIS A FUNCTION
        % Obtain eccentricity and polar angle
        [TH,R]           = cart2pol(dt.synth.x0, dt.synth.y0);
        dt.synth.angle   = rad2deg(TH);
        dt.synth.eccen   = R;
        dt.synth.aspect  = dt.synth.sMaj ./ dt.synth.sMin;
        
        [TH,R]           = cart2pol(dt.(tool).x0, dt.(tool).y0);
        dt.(tool).angle  = rad2deg(TH);
        dt.(tool).eccen  = R;
        dt.(tool).aspect = dt.(tool).sMaj  ./ dt.(tool).sMin;
        
        % Check that we are getting the values we want
        xvalues = unique(dt.synth.eccen);
        isclose(linspace(1,9,8)',xvalues,'tolerance',0.001);
        
        % Filter all that we can filter
        % Noise levels
        dt = dt(dt.noiseLevel==nlvl,:);
        % Assert and remove the rest options
        nls=unique(dt.noiseLevel);assert(nls==nlvl);
        % We want to use just its own HRF, remove the vista one
        dt = dt(dt.HRFtype==string(useHRF),:);
        
        % Select a size, lets take the smallest one for now
        % Create a subtable of eight values for circles: aspect = 1
        dtcirc  = dt(dt.synth.aspect==1,:);
        dtcirc  = dtcirc(dtcirc.synth.sMaj==checksize,:);
        assert(unique(dtcirc.synth.sMin)==checksize)
        % Create a subtable of eight values for ellipses: aspect = 2
        dtellip = dt(dt.synth.aspect==2,:);
        dtellip = dtellip(dtellip.synth.sMaj == ellipsizes{ns}(1),:);
        assert(unique(dtellip.synth.sMin)    == ellipsizes{ns}(2));
        
        % Obtain eccen vals, this is going to be the x axis
        if isequal(dtcirc.synth.eccen, dtellip.synth.eccen)
            eccenvals = dtcirc.synth.eccen;
            Cs        = 0.65*distinguishable_colors(1+length(eccenvals),'w');
        else
            error('Check the unique eccentricity values')
        end
        
        ystart=zeros(size(eccenvals));
        ystop=8*ones(size(eccenvals));
        plot([eccenvals.';eccenvals.'],[ystart.';ystop.'],'LineWidth',.7,'LineStyle','-.','Color','k')
        hold on
        
        % If we want to show only the aspect ratio
        a = plot(eccenvals, dtcirc.(tool).aspect,'-kx');hold on
        b = plot(eccenvals, dtellip.(tool).aspect,'--ko');
        % If we want to show in the same plot distance from GT in eccentricity
        % a = plot(dtcirc.(tool).eccen, dtcirc.(tool).aspect,'-kx');hold on
        % b = plot(dtellip.(tool).eccen, dtellip.(tool).aspect,'--ko');
        
        % Apply percentiles and plot individually
        title(strrep(sprintf('%s_TR-%i_Dur-%is_Size-%0.1g',tool,tr,duration,checksize),'_','\_'))
        xlabel('Eccentricity (deg)')
        ylabel('pRF aspect ratio')
        ylim([0,8]);
        set(gca, 'FontSize', 16)
        legend([a,b],{['G.T. aspect ratio = 1 (' ...
            num2str(checksize) 'deg /' num2str(checksize) 'deg)'], ...
            ['G.T. aspect ratio = 2 (' ...
            num2str(ellipsizes{ns}(1)) 'deg/' num2str(ellipsizes{ns}(2)) 'deg)']})
        
        % STORE VALUES FOR LATER
        CT1tool   = [CT1tool; dtcirc];
        CT2tool   = [CT2tool; dtellip];
    end
    % STORE VALUES FOR LATER
    CT1.(tool) = CT1tool;
    CT2.(tool) = CT2tool;
end

fname = fullfile(saveTo, strcat(fnameRoot,['.' ext]));
saveas(gcf,fname,ext);
fprintf('\nSaved %s\n', fname)



%% STATS SECTION
% {
% CIRCLE
if isequal(CT1.afni6.synth, CT1.vista6.synth)  
    synth   = CT1.afni6.synth;
    afnitt  = CT1.afni6.afni6;
    vistatt = CT1.vista6.vista6;
else
    error('Synth values are not equal'); 
end

% Calculate Euclidian distances for the centers
afniCenterDist1  = sqrt((synth.x0 - afnitt.x0).^2 + (synth.y0 - afnitt.y0).^2);
vistaCenterDist1 = sqrt((synth.x0 - vistatt.x0).^2 + (synth.y0 - vistatt.y0).^2);

% Differences in aspect ratios
afniAspect1  = abs(synth.aspect - afnitt.aspect);
vistaAspect1 = abs(synth.aspect - vistatt.aspect);

fprintf('\nASPECT RATIO GROUND TRUTH = 1')
fprintf('\npRF Center')
fprintf('\n[Figure 2: mean(std) vs G.T.] AFNI6: %.2g(%.2g), VISTA6: %.2g(%.2g)',...)
    mean(afniCenterDist1), std(afniCenterDist1), ...
    mean(vistaCenterDist1), std(vistaCenterDist1))
[H,P,CI,STATS] = ttest(afniCenterDist1, vistaCenterDist1);
fprintf('\n[Figure 2: paired t-test AFNIvsVISTA] t: %.2g, p: %.2g',STATS.tstat,P)

fprintf('\npRF Aspect Ratio')
fprintf('\n[Figure 2: mean(std) vs G.T.] AFNI6: %.2g(%.2g), VISTA6: %.2g(%.2g)',...)
    mean(afniAspect1), std(afniAspect1), ...
    mean(vistaAspect1), std(vistaAspect1))
[H,P,CI,STATS] = ttest(afniAspect1, vistaAspect1);
fprintf('\n[Figure 2: paired t-test AFNIvsVISTA] t: %.2g, p: %.2g\n',STATS.tstat,P)

% ELLIPSE
if isequal(CT2.afni6.synth, CT2.vista6.synth)  
    synth   = CT2.afni6.synth;
    afnitt  = CT2.afni6.afni6;
    vistatt = CT2.vista6.vista6;
else
    error('Synth values are not equal'); 
end

% Calculate Euclidian distances for the centers
afniCenterDist2  = sqrt((synth.x0 - afnitt.x0).^2 + (synth.y0 - afnitt.y0).^2);
vistaCenterDist2 = sqrt((synth.x0 - vistatt.x0).^2 + (synth.y0 - vistatt.y0).^2);

% Differences in aspect ratios
afniAspect2  = abs(synth.aspect - afnitt.aspect);
vistaAspect2 = abs(synth.aspect - vistatt.aspect);
fprintf('\nASPECT RATIO GROUND TRUTH = 2')
fprintf('\npRF Center')
fprintf('\n[Figure 2: mean(std) vs G.T.] AFNI6: %.2g(%.2g), VISTA6: %.2g(%.2g)',...)
    mean(afniCenterDist2), std(afniCenterDist2), ...
    mean(vistaCenterDist2), std(vistaCenterDist2))
[H,P,CI,STATS] = ttest(afniCenterDist2,vistaCenterDist2);
fprintf('\n[Figure 2: paired t-test AFNIvsVISTA] t: %.2g, p: %.2g\n',STATS.tstat,P)

fprintf('\npRF Aspect Ratio')
fprintf('\n[Figure 2: mean(std) vs G.T.] AFNI6: %.2g(%.2g), VISTA6: %.2g(%.2g)',...)
    mean(afniAspect2), std(afniAspect2), ...
    mean(vistaAspect2), std(vistaAspect2))
[H,P,CI,STATS] = ttest(afniAspect2,vistaAspect2);
fprintf('\n[Figure 2: paired t-test AFNIvsVISTA] t: %.2g, p: %.2g\n',STATS.tstat,P)
%}

end

