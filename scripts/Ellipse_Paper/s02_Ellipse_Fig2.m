function s02_Ellipse_Fig2
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
for nt=1:length(tools)
    tool   = tools{nt};
    useHRF = useHRFs{nt};
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
        % Obtain eccen  vals, this is going to be the x axis
        eccenvals = unique(dt.synth.eccen);
        Cs         = 0.65*distinguishable_colors(1+length(eccenvals),'w');
        
        % Select a size, lets take the smalles one for now
        dtcirc  = dt(dt.synth.aspect==1,:);
        dtcirc  = dtcirc(dtcirc.synth.sMaj==checksize,:);
        assert(unique(dtcirc.synth.sMin)==checksize)
        aspect1   = unique(dtcirc.(tool).aspect);
        
        dtellip = dt(dt.synth.aspect==2,:);
        dtellip = dtellip(dtellip.synth.sMaj == ellipsizes{ns}(1),:);
        assert(unique(dtellip.synth.sMin)    == ellipsizes{ns}(2));
        aspect2   = unique(dtellip.(tool).aspect);
        
        
        ystart=zeros(size(eccenvals));
        ystop=8*ones(size(eccenvals));
        plot([eccenvals.';eccenvals.'],[ystart.';ystop.'],'LineWidth',.7,'LineStyle','-.','Color','k')
        hold on
        
        a = plot(eccenvals, aspect1,'-kx');
        b = plot(eccenvals, aspect2,'--ko');
        % Add dashed lines with GT
        % plot([0,max(eccenvals)],[2,2],'LineWidth',1.5,'LineStyle','--','Color','c');
        % plot([0,max(eccenvals)],[1,1],'LineWidth',1.5,'LineStyle','--','Color',0.75*[0 1 0])
        
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
    end
end

saveas(gcf,fullfile(saveTo, strcat(fnameRoot,['.' ext])),ext);

end

