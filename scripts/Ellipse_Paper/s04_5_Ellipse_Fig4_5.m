function s04_5_Ellipse_Fig4_5
% Make Figures 4 and 5
%
% See also
%  s00_MainFiguresScript

%% Plotting parameters
ext  = 'png'; % Could be svg
saveTo = fullfile(pmRootPath,'local','figures');  % Folder path
if ~exist(saveTo,'dir'), mkdir(saveTo); end

centerPerc = 50;  % Center percentage that we plot

%% ECCENTRICITY TR=2
sub = 'ellipse'; ses = 'eccsv2';
p = fullfile(pmRootPath,'local',sub,'BIDS','derivatives','prfreport',['sub-' sub],['ses-' ses]);
f = ['sub-' sub '_ses-' ses '-prf_acq-normal_run-01_bold.mat'];

A2 = load(fullfile(p,f));
% It seems that AFNI-s theta was not correctly corrected. This has been fixed now.
% It only affects to results in sub-ellipse/ses-*v2
% if strcmp(tool,'afni6');A.compTable.afni6.Th = A.compTable.afni6.Th + deg2rad(90);end
% Add the SNR values (this will come from prfreport in the future)
sub = 'ellipse'; ses = 'eccsv2SNR';
p = fullfile(pmRootPath,'local',sub,'BIDS','derivatives','prfsynth',['sub-' sub],['ses-' ses]);
f = ['sub-' sub '_ses-' ses '_task-prf_acq-normal_run-01_bold.json'];
B = struct2table(jsondecode(fileread(fullfile(p,f))));
A2.compTable.SNR = B.SNR;

% SAME HRF; RATIO 1 and 2
fnameBegin = 'Fig4-5_EccSimHRFokTR2';
% ext        = 'svg';
nlvls      = {"mid","low"};
% centerPerc = 50;
eccenInGT  = true;
checksizes = [0.5,1,2,3];
ellipsizes = {[1,0.5],[2,1],[4,2],[6,3]};
xlims      = [0,10];
ylims      = [0,10];
tools      = {'afni6'   , 'vista6'};
useHRFs    = {'afni_spm', 'vista_twogammas'};
duration   = 400;
tr         = 2;
N          = 100;  % Repeated values
nrow = 2; ncol=4;
for nn = 1:length(nlvls)
    nlvl = nlvls{nn};
    
    % Filter noise values
    DT   = A2.compTable(A2.compTable.noiseLevel==nlvl,:);
    
    % Create main plot with the ground truth lines
    fnameEnd = sprintf('TR-%i_Dur-%is_Noise-%s_C.I.-%i',...
        tr,duration,nlvl,centerPerc);
    fnameRoot = strcat(fnameBegin,'-', fnameEnd);
    disp(fnameRoot)
    kk = mrvNewGraphWin(fnameRoot);
    % Fig size is relative to the screen used. This is for laptop at 1900x1200
    set(kk,'Position',[0.007 0.62  1  0.5]);
    np=0;
    for nt=1:length(tools)
        tool   = tools{nt};
        useHRF = useHRFs{nt};
        dt = DT;
        % MAKE THIS A FUNCTION
        % Obtain eccentricity and polar angle
        [TH,R]         = cart2pol(dt.synth.x0, dt.synth.y0);
        dt.synth.angle = rad2deg(TH);
        dt.synth.eccen = R;
        dt.synth.aspect= dt.synth.sMaj ./ dt.synth.sMin;
        
        [TH,R]           = cart2pol(dt.(tool).x0, dt.(tool).y0);
        dt.(tool).angle  = rad2deg(TH);
        dt.(tool).eccen  = R;
        dt.(tool).aspect = dt.(tool).sMaj  ./ dt.(tool).sMin;
        
        for ns=1:length(checksizes)
            np=np+1;
            subplot(nrow,ncol,np)
            checksize  = checksizes(ns);
            ellipsize  = ellipsizes{ns};
            
            
            % Check that we are getting the values we want
            xvalues = unique(dt.synth.eccen);
            isclose(linspace(1,9,8)',xvalues,'tolerance',0.001);
            
            % Filter all that we can filter
            % Assert and remove the rest options
            nls=unique(dt.noiseLevel);assert(nls==nlvl);
            % Check percentage is 100 based
            if centerPerc < 1; centerPerc = centerPerc*100; end
            % Define the required confidence intervals as two percentiles
            twoTailedRange = (100 - centerPerc) / 2;
            % We want to use just its own HRF, remove the vista one
            dt = dt(dt.HRFtype==string(useHRF),:);
            
            
            
            % Aspect ratio: 1
            dtcirc = dt(dt.synth.aspect==1,:);
            nls=unique(dtcirc.synth.aspect);assert(nls==1);
            % Select the size
            dtcirc = dtcirc(dtcirc.synth.sMaj==checksize,:);
            assert(unique(dtcirc.synth.sMin)==checksize)
            SNRcirc     = dtcirc.SNR;
            meanSNRcirc = mean(SNRcirc);
            stdSNRcirc  = std(SNRcirc);
            
            % Aspect ratio: 2
            dtellip = dt(dt.synth.aspect==2,:);
            nls=unique(dtellip.synth.aspect);assert(nls==2);
            % Select the size
            dtellip = dtellip(dtellip.synth.sMaj==ellipsize(1),:);
            assert(unique(dtellip.synth.sMin)==ellipsize(2))
            SNRellip     = dtellip.SNR;
            meanSNRellip = mean(SNRellip);
            stdSNRellip  = std(SNRellip);
            
            % Obtain eccen  vals, this is going to be the x axis
            eccenvals = unique(dt.synth.eccen);
            
            ystart=zeros(size(eccenvals));
            ystop=8*ones(size(eccenvals));
            plot([eccenvals.';eccenvals.'],[ystart.';ystop.'], ...
                'LineWidth',.7,'LineStyle','-.','Color','k')
            hold on
            % plot([0,max(eccenvals)],[1,1],'LineWidth',1.5,'LineStyle','--','Color',0.75*[0 1 0])
            % plot([0,max(eccenvals)],[2,2],'LineWidth',1.5,'LineStyle','--','Color','c')
            % Cs              = 0.65*distinguishable_colors(1+length(eccenvals),'w');
            
            % Apply percentiles and plot individually
            for ne=1:length(eccenvals)
                % C           = Cs(ne,:);
                ecc         = eccenvals(ne);
                aspectcirc  = dtcirc.(tool).aspect(dtcirc.synth.eccen==ecc);
                aspectellip = dtellip.(tool).aspect(dtellip.synth.eccen==ecc);
                realeccencirc   = dtcirc.(tool).eccen(dtcirc.synth.eccen==ecc);
                realeccenellip  = dtellip.(tool).eccen(dtellip.synth.eccen==ecc);
                Bcirc           = prctile(aspectcirc, [twoTailedRange, 100 - twoTailedRange]);
                Bellip          = prctile(aspectellip, [twoTailedRange, 100 - twoTailedRange]);
                inRangecirc     = aspectcirc>=Bcirc(1) & aspectcirc<=Bcirc(2);
                inRangeellip    = aspectellip>=Bellip(1) & aspectellip<=Bellip(2);
                % Apply
                aspectcicirc    = aspectcirc(inRangecirc);
                realeccencicirc = realeccencirc(inRangecirc);
                
                aspectciellip    = aspectellip(inRangeellip);
                realeccenciellip = realeccenellip(inRangeellip);
                
                % Medians
                aspectmedcirc   = median(aspectcicirc);
                aspectmincirc   = min(aspectcicirc);
                aspectmaxcirc   = max(aspectcicirc);
                realeccenmedcirc= median(realeccencicirc);
                realeccenmincirc= min(realeccencicirc);
                realeccenmaxcirc= max(realeccencicirc);
                
                aspectmedellip   = median(aspectciellip);
                aspectminellip   = min(aspectciellip);
                aspectmaxellip   = max(aspectciellip);
                realeccenmedellip= median(realeccenciellip);
                realeccenminellip= min(realeccenciellip);
                realeccenmaxellip= max(realeccenciellip);
                
                
                
                % Plot it
                if eccenInGT
                    as = scatter(ecc,aspectmedcirc,80,'k','filled');
                    a  = plot(ecc * [1,1],...
                        [aspectmincirc  , aspectmaxcirc], ...
                        'Color','k','LineStyle','-','LineWidth',3);  % 0.75*[0 1 0]
                    
                    bs = scatter(ecc+.15,aspectmedellip,80,'k^','filled');
                    b  = plot((ecc+.15) * [1,1],...
                        [aspectminellip  , aspectmaxellip], ...
                        'Color','k','LineStyle',':','LineWidth',2);
                else
                    scatter(realeccenmed,aspectmed,60,0.75*[0 1 0],'filled')
                    hax = plot([realeccenmin, realeccenmax],...
                        aspectmed*[1,1], ...
                        'Color',0.75*[0 1 0],'LineStyle','-','LineWidth',2); % 'Color','k',
                    vax = plot(realeccenmed * [1,1],...
                        [aspectmin  , aspectmax], ...
                        'Color',0.75*[0 1 0],'LineStyle','-','LineWidth',2); %
                end
            end
            % SNR will be calculated at the level of the graph
            text(1.1*xlims(1),1.1*ylims(1), ...
                sprintf('SNRcirc:%.2g(±%.2g) | SNRellip:%.2g(±%.2g)', ...
                meanSNRcirc, stdSNRcirc,meanSNRellip, stdSNRellip), ...
                'FontWeight','bold','FontSize',12)
            legend([as,bs],...
                {sprintf('G.T. Aspect = 1(%g deg/%g deg)',checksize,checksize), ...
                sprintf('G.T. Aspect = 2(%g deg/%g deg)',ellipsize(1),ellipsize(2))})
            title(strrep(sprintf('%s_TR-%i_Dur-%is_Noise-%s_C.I.-%i_size-%0.1g',...
                tool,tr,duration,nlvl,centerPerc,checksize),'_','\_'))
            
            xlabel('Eccentricity')
            ylabel('pRF aspect ratio')
            ylim([0,8]);
            set(gca, 'FontSize', 16)
        end
    end
    saveas(gcf,fullfile(saveTo, strcat(fnameRoot,['.' ext])),ext);
end


%% ECCENTRICITY TR=1
sub = 'ellipse'; ses = 'eccsv2TR1';
p = fullfile(pmRootPath,'local',sub,'BIDS','derivatives','prfreport',['sub-' sub],['ses-' ses]);
f = ['sub-' sub '_ses-' ses '-prf_acq-normal_run-01_bold.mat'];
A1 = load(fullfile(p,f));

% SAME HRF; RATIO 1 and 2
fnameBegin = 'Fig4-5_EccSimHRFokTR1';
% ext        = 'svg';
nlvls      = {"mid","low"};
% centerPerc = 50;
eccenInGT  = true;
checksizes = [0.5,1,2,3];
ellipsizes = {[1,0.5],[2,1],[4,2],[6,3]};
xlims       = [0,10];
ylims       = [0,10];
tools      = {'afni6'   , 'vista6'};
useHRFs    = {'afni_spm', 'vista_twogammas'};
duration   = 400;
tr         = 1;
nrow = 2; ncol=4;
for nlvl = nlvls
    nlvl = nlvl{:};
    % Create main plot with the ground truth lines
    fnameEnd = sprintf('TR-%i_Dur-%is_Noise-%s_C.I.-%i',...
        tr,duration,nlvl,centerPerc);
    fnameRoot = strcat(fnameBegin,'-', fnameEnd);
    disp(fnameRoot)
    kk = mrvNewGraphWin(fnameRoot);
    % Fig size is relative to the screen used. This is for laptop at 1900x1200
    set(kk,'Position',[0.007 0.62  1  0.5]);
    np=0;
    for nt=1:length(tools)
        tool   = tools{nt};
        useHRF = useHRFs{nt};
        for ns=1:length(checksizes)
            np=np+1;
            subplot(nrow,ncol,np)
            checksize  = checksizes(ns);
            ellipsize  = ellipsizes{ns};
            dt         = A1.compTable;
            % MAKE THIS A FUNCTION
            % Obtain eccentricity and polar angle
            [TH,R]         = cart2pol(dt.synth.x0, dt.synth.y0);
            dt.synth.angle = rad2deg(TH);
            dt.synth.eccen = R;
            dt.synth.aspect= dt.synth.sMaj ./ dt.synth.sMin;
            
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
            % Check percentage is 100 based
            if centerPerc < 1; centerPerc = centerPerc*100; end
            % Define the required confidence intervals as two percentiles
            twoTailedRange = (100 - centerPerc) / 2;
            % We want to use just its own HRF, remove the vista one
            dt = dt(dt.HRFtype==string(useHRF),:);
            
            % Aspect ratio: 1
            dtcirc = dt(dt.synth.aspect==1,:);
            nls=unique(dtcirc.synth.aspect);assert(nls==1);
            % Select the size
            dtcirc = dtcirc(dtcirc.synth.sMaj==checksize,:);
            assert(unique(dtcirc.synth.sMin)==checksize)
            SNRcirc     = dtcirc.SNR;
            meanSNRcirc = mean(SNRcirc);
            stdSNRcirc  = std(SNRcirc);
            
            % Aspect ratio: 2
            dtellip = dt(dt.synth.aspect==2,:);
            nls=unique(dtellip.synth.aspect);assert(nls==2);
            % Select the size
            dtellip = dtellip(dtellip.synth.sMaj==ellipsize(1),:);
            assert(unique(dtellip.synth.sMin)==ellipsize(2))
            SNRellip     = dtellip.SNR;
            meanSNRellip = mean(SNRellip);
            stdSNRellip  = std(SNRellip);
            
            % Obtain eccen  vals, this is going to be the x axis
            eccenvals = unique(dt.synth.eccen);
            
            ystart=zeros(size(eccenvals));
            ystop=8*ones(size(eccenvals));
            plot([eccenvals.';eccenvals.'],[ystart.';ystop.'], ...
                'LineWidth',.7,'LineStyle','-.','Color','k')
            hold on
            % plot([0,max(eccenvals)],[1,1],'LineWidth',1.5,'LineStyle','--','Color',0.75*[0 1 0])
            % plot([0,max(eccenvals)],[2,2],'LineWidth',1.5,'LineStyle','--','Color','c')
            % Cs              = 0.65*distinguishable_colors(1+length(eccenvals),'w');
            
            % Apply percentiles and plot individually
            for ne=1:length(eccenvals)
                % C           = Cs(ne,:);
                ecc         = eccenvals(ne);
                aspectcirc  = dtcirc.(tool).aspect(dtcirc.synth.eccen==ecc);
                aspectellip = dtellip.(tool).aspect(dtellip.synth.eccen==ecc);
                realeccencirc   = dtcirc.(tool).eccen(dtcirc.synth.eccen==ecc);
                realeccenellip  = dtellip.(tool).eccen(dtellip.synth.eccen==ecc);
                Bcirc           = prctile(aspectcirc, [twoTailedRange, 100 - twoTailedRange]);
                Bellip          = prctile(aspectellip, [twoTailedRange, 100 - twoTailedRange]);
                inRangecirc     = aspectcirc>=Bcirc(1) & aspectcirc<=Bcirc(2);
                inRangeellip    = aspectellip>=Bellip(1) & aspectellip<=Bellip(2);
                % Apply
                aspectcicirc    = aspectcirc(inRangecirc);
                realeccencicirc = realeccencirc(inRangecirc);
                
                aspectciellip    = aspectellip(inRangeellip);
                realeccenciellip = realeccenellip(inRangeellip);
                
                % Medians
                aspectmedcirc   = median(aspectcicirc);
                aspectmincirc   = min(aspectcicirc);
                aspectmaxcirc   = max(aspectcicirc);
                realeccenmedcirc= median(realeccencicirc);
                realeccenmincirc= min(realeccencicirc);
                realeccenmaxcirc= max(realeccencicirc);
                
                aspectmedellip   = median(aspectciellip);
                aspectminellip   = min(aspectciellip);
                aspectmaxellip   = max(aspectciellip);
                realeccenmedellip= median(realeccenciellip);
                realeccenminellip= min(realeccenciellip);
                realeccenmaxellip= max(realeccenciellip);
                
                % Plot it
                if eccenInGT
                    as = scatter(ecc,aspectmedcirc,80,'k','filled');
                    a  = plot(ecc * [1,1],...
                        [aspectmincirc  , aspectmaxcirc], ...
                        'Color','k','LineStyle','-','LineWidth',3);  % 0.75*[0 1 0]
                    
                    bs = scatter(ecc+.15,aspectmedellip,80,'k^','filled');
                    b  = plot((ecc+.15) * [1,1],...
                        [aspectminellip  , aspectmaxellip], ...
                        'Color','k','LineStyle',':','LineWidth',2);
                else
                    scatter(realeccenmed,aspectmed,60,0.75*[0 1 0],'filled')
                    hax = plot([realeccenmin, realeccenmax],...
                        aspectmed*[1,1], ...
                        'Color',0.75*[0 1 0],'LineStyle','-','LineWidth',2); % 'Color','k',
                    vax = plot(realeccenmed * [1,1],...
                        [aspectmin  , aspectmax], ...
                        'Color',0.75*[0 1 0],'LineStyle','-','LineWidth',2); %
                end
            end
            % SNR will be calculated at the level of the graph
            text(1.1*xlims(1),1.1*ylims(1), ...
                sprintf('SNRcirc:%.2g(±%.2g) | SNRellip:%.2g(±%.2g)', ...
                meanSNRcirc, stdSNRcirc,meanSNRellip, stdSNRellip), ...
                'FontWeight','bold','FontSize',12)
            legend([as,bs],...
                {sprintf('G.T. Aspect = 1(%g deg/%g deg)',checksize,checksize), ...
                sprintf('G.T. Aspect = 2(%g deg/%g deg)',ellipsize(1),ellipsize(2))})
            title(strrep(sprintf('%s_TR-%i_Dur-%is_Noise-%s_C.I.-%i_size-%0.1g',...
                tool,tr,duration,nlvl,centerPerc,checksize),'_','\_'))
            
            xlabel('Eccentricity')
            ylabel('pRF aspect ratio')
            ylim([0,8]);
            set(gca, 'FontSize', 16)
        end
    end
    saveas(gcf,fullfile(saveTo, strcat(fnameRoot,['.' ext])),ext);
end

%% SAVE A1 and A2 for CUMSUM

dt = [A1.compTable; A2.compTable];
% Calculate the overal aspect ratio
tools = {'synth','afni6','vista6'};
for nt=1:length(tools)
    dt.(tools{nt}).aspect = dt.(tools{nt}).sMaj ./ dt.(tools{nt}).sMin;
    [TH,R] = cart2pol(dt.(tools{nt}).x0, dt.(tools{nt}).y0);
    dt.(tools{nt}).angle = rad2deg(TH);
    dt.(tools{nt}).eccen = R;
    dt.(tools{nt}).area  = pmEllipseArea(2*dt.(tools{nt}).sMaj, 2*dt.(tools{nt}).sMin);
end
% GT aspect ratio is always one
dtaspect2 = dt(dt.synth.aspect==2,:);
dt        = dt(dt.synth.aspect==1,:);

save(fullfile(pmRootPath,'local','A1A2dt.mat'),'dt')
save(fullfile(pmRootPath,'local','A1A2dtaspect2.mat'),'dtaspect2')


%% SIZES TR=2
sub = 'ellipse'; ses = 'sizesv2';
p = fullfile(pmRootPath,'local',sub,'BIDS','derivatives','prfreport',['sub-' sub],['ses-' ses]);
f = ['sub-' sub '_ses-' ses '-prf_acq-normal_run-01_bold.mat'];
load(fullfile(p,f))

% Add the SNR values (this will come from prfreport in the future)
sub = 'ellipse'; ses = 'sizesv2SNR';
p = fullfile(pmRootPath,'local',sub,'BIDS','derivatives','prfsynth',['sub-' sub],['ses-' ses]);
f = ['sub-' sub '_ses-' ses '_task-prf_acq-normal_run-01_bold.json'];
B = struct2table(jsondecode(fileread(fullfile(p,f))));
compTable.SNR = B.SNR;
B2.compTable = compTable;



fnameBegin = 'Fig4-5_SizeSimTR2';
% ext        = 'svg';
nlvls       = {"mid","low"};
% centerPerc = 50;
eccenInGT  = true;
tools      = {'vista6'          , 'afni6'};  % 'vista6' 'afni6' 'vista4' 'afni4'
useHRFs    = {'vista_twogammas' , 'afni_spm'};
duration   = 400;
tr         = 2;
% tool       = 'afni6';
% useHRF     = 'afni_spm';
for nlvl = nlvls
    nlvl = nlvl{:};
    for nt=1:length(tools)
        tool   = tools{nt};
        useHRF = useHRFs{nt};
        dt     = compTable;
        % MAKE THIS A FUNCTION
        % Obtain eccentricity and polar angle
        [TH,R]         = cart2pol(dt.synth.x0, dt.synth.y0);
        dt.synth.angle = rad2deg(TH);
        dt.synth.eccen = R;
        dt.synth.aspect= dt.synth.sMaj ./ dt.synth.sMin;
        
        [TH,R]           = cart2pol(dt.(tool).x0, dt.(tool).y0);
        dt.(tool).angle  = rad2deg(TH);
        dt.(tool).eccen  = R;
        dt.(tool).aspect = dt.(tool).sMaj  ./ dt.(tool).sMin;
        
        % Filter all that we can filter
        % Noise levels
        dt = dt(dt.noiseLevel==nlvl,:);
        % Assert and remove the rest options
        nls=unique(dt.noiseLevel);assert(nls==nlvl);
        % Aspect ratio: start with synthesized aspect ratio = 1
        dt = dt(dt.synth.aspect==1,:);
        nls=unique(dt.synth.aspect);assert(nls==1);
        
        % Check percentage is 100 based
        if centerPerc < 1; centerPerc = centerPerc*100; end
        % Define the required confidence intervals as two percentiles
        twoTailedRange = (100 - centerPerc) / 2;
        
        % We want to use just its own HRF, remove the other one
        dt = dt(dt.HRFtype==string(useHRF),:);
        
        % Obtain size  vals, this is going to be the x axis
        sizes = unique(dt.synth.sMaj);
        
        
        % Create main plot with the ground truth lines
        fnameEnd = sprintf('%s_TR-%i_Dur-%is_Noise-%s_C.I.-%i',...
            tool,tr,duration,nlvl,centerPerc);
        fnameRoot = strcat(fnameBegin,'-', fnameEnd);
        disp(fnameRoot)
        kk = mrvNewGraphWin(fnameRoot);
        % Fig size is relative to the screen used. This is for laptop at 1900x1200
        set(kk,'Position',[0.007 0.62  0.4  0.4]);
        ystart=ones(size(sizes));
        ystop=6*ones(size(sizes));
        plot([sizes.';sizes.'],[ystart.';ystop.'], ...
            'LineWidth',.7,'LineStyle','-.','Color','k')
        hold on
        % plot([0,6],[1,1],'LineWidth',1.5,'LineStyle','--','Color','k') % 0.75*[0 1 0])
        % Cs              = 0.65*distinguishable_colors(1+length(sizes),'w');
        % Apply percentiles and plot individually
        for ne=1:length(sizes)
            ecc          = sizes(ne);
            aspect       = dt.(tool).aspect(dt.synth.sMaj==ecc);
            realeccen    = dt.(tool).eccen(dt.synth.sMaj==ecc);
            B            = prctile(aspect, [twoTailedRange, 100 - twoTailedRange]);
            inRange      = aspect>=B(1) & aspect<=B(2);
            snrperecc    = dt.SNR(dt.synth.sMaj==ecc);
            % Apply
            aspectci     = aspect(inRange);
            realeccenci  = realeccen(inRange);
            % Medians
            aspectmed    = median(aspectci);
            aspectmin    = min(aspectci);
            aspectmax    = max(aspectci);
            
            realeccenmed = median(realeccenci);
            realeccenmin = min(realeccenci);
            realeccenmax = max(realeccenci);
            
            snrpereccmed = median(snrperecc);
            
            % Plot it
            if eccenInGT
                scatter(ecc,aspectmed,80,'k','filled')
                vax = plot(ecc * [1,1],...
                    [aspectmin  , aspectmax], ...
                    'Color','k','LineStyle','-','LineWidth',3); %
                
                text(ecc, 0, sprintf(' %0.2gdB',snrpereccmed),...
                    'FontSize',12,'Rotation',90)
            else
                scatter(realeccenmed,aspectmed,60,C,'filled')
                hax = plot([realeccenmin, realeccenmax],...
                    aspectmed*[1,1], ...
                    'Color','k','LineStyle','-','LineWidth',2); % 'Color','k',
                vax = plot(realeccenmed * [1,1],...
                    [aspectmin  , aspectmax], ...
                    'Color','k','LineStyle','-','LineWidth',2); %
            end
        end
        
        title(strrep(fnameRoot,'_','\_'))
        xlabel('Radius size (dashed=ground truth)')
        ylabel('pRF aspect ratio (ground truth=1)')
        ylim([0,8]);
        set(gca, 'FontSize', 16)
        saveas(gcf,fullfile(saveTo, strcat(fnameRoot,['.' ext])),ext);
    end
end

%% SIZES TR=1

sub = 'ellipse'; ses = 'sizesv2TR1';
p = fullfile(pmRootPath,'local',sub,'BIDS','derivatives','prfreport',['sub-' sub],['ses-' ses]);
f = ['sub-' sub '_ses-' ses '-prf_acq-normal_run-01_bold.mat'];
load(fullfile(p,f))
% Add the SNR values (this will come from prfreport in the future)
sub = 'ellipse'; ses = 'sizesv2SNR';
p = fullfile(pmRootPath,'local',sub,'BIDS','derivatives','prfsynth',['sub-' sub],['ses-' ses]);
f = ['sub-' sub '_ses-' ses '_task-prf_acq-normal_run-01_bold.json'];
B = struct2table(jsondecode(fileread(fullfile(p,f))));
compTable.SNR = B.SNR;
B1.compTable = compTable;


fnameBegin = 'Fig4-5_SizeSimTR1';
% ext        = 'svg';
nlvls       = {"mid","low"};
% centerPerc = 50;
eccenInGT  = true;
tools      = {'vista6'          , 'afni6'};  % 'vista6' 'afni6' 'vista4' 'afni4'
useHRFs    = {'vista_twogammas' , 'afni_spm'};
duration   = 400;
tr         = 1;
% tool       = 'afni6';
% useHRF     = 'afni_spm';
for nlvl = nlvls
    nlvl = nlvl{:};
    for nt=1:length(tools)
        tool   = tools{nt};
        useHRF = useHRFs{nt};
        dt     = compTable;
        % MAKE THIS A FUNCTION
        % Obtain eccentricity and polar angle
        [TH,R]         = cart2pol(dt.synth.x0, dt.synth.y0);
        dt.synth.angle = rad2deg(TH);
        dt.synth.eccen = R;
        dt.synth.aspect= dt.synth.sMaj ./ dt.synth.sMin;
        
        [TH,R]           = cart2pol(dt.(tool).x0, dt.(tool).y0);
        dt.(tool).angle  = rad2deg(TH);
        dt.(tool).eccen  = R;
        dt.(tool).aspect = dt.(tool).sMaj  ./ dt.(tool).sMin;
        
        % Filter all that we can filter
        % Noise levels
        dt = dt(dt.noiseLevel==nlvl,:);
        % Assert and remove the rest options
        nls=unique(dt.noiseLevel);assert(nls==nlvl);
        % Aspect ratio: start with synthesized aspect ratio = 1
        dt = dt(dt.synth.aspect==1,:);
        nls=unique(dt.synth.aspect);assert(nls==1);
        
        % Check percentage is 100 based
        if centerPerc < 1; centerPerc = centerPerc*100; end
        % Define the required confidence intervals as two percentiles
        twoTailedRange = (100 - centerPerc) / 2;
        
        % We want to use just its own HRF, remove the other one
        dt = dt(dt.HRFtype==string(useHRF),:);
        
        % Obtain size  vals, this is going to be the x axis
        sizes = unique(dt.synth.sMaj);
        
        
        % Create main plot with the ground truth lines
        fnameEnd = sprintf('%s_TR-%i_Dur-%is_Noise-%s_C.I.-%i',...
            tool,tr,duration,nlvl,centerPerc);
        fnameRoot = strcat(fnameBegin,'-', fnameEnd);
        disp(fnameRoot)
        kk = mrvNewGraphWin(fnameRoot);
        % Fig size is relative to the screen used. This is for laptop at 1900x1200
        set(kk,'Position',[0.007 0.62  0.4  0.4]);
        ystart=ones(size(sizes));
        ystop=6*ones(size(sizes));
        plot([sizes.';sizes.'],[ystart.';ystop.'], ...
            'LineWidth',.7,'LineStyle','-.','Color','k')
        hold on
        % plot([0,6],[1,1],'LineWidth',1.5,'LineStyle','--','Color','k') % 0.75*[0 1 0])
        % Cs              = 0.65*distinguishable_colors(1+length(sizes),'w');
        % Apply percentiles and plot individually
        for ne=1:length(sizes)
            ecc          = sizes(ne);
            aspect       = dt.(tool).aspect(dt.synth.sMaj==ecc);
            realeccen    = dt.(tool).eccen(dt.synth.sMaj==ecc);
            B            = prctile(aspect, [twoTailedRange, 100 - twoTailedRange]);
            inRange      = aspect>=B(1) & aspect<=B(2);
            snrperecc    = dt.SNR(dt.synth.sMaj==ecc);
            % Apply
            aspectci     = aspect(inRange);
            realeccenci  = realeccen(inRange);
            % Medians
            aspectmed    = median(aspectci);
            aspectmin    = min(aspectci);
            aspectmax    = max(aspectci);
            
            realeccenmed = median(realeccenci);
            realeccenmin = min(realeccenci);
            realeccenmax = max(realeccenci);
            
            snrpereccmed = median(snrperecc);
            
            % Plot it
            if eccenInGT
                scatter(ecc,aspectmed,80,'k','filled')
                vax = plot(ecc * [1,1],...
                    [aspectmin  , aspectmax], ...
                    'Color','k','LineStyle','-','LineWidth',3); %
                
                text(ecc, 0, sprintf('%0.2gdB',snrpereccmed),...
                    'FontSize',12,'Rotation',90)
            else
                scatter(realeccenmed,aspectmed,60,C,'filled')
                hax = plot([realeccenmin, realeccenmax],...
                    aspectmed*[1,1], ...
                    'Color','k','LineStyle','-','LineWidth',2); % 'Color','k',
                vax = plot(realeccenmed * [1,1],...
                    [aspectmin  , aspectmax], ...
                    'Color','k','LineStyle','-','LineWidth',2); %
            end
        end
        
        title(strrep(fnameRoot,'_','\_'))
        xlabel('Radius size (dashed=ground truth)')
        ylabel('pRF aspect ratio (ground truth=1)')
        ylim([0,8]);
        set(gca, 'FontSize', 16)
        saveas(gcf,fullfile(saveTo, strcat(fnameRoot,['.' ext])),ext);
    end
end



%% SAVE B1 and B2 for CUMSUM plot at the end
dt = [B1.compTable; B2.compTable];
% Calculate the overal aspect ratio
tools = {'synth','afni6','vista6'};
for nt=1:length(tools)
    dt.(tools{nt}).aspect = dt.(tools{nt}).sMaj ./ dt.(tools{nt}).sMin;
    [TH,R] = cart2pol(dt.(tools{nt}).x0, dt.(tools{nt}).y0);
    dt.(tools{nt}).angle = rad2deg(TH);
    dt.(tools{nt}).eccen = R;
    dt.(tools{nt}).area  = pmEllipseArea(2*dt.(tools{nt}).sMaj, 2*dt.(tools{nt}).sMin);
end
% GT aspect ratio is always one
dt = dt(dt.synth.aspect==1,:);

save(fullfile(pmRootPath,'local','B1B2dt.mat'),'dt')



%% Combined histograms: aspect1 (run twice, tr=1, tr=2)

% Plot the histograms for afni and vista with low and mid
% Read the synthetic data

A1A2 = load(fullfile(pmRootPath,'local','A1A2dt.mat'));
A1A2 = A1A2.dt;
% Apply the same restrictions as above, to the synth table
%{
A1A2 = A1A2(A1A2.synth.sMaj > sMajMIN & ...
            A1A2.synth.sMin > sMinMIN & ...
            A1A2.synth.sMaj < sMajMAX & ...
            A1A2.synth.eccen > eccenMIN & ...
            A1A2.synth.eccen < eccenMAX,:);
%}
B1B2 = load(fullfile(pmRootPath,'local','B1B2dt.mat'));
B1B2 = B1B2.dt;
% Apply the same restrictions as above, to the synth table
%{
B1B2 = B1B2(B1B2.synth.sMaj > sMajMIN & ...
            B1B2.synth.sMin > sMinMIN & ...
            B1B2.synth.sMaj < sMajMAX & ...
            B1B2.synth.eccen > eccenMIN & ...
            B1B2.synth.eccen < eccenMAX,:);
%}

A1A2 = [A1A2;B1B2];  % synth.aspect is already = 1
A1A2.afni6.TR  = A1A2.TR;
A1A2.vista6.TR = A1A2.TR;
sMajMIN     = 1;
sMinMIN     = 1;
sMajMAX     = 4;
eccenMIN    = 2;
eccenMAX    = 6;
tr          = 1;

histbins    = 100;

afnilow   = A1A2.afni6(A1A2.noiseLevel=="low" & A1A2.HRFtype=="afni_spm" & ...
    A1A2.synth.sMaj > sMajMIN & ...
    A1A2.synth.sMin > sMinMIN & ...
    A1A2.synth.sMaj < sMajMAX & ...
    A1A2.synth.eccen > eccenMIN & ...
    A1A2.synth.eccen < eccenMAX,:);
B=prctile(afnilow.aspect,[5,95]);inR=afnilow.aspect>=B(1) & afnilow.aspect<=B(2);
afnilow=afnilow(inR,:);
afnimid   = A1A2.afni6(A1A2.noiseLevel=="mid" & A1A2.HRFtype=="afni_spm" & ...
    A1A2.synth.sMaj > sMajMIN & ...
    A1A2.synth.sMin > sMinMIN & ...
    A1A2.synth.sMaj < sMajMAX & ...
    A1A2.synth.eccen > eccenMIN & ...
    A1A2.synth.eccen < eccenMAX,:);
B=prctile(afnimid.aspect,[5,95]);inR=afnimid.aspect>=B(1) & afnimid.aspect<=B(2);
afnimid=afnimid(inR,:);
vistalow  = A1A2.vista6(A1A2.noiseLevel=="low" & A1A2.HRFtype=="vista_twogammas" & ...
    A1A2.synth.sMaj > sMajMIN & ...
    A1A2.synth.sMin > sMinMIN & ...
    A1A2.synth.sMaj < sMajMAX & ...
    A1A2.synth.eccen > eccenMIN & ...
    A1A2.synth.eccen < eccenMAX,:);
B=prctile(vistalow.aspect,[5,95]);inR=vistalow.aspect>=B(1) & vistalow.aspect<=B(2);
vistalow=vistalow(inR,:);
vistamid  = A1A2.vista6(A1A2.noiseLevel=="mid" & A1A2.HRFtype=="vista_twogammas" & ...
    A1A2.synth.sMaj > sMajMIN & ...
    A1A2.synth.sMin > sMinMIN & ...
    A1A2.synth.sMaj < sMajMAX & ...
    A1A2.synth.eccen > eccenMIN & ...
    A1A2.synth.eccen < eccenMAX,:);
B=prctile(vistamid.aspect,[5,95]);inR=vistamid.aspect>=B(1) & vistamid.aspect<=B(2);
vistamid=vistamid(inR,:);

fnameRoot = sprintf('Fig4-5_Histograms_Synth_Aspect-1_TR-%i',tr);
disp(fnameRoot)
% saveToType = 'svg';
kk = mrvNewGraphWin(fnameRoot);
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  1  1]);

subplot(2,2,1)
aspect = afnilow.aspect(afnilow.TR==tr);
h = histogram(aspect, histbins,'Normalization','probability'); hold on
set(h,'LineWidth',2,'EdgeColor','k','FaceAlpha',1,'FaceColor','k');hold on
medaspect = median(aspect);
a=plot(medaspect*[1,1],[0,max(h.Values)],'r-','LineWidth',1);
title(sprintf('Afni low noise, TR=%g',tr))
xlabel('Aspect Ratio (GT Aspect = 1)')
set(gca,'FontName', 'Arial','FontSize',16)

subplot(2,2,2)
aspect = afnimid.aspect(afnimid.TR==tr);
h = histogram(aspect,histbins,'Normalization','probability'); hold on
set(h,'LineWidth',2,'EdgeColor','k','FaceAlpha',1,'FaceColor','k');hold on
medaspect = median(aspect);
a=plot(medaspect*[1,1],[0,max(h.Values)],'r-','LineWidth',1);
title(sprintf('Afni mid noise, TR=%g',tr))
xlabel('Aspect Ratio (GT Aspect = 1)')
set(gca,'FontName', 'Arial','FontSize',16)

subplot(2,2,3)
aspect = vistalow.aspect(vistalow.TR==tr);
h = histogram(aspect,histbins,'Normalization','probability'); hold on
set(h,'LineWidth',2,'EdgeColor','k','FaceAlpha',1,'FaceColor','k');hold on
medaspect = median(aspect);
a=plot(medaspect*[1,1],[0,max(h.Values)],'r-','LineWidth',1);
title(sprintf('mrVista low noise, TR=%g',tr))
xlabel('Aspect Ratio (GT Aspect = 1)')
set(gca,'FontName', 'Arial','FontSize',16)

subplot(2,2,4)
aspect = vistamid.aspect(vistamid.TR==tr);
h = histogram(aspect,histbins,'Normalization','probability'); hold on
set(h,'LineWidth',2,'EdgeColor','k','FaceAlpha',1,'FaceColor','k');hold on
medaspect = median(aspect);
a=plot(medaspect*[1,1],[0,max(h.Values)],'r-','LineWidth',1);
title(sprintf('mrVista mid noise, TR=%g',tr))
xlabel('Aspect Ratio (GT Aspect = 1)')
set(gca,'FontName', 'Arial','FontSize',16)


saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.',ext)),ext);


%% Combined histograms: aspect2 (run twice, tr=1, tr=2)
% Plot the histograms for afni and vista with low and mid
% Read the synthetic data

A1A2 = load(fullfile(pmRootPath,'local','A1A2dtaspect2.mat'));
A1A2 = A1A2.dtaspect2;

A1A2.afni6.TR  = A1A2.TR;
A1A2.vista6.TR = A1A2.TR;
sMajMIN     = 1;
sMajMAX     = 4;
eccenMIN    = 2;
eccenMAX    = 6;
aspectMAX   = 5;
tr          = 1;

histbins    = 50;

afnilow   = A1A2.afni6(A1A2.noiseLevel=="low" & A1A2.HRFtype=="afni_spm" & ...
    A1A2.synth.sMaj > sMajMIN & ...
    A1A2.synth.sMaj < sMajMAX & ...
    A1A2.afni6.aspect < aspectMAX & ...
    A1A2.synth.eccen  > eccenMIN & ...
    A1A2.synth.eccen  < eccenMAX,:);
% B=prctile(afnilow.aspect,[5,95]);inR=afnilow.aspect>=B(1) & afnilow.aspect<=B(2);
% afnilow=afnilow(inR,:);
afnimid   = A1A2.afni6(A1A2.noiseLevel=="mid" & A1A2.HRFtype=="afni_spm" & ...
    A1A2.synth.sMaj > sMajMIN & ...
    A1A2.synth.sMaj < sMajMAX & ...
    A1A2.afni6.aspect < aspectMAX & ...
    A1A2.synth.eccen > eccenMIN & ...
    A1A2.synth.eccen < eccenMAX,:);
% B=prctile(afnimid.aspect,[5,95]);inR=afnimid.aspect>=B(1) & afnimid.aspect<=B(2);
% afnimid=afnimid(inR,:);
vistalow  = A1A2.vista6(A1A2.noiseLevel=="low" & A1A2.HRFtype=="vista_twogammas" & ...
    A1A2.synth.sMaj > sMajMIN & ...
    A1A2.synth.sMaj < sMajMAX & ...
    A1A2.vista6.aspect < aspectMAX & ...
    A1A2.synth.eccen > eccenMIN & ...
    A1A2.synth.eccen < eccenMAX,:);
% B=prctile(vistalow.aspect,[5,95]);inR=vistalow.aspect>=B(1) & vistalow.aspect<=B(2);
% vistalow=vistalow(inR,:);
vistamid  = A1A2.vista6(A1A2.noiseLevel=="mid" & A1A2.HRFtype=="vista_twogammas" & ...
    A1A2.synth.sMaj > sMajMIN & ...
    A1A2.synth.sMaj < sMajMAX & ...
    A1A2.vista6.aspect < aspectMAX & ...
    A1A2.synth.eccen > eccenMIN & ...
    A1A2.synth.eccen < eccenMAX,:);
% B=prctile(vistamid.aspect,[5,95]);inR=vistamid.aspect>=B(1) & vistamid.aspect<=B(2);
% vistamid=vistamid(inR,:);

fnameRoot = sprintf('Fig4-5_Histograms_Synth_Aspect-2_TR-%i',tr);
disp(fnameRoot)
% saveToType = 'svg';
kk = mrvNewGraphWin(fnameRoot);
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  1  1]);

subplot(2,2,1)
aspect = afnilow.aspect(afnilow.TR==tr);
h = histogram(aspect, histbins,'Normalization','probability'); hold on
set(h,'LineWidth',2,'EdgeColor','k','FaceAlpha',1,'FaceColor','k');hold on
medaspect = median(aspect);
a=plot(medaspect*[1,1],[0,max(h.Values)],'r-','LineWidth',1);
title(sprintf('Afni low noise, TR=%g',tr))
xlabel('Aspect Ratio (GT Aspect = 2)')
xlim([1,5])
set(gca,'FontName', 'Arial','FontSize',16)

subplot(2,2,2)
aspect = afnimid.aspect(afnimid.TR==tr);
h = histogram(aspect,histbins,'Normalization','probability'); hold on
set(h,'LineWidth',2,'EdgeColor','k','FaceAlpha',1,'FaceColor','k');hold on
medaspect = median(aspect);
a=plot(medaspect*[1,1],[0,max(h.Values)],'r-','LineWidth',1);
title(sprintf('Afni mid noise, TR=%g',tr))
xlabel('Aspect Ratio (GT Aspect = 2)')
xlim([1,5])
set(gca,'FontName', 'Arial','FontSize',16)

subplot(2,2,3)
aspect = vistalow.aspect(vistalow.TR==tr);
h = histogram(aspect,histbins,'Normalization','probability'); hold on
set(h,'LineWidth',2,'EdgeColor','k','FaceAlpha',1,'FaceColor','k');hold on
medaspect = median(aspect);
a=plot(medaspect*[1,1],[0,max(h.Values)],'r-','LineWidth',1);
title(sprintf('mrVista low noise, TR=%g',tr))
xlabel('Aspect Ratio (GT Aspect = 2)')
xlim([1,5])
set(gca,'FontName', 'Arial','FontSize',16)

subplot(2,2,4)
aspect = vistamid.aspect(vistamid.TR==tr);
h = histogram(aspect,histbins,'Normalization','probability'); hold on
set(h,'LineWidth',2,'EdgeColor','k','FaceAlpha',1,'FaceColor','k');hold on
medaspect = median(aspect);
a=plot(medaspect*[1,1],[0,max(h.Values)],'r-','LineWidth',1);
title(sprintf('mrVista mid noise, TR=%g',tr))
xlabel('Aspect Ratio (GT Aspect = 2)')
xlim([1,5])
set(gca,'FontName', 'Arial','FontSize',16)


saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.',ext)),ext);

end

