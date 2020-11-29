function s03_Ellipse_Fig3(saveTo, saveToType)
if ~isfolder(saveTo); mkdir(saveTo); end
% saveTo = '~/gDrive/STANFORD/PROJECTS/2019_PRF_Validation_methods_(Gari)/__PUBLISH__/ELLIPTICAL/Figures/RAW';

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


% RATIO 1
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
% saveToType  = 'svg';

numanalysis = length(tools);

useHRFs = {};
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
    mm        = mrvNewGraphWin(fnameRoot,[]);  % add off to run it in the server or Docker container
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
        'newWin'    , newWin ,'saveTo'     ,'','saveToType',saveToType)
    % Hist
    subplot(numanalysis,numanalysis,2)
    tt = A.compTable.afni6(A.compTable.noiseLevel==string(nslvl{:}) & ...
                           A.compTable.HRFtype==string(useHRF) & ...
                           A.compTable.synth.sMaj==2 & ...
                           A.compTable.synth.sMin==2 & ...
                           A.compTable.synth.x0==3.1315 & ...
                           A.compTable.synth.y0==3.1315,:);
    aspect    = tt.sMaj  ./ tt.sMin;
    h = histogram(aspect,15);hold on
    medaspect = median(aspect);
    plot(medaspect*[1,1],[0,18],'r-','LineWidth',1)
    set(h,'LineWidth',2,'EdgeColor','k','FaceAlpha',1,'FaceColor','k');hold on
    title('AFNI Elliptical')
    xlabel('Aspect Ratio')    
    set(gca,'FontName', 'Arial','FontSize',16)
    
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
        'newWin'    , newWin ,'saveTo'     ,'','saveToType',saveToType)
    % Hist
    subplot(numanalysis,numanalysis,4)
    tt = A.compTable.vista6(A.compTable.noiseLevel==string(nslvl{:}) & ...
                            A.compTable.HRFtype==string(useHRF) & ...
                            A.compTable.synth.sMaj==2 & ...
                            A.compTable.synth.sMin==2 & ...
                            A.compTable.synth.x0==3.1315 & ...
                            A.compTable.synth.y0==3.1315,:);
    aspect    = tt.sMaj  ./ tt.sMin;
    h = histogram(aspect,15); hold on
    set(h,'LineWidth',2,'EdgeColor','k','FaceAlpha',1,'FaceColor','k');hold on
    medaspect = median(aspect);
    plot(medaspect*[1,1],[0,18],'r-','LineWidth',1)
    title('mrVista Elliptical')
    xlabel('Aspect Ratio')
    set(gca,'FontName', 'Arial','FontSize',16)
    
    saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.',saveToType)),saveToType);
end












% RATIO 2
%{
% Generic values coming from the config.json
onlyCenters = false;
userfsize   = 2;
userfsizemin= 1;
location    = [3.1315,3.1315];  % [3,3]; % 
useHRF      = {};
centerPerc  = 90;
lineStyle   = '-';
lineWidth   = 0.7;
fontsize    = 14;
noiselevel  = {'low','mid'};
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
saveToType  = 'png';
numanalysis = length(tools);
useHRFs     = {};



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
    fnameRoot = sprintf('R2_CloudPlots4x4_Noise-%s_Size-%i', nslvl{:}, userfsize);
    mm        = mrvNewGraphWin(fnameRoot,[]);  % add off to run it in the server or Docker container
    set(mm,'Units','centimeters','Position',[0 0 10*numanalysis 10*numanalysis]);
    np      = 0;
    for tool = tools 
        for useHRF = useHRFs
            np=np+1;
            subplot(numanalysis,numanalysis,np)
            pmCloudOfResults(A.compTable   , tool ,'onlyCenters',onlyCenters ,...
                'userfsize' , userfsize, 'userfsizemin',userfsizemin,...
                'centerPerc', centerPerc    ,'useHRF'     ,useHRF{:},...
                'lineStyle' , lineStyle, ...
                'lineWidth' , lineWidth     ,'noiselevel' ,nslvl{:} , ...
                'useellipse', useellipse, 'addsnr',true,...
                'location'  , location,...
                'addtext'   , addtext, 'adddice',false,'addsnr',true,...
                'color'     , color, 'xlims',xlims,'ylims',ylims,'fontsize', fontsize, ...
                'xtick'     , xtick,'ytick',ytick, 'addcibar', addcibar,'addcihist', addcihist,  ...
                'newWin'    , newWin ,'saveTo'     ,'','saveToType',saveToType)
        end;end
    set(gca,'FontName', 'Arial')
    saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.',saveToType)),saveToType);
end
%}

end

