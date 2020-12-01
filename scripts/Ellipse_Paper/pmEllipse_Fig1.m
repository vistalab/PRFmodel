function s01_Ellipse_Fig1
% Create Figure 1
%
% See also
%  s00_MainFiguresScript
%

%%
ext  = 'png'; % Could be svg
saveTo = fullfile(pmRootPath,'local','figures');  % Folder path
if ~exist(saveTo,'dir'), mkdir(saveTo); end


%% Load the noiseless plots: accuracy

% The file with the fits produced from the full calculation.
sub = 'ellipse'; ses = 'noiselesssimplev2';
p = fullfile(pmRootPath,'local',sub,'BIDS','derivatives','prfreport',['sub-' sub],['ses-' ses]);
if ~isfolder(p); mkdir(p); end
f = ['sub-' sub '_ses-' ses '-prf_acq-normal_run-01_bold.mat'];

theFitFile = fullfile(p,f);
if isfile(theFitFile)
    % If the file already exists, just load it
    load(theFitFile,'compTable');
else
    disp('Calculating ground-truth, noise-free, data.')
    afnicompTable  = pmNoiseFreeTests('afni6' , 'ellipse', true, 'usenifti',true);
    vistacompTable = pmNoiseFreeTests('vista6', 'ellipse', true, 'usenifti',true);
    
    disp('Saveing ground-truth, noise-free calculations.')
    compTable        = afnicompTable;
    compTable.vista6 = vistacompTable.vista6;
    if ~exist(p,'dir'); mkdir(p); end
    save(theFitFile, 'compTable')
end


%% RATIO 1
fnameRoot = 'Fig1-A_ELLIP_NoiselessCloudPoints4ratios_RATIO1'; % ext = 'svg';
kk = mrvNewGraphWin(fnameRoot);
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  0.4  0.3]);
nrows  = 2; ncols = 4;
ratios = [0.5,1,2,3];

% Apply params to all
nslvl  = 'none';
addcihist = false;

% Plot each tool separately

% plot AFNI
for nr = 1:length(ratios)
    subplot(nrows,ncols,nr)
    r      = ratios(nr);
    tools  = {'afni6'};
    useHRF = 'afni_spm';
    switch r
        case 0.5,   sMin=0.5; sMaj=0.5; useellipse=true;
        case 1  ,   sMin=1  ; sMaj=1  ; useellipse=true;
        case 2  ,   sMin=2  ; sMaj=2  ; useellipse=true;
        case 3  ,   sMin=3  ; sMaj=3  ; useellipse=true;
        otherwise, error('Ratio %i not contemplated',r)
    end
    
    pmCloudOfResults(compTable, tools ,'onlyCenters',false , ...
        'userfsize' , sMaj, 'userfsizemin' , sMin, 'useellipse',useellipse, ...
        'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
        'lineWidth' , 1     ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
        'centerDistr', false,'synthbluelinewidth',1.5,...
        'xlims',[0, 6],'ylims',[0, 6], 'xtick',[0,1,2,3,4,5,6],'ytick',[0,1,2,3,4,5,6], ...
        'newWin'    , false ,'saveTo'     ,'','saveToType','svg')
end
% Plot mrVista
for nr = 1:length(ratios)
    subplot(nrows,ncols,nr+length(ratios))
    r      = ratios(nr);
    tools  = {'vista6'};
    useHRF = 'afni_spm';  % it is really vista_twogammas, but not in the table
    switch r
        case 0.5,   sMin=0.5; sMaj=0.5; useellipse=true;
        case 1  ,   sMin=1  ; sMaj=1  ; useellipse=true;
        case 2  ,   sMin=2  ; sMaj=2  ; useellipse=true;
        case 3  ,   sMin=3  ; sMaj=3  ; useellipse=true;
        otherwise, error('Ratio %i not contemplated',r)
    end
    
    pmCloudOfResults(compTable   , tools ,'onlyCenters',false , ...
        'userfsize'  , sMaj, 'userfsizemin' , sMin, 'useellipse',useellipse, ...
        'centerPerc' , 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
        'lineWidth'  , 1     ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
        'centerDistr', false,'synthbluelinewidth',1.5,...
        'xlims',[0, 6],'ylims',[0, 6], 'xtick',[0,1,2,3,4,5,6],'ytick',[0,1,2,3,4,5,6], ...
        'newWin'    , false ,'saveTo'     ,'','saveToType','svg')
end
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,['.' ext])),ext);

% {
% RATIO others
fnameRoot = 'Fig1-B_ELLIP_NoiselessCloudPoints4ratios_RATIOrest'; % ext = 'svg';
kk = mrvNewGraphWin(fnameRoot);
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  0.4  0.3]);
nrows  = 2; ncols = 4;
ratios = [1.5,2,3,4];

% Apply params to all
nslvl  = 'none';
addcihist = false;

% Plot each tool separately
% Plot mrVista
for nr = 1:length(ratios)
    subplot(nrows,ncols,nr+length(ratios))
    r      = ratios(nr);
    tools  = {'vista6'};
    useHRF = 'afni_spm';  % see above
    switch r
        case 1.5, sMin=2; sMaj=3; useellipse=true;
        case 2,   sMin=1; sMaj=2; useellipse=true;
        case 3,   sMin=1; sMaj=3; useellipse=true;
        case 4,   sMin=.5; sMaj=2; useellipse=true;
        otherwise, error('Ratio %i not contemplated',r)
    end
    pmCloudOfResults(compTable   , tools ,'onlyCenters',false , ...
        'userfsize'  , sMaj, 'userfsizemin' , sMin, 'useellipse',useellipse, ...
        'centerPerc' , 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
        'lineWidth'  , 1     ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
        'centerDistr', false,'synthbluelinewidth',1.5,...
        'xlims',[0, 6],'ylims',[0, 6], 'xtick',[0,1,2,3,4,5,6],'ytick',[0,1,2,3,4,5,6], ...
        'newWin'    , false ,'saveTo'     ,'','saveToType','svg')
end
% plot afni
for nr = 1:length(ratios)
    subplot(nrows,ncols,nr)
    r      = ratios(nr);
    tools  = {'afni6'};
    useHRF = 'afni_spm';
    switch r
        case 1.5, sMin=2; sMaj=3; useellipse=true;
        case 2,   sMin=1; sMaj=2; useellipse=true;
        case 3,   sMin=1; sMaj=3; useellipse=true;
        case 4,   sMin=.5; sMaj=2; useellipse=true;
        otherwise, error('Ratio %i not contemplated',r)
    end
    A = compTable;
    A.afni6.Th = A.afni6.Th + deg2rad(90);
    pmCloudOfResults(A, tools ,'onlyCenters',false , ...
        'userfsize' , sMaj, 'userfsizemin' , sMin, 'useellipse',useellipse, ...
        'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
        'lineWidth' , 1     ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
        'centerDistr', false,'synthbluelinewidth',1.5,...
        'xlims',[0, 6],'ylims',[0, 6], 'xtick',[0,1,2,3,4,5,6],'ytick',[0,1,2,3,4,5,6], ...
        'newWin'    , false ,'saveTo'     ,'','saveToType','svg')
end

saveas(gcf,fullfile(saveTo, strcat(fnameRoot,['.' ext])),ext);

end

