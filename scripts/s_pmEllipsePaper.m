% {
clear all; close all; clc
p = '/Users/glerma/gDrive/STANFORD/PROJECTS/2019_PRF_Validation_methods_(Gari)/__PUBLISH__/ELLIPTICAL';
f = 'sub-ellipse_ses-sess02-prf_acq-normal_run-01_bold.mat';
load(fullfile(p,f))

saveTo = '~/gDrive/STANFORD/PROJECTS/2019_PRF_Validation_methods_(Gari)/__PUBLISH__/ELLIPTICAL/Figures/RAW';
%}

%% NOTES FOR THE ANALYSIS
% 1. What was the stimulus that Baker used? You might want to try your
%    simulations with the same stimulus sequence. 
% 2. Is there a bias in the
%    distribution of angles of the major axis, or is the distribution pretty flat
%    (as a function of angle)? 
% 3. If you use an ellipse as ground truth, and add
%    low or moderate noise. how much elongation do you need in order to reliably
%    recover the correct angle? For example, if the ratio of major axis to minor
%    axis is 1.1, can you correctly recover the alpha? What about 1.2? 2? This
%    would be one way of quantifying the usefulness of the elliptical model--if it
%    can only accurately characterize pRF which are implausibly elongated, then
%       it's not a very useful model.

% Silson 2018 J Neuro uses
%    - Bar sweeps: 38 sec
%        This means that the stim duration needs to be:
%        pm.Stimulus.durationSecs = 400 secs
%    - TR = 2 sec
%    - Eccentricities:
%        vals = linspace(1,9,8);
%        usevals = sqrt((vals.^2)/2)
%        usevals = [0.7071,1.5152,2.3234,3.1315,3.9396,4.7477,5.5558,6.3640]; 
% Use as a input
    



kk = mrvNewGraphWin('NoiselessCloudPoints','wide');
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  0.2  0.3]);
% subplot(1,4,1)
tools  = {'afni'};
useHRF = 'afni_spm';
nslvl  = 'mid';
pmCloudOfResults(compTable   , tools ,'onlyCenters', false ,'userfsize' , 2, 'centerdistr',false,...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'useellipse',true, 'lineWidth' , .7     ,'noiselevel' ,nslvl , 'addtext',true, ...
                 'color', [0.5,0.5,0.5], 'xlims',[-3, 6],'ylims',[-3,6],...
                 'addcihist', false, 'xtick',[-2:6],'ytick',[-2:6],...
                 'location', [3,0], ... % 'all', ... % [3,3], ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')
             
             
             
%% NOTES ON THETA

% STANDARDIZE
% In our tool, Theta is going to be the anticlockwise angle of sigmaMajor, from
% X, like in tha mathematical definition:
%  https://en.wikipedia.org/wiki/Polar_coordinate_system#/media/File:Examples_of_Polar_Coordinates.svg
% The important thing here is that it's the angle of sigmaMajor, or sigma_y in AFNI
% (see below)


% AFNI
% This is how we read it:
% pmEstimates.Theta = results(:,6);
% This is the explanation of Theta by AFNI
%                   Given stimulus images over time s(x,y,t), find x0, y0, sigma, R and
%                   theta values that produce a best fit of the model to the
%                   data.  Here x0, y0 are taken to be the center of the
%                   population receptive field, sigma is the minor width of it
%                   (sigma_x, below), sigrat R is the ratio (sigma_y / sigma_x),
%                   and theta is the rotation from the y-direction major axis
%                   (so zero is in the positive y-direction).
% 
%                   We assume sigma_y >= sigma_x and refer to sigrat >= 1, since that
%                   sufficiently represents all possibilities.  The reciprocol would
%                   come from the negative complimentary angle, and would therefore be a
%                   redundant solution.
% 
%                    parameter domains:
%                      x,y        : [-1,1], scaled by the mask, itself
%                      sigma      : (0,1], where 1 means the mask radius
%                      R (sigrat) : [1,inf), since sigma defines the smaller size
%                      theta      : [-PI/2, PI/2), since rotation by PI has no effect

%% Create the fig for sup mat asked by jon
A = load('/Users/glerma/gDrive/STANFORD/PROJECTS/2019_PRF_Validation_methods_(Gari)/__PUBLISH__/PAPER_SUBMISSION01/Figures/RAW/sub-paper_ses-sess01-prf_acq-normal_run-01_bold.mat');
kk = mrvNewGraphWin('NoiselessCloudPoints3sizes','wide');
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  0.8  0.9]);

% Apply params to all
nslvl  = 'low';
addcihist = true;

% Plot individuals
subplot(3,4,1)
tools  = {'vista'};
useHRF = 'vista_twogammas';
pmCloudOfResults(A.compTable   , tools ,'onlyCenters',false ,'userfsize' , 0.25, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 2     ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
                ... 'xlims',[2.75, 3.25],'ylims',[2.75, 3.25], 'xtick',[2.5:.125:3.5],'ytick',[2.5:.125:3.5], ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')

subplot(3,4,2)
tools  = {'afni'};
useHRF = 'afni_spm';
% nslvl  = 'none';
pmCloudOfResults(A.compTable   , tools ,'onlyCenters',false ,'userfsize' , 0.25, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 2     ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
                 ... 'xlims',[2.75, 3.25],'ylims',[2.75, 3.25], 'xtick',[2.5:.125:3.5],'ytick',[2.5:.125:3.5], ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')

subplot(3,4,3)
tools  = {'popeye'};
useHRF = 'popeye_twogammas';
% nslvl  = 'none';
pmCloudOfResults(A.compTable   , tools ,'onlyCenters',false ,'userfsize' , 0.25, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 2     ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
                 ... 'xlims',[2.75, 3.25],'ylims',[2.75, 3.25], 'xtick',[2.5:.125:3.5],'ytick',[2.5:.125:3.5], ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')

subplot(3,4,4)
tools  = {'aprf'};
useHRF = 'canonical';
% nslvl  = 'none';
pmCloudOfResults(A.compTable   , tools ,'onlyCenters',false ,'userfsize' , 0.25, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 1.5   ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
                 ... 'xlims',[2.75, 3.25],'ylims',[2.75, 3.25], 'xtick',[2.5:.125:3.5],'ytick',[2.5:.125:3.5], ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')

subplot(3,4,5)
tools  = {'vista'};
useHRF = 'vista_twogammas';
% nslvl  = 'none';
pmCloudOfResults(A.compTable   , tools ,'onlyCenters',false ,'userfsize' , 1, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 2     ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
                ...  'xlims',[2, 4],'ylims',[2, 4], 'xtick',[2:.5:4],'ytick',[2:.5:4], ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')

subplot(3,4,6)
tools  = {'afni'};
useHRF = 'afni_spm';
% nslvl  = 'none';
pmCloudOfResults(A.compTable   , tools ,'onlyCenters',false ,'userfsize' , 1, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 2     ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
                ...  'xlims',[2, 4],'ylims',[2, 4], 'xtick',[2:.5:4],'ytick',[2:.5:4], ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')

subplot(3,4,7)
tools  = {'popeye'};
useHRF = 'popeye_twogammas';
% nslvl  = 'none';
pmCloudOfResults(A.compTable   , tools ,'onlyCenters',false ,'userfsize' , 1, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 2     ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
                ...  'xlims',[2, 4],'ylims',[2, 4], 'xtick',[2:.5:4],'ytick',[2:.5:4], ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')

subplot(3,4,8)
tools  = {'aprf'};
useHRF = 'canonical';
% nslvl  = 'none';
pmCloudOfResults(A.compTable   , tools ,'onlyCenters',false ,'userfsize' , 1, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 1.5   ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
                ...  'xlims',[2, 4],'ylims',[2, 4], 'xtick',[2:.5:4],'ytick',[2:.5:4], ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')
             
             
subplot(3,4,9)
tools  = {'vista'};
useHRF = 'vista_twogammas';
% nslvl  = 'none';
pmCloudOfResults(A.compTable   , tools ,'onlyCenters',false ,'userfsize' , 2, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 2     ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')

subplot(3,4,10)
tools  = {'afni'};
useHRF = 'afni_spm';
% nslvl  = 'none';
pmCloudOfResults(A.compTable   , tools ,'onlyCenters',false ,'userfsize' , 2, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 2     ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')

subplot(3,4,11)
tools  = {'popeye'};
useHRF = 'popeye_twogammas';
% nslvl  = 'none';
pmCloudOfResults(A.compTable   , tools ,'onlyCenters',false ,'userfsize' , 2, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 2     ,'noiselevel' ,nslvl ,'addcihist', addcihist, ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')

subplot(3,4,12)
tools  = {'aprf'};
useHRF = 'canonical';
% nslvl  = 'none';
pmCloudOfResults(A.compTable   , tools ,'onlyCenters',false ,'userfsize' , 2, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 1.5   ,'noiselevel' ,nslvl ,'addcihist', addcihist, ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')             
             
             
             
fnameRoot = 'Noisefree_accuracy_3sizes_lownoise';
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.svg')),'svg');



%% Noiseless plots: accuracy
% Calculate data first:
[afnicompTable, afnitSeries, afniresults] = pmNoiseFreeTests('afni6','ellipse',true);
[vistacompTable, vistatSeries, vistaresults] = pmNoiseFreeTests('vista6','ellipse',true);



kk = mrvNewGraphWin('ELLIP_NoiselessCloudPoints4ratios','wide');
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  0.4  0.3]);
nrows  = 2; ncols = 4;
ratios = [1,1.5,2,3];

% Apply params to all
nslvl  = 'none';
addcihist = true;

% Plot each tool separately
% Plot mrVista
for nr = 1:length(ratios)
    subplot(nrows,ncols,nr)
    r      = ratios(nr);
    tools  = {'vista6'};
    useHRF = 'vista_twogammas';
    switch r
        case 1,   sMin=2; sMaj=2; useellipse=true;
        case 1.5, sMin=2; sMaj=3; useellipse=true;
        case 2,   sMin=1; sMaj=2; useellipse=true;
        case 3,   sMin=1; sMaj=3; useellipse=true;
        otherwise, error('Ratio %i not contemplated',r)
    end
    pmCloudOfResults(vistacompTable   , tools ,'onlyCenters',false , ...
                     'userfsize' , sMaj, 'userfsizemin' , sMin, 'useellipse',useellipse, ...
                     'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                     'lineWidth' , 2     ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
                    ... 'xlims',[2.75, 3.25],'ylims',[2.75, 3.25], 'xtick',[2.5:.125:3.5],'ytick',[2.5:.125:3.5], ...
                     'newWin'    , false ,'saveTo'     ,'','saveToType','svg')
end
% plot afni
for nr = 1:length(ratios)
    subplot(nrows,ncols,nr+length(ratios))
    r      = ratios(nr);
    tools  = {'afni6'};
    useHRF = 'afni_spm';
    switch r
        case 1,   sMin=2; sMaj=2; useellipse=true;
        case 1.5, sMin=2; sMaj=3; useellipse=true;
        case 2,   sMin=1; sMaj=2; useellipse=true;
        case 3,   sMin=1; sMaj=3; useellipse=true;
        otherwise, error('Ratio %i not contemplated',r)
    end
    pmCloudOfResults(afnicompTable, tools ,'onlyCenters',false , ...
        'userfsize' , sMaj, 'userfsizemin' , sMin, 'useellipse',useellipse, ...
        'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
        'lineWidth' , 2     ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
        ... 'xlims',[2.75, 3.25],'ylims',[2.75, 3.25], 'xtick',[2.5:.125:3.5],'ytick',[2.5:.125:3.5], ...
        'newWin'    , false ,'saveTo'     ,'','saveToType','svg')
end

fnameRoot = 'ELLIP_NoiselessCloudPoints4ratios';
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.svg')),'svg');             
             
             
%% HRF EFFECT PLOTS
    COMBINE_PARAMETERS                       = struct();
    COMBINE_PARAMETERS.RF.Centerx0           = [3];
    COMBINE_PARAMETERS.RF.Centery0           = [3];  
    COMBINE_PARAMETERS.RF.sigmaMajor         = [2];  
    COMBINE_PARAMETERS.RF.sigmaMinor         = 'same';
    COMBINE_PARAMETERS.TR                    = 1;

    HRF                                      = struct();
    HRF(1).Type                              = 'boynton';  
    HRF(1).normalize                         = 'height'; 
    HRF(1).params.n = 3;
    HRF(1).params.tau = 1.08;
    HRF(1).params.delay = 2.05;
        
    HRF(2).Type                              = 'boynton';
    HRF(2).normalize                         = 'height'; 
    HRF(2).params.n = 3;
    HRF(2).params.tau = 1.38;
    HRF(2).params.delay = 2;
    
    HRF(3).Type                              = 'boynton';
    HRF(3).normalize                         = 'height'; 
    HRF(3).params.n = 3;
    HRF(3).params.tau = 1.68;
    HRF(3).params.delay = 1.75;
    
    HRF(4).Type                              = 'boynton';
    HRF(4).normalize                         = 'height'; 
    HRF(4).params.n = 3;
    HRF(4).params.tau = 1.935;
    HRF(4).params.delay = 1.65;
    
    COMBINE_PARAMETERS.HRF                   = HRF;
        Noise                                = struct();
        Noise(1).seed                        = 'none';
    COMBINE_PARAMETERS.Noise                 = Noise;
    synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS, 'repeats',1);
    synthDT = pmForwardModelCalculate(synthDT);
    sDT = synthDT;
    
    %% Solve it
    boyntonresultsvista = pmModelFit(sDT, 'vista','model','one oval gaussian');
    % boyntonresultsafni = pmModelFit(sDT, 'afni6');
    
    %% Create comptTable
    paramDefaults = {'Centerx0','Centery0','Theta','sigmaMinor','sigmaMajor'};
    boyntoncompTable  = pmResultsCompare(sDT, {'aprf'}, {boyntonresultsvista}, ...
        'params', paramDefaults, ...
        'shorten names',true, ...
        'dotSeries', false);
    
    %% Plot it
    hh = mrvNewGraphWin('HRF comparison');
    set(hh,'Position',[0.007 0.62  0.8  0.8]);

    Cs  = 0.65 * distinguishable_colors(6,'w');
    
    % Create the fit plots with the ground truth
    tools  = {'aprf'}; nslvl  = 'none';
    HRFs   = {'boynton','boynton','boynton','boynton'};
    for ii=1:height(boyntoncompTable)
        subplot(2,5,ii)
        useHRF = HRFs{ii};
        ttable = boyntoncompTable(ii,:);
        pmCloudOfResults(ttable   , tools ,'onlyCenters',false ,'userfsize' , 2, ...
            'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-','color',Cs(ii+1,:), ...
            'lineWidth' , 2     ,'noiselevel' ,nslvl , ...
            'useellipse', true, ...
            'newWin'    , false ,'saveTo'     ,'','saveToType','svg')
    end
    
    subplot(2,5,[6:10])
    a   = [];
    leg = [];
    for ii = 1:height(sDT)
        thispm = sDT.pm(ii);
        if ii~=5
            line='-';
            % thisleg = {sprintf('Boynton %i (width=%1d)',ii,thispm.HRF.width)};
            thisleg = {sprintf('Boynton %i',ii)};
        else
            % thisleg = {sprintf('aPRF canonical (width=%1d)',thispm.HRF.width)};
            thisleg = {sprintf('aPRF canonical')};
            line='-.';
        end
        a = [a;thispm.HRF.plot('window',false,'dots',false,'addwidth',false,...
            'xlims',[0,20],'color',Cs(ii+1,:),'line',line)];
        leg = [leg;thisleg];
    end
    legend(a,{'Boynton 1','Boynton 2','Boynton 3','Boynton 4','aprf\_canonical'})
    legend(a,leg)
    title('Boynton HRFs modulated in width and canonical aprf')
    xticks([0:20])
    
    fnameRoot = 'ELLIP_HRF_and_width_VISTA_';
    % fnameRoot = 'ELLIP_HRF_and_width_AFNI_';
    saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.svg')),'svg');

             
 %% another plot            
             
kk = mrvNewGraphWin('NoiselessCloudPoints','wide');
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  0.2  0.3]);
% subplot(1,4,1)
tools  = {'afni'};
useHRF = 'afni_spm';
nslvl  = 'mid';
pmCloudOfResults(compTable   , tools ,'onlyCenters', false ,'userfsize' , 2, 'centerdistr',false,...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'useellipse',true, 'lineWidth' , .7     ,'noiselevel' ,nslvl , 'addtext',true, ...
                 'color', [0.5,0.5,0.5], 'xlims',[-3, 6],'ylims',[-3,6],...
                 'addcihist', false, 'xtick',[-2:6],'ytick',[-2:6],...
                 'location', [3,0], ... % 'all', ... % [3,3], ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')





















