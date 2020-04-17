%% TODO: SAVE ALL IMPORTANT FILES, not in local


% {
clear all; close all; clc
% p = '/Users/glerma/gDrive/STANFORD/PROJECTS/2019_PRF_Validation_methods_(Gari)/__PUBLISH__/ELLIPTICAL';
% f = 'sub-ellipse_ses-sess02-prf_acq-normal_run-01_bold.mat';


saveTo = '~/gDrive/STANFORD/PROJECTS/2019_PRF_Validation_methods_(Gari)/__PUBLISH__/ELLIPTICAL/Figures/RAW';
%}

%% NOTES FOR THE ANALYSIS
% 1. DONE (similar stimulus, see below): What was the stimulus that Baker used? You might want to try your
%    simulations with the same stimulus sequence. 
% 2. DONE: Is there a bias in the
%    distribution of angles of the major axis, or is the distribution pretty flat
%    (as a function of angle)? 
% 3. TODO: If you use an ellipse as ground truth, and add
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
tools  = {'afni4'};
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

%% HRF EFFECT PLOTS
    COMBINE_PARAMETERS                       = struct();
    COMBINE_PARAMETERS.RF.Centerx0           = [3];
    COMBINE_PARAMETERS.RF.Centery0           = [3];  
    COMBINE_PARAMETERS.RF.sigmaMajor         = [2];  
    COMBINE_PARAMETERS.RF.sigmaMinor         = 'same';
    COMBINE_PARAMETERS.TR                    = 1;
    COMBINE_PARAMETERS.Stimulus.durationSecs = 200;

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
        Noise(1).seed                        = 'random';
        Noise(1).voxel                       = 'mid';
    COMBINE_PARAMETERS.Noise                 = Noise;
    synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS, 'repeats',1);
    synthDT = pmForwardModelCalculate(synthDT);
    sDT = synthDT;
    
    % Solve it
    % boyntonresultsvista = pmModelFit(sDT, 'vista','model','one oval gaussian');
    options.afni.model      = 'afni6';
    options.afni.hrf        = 'spm';
    boyntonresultsafni = pmModelFit(sDT, 'afni', 'options', options);
    
    % Create comptTable
    paramDefaults = {'Centerx0','Centery0','Theta','sigmaMinor','sigmaMajor'};
    boyntoncompTable  = pmResultsCompare(sDT, {'aprf'}, {boyntonresultsafni}, ...
        'params', paramDefaults, ...
        'shorten names',true, ...
        'dotSeries', false);
    
    % Plot it
    hh = mrvNewGraphWin('HRF comparison');
    set(hh,'Position',[0.007 0.62  0.8  0.8]);

    Cs  = 0.65 * distinguishable_colors(6,'w');
    
    % Create the fit plots with the ground truth
    tools  = {'aprf'}; nslvl  = 'mid';
    HRFs   = {'boynton','boynton','boynton','boynton'};
    for ii=1:height(boyntoncompTable)
        subplot(2,5,ii)
        useHRF = HRFs{ii};
        ttable = boyntoncompTable(ii,:);
        pmCloudOfResults(ttable   , tools ,'onlyCenters',false ,'userfsize' , 2, ...
            'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-','color',Cs(ii+1,:), ...
            'lineWidth' , 2     ,'noiselevel' ,nslvl , ...
            'useellipse', true, ...
            'xlims',[0, 7],'ylims',[0, 7], 'xtick',[1:1:6],'ytick',[1:1:6], ...
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

 % another plot            
             
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
                          
%%
% [kkvistacompTable, kkvistatSeries, kkvistaresults] = pmNoiseFreeTests('vista6','ellipse',true);             
% [kkafnicompTable, kkafnitSeries, kkafniresults] = pmNoiseFreeTests('afni6','ellipse',true);

%% CircularVsElliptical
if CircularVsElliptical
        COMBINE_PARAMETERS                       = struct();
        COMBINE_PARAMETERS.RF.Centerx0           = [3];
        COMBINE_PARAMETERS.RF.Centery0           = [3];  
        COMBINE_PARAMETERS.RF.sigmaMajor         = [2];  
        COMBINE_PARAMETERS.RF.sigmaMinor         = 'same';
        COMBINE_PARAMETERS.TR                    = 1.5;

        HRF                                      = struct();
        HRF(1).Type                              = 'vista_twogammas'; 
        HRF(1).normalize                         = 'height'; 
        HRF(2).Type                              = 'afni_spm'; 
        HRF(2).normalize                         = 'height'; 

        COMBINE_PARAMETERS.HRF                   = HRF;
            Noise                                = struct();
            Noise(1).voxel                       = 'low';
            Noise(1).seed                        = 'none';
            Noise(1).jitter                      = [0, 0];  % [0.1, 0.1];
        COMBINE_PARAMETERS.Noise                 = Noise;

        % This is the same one as before, but now we want to do the slow stimuli version
        % by Jon's suggestion
        COMBINE_PARAMETERS.Stimulus.durationSecs = 300;

        synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS, 'repeats',10);
        synthDT = pmForwardModelCalculate(synthDT);
        sDT = synthDT;

        %% Solve it with ellipticals

        mrvista_results = pmModelFit(sDT, 'vista','model','onegaussian');
        afni4_results   = pmModelFit(sDT, 'afni_4');
        afni6_results   = pmModelFit(sDT, 'afni_6');
        mrvistaoval_results = pmModelFit(sDT, 'vista','model','oneovalgaussian');
        

        %% Create comptTable
        paramDefaults = {'Centerx0','Centery0','Theta','sigmaMinor','sigmaMajor'};
        compTable  = pmResultsCompare(sDT, {'vista','afni_4','vistaoval','afni_6'}, ...
                                           {mrvista_results, afni4_results, mrvistaoval_results, afni6_results}, ...
            'params', paramDefaults, ...
            'shorten names',true, ...
            'dotSeries', false);

        %% Plot it
        hh = mrvNewGraphWin('circularElliptical');
        set(hh,'Position',[0.007 0.62  0.4  0.8]);
        nrows=2; ncols=2;
        Cs  = 0.65 * distinguishable_colors(6,'w');
        xlims = [-2,8];
        ylims = xlims;
        nlvl  = 'none';

        subplot(nrows,ncols,1)
        pmCloudOfResults(compTable   , {'vista'} ,'onlyCenters',false ,'userfsize' , 2, ...
                     'centerPerc', 90    ,'useHRF'     ,'vista_twogammas' ,'lineStyle' , '-', ...
                     'lineWidth' , .7     ,'noiselevel' ,nlvl , 'addtext',true, ...
                     'color', [0.5,0.5,0.5], 'xlims',xlims,'ylims',ylims,...
                     'xtick',[1:5],'ytick',[1:5], 'addcibar', false, 'useEllipse', false, ...
                     'newWin'    , false ,'saveTo'     ,'','saveToType','svg')
        title('vista circular')

        subplot(nrows,ncols,2)
        pmCloudOfResults(compTable   , {'afni_4'} ,'onlyCenters',false ,'userfsize' , 2, ...
                     'centerPerc', 90    ,'useHRF'     ,'afni_spm' ,'lineStyle' , '-', ...
                     'lineWidth' , .7     ,'noiselevel' ,nlvl , 'addtext',true, ...
                     'color', [0.5,0.5,0.5], 'xlims',xlims,'ylims',ylims,...
                     'xtick',[1:5],'ytick',[1:5], 'addcibar', false, 'useEllipse', false, ...
                     'newWin'    , false ,'saveTo'     ,'','saveToType','svg')
        title('afni circular')

        subplot(nrows,ncols,3)
        pmCloudOfResults(compTable   , {'vistaoval'} ,'onlyCenters',false ,'userfsize' , 2, ...
                     'centerPerc', 90    ,'useHRF'     ,'vista_twogammas' ,'lineStyle' , '-', ...
                     'lineWidth' , .7     ,'noiselevel' ,nlvl , 'addtext',true, ...
                     'color', [0.5,0.5,0.5], 'xlims',xlims,'ylims',ylims,...
                     'xtick',[1:5],'ytick',[1:5], 'addcibar', false, 'useEllipse', true, ...
                     'newWin'    , false ,'saveTo'     ,'','saveToType','svg')
        title('vista elliptical')

        subplot(nrows,ncols,4)
        pmCloudOfResults(compTable   , {'afni_6'} ,'onlyCenters',false ,'userfsize' , 2, ...
                     'centerPerc', 90    ,'useHRF'     ,'afni_spm' ,'lineStyle' , '-', ...
                     'lineWidth' , .7     ,'noiselevel' ,nlvl , 'addtext',true, ...
                     'color', [0.5,0.5,0.5], 'xlims',xlims,'ylims',ylims,...
                     'xtick',[1:5],'ytick',[1:5], 'addcibar', false, 'useEllipse', true, ...
                     'newWin'    , false ,'saveTo'     ,'','saveToType','svg')
        title('afni elliptical')
    
    
    fnameRoot = 'CircularElliptical';
    saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.svg')),'svg');
    
    
end

%% FIG 1: Noiseless plots: accuracy
sub = 'ellipse'; ses = 'noiselesssimplev2';
p = ['/Users/glerma/toolboxes/PRFmodel/local/' sub '/BIDS/derivatives/prfreport/sub-' sub '/ses-' ses];
if ~isdir(p); mkdir(p); end
f = ['sub-' sub '_ses-' ses '-prf_acq-normal_run-01_bold.mat'];
if isfile(fullfile(p,f))
    load(fullfile(p,f))
else
    % Calculate data
    [afnicompTable , afnitSeries , afniresults]  = pmNoiseFreeTests('afni6' , 'ellipse', true);
    [vistacompTable, vistatSeries, vistaresults] = pmNoiseFreeTests('vista6', 'ellipse', true);
    % Save it so that we don't need to generate every time
    compTable        = afnicompTable;
    compTable.vista6 = vistacompTable.vista6;
    save(fullfile(p,f), 'compTable')
end


% RATIO 1
fnameRoot = 'ELLIP_NoiselessCloudPoints4ratios_RATIO1'; ext = 'png';
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
    pmCloudOfResults(afnicompTable, tools ,'onlyCenters',false , ...
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
    useHRF = 'vista_twogammas';
    switch r
        case 0.5,   sMin=0.5; sMaj=0.5; useellipse=true;
        case 1  ,   sMin=1  ; sMaj=1  ; useellipse=true;
        case 2  ,   sMin=2  ; sMaj=2  ; useellipse=true;
        case 3  ,   sMin=3  ; sMaj=3  ; useellipse=true;
        otherwise, error('Ratio %i not contemplated',r)
    end
    pmCloudOfResults(vistacompTable   , tools ,'onlyCenters',false , ...
                     'userfsize'  , sMaj, 'userfsizemin' , sMin, 'useellipse',useellipse, ...
                     'centerPerc' , 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                     'lineWidth'  , 1     ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
                     'centerDistr', false,'synthbluelinewidth',1.5,...
                     'xlims',[0, 6],'ylims',[0, 6], 'xtick',[0,1,2,3,4,5,6],'ytick',[0,1,2,3,4,5,6], ...
                     'newWin'    , false ,'saveTo'     ,'','saveToType','svg')
end
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,['.' ext])),ext);    















% RATIO others
fnameRoot = 'ELLIP_NoiselessCloudPoints4ratios_RATIOrest'; ext = 'png';
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
    useHRF = 'vista_twogammas';
    switch r
        case 1.5, sMin=2; sMaj=3; useellipse=true;
        case 2,   sMin=1; sMaj=2; useellipse=true;
        case 3,   sMin=1; sMaj=3; useellipse=true;
        case 4,   sMin=.5; sMaj=2; useellipse=true;
        otherwise, error('Ratio %i not contemplated',r)
    end
    pmCloudOfResults(vistacompTable   , tools ,'onlyCenters',false , ...
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
    A = afnicompTable;
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

%% FIG 2: Noiseless: eccentricity plots, but noiseless
sub = 'ellipse'; ses = 'noiselessv2';
p = ['/Users/glerma/toolboxes/PRFmodel/local/' sub '/BIDS/derivatives/prfreport/sub-' sub '/ses-' ses];
f = ['sub-' sub '_ses-' ses '-prf_acq-normal_run-01_bold.mat'];
if isfile(fullfile(p,f))
    load(fullfile(p,f))
else
    % Calculate data
    [afnicompTable , afnitSeries , afniresults]  = pmNoiseFreeTests('afni6' ,'eccen',true);
    [vistacompTable, vistatSeries, vistaresults] = pmNoiseFreeTests('vista6','eccen',true);
    % Save it so that we don't need to generate every time
    compTable        = afnicompTable;
    compTable.vista6 = vistacompTable.vista6;
    save(fullfile(p,f), 'compTable')
end



% INDIVIDUAL PLOTS
%{
locs      = unique(vistacompTable.synth.x0);
for nl = 1:length(locs)
    location = locations(nl,:);
    fnameRoot = sprintf('ELLIP_NoiselessCloudPoints4ratios_ECCEN-%1.2f',location(1)); ext = 'png';
    kk = mrvNewGraphWin(fnameRoot);
    % Fig size is relative to the screen used. This is for laptop at 1900x1200
    set(kk,'Position',[0.007 0.62  0.4  0.3]);
    nrows  = 2; ncols = 4;
    ratios = [0.5,1,2,3];

    % Apply params to all
    nslvl  = 'none';
    addcihist = false;
    % plot AFNI
    for nr = 1:length(ratios)
        subplot(nrows,ncols,nr+length(ratios))
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
        pmCloudOfResults(afnicompTable, tools ,'onlyCenters',false , ...
            'userfsize' , sMaj, 'userfsizemin' , sMin, 'useellipse',useellipse, ...
            'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
            'lineWidth' , 1     ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
            'centerDistr', false,'synthbluelinewidth',1.5,...
            'location', location, ...
            'xlims',[location(1)-3.5, location(1)+3.5],...
            'ylims',[location(1)-3.5, location(1)+3.5], ... % 'xtick',[0,1,2,3,4,5,6],'ytick',[0,1,2,3,4,5,6], ...
            'newWin'    , false ,'saveTo'     ,'','saveToType','svg')
    end
    % Plot each tool separately
    % Plot mrVista
    for nr = 1:length(ratios)
        subplot(nrows,ncols,nr)
        r      = ratios(nr);
        tools  = {'vista6'};
        useHRF = 'vista_twogammas';
        switch r
            case 0.5,   sMin=0.5; sMaj=0.5; useellipse=true;
            case 1  ,   sMin=1  ; sMaj=1  ; useellipse=true;
            case 2  ,   sMin=2  ; sMaj=2  ; useellipse=true;
            case 3  ,   sMin=3  ; sMaj=3  ; useellipse=true;
            otherwise, error('Ratio %i not contemplated',r)
        end
        pmCloudOfResults(vistacompTable   , tools ,'onlyCenters',false , ...
                         'userfsize'  , sMaj, 'userfsizemin' , sMin, 'useellipse',useellipse, ...
                         'centerPerc' , 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                         'lineWidth'  , 1     ,'noiselevel' ,nslvl , 'addcihist', addcihist,...
                         'centerDistr', false,'synthbluelinewidth',1.5,...
                         'location', location, ...
                         'xlims',[location(1)-3.5, location(1)+3.5], ...
                         'ylims',[location(1)-3.5, location(1)+3.5], ... % 'xtick',[0,1,2,3,4,5,6],'ytick',[0,1,2,3,4,5,6], ...
                         'newWin'    , false ,'saveTo'     ,'','saveToType','svg')
    end
    saveas(gcf,fullfile(saveTo, strcat(fnameRoot,['.' ext])),ext);
end   
%}



% SUMMARY PLOT: RATIO 1 and 2
% As a summary, plot the eccen vs aspect ratio
fnameBegin = 'NoiselessEccSimRatio1and2';
ext        = 'png';
nlvl       = "none";
eccenInGT  = true;
checksizes = [0.5    ,     1,       2,    3];
ellipsizes = {[1,0.5], [2,1], [3,1.5], [4,2]};
tools      = {'afni6'          , 'vista6'};  % 'vista6' 'afni6' 'vista4' 'afni4'
% for vista is vista_twogammas, but only one value in synth, see B table
useHRFs    = {'afni_spm'       , 'afni_spm' };
duration   = 400;
tr         = 2;
nrow = 2, ncol = 4;
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
        xlabel('Eccentricity')
        ylabel('pRF aspect ratio')
        ylim([0,8]);
        set(gca, 'FontSize', 16)
        legend([a,b],{'Ground truth aspect ratio = 1', ...
                     ['Ground truth aspect ratio = 2 (' ...
                     num2str(ellipsizes{ns}(1)) '/' num2str(ellipsizes{ns}(2)) ')']})
    end
end
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,['.' ext])),ext);    

%% FIGURE 7 equivalent with Vista6 and Afni6
sub = 'ellipse'; ses = 'eccsv2';
p = ['/Users/glerma/toolboxes/PRFmodel/local/' sub '/BIDS/derivatives/prfreport/sub-' sub '/ses-' ses];
f = ['sub-' sub '_ses-' ses '-prf_acq-normal_run-01_bold.mat'];
A = load(fullfile(p,f))
% It seems that AFNI-s theta was not correctly corrected. This has been fixed now. 
% It only affects to results in sub-ellipse/ses-*v2
if strcmp(tool,'afni6');A.compTable.afni6.Th = A.compTable.afni6.Th + deg2rad(90);end
% Add the SNR values (this will come from prfreport in the future)
sub = 'ellipse'; ses = 'eccsv2SNR';
p = ['/Users/glerma/toolboxes/PRFmodel/local/' sub '/BIDS/derivatives/prfsynth/sub-' sub '/ses-' ses];
f = ['sub-' sub '_ses-' ses '_task-prf_acq-normal_run-01_bold.json'];
B = struct2table(jsondecode(fileread(fullfile(p,f))));
A.compTable.SNR = B.SNR;


tools   = {'vista6','afni6'};
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
    fnameRoot = sprintf('R1_CloudPlots4x4_Noise-%s_Size-%i', nslvl{:}, userfsize);
    mm        = mrvNewGraphWin(fnameRoot,[]);  % add off to run it in the server or Docker container
    set(mm,'Units','centimeters','Position',[0 0 10*numanalysis 10*numanalysis]);
    np      = 0;
    for tool = tools; for useHRF = useHRFs
            np=np+1;
            subplot(numanalysis,numanalysis,np)
            pmCloudOfResults(A.compTable   , tool ,'onlyCenters',onlyCenters ,...
                'userfsize' , userfsize, ...
                'centerPerc', centerPerc    ,'useHRF'     ,useHRF{:},...
                'lineStyle' , lineStyle, ...
                'lineWidth' , lineWidth     ,'noiselevel' ,nslvl{:} , ...
                'useellipse', useellipse, 'addsnr',true,...
                'location',location,...
                'addtext',addtext, 'adddice',false,'addsnr',true,...
                'color', color, 'xlims',xlims,'ylims',ylims,'fontsize', fontsize, ...
                'xtick',xtick,'ytick',ytick, 'addcibar', addcibar,'addcihist', addcihist,  ...
                'newWin'    , newWin ,'saveTo'     ,'','saveToType',saveToType)
        end;end
    set(gca,'FontName', 'Arial')
    saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.',saveToType)),saveToType);
end












% RATIO 2
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

%% Silson 2018 plot: ECCEN vs ASPECT
sub = 'ellipse'; ses = 'eccsv2';
p = ['/Users/glerma/toolboxes/PRFmodel/local/' sub '/BIDS/derivatives/prfreport/sub-' sub '/ses-' ses];
f = ['sub-' sub '_ses-' ses '-prf_acq-normal_run-01_bold.mat'];
A = load(fullfile(p,f))
% It seems that AFNI-s theta was not correctly corrected. This has been fixed now. 
% It only affects to results in sub-ellipse/ses-*v2
if strcmp(tool,'afni6');A.compTable.afni6.Th = A.compTable.afni6.Th + deg2rad(90);end
% Add the SNR values (this will come from prfreport in the future)
sub = 'ellipse'; ses = 'eccsv2SNR';
p = ['/Users/glerma/toolboxes/PRFmodel/local/' sub '/BIDS/derivatives/prfsynth/sub-' sub '/ses-' ses];
f = ['sub-' sub '_ses-' ses '_task-prf_acq-normal_run-01_bold.json'];
B = struct2table(jsondecode(fileread(fullfile(p,f))));
A.compTable.SNR = B.SNR;

% SAME HRF; RATIO 1 and 2
fnameBegin = 'SilsonEccSimHRFok';
ext        = 'png';
nlvls      = {"mid","low"};
centerPerc = 50;
eccenInGT  = true;
checksizes = [0.5,1,2,3];
ellipsizes = {[1,0.5],[2,1],[4,2],[6,3]};
tools      = {'vista6'          , 'afni6'};  % 'vista6' 'afni6' 'vista4' 'afni4'
useHRFs    = {'vista_twogammas' ,'afni_spm'};
duration   = 400;
tr         = 2;
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
            dt         = A.compTable;
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
            plot([0,max(eccenvals)],[1,1],'LineWidth',1.5,'LineStyle','--','Color',0.75*[0 1 0])
            plot([0,max(eccenvals)],[2,2],'LineWidth',1.5,'LineStyle','--','Color','c')
            Cs              = 0.65*distinguishable_colors(1+length(eccenvals),'w');

            % Apply percentiles and plot individually
            for ne=1:length(eccenvals)
                C           = Cs(ne,:);
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
                    a = scatter(ecc,aspectmedcirc,80,0.75*[0 1 0],'filled');
                    vax = plot(ecc * [1,1],...
                        [aspectmincirc  , aspectmaxcirc], ...
                        'Color',0.75*[0 1 0],'LineStyle','-','LineWidth',3);
                    
                    b = scatter(ecc+.15,aspectmedellip,80,'c','filled');
                    vax = plot((ecc+.15) * [1,1],...
                        [aspectminellip  , aspectmaxellip], ...
                        'Color','c','LineStyle','-','LineWidth',3);
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
            legend([a,b],{'G.T. Aspect = 1', 'G.T. Aspect = 2'})
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

% THE OTHER HRF; RATIO 1 and 2
fnameBegin = 'SilsonEccSimHRFinv';
ext        = 'png';
nlvls      = {"mid","low"};
centerPerc = 50;
eccenInGT  = true;
checksizes = [0.5,1,2,3];
ellipsizes = {[1,0.5],[2,1],[4,2],[6,3]};
tools      = {'vista6'          , 'afni6'};  % 'vista6' 'afni6' 'vista4' 'afni4'
useHRFs    = {'afni_spm','vista_twogammas'};
duration   = 400;
tr         = 2;
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
            dt         = A.compTable;
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
            plot([0,max(eccenvals)],[1,1],'LineWidth',1.5,'LineStyle','--','Color',0.75*[0 1 0])
            plot([0,max(eccenvals)],[2,2],'LineWidth',1.5,'LineStyle','--','Color','c')
            Cs              = 0.65*distinguishable_colors(1+length(eccenvals),'w');

            % Apply percentiles and plot individually
            for ne=1:length(eccenvals)
                C           = Cs(ne,:);
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
                    a = scatter(ecc,aspectmedcirc,80,0.75*[0 1 0],'filled');
                    vax = plot(ecc * [1,1],...
                        [aspectmincirc  , aspectmaxcirc], ...
                        'Color',0.75*[0 1 0],'LineStyle','-','LineWidth',3);
                    
                    b = scatter(ecc+.15,aspectmedellip,80,'c','filled');
                    vax = plot((ecc+.15) * [1,1],...
                        [aspectminellip  , aspectmaxellip], ...
                        'Color','c','LineStyle','-','LineWidth',3);
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
            legend([a,b],{'G.T. Aspect = 1', 'G.T. Aspect = 2'})
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

%% Silson 2018 plot: SIZE vs ASPECT
sub = 'ellipse'; ses = 'sizesv2';
p = ['/Users/glerma/toolboxes/PRFmodel/local/' sub '/BIDS/derivatives/prfreport/sub-' sub '/ses-' ses];
f = ['sub-' sub '_ses-' ses '-prf_acq-normal_run-01_bold.mat'];
load(fullfile(p,f))



fnameBegin = 'SilsonSizeSim';
ext        = 'png';
nlvls       = {"mid","low"};
centerPerc = 50;
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
    ystart=zeros(size(sizes));
    ystop=6*ones(size(sizes));
    plot([sizes.';sizes.'],[ystart.';ystop.'], ...
        'LineWidth',.7,'LineStyle','-.','Color','k')
    hold on
    plot([0,6],[1,1],'LineWidth',1.5,'LineStyle','--','Color',0.75*[0 1 0])
    Cs              = 0.65*distinguishable_colors(1+length(sizes),'w');
    % Apply percentiles and plot individually
    for ne=1:length(sizes)
        C           = Cs(ne,:);
        ecc         = sizes(ne);
        aspect      = dt.(tool).aspect(dt.synth.sMaj==ecc);
        realeccen   = dt.(tool).eccen(dt.synth.sMaj==ecc);
        B           = prctile(aspect, [twoTailedRange, 100 - twoTailedRange]);
        inRange     = aspect>=B(1) & aspect<=B(2);
        % Apply
        aspectci    = aspect(inRange);
        realeccenci = realeccen(inRange);
        % Medians
        aspectmed   = median(aspectci);
        aspectmin   = min(aspectci);
        aspectmax   = max(aspectci);
        realeccenmed= median(realeccenci);
        realeccenmin= min(realeccenci);
        realeccenmax= max(realeccenci);

        % Plot it

        if eccenInGT
            scatter(ecc,aspectmed,80,C,'filled')
            vax = plot(ecc * [1,1],...
                [aspectmin  , aspectmax], ...
                'Color',C,'LineStyle','-','LineWidth',3); %
        else
            scatter(realeccenmed,aspectmed,60,C,'filled')
            hax = plot([realeccenmin, realeccenmax],...
                aspectmed*[1,1], ...
                'Color',C,'LineStyle','-','LineWidth',2); % 'Color','k',
            vax = plot(realeccenmed * [1,1],...
                [aspectmin  , aspectmax], ...
                'Color',C,'LineStyle','-','LineWidth',2); %
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

%% Silson 2018 plot: SIZEfit vs SIZEsynth
sub = 'ellipse'; ses = 'sizesv2';
p = ['/Users/glerma/toolboxes/PRFmodel/local/' sub '/BIDS/derivatives/prfreport/sub-' sub '/ses-' ses];
f = ['sub-' sub '_ses-' ses '-prf_acq-normal_run-01_bold.mat'];
load(fullfile(p,f))



fnameBegin = 'SilsonSizeAndSize';
ext        = 'png';
nlvls       = {"mid","low"};
centerPerc = 90;
eccenInGT  = true;
tools      = {'vista6'          , 'afni6'};  % 'vista6' 'afni6' 'vista4' 'afni4'
useHRFs    = {'vista_twogammas' , 'afni_spm'};
duration   = 400;
tr         = 2;
for nlvl = nlvls
    nlvl = nlvl{:};
for nt=1:length(tools)
    tool   = tools{nt};
    useHRF = useHRFs{nt};
    dt     = compTable;
    % We want to use just its own HRF, remove the other one
    dt = dt(dt.HRFtype==string(useHRF),:);
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
    ystart=zeros(size(sizes));
    ystop=6*ones(size(sizes));
    plot([sizes.';sizes.'],[ystart.';ystop.'], ...
        'LineWidth',.7,'LineStyle','-.','Color','k');
    hold on
    p = identityLine(gca);
    set(p,'LineWidth',2,'LineStyle','--','Color',0.75*[0 1 0]);
    Cs              = 0.65*distinguishable_colors(1+length(sizes),'w');
    % Apply percentiles and plot individually
    for ne=1:length(sizes)
        C           = Cs(ne,:);
        sz          = sizes(ne);
        smajs       = dt.(tool).sMaj(dt.synth.sMaj==sz);
        smins       = dt.(tool).sMin(dt.synth.sMin==sz);
        % aspect      = dt.(tool).aspect(dt.synth.sMaj==sz);
        % realeccen   = dt.(tool).eccen(dt.synth.sMaj==sz);
        B           = prctile(smajs, [twoTailedRange, 100 - twoTailedRange]);
        inRange     = smajs>=B(1) & smajs<=B(2);
        % Apply
        smajsci      = smajs(inRange);
        sminsci      = smins(inRange);
        % Medians
        smajmed      = median(smajsci);
        smajmin   = min(smajsci);
        smajmax   = max(smajsci);
        sminmed   = median(sminsci);
        sminmin   = min(sminsci);
        sminmax   = max(sminsci);

        
        % Plot it

        if eccenInGT
            % smaj
            scatter(sz,smajmed,80,'b','filled')
            vaxmaj = plot(sz * [1,1],...
                [smajmin  , smajmax], ...
                'Color','b','LineStyle','-','LineWidth',3); %
            % smin
            scatter(sz,sminmed,80,'r','filled')
            vaxmin = plot(sz * [1,1],...
                [sminmin  , sminmax], ...
                'Color','r','LineStyle','-','LineWidth',3); %
        else
            scatter(realeccenmed,aspectmed,60,C,'filled')
            hax = plot([realeccenmin, realeccenmax],...
                aspectmed*[1,1], ...
                'Color',C,'LineStyle','-','LineWidth',2); % 'Color','k',
            vax = plot(realeccenmed * [1,1],...
                [aspectmin  , aspectmax], ...
                'Color',C,'LineStyle','-','LineWidth',2); %
        end
    end
    title(strrep(fnameRoot,'_','\_'))
    xlabel('Radius size GT. ')
    ylabel('Radius size fitted.')
    ylim([0,10]); 
    set(gca, 'FontSize', 16)
    legend([vaxmaj,vaxmin],{'sigMajor','sigMinor'}, 'location','best')
    saveas(gcf,fullfile(saveTo, strcat(fnameRoot,['.' ext])),ext);    
end;end

% Unique eccentricities are xvalues
vAspect = zeros(100,length(xvalues));
aAspect = zeros(100,length(xvalues));

for xx= 1:length(xvalues)
    
    lst = (vista_y(:,1) == xvalues(xx));
    vAspect(:,xx) = vista_y(lst,2);
    
    lst = (afni_y(:,1) == xvalues(xx));
    aAspect(:,xx) = afni_y(lst,2);
end

 mrvNewGraphWin; vHist = histogram(vAspect(:,1)); 
 vHist.BinEdges = [1:0.5:5];
 title('Vista')
 mrvNewGraphWin; aHist = histogram(aAspect(:,1)); 
 aHist.BinEdges = [1:0.5:5]; title('Afni')
 
 median(vAspect)
 median(aAspect)
 xvalues
 
 
 median(afni_y(:,2))
 median(vista_y(:,2))


mrvNewGraphWin;plot(xvalues,median(aAspect),'bo')
set(gca,'ylim',[1,5])

mrvNewGraphWin;plot(xvalues,median(vAspect),'bo')
set(gca,'ylim',[1,5])

%% Polar angles and Thetas
% DATA
sub = 'ellipse'; ses = 'sess04';
p = ['/Users/glerma/toolboxes/PRFmodel/local/' sub '/BIDS/derivatives/prfreport/sub-' sub '/ses-' ses '/v2'];
f = ['sub-' sub '_ses-' ses '-prf_acq-normal_run-01_bold.mat'];
SS = load(fullfile(p,f));
A = SS.compTable;
% Silson 2018 plot: THETA vs ANGLE

% VARIABLES
fnameBegin     = 'SilsonSimTheta';
nlvls          = {"mid","low"};
tools          = {'vista6','afni6'};
centerPerc     = 50;
eccenInGT      = true;
checksizes     = [0.5,1,2,3];
gtaspectratio  = 1;
eccentricities = linspace(1,9,8);
eccentricities = round(eccentricity,2)';

for nlvl = nlvls
    nlvl = nlvl{:};
    for tool = tools
        tool = tool{:};
        for nr=1:checksizes
            rfsize = checksizes(nr); %   end;end;end
            % Use only this noise:
            dt = A(A.noiseLevel==nlvl,:);
            % Use only this size. We want ratio==1, so both the same:
            dt = dt(dt.synth.sMaj==rfsize,:);
            dt = dt(dt.synth.sMin==rfsize,:);
            
            
            % MAKE THIS A FUNCTION 
            % Obtain eccentricity and polar angle
            [TH,R]         = cart2pol(dt.synth.x0, dt.synth.y0);
            dt.synth.angle = round(rad2deg(TH),2);
            dt.synth.eccen = round(R,2);
            dt.synth.aspect= round(dt.synth.sMaj ./ dt.synth.sMin);

            % Do it for the specific tool
            [TH,R]         = cart2pol(dt.(tool).x0, dt.(tool).y0);
            dt.(tool).angle  = rad2deg(TH);
            dt.(tool).eccen  = R;
            dt.(tool).aspect = dt.(tool).sMaj  ./ dt.(tool).sMin;

            % Filter all that we can filter
            % Eccentricity
            % dt = dt(dt.synth.eccen==eccentricity,:);
            
            % Aspect ratio: start with synthesized aspect ratio = 1
            % dt = dt(dt.synth.aspect==gtaspectratio,:);
            
            % Check percentage is 100 based
            if centerPerc < 1; centerPerc = centerPerc*100; end
            % Define the required confidence intervals as two percentiles
            twoTailedRange = (100 - centerPerc) / 2;


            % Create main plot with the ground truth lines
            fnameEnd  = sprintf('%s_TR-2_Dur-400s_Noise-%s_CI-%i_GTsize-%ideg',tool,nlvl,centerPerc,unique(dt.synth.sMaj));
            fnameRoot = strcat(fnameBegin,'-', fnameEnd); 
            disp(fnameRoot)
            kk = mrvNewGraphWin(fnameRoot);
            % Fig size is relative to the screen used. This is for laptop at 1900x1200
            set(kk,'Position',[0.007 0.62  0.4  0.4]);
            centers = unique(dt.synth{:,{'x0','y0'}}, 'rows');
            radii   = unique(dt.synth.sMaj)/2; % Viscircles needs radius and sigma-s are diameters
            % Create color vector
            Cs = 0.65*distinguishable_colors(1+size(centers,1),'w');

            subplot(2,2,2)
            polarhistogram(dt.(tool).Th, 90, 'DisplayStyle','bar')
            thetalim([-90, 90])
            title('Absolute theta values')

            subplot(2,2,4)
            thdeg         = rad2deg(dt.(tool).Th);
            relativetheta = abs(dt.(tool).angle - thdeg);
            polarhistogram(deg2rad(relativetheta), 90, 'DisplayStyle','bar')
            % thetalim([-90, 90])
            title('Relative theta (|prf angle - theta)| ')

            subplot(2,2,[1,3])
            % Plot GT centers
            scatter(centers(:,1), centers(:,2),40,Cs(1,:),'filled')
            hold on
            % Plot GT circles
            viscircles(centers,radii*ones(size(centers,1),1),'LineWidth',2,'LineStyle','--','Color',Cs(1,:));
            axis equal; grid
            xlim([0,5]); ylim([-5,5])
            ax = gca;
            ax.XAxisLocation = 'origin';
            ax.YAxisLocation = 'origin';
            set(ax,'FontSize',16);


            % INTRODUCE THIS PLOT TOO
            [y,x] = ksdensity(dt.(tool).aspect);
            plot(x,y,'k-');
            xlabel('Aspect ratio')



            % Apply percentiles and plot individually
            anglevals = unique(dt.synth.angle);
            for ne=1:length(anglevals)
                C            = Cs(ne+1,:);
                ang          = anglevals(ne);
                ecc          = unique(dt.synth.eccen);
                aspect       = dt.(tool).aspect(dt.synth.angle==ang);
                realeccen    = dt.(tool).eccen(dt.synth.angle==ang);
                realangle    = dt.(tool).angle(dt.synth.angle==ang);
                x0           = dt.(tool).x0(dt.synth.angle==ang);
                y0           = dt.(tool).y0(dt.synth.angle==ang);
                smin         = dt.(tool).sMin(dt.synth.angle==ang);
                smaj         = dt.(tool).sMaj(dt.synth.angle==ang);
                % DECIDE THIS, control CI with aspect?
                B            = prctile(aspect, [twoTailedRange, 100 - twoTailedRange]);
                inRange      = aspect>=B(1) & aspect<=B(2);
                % Apply
                aspectci     = aspect(inRange);
                realeccenci  = realeccen(inRange);
                realangleci  = realangle(inRange);
                realx0ci     = x0(inRange);
                realy0ci     = y0(inRange);
                realsminci   = smin(inRange);
                realsmajci   = smaj(inRange);
                % Medians
                aspectmed    = median(aspectci);
                aspectmin    = min(aspectci);
                aspectmax    = max(aspectci);

                realeccenmed = median(realeccenci);
                realeccenmin = min(realeccenci);
                realeccenmax = max(realeccenci);

                realanglemed = median(realangleci);
                realanglemin = min(realangleci);
                realanglemax = max(realangleci);

                realx0med    = median(realx0ci);
                realx0min    = min(realx0ci);
                realx0max    = max(realx0ci);

                realy0med    = median(realy0ci);
                realy0min    = min(realy0ci);
                realy0max    = max(realy0ci);

                reasmin0med  = median(realsminci);
                realsminmin  = min(realsminci);
                realsminmax  = max(realsminci);

                realsmajmed  = median(realsmajci);
                realsmajmin  = min(realsmajci);
                realsmajmax  = max(realsmajci);
                % Plot it

                if eccenInGT
                    scatter(realx0med, realy0med, 80,C,'filled')
                    for ne = 1:length(aspectci)
                        h = drawellipse(realx0ci(ne),realy0ci(ne),realangleci(ne),realsmajci(ne)/2,realsminci(ne)/2);
                        set(h,'LineWidth',.7,'LineStyle','-','Color',[.5 .5 .5]);
                        hold on
                    end
                    % vax = plot([realx0min  , realx0max],...
                    %            [realy0min  , realy0max], ...
                    %             'Color',C,'LineStyle','-','LineWidth',2); %
                    fitEllipse(realx0ci, realy0ci, [0,0,0],centerPerc)
                else
                    scatter(realeccenmed,aspectmed,60,C,'filled')
                    hax = plot([realeccenmin, realeccenmax],...
                        aspectmed*[1,1], ...
                        'Color',C,'LineStyle','-','LineWidth',2); % 'Color','k',
                    vax = plot(realeccenmed * [1,1],...
                        [aspectmin  , aspectmax], ...
                        'Color',C,'LineStyle','-','LineWidth',2); %
                end
            end
            title(strrep(fnameRoot,'_','\_'))
            xlabel('x deg ')
            ylabel('y deg')
            set(gca, 'FontSize', 16)

end;end;end
