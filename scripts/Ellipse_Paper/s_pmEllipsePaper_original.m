%% TODO: SAVE ALL IMPORTANT FILES, not in local


% {
clear all; close all; clc
% p = '/Users/glerma/gDrive/STANFORD/PROJECTS/2019_PRF_Validation_methods_(Gari)/__PUBLISH__/ELLIPTICAL';
% f = 'sub-ellipse_ses-sess02-prf_acq-normal_run-01_bold.mat';
tbUse prfmodel

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
clear all; close all; clc
saveTo = '~/gDrive/STANFORD/PROJECTS/2019_PRF_Validation_methods_(Gari)/__PUBLISH__/ELLIPTICAL/Figures/RAW';

sub = 'ellipse'; ses = 'noiselesssimplev2';
p = ['/Users/glerma/toolboxes/PRFmodel/local/' sub '/BIDS/derivatives/prfreport/sub-' sub '/ses-' ses];
if ~isdir(p); mkdir(p); end
f = ['sub-' sub '_ses-' ses '-prf_acq-normal_run-01_bold.mat'];
if isfile(fullfile(p,f))
    load(fullfile(p,f));
else
    % Calculate data
    [afnicompTable , afnitSeries , afniresults]  = pmNoiseFreeTests('afni6' , 'ellipse', true, 'usenifti',true);
    [vistacompTable, vistatSeries, vistaresults] = pmNoiseFreeTests('vista6', 'ellipse', true, 'usenifti',true);
    % Save it so that we don't need to generate every time
    compTable        = afnicompTable;
    compTable.vista6 = vistacompTable.vista6;
    save(fullfile(p,f), 'compTable')
end


% RATIO 1
fnameRoot = 'ELLIP_NoiselessCloudPoints4ratios_RATIO1'; ext = 'svg';
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
fnameRoot = 'ELLIP_NoiselessCloudPoints4ratios_RATIOrest'; ext = 'svg';
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
%}

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
locations = [locs,locs];
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
ext        = 'svg';
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

%% FIG3: PLOS Comp Bio Figure 7 equivalent with Vista6 and Afni6, and histograms
clear all; close all; clc
saveTo = '~/gDrive/STANFORD/PROJECTS/2019_PRF_Validation_methods_(Gari)/__PUBLISH__/ELLIPTICAL/Figures/RAW';

sub = 'ellipse'; ses = 'eccsv2';
p = ['/Users/glerma/toolboxes/PRFmodel/local/' sub '/BIDS/derivatives/prfreport/sub-' sub '/ses-' ses];
f = ['sub-' sub '_ses-' ses '-prf_acq-normal_run-01_bold.mat'];
A = load(fullfile(p,f))
% It seems that AFNI-s theta was not correctly corrected. This has been fixed now. 
% It only affects to results in sub-ellipse/ses-*v2
% if strcmp(tool,'afni6');A.compTable.afni6.Th = A.compTable.afni6.Th + deg2rad(90);end
% Add the SNR values (this will come from prfreport in the future)
sub = 'ellipse'; ses = 'eccsv2SNR';
p = ['/Users/glerma/toolboxes/PRFmodel/local/' sub '/BIDS/derivatives/prfsynth/sub-' sub '/ses-' ses];
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
saveToType  = 'svg';

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

%% Silson 2018 plot: TR=2 ECCEN vs ASPECT
clear all; close all; clc
saveTo = '~/gDrive/STANFORD/PROJECTS/2019_PRF_Validation_methods_(Gari)/__PUBLISH__/ELLIPTICAL/Figures/RAW';

sub = 'ellipse'; ses = 'eccsv2';
p = ['/Users/glerma/toolboxes/PRFmodel/local/' sub '/BIDS/derivatives/prfreport/sub-' sub '/ses-' ses];
f = ['sub-' sub '_ses-' ses '-prf_acq-normal_run-01_bold.mat'];
A2 = load(fullfile(p,f))
% It seems that AFNI-s theta was not correctly corrected. This has been fixed now. 
% It only affects to results in sub-ellipse/ses-*v2
% if strcmp(tool,'afni6');A.compTable.afni6.Th = A.compTable.afni6.Th + deg2rad(90);end
% Add the SNR values (this will come from prfreport in the future)
sub = 'ellipse'; ses = 'eccsv2SNR';
p = ['/Users/glerma/toolboxes/PRFmodel/local/' sub '/BIDS/derivatives/prfsynth/sub-' sub '/ses-' ses];
f = ['sub-' sub '_ses-' ses '_task-prf_acq-normal_run-01_bold.json'];
B = struct2table(jsondecode(fileread(fullfile(p,f))));
A2.compTable.SNR = B.SNR;

% SAME HRF; RATIO 1 and 2
fnameBegin = 'SilsonEccSimHRFok';
ext        = 'svg';
nlvls      = {"mid","low"};
centerPerc = 50;
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











% THE OTHER HRF; RATIO 1 and 2
%{
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
%}

%% Silson 2018 plot: TR=1 ECCEN vs ASPECT
% clear all; close all; clc
saveTo = '~/gDrive/STANFORD/PROJECTS/2019_PRF_Validation_methods_(Gari)/__PUBLISH__/ELLIPTICAL/Figures/RAW';


sub = 'ellipse'; ses = 'eccsv2TR1';
p = ['/Users/glerma/toolboxes/PRFmodel/local/' sub '/BIDS/derivatives/prfreport/sub-' sub '/ses-' ses];
f = ['sub-' sub '_ses-' ses '-prf_acq-normal_run-01_bold.mat'];
A1 = load(fullfile(p,f))

% SAME HRF; RATIO 1 and 2
fnameBegin = 'SilsonEccSimHRFokTR1';
ext        = 'svg';
nlvls      = {"mid","low"};
centerPerc = 50;
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

% THE OTHER HRF; RATIO 1 and 2
%{
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
%}

%% SAVE A1 and A2 for CUMSUM plot at the end
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

%% Silson 2018 plot: TR=2 SIZE vs ASPECT
clear all; close all; clc
saveTo = '~/gDrive/STANFORD/PROJECTS/2019_PRF_Validation_methods_(Gari)/__PUBLISH__/ELLIPTICAL/Figures/RAW';

sub = 'ellipse'; ses = 'sizesv2';
p = ['/Users/glerma/toolboxes/PRFmodel/local/' sub '/BIDS/derivatives/prfreport/sub-' sub '/ses-' ses];
f = ['sub-' sub '_ses-' ses '-prf_acq-normal_run-01_bold.mat'];
load(fullfile(p,f))
% Add the SNR values (this will come from prfreport in the future)
sub = 'ellipse'; ses = 'sizesv2SNR';
p = ['/Users/glerma/toolboxes/PRFmodel/local/' sub '/BIDS/derivatives/prfsynth/sub-' sub '/ses-' ses];
f = ['sub-' sub '_ses-' ses '_task-prf_acq-normal_run-01_bold.json'];
B = struct2table(jsondecode(fileread(fullfile(p,f))));
compTable.SNR = B.SNR;
B2.compTable = compTable;



fnameBegin = 'SilsonSizeSim';
ext        = 'svg';
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

%% Silson 2018 plot: TR=1 SIZE vs ASPECT
% clear all; close all; clc
saveTo = '~/gDrive/STANFORD/PROJECTS/2019_PRF_Validation_methods_(Gari)/__PUBLISH__/ELLIPTICAL/Figures/RAW';

sub = 'ellipse'; ses = 'sizesv2TR1';
p = ['/Users/glerma/toolboxes/PRFmodel/local/' sub '/BIDS/derivatives/prfreport/sub-' sub '/ses-' ses];
f = ['sub-' sub '_ses-' ses '-prf_acq-normal_run-01_bold.mat'];
load(fullfile(p,f))
% Add the SNR values (this will come from prfreport in the future)
sub = 'ellipse'; ses = 'sizesv2SNR';
p = ['/Users/glerma/toolboxes/PRFmodel/local/' sub '/BIDS/derivatives/prfsynth/sub-' sub '/ses-' ses];
f = ['sub-' sub '_ses-' ses '_task-prf_acq-normal_run-01_bold.json'];
B = struct2table(jsondecode(fileread(fullfile(p,f))));
compTable.SNR = B.SNR;
B1.compTable = compTable;


fnameBegin = 'SilsonSizeSimTR1';
ext        = 'svg';
nlvls       = {"mid","low"};
centerPerc = 50;
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
clear all; close all; clc
saveTo = '~/gDrive/STANFORD/PROJECTS/2019_PRF_Validation_methods_(Gari)/__PUBLISH__/ELLIPTICAL/Figures/RAW';

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

fnameRoot = sprintf('Histograms_Synth_Aspect-1_TR-%i',tr);
disp(fnameRoot)
saveToType = 'svg';
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


saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.',saveToType)),saveToType);

%% Combined histograms: aspect2 (run twice, tr=1, tr=2)
clear all; close all; clc
saveTo = '~/gDrive/STANFORD/PROJECTS/2019_PRF_Validation_methods_(Gari)/__PUBLISH__/ELLIPTICAL/Figures/RAW';

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

fnameRoot = sprintf('Histograms_Synth_Aspect-2_TR-%i',tr);
disp(fnameRoot)
saveToType = 'svg';
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


saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.',saveToType)),saveToType);

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
% We are not doing this, if it is not elliptical, then do not use it. 


% DATA
proj = 'ellipse'; sub = 'ellipse'; ses = 'thetasv2';
p    = ['/Users/glerma/toolboxes/PRFmodel/local/' ...
         sub '/BIDS/derivatives/prfreport/sub-' sub '/ses-' ses];
f    = ['sub-' sub '_ses-' ses '-prf_acq-normal_run-01_bold.mat'];
SS   = load(fullfile(p,f));
A    = SS.compTable;
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

%% Noisy circles
% NoiseCircleTest (Simulation:  Take a circle and make two noisy estimates of
% its diameter.  (D1 + noise)/(D2 + noise).  This ratio will always be centered
% on 1.  But now, use the same data to estimate
% (max(D1+noise,D2+noise))/min((D1+noise,D2+noise).  This will always be > 1.
% How much greater?  If your estimate is < (1.2?  1.3?) then the data are
% consistent with a circle.

% Noisy Circle Tests.m
% We have some results in mrVista with simulated data that are not perfect but
% that can be explained by the very nature of the calculation. 

noiseFactor  = .7; % Std. Dev.
radiuses     = [0.25:0.25:6]';  % degrees
% freeRatiosMn = zeros(size(radiuses)); 
% freeRatiosSt = zeros(size(radiuses)); 
limRatiosMed  = zeros(size(radiuses)); 
limRatiosMin  = zeros(size(radiuses)); 
limRatiosMax  = zeros(size(radiuses)); 
for nr=1:length(radiuses)
    radius = radiuses(nr);
    % Generate the random two radiuses
    rng(44444,'twister')
    n1    = noiseFactor * randn(1000,1);
    R1    = radius + n1;
    R1ind = R1 > 0.01;
    
    rng(54321,'twister')
    n2    = noiseFactor * randn(1000,1);
    R2    = radius + n2;
    R2ind = R2 > 0.01;
    
    % Combine all positive radius (AND)
    Rind = R1ind & R2ind;
    R1   = R1(Rind);
    R2   = R2(Rind);

    % Calculate the ratio that is always positive
    limRatio         = max([R1,R2],[],2) ./ min([R1,R2],[],2);
    
    % Calculate median, min and max
    limRatiosMed(nr)  = median(limRatio);
    limRatiosMin(nr)  = min(limRatio);
    limRatiosMax(nr)  = max(limRatio);
end
    
ext       = 'svg';
fnameRoot = 'RadiusSims';
disp(fnameRoot)
kk = mrvNewGraphWin(fnameRoot);
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  0.4  0.4]);

ystart =     ones(size(radiuses));
ystop  = 6 * ones(size(radiuses));
plot([radiuses.';radiuses.'],[ystart.';ystop.'], ...
    'LineWidth',.7,'LineStyle','-.','Color','k')
hold on
% plot([0,6],[1,1],'LineWidth',1.5,'LineStyle','--','Color','k') % 0.75*[0 1 0])
% Cs              = 0.65*distinguishable_colors(1+length(sizes),'w');
% Apply percentiles and plot individually
for ne=1:length(radiuses)
    rad          = radiuses(ne);
    limRatioMed  = limRatiosMed(ne);
    limRatioMin  = limRatiosMin(ne);
    limRatioMax  = limRatiosMax(ne);
        
    % Plot it
    b   = scatter(rad,limRatioMed,80,'ko','filled');
    vax = plot(rad * [1,1],[limRatioMin, limRatioMax], ...
             'Color','k','LineStyle',':','LineWidth',3);    
end

% legend([a,b], {'Free Ratios','Limited Ratios'})
title(strrep(fnameRoot,'_','\_'))
xlabel('Radius')
ylabel('Aspect Ratio')
ylim([0,8]);
set(gca, 'FontSize', 16)
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,['.' ext])),ext);

%% READ: Silson 2018 plot: Real Data 7T

% Filter data:
% Size: both radius at least 1 deg
% Eccen: bigger than 2deg, less than 6

% TASKS:
% Download new data
% 3 subject mixed
% left and right mixed
% - plot size vs aspect as well
% - make range a little bit biiger to find errors predicted by synthetic
% - spearate V1v-d from V2v-d
% plot histograms of aspect ratios
% think about the limitations we add to the data?
% V1-V2 should be around 1-2 deg, check kendrick kay and noah
% filter by R2, bigger than 45%
% check radiality: correlation plot between angle and Theta, if it is organized
% around zero, there is a relationship (the difference is zero)

clear all; close all;
saveTo = '~/gDrive/STANFORD/PROJECTS/2019_PRF_Validation_methods_(Gari)/__PUBLISH__/ELLIPTICAL/Figures/RAW';

proj   = 'realdata';
% tools  = {'afni6','afni4','vista6','vista4'};
tools  = {'vista6','vista4'};
subs   = {'115017','164131','536647'}; 
ses    = '01';
run    = '01';
% Read labels from FS
setenv('SUBJECTS_DIR','/Applications/freesurfer/subjects')
hemis  = {'lh','rh'};
% The mask labels where copied to the fsaverage folder, but they are in
% prfmodel/scripts/realDataAnalysis/
labels = {'V1_exvivo','V2_exvivo', ...
          'V1_exvivo.thresh','V2_exvivo.thresh', 'dorsalV1V2V3Mask'};
kklab = struct();
for nh=[1,2]; for nl=1:length(labels)
    tmp = read_label('fsaverage',[hemis{nh} '.' labels{nl}]);
    kklab.(hemis{nh}).(strrep(labels{nl},'.','_')) = tmp(:,1) + 1;  % fs is a zero based index
    
end;end
lhben   = MRIread(fullfile(pmRootPath,'scripts/readDataAnalysis/fsaverageLabels','lh.benson14_varea.v4_0.mgz'));
rhben   = MRIread(fullfile(pmRootPath,'scripts/readDataAnalysis/fsaverageLabels','rh.benson14_varea.v4_0.mgz'));
ltmpmgh = lhben;
rtmpmgh = rhben;
kklab.lh.bensonV1V2V3 = lhben.vol;
kklab.rh.bensonV1V2V3 = rhben.vol;
% Now we need to obtain the final labels, separating dorsal and ventral as well
kklab.lh.V1  = find(kklab.lh.bensonV1V2V3==1);
kklab.lh.V2  = find(kklab.lh.bensonV1V2V3==2);
kklab.lh.V3  = find(kklab.lh.bensonV1V2V3==3);
kklab.rh.V1  = find(kklab.rh.bensonV1V2V3==1);
kklab.rh.V2  = find(kklab.rh.bensonV1V2V3==2);
kklab.rh.V3  = find(kklab.rh.bensonV1V2V3==3);
% Separated dorsal
il1=ismember(kklab.lh.V1,kklab.lh.dorsalV1V2V3Mask); kklab.lh.V1d = kklab.lh.V1(il1);
il2=ismember(kklab.lh.V2,kklab.lh.dorsalV1V2V3Mask); kklab.lh.V2d = kklab.lh.V2(il2);
il3=ismember(kklab.lh.V3,kklab.lh.dorsalV1V2V3Mask); kklab.lh.V3d = kklab.lh.V3(il3);
ir1=ismember(kklab.rh.V1,kklab.lh.dorsalV1V2V3Mask); kklab.rh.V1d = kklab.rh.V1(ir1);
ir2=ismember(kklab.rh.V2,kklab.lh.dorsalV1V2V3Mask); kklab.rh.V2d = kklab.rh.V2(ir2);
ir3=ismember(kklab.rh.V3,kklab.lh.dorsalV1V2V3Mask); kklab.rh.V3d = kklab.rh.V3(ir3);
% Separated ventral
kklab.lh.V1v = kklab.lh.V1(~il1);
kklab.lh.V2v = kklab.lh.V2(~il2);
kklab.lh.V3v = kklab.lh.V3(~il3);
kklab.rh.V1v = kklab.rh.V1(~ir1);
kklab.rh.V2v = kklab.rh.V2(~ir2);
kklab.rh.V3v = kklab.rh.V3(~ir3);

% Read the data, but we need to get one with all the data together, and another
% one with regions separated. 

% TODO: WORST CODE EVER BELOW
% JUST ADD TEH WHOLE fsaverage and then filter it...    

compTable   = struct();
bylabel     = struct();
bylabelsums = struct();
for nt=1:length(tools)
    tool  = tools{nt};
    disp(tool)   
    subestimates   = table();
    bylabel.(tool) = struct();
    % Do the links here to have life easier afterwards
    bylabelsums.(tool).V1    =  table();
    bylabelsums.(tool).V2    =  table();
    bylabelsums.(tool).V3    =  table();
    % lh-rh-Dorsal
    bylabelsums.(tool).V1d   =  table();
    bylabelsums.(tool).V2d   =  table();
    bylabelsums.(tool).V3d   =  table();
    % lh-rh-Ventral
    bylabelsums.(tool).V1v   =  table();
    bylabelsums.(tool).V2v   =  table();
    bylabelsums.(tool).V3v   =  table();
    for ns=1:length(subs)
        sub   = subs{ns};
        disp(sub)
        p     = ['/Users/glerma/toolboxes/PRFmodel/local/' proj '/BIDS/derivatives/prfanalyze-' tool '/sub-' sub '/ses-' ses];
        bylabel.(tool).(['s' sub]) = struct();
        % Copy the code from prfReportWrapper to create the tables
        %% Read the results back to Matlab
        % Read all the nifti-s for this type of tool 
        cd(p)
        niftis = dir(['*run-' run '*.nii.gz']);
        if length(niftis) < 6
            error('Not enough nifti result files in output, at least there should be x0,y0,sigmamajor,sigmaminor,theta and modelpred')
        end
        filesRead = 0;

        pmEstimates = table();
        bylabel.(tool).(['s' sub]).lhV1_exvivo =  table();
        bylabel.(tool).(['s' sub]).lhV2_exvivo =  table();
        bylabel.(tool).(['s' sub]).rhV1_exvivo =  table();
        bylabel.(tool).(['s' sub]).rhV2_exvivo =  table();
        
        % Benson V1-V2-V3
        bylabel.(tool).(['s' sub]).lhV1  =  table();
        bylabel.(tool).(['s' sub]).lhV2  =  table();
        bylabel.(tool).(['s' sub]).lhV3  =  table();
        bylabel.(tool).(['s' sub]).rhV1  =  table();
        bylabel.(tool).(['s' sub]).rhV2  =  table();
        bylabel.(tool).(['s' sub]).rhV3  =  table();
        % Dorsal
        bylabel.(tool).(['s' sub]).lhV1d =  table();
        bylabel.(tool).(['s' sub]).lhV2d =  table();
        bylabel.(tool).(['s' sub]).lhV3d =  table();
        bylabel.(tool).(['s' sub]).rhV1d =  table();
        bylabel.(tool).(['s' sub]).rhV2d =  table();
        bylabel.(tool).(['s' sub]).rhV3d =  table();
        % Ventral
        bylabel.(tool).(['s' sub]).lhV1v =  table();
        bylabel.(tool).(['s' sub]).lhV2v =  table();
        bylabel.(tool).(['s' sub]).lhV3v =  table();
        bylabel.(tool).(['s' sub]).rhV1v =  table();
        bylabel.(tool).(['s' sub]).rhV2v =  table();
        bylabel.(tool).(['s' sub]).rhV3v =  table();
        % Do the links here to have life easier afterwards
        bylabel.(tool).(['s' sub]).V1    =  table();
        bylabel.(tool).(['s' sub]).V2    =  table();
        bylabel.(tool).(['s' sub]).V3    =  table();
        % lh-rh-Dorsal
        bylabel.(tool).(['s' sub]).V1d   =  table();
        bylabel.(tool).(['s' sub]).V2d   =  table();
        bylabel.(tool).(['s' sub]).V3d   =  table();
        % lh-rh-Ventral
        bylabel.(tool).(['s' sub]).V1v   =  table();
        bylabel.(tool).(['s' sub]).V2v   =  table();
        bylabel.(tool).(['s' sub]).V3v   =  table();
        
        for nn=1:length(niftis)
            fname = niftis(nn).name
            % Assume always it will be .nii.gz and that the . has not been used in the filename
            trunkname = split(fname,'.');
            if length(trunkname) ~= 3; error('%s should be a filename with no dots and .nii.gz',fname);end	
            trunkname = trunkname{1};
            resname   = split(trunkname, '_');
            resname   = resname{end};
            % read the nifti and squeeze the result matrix
            fprintf('Attempting to read %s\n',fullfile(pwd,fname))
            tmp       = niftiRead(fname);
            data      = squeeze(tmp.data);	
            % Data can be separated in different labels
            lhtmpfs   = zeros(1,163842);
            rhtmpfs   = zeros(1,163842);
            % Assign our results to a real brain
            lV1l = length(kklab.lh.V1_exvivo);
            lV2l = length(kklab.lh.V2_exvivo);
            rV1l = length(kklab.rh.V1_exvivo);
            rV2l = length(kklab.rh.V2_exvivo);
            lhtmpfs(kklab.lh.V1_exvivo) = data(1:lV1l,1)';
            lhtmpfs(kklab.lh.V2_exvivo) = data(lV1l+1 : lV1l+lV2l,1)';
            rhtmpfs(kklab.rh.V1_exvivo) = data(lV1l+lV2l+1:lV1l+lV2l+rV1l,1)';
            rhtmpfs(kklab.rh.V2_exvivo) = data(lV1l+lV2l+rV1l+1:lV1l+lV2l+rV1l+rV2l,1)';
            % Now is better if we use Noah's separation, this is in kklab.lh.bensonV1V2V3
            %{
            lV1fnamemgz = fullfile(p,[trunkname  '_lh_V1.mgz']);
            lV2fnamemgz = fullfile(p,[trunkname  '_lh_V2.mgz']);
            rV1fnamemgz = fullfile(p,[trunkname  '_rh_V1.mgz']);
            rV2fnamemgz = fullfile(p,[trunkname  '_rh_V2.mgz']);
            lV1file=ltmpmgh;lV1file.vol=lhtmpfs;MRIwrite(lV1file,lV1fnamemgz);
            lV2file=ltmpmgh;lV2file.vol=lhtmpfs;MRIwrite(lV2file,lV2fnamemgz);
            rV1file=rtmpmgh;rV1file.vol=rhtmpfs;MRIwrite(rV1file,rV1fnamemgz);
            rV2file=rtmpmgh;rV2file.vol=rhtmpfs;MRIwrite(rV2file,rV2fnamemgz);
            %}
            % asign it to the table
            switch resname
                case 'sigmaminor', resname = 'sMin';
                case 'sigmamajor', resname = 'sMaj';
                case 'theta'     , resname = 'Th';
                case 'centerx0'  , resname = 'x0';
                case 'centery0'  , resname = 'y0';
            end
            pmEstimates.(resname) = data;	        
            
            % Now assign it to a structute with all values separated
            bylabel.(tool).(['s' sub]).lhV1_exvivo.(resname) = data(1:lV1l,1);
            bylabel.(tool).(['s' sub]).lhV2_exvivo.(resname) = data(lV1l+1 : lV1l+lV2l,1);
            bylabel.(tool).(['s' sub]).rhV1_exvivo.(resname) = data(lV1l+lV2l+1:lV1l+lV2l+rV1l,1);
            bylabel.(tool).(['s' sub]).rhV2_exvivo.(resname) = data(lV1l+lV2l+rV1l+1:lV1l+lV2l+rV1l+rV2l,1);
            % Benson V1-V2-V3
            bylabel.(tool).(['s' sub]).lhV1.(resname) = lhtmpfs(kklab.lh.V1)';
            bylabel.(tool).(['s' sub]).lhV2.(resname) = lhtmpfs(kklab.lh.V2)';
            bylabel.(tool).(['s' sub]).lhV3.(resname) = lhtmpfs(kklab.lh.V3)';
            bylabel.(tool).(['s' sub]).rhV1.(resname) = rhtmpfs(kklab.rh.V1)';
            bylabel.(tool).(['s' sub]).rhV2.(resname) = rhtmpfs(kklab.rh.V2)';
            bylabel.(tool).(['s' sub]).rhV3.(resname) = rhtmpfs(kklab.rh.V3)';
            % Dorsal
            bylabel.(tool).(['s' sub]).lhV1d.(resname) = lhtmpfs(kklab.lh.V1d)';
            bylabel.(tool).(['s' sub]).lhV2d.(resname) = lhtmpfs(kklab.lh.V2d)';
            bylabel.(tool).(['s' sub]).lhV3d.(resname) = lhtmpfs(kklab.lh.V3d)';
            bylabel.(tool).(['s' sub]).rhV1d.(resname) = rhtmpfs(kklab.rh.V1d)';
            bylabel.(tool).(['s' sub]).rhV2d.(resname) = rhtmpfs(kklab.rh.V2d)';
            bylabel.(tool).(['s' sub]).rhV3d.(resname) = rhtmpfs(kklab.rh.V3d)';
            % Ventral
            bylabel.(tool).(['s' sub]).lhV1v.(resname) = lhtmpfs(kklab.lh.V1v)';
            bylabel.(tool).(['s' sub]).lhV2v.(resname) = lhtmpfs(kklab.lh.V2v)';
            bylabel.(tool).(['s' sub]).lhV3v.(resname) = lhtmpfs(kklab.lh.V3v)';
            bylabel.(tool).(['s' sub]).rhV1v.(resname) = rhtmpfs(kklab.rh.V1v)';
            bylabel.(tool).(['s' sub]).rhV2v.(resname) = rhtmpfs(kklab.rh.V2v)';
            bylabel.(tool).(['s' sub]).rhV3v.(resname) = rhtmpfs(kklab.rh.V3v)';
            % Do the links here to have life easier afterwards
            bylabel.(tool).(['s' sub]).V1.(resname) = [lhtmpfs(kklab.lh.V1),rhtmpfs(kklab.rh.V1)]';
            bylabel.(tool).(['s' sub]).V2.(resname) = [lhtmpfs(kklab.lh.V2),rhtmpfs(kklab.rh.V2)]';
            bylabel.(tool).(['s' sub]).V3.(resname) = [lhtmpfs(kklab.lh.V3),rhtmpfs(kklab.rh.V3)]';
            % lh-rh-Dorsal
            bylabel.(tool).(['s' sub]).V1d.(resname) = [lhtmpfs(kklab.lh.V1d),rhtmpfs(kklab.rh.V1d)]';
            bylabel.(tool).(['s' sub]).V2d.(resname) = [lhtmpfs(kklab.lh.V2d),rhtmpfs(kklab.rh.V2d)]';
            bylabel.(tool).(['s' sub]).V3d.(resname) = [lhtmpfs(kklab.lh.V3d),rhtmpfs(kklab.rh.V3d)]';
            % lh-rh-Ventral
            bylabel.(tool).(['s' sub]).V1v.(resname) = [lhtmpfs(kklab.lh.V1v),rhtmpfs(kklab.rh.V1v)]';
            bylabel.(tool).(['s' sub]).V2v.(resname) = [lhtmpfs(kklab.lh.V2v),rhtmpfs(kklab.rh.V2v)]';
            bylabel.(tool).(['s' sub]).V3v.(resname) = [lhtmpfs(kklab.lh.V3v),rhtmpfs(kklab.rh.V3v)]';
        end
        % Concatenate different subjects
        subestimates = [subestimates;pmEstimates];
        
        bylabelsums.(tool).V1    =  [bylabelsums.(tool).V1; bylabel.(tool).(['s' sub]).V1];
        bylabelsums.(tool).V2    =  [bylabelsums.(tool).V2; bylabel.(tool).(['s' sub]).V2];
        bylabelsums.(tool).V3    =  [bylabelsums.(tool).V3; bylabel.(tool).(['s' sub]).V3];
        % lh-rh-Dorsal
        bylabelsums.(tool).V1d   =  [bylabelsums.(tool).V1d; bylabel.(tool).(['s' sub]).V1d];
        bylabelsums.(tool).V2d   =  [bylabelsums.(tool).V2d; bylabel.(tool).(['s' sub]).V2d];
        bylabelsums.(tool).V3d   =  [bylabelsums.(tool).V3d; bylabel.(tool).(['s' sub]).V3d];
        % lh-rh-Ventral
        bylabelsums.(tool).V1v   =  [bylabelsums.(tool).V1v; bylabel.(tool).(['s' sub]).V1v];
        bylabelsums.(tool).V2v   =  [bylabelsums.(tool).V2v; bylabel.(tool).(['s' sub]).V2v];
        bylabelsums.(tool).V3v   =  [bylabelsums.(tool).V3v; bylabel.(tool).(['s' sub]).V3v];
    end
    compTable.(tool) = subestimates;
    compTable.(tool).R2 = calccod(compTable.(tool).modelpred', compTable.(tool).testdata')';
end

nonfilteredbuylabelsums = bylabelsums;

%% PLOT 6 (A,B,C,D): Silson 2018 plot: Real Data 7T, ONLY VISTA
% restart every time
bylabelsums = nonfilteredbuylabelsums;

doSave     = true;
ext        = 'svg';
centerPerc = 90;
eccenInGT  = true;
xlims      = [0,10];
ylims      = [0,10];
tools      = {'vista6'}; 
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
sMajMIN    = 1 ;%.5; % 1;
sMinMIN    = 1 ;%.5; % .75;
sMajMAX    = 3% 1.5% 3; % 4;
eccenMIN   = 2;%1% 2;
eccenMAX   = 6% 6;
minR2      = 0.25;
% How many bins
NeccenBins = 6;
NareaBins  = NeccenBins;
%close all
tools  = {'vista4','vista6'};
subs   = {'115017','164131','536647'};

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
sub = 'ellipse'; ses = 'tr1dur300v2';
p = ['/Users/glerma/toolboxes/PRFmodel/local/' sub '/BIDS/derivatives/prfreport/sub-' sub '/ses-' ses];
f = ['sub-' sub '_ses-' ses '-prf_acq-normal_run-01_bold.mat'];
tools = {'synth','vista4','vista6'};

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
%{
A1A2 = A1A2(A1A2.vista6.sMaj >= 0.5 & ...
             A1A2.vista6.sMaj < sMajMAX & ...
             A1A2.vista6.eccen > eccenMIN & ...
             A1A2.vista6.eccen < eccenMAX & ...
         A1A2.HRFtype=="vista_twogammas", :);
%}
% {
A1A2 = A1A2(A1A2.vista6.sMaj >= sMajMIN  & ...
             A1A2.vista6.sMaj <= sMajMAX & ...
             A1A2.vista6.eccen >= eccenMIN & ...
             A1A2.vista6.eccen <= eccenMAX , :); % & ...
         % A1A2.HRFtype=="vista_twogammas", :);
%}
%{
A1A2 = A1A2(A1A2.HRFtype=="vista_twogammas", :);
%}
%{
A1A2 = A1A2(A1A2.HRFtype=="afni_spm", :);
%}
unique(A1A2.synth.sMaj)
unique(A1A2.synth.sMin)
unique(A1A2.synth.eccen)
unique(A1A2.HRFtype)
unique(A1A2.noiseLevel)

aspect1  = A1A2.vista6.aspect(A1A2.noiseLevel=="low");
B1=prctile(aspect1, [5, 95]);inRange1 = aspect1 >= B1(1) & aspect1 <= B1(2);
aspect1  = aspect1(inRange1);
sprintf('Low noise: Min aspect ratio for vista 6 is %g and max is %g', min(aspect1),max(aspect1))


aspect2  = A1A2.vista6.aspect(A1A2.noiseLevel=="mid");
B2=prctile(aspect2, [5, 95]);inRange2 = aspect2 >= B2(1) & aspect2 <= B2(2);
aspect2  = aspect2(inRange2);
sprintf('Mid noise: Min aspect ratio for vista 6 is %g and max is %g', min(aspect2),max(aspect2))

mediansyntheticdatalow = median(aspect1);
mediansyntheticdatamid = median(aspect2);

% Discretize by label, to bin the eccentricities
A1A2.vista6.Y = zeros(size(A1A2.vista6.aspect));
A1A2.vista6.Y(A1A2.noiseLevel=="mid") = discretize(A1A2.vista6.eccen(A1A2.noiseLevel=="mid"),eccenvalues); 
A1A2.vista6.Y(A1A2.noiseLevel=="low") = discretize(A1A2.vista6.eccen(A1A2.noiseLevel=="low"),eccenvalues); 
A1A2low = A1A2(A1A2.noiseLevel=="low", :);
A1A2mid = A1A2(A1A2.noiseLevel=="mid", :);
% Apply percentiles and plot individually
% Create the vectors and then plot all together
vistaMedLowEcc = zeros(1,length(eccenvalues)-1);
vistaMedMidEcc = zeros(1,length(eccenvalues)-1);
for ne=1:(length(eccenvalues)-1)
    aspecclow = A1A2low.vista6.aspect(A1A2low.vista6.Y==ne);
    aspeccmid = A1A2mid.vista6.aspect(A1A2mid.vista6.Y==ne);
    % Median
    vistaMedLowEcc(ne) = median(aspecclow);
    vistaMedMidEcc(ne) = median(aspeccmid);
end

% prepare data
tool = 'vista6';
for nl  = 1:length(useLabels)
    lab = useLabels{nl};
    % Discretize by label, to bin the eccentricities
    bylabelsums.(tool).(lab).Y = discretize(bylabelsums.(tool).(lab).eccen,eccenvalues); 
    % Apply percentiles and plot individually
    % Create the vectors and then plot all together
    aspectmedecc.(lab) = zeros(1,length(eccenvalues)-1);
    aspectminecc.(lab) = zeros(1,length(eccenvalues)-1);
    aspectmaxecc.(lab) = zeros(1,length(eccenvalues)-1);
    
    for ne=1:(length(eccenvalues)-1)
        % ECC - ASPECT
        aspecc = bylabelsums.(tool).(lab).aspect(bylabelsums.(tool).(lab).Y==ne);
        % Median and std
        if isempty(aspecc)
            aspectmedecc.(lab)(ne) = 0;
            aspectminecc.(lab)(ne) = 0;
            aspectmaxecc.(lab)(ne) = 0;
        else
            aspectmedecc.(lab)(ne) = median(aspecc);
            aspectminecc.(lab)(ne) = min(aspecc);  % They use SEM, check
            aspectmaxecc.(lab)(ne) = max(aspecc);
        end
    end
end









% PLOT 6A
fnameBegin = 'RealData_SilsonEcc&Size';
% Create main plot with the ground truth lines
fnameEnd = sprintf('TR-%i_Dur-%is_C.I.-%i',tr,duration,centerPerc);
fnameRoot = strcat(fnameBegin,'-', fnameEnd);
disp(fnameRoot)
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
ylim([1,2])
xticks(Emidpoints)
set(gca, 'FontSize', 16) 


% SAVE 6A
if doSave;saveas(gcf,fullfile(saveTo, strcat(fnameRoot,['.' ext])),ext);  end



% PLOT 6B
aspects = [];
fnameBegin = 'RealData_AspectHistogram_Separated';
% Create main plot with the ground truth lines
fnameEnd = sprintf('TR-%i_Dur-%is_C.I.-%i',tr,duration,centerPerc);
fnameRoot = strcat(fnameBegin,'-', fnameEnd);
disp(fnameRoot)
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
    aspects = [aspects;aspectvista];
    
    % Plot it
    h = histogram(aspectvista,20,'Normalization','probability');
    set(h,'LineWidth',2,'EdgeColor',[.5 .5 .5],'EdgeAlpha',0,'FaceAlpha',1,'FaceColor',[.5 .5 .5]);hold on
    plot(median(aspectvista)*[1,1],[0,max(h.Values)],'r-')    
    
    % Add the low noise and mid noise lines now

    
    xlim([1,2.8])
    ylim([0,0.35])
    % legend({lab,'median','Synth Low Noise','Synth Mid Noise'});
    title(lab)
    xlabel('Aspect Ratio')
end
% SAVE PLOT 6B
if doSave;saveas(gcf,fullfile(saveTo, strcat(fnameRoot,['.' ext])),ext);  end






% PLOT 6C
fnameBegin = 'RealData_AspectHistogram_Combined';
% Create main plot with the ground truth lines
fnameEnd = sprintf('TR-%i_Dur-%is_C.I.-%i',tr,duration,centerPerc);
fnameRoot = strcat(fnameBegin,'-', fnameEnd);
disp(fnameRoot)
kk = mrvNewGraphWin(fnameRoot);
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  .5  .5]);

binWidth = 0.05;

aspect2  = aspect2(aspect2 < 4);
% hmid = histogram(aspect2,'DisplayStyle','stairs','BinWidth',h.BinWidth,'Normalization','probability');
hmid = histogram(aspect2,'BinWidth',binWidth,'Normalization','probability');
% set(hmid,'LineWidth',2,'EdgeColor',[.5 .5 .5 ],'LineStyle','-','EdgeAlpha',.5,'FaceAlpha',.5,'FaceColor',[.5 .5 .5 ]);
set(hmid,'LineWidth',2,'EdgeColor','k','FaceAlpha',1,'FaceColor','k');hold on
blow = plot(median(aspect2)*[1,1],[0,.1],'LineWidth',2,'Color','k','LineStyle','-'); 


h = histogram(aspects,35,'Normalization','probability','BinWidth',binWidth);hold on
% set(h,'LineWidth',2,'EdgeColor','k','FaceAlpha',1,'FaceColor','k');
set(h,'LineWidth',2,'EdgeColor',[.5 .5 .5 ],'LineStyle','-','EdgeAlpha',0,'FaceAlpha',.75,'FaceColor',[.5 .5 .5 ]);
a = plot(median(aspects)*[1,1],[0,.1],'Color',[.5 .5 .5 ]);

tool = 'vista6';
% Add the low noise and mid noise lines now
%  aspect1  = aspect1(aspect1 < 4);
% hlow = histogram(aspect1,'DisplayStyle','stairs','BinWidth',h.BinWidth,'Normalization','probability');
% hlow = histogram(aspect1,'BinWidth',h.BinWidth,'Normalization','probability');
% set(hlow,'LineWidth',2,'EdgeColor',[1 .5 .5 ],'LineStyle','-','EdgeAlpha',0,'FaceAlpha',.5,'FaceColor',[1 .5 .5 ]);
% alow = plot(median(aspect1)*[1,1],[0,.1],'LineWidth',2,'Color',[1 .5 .5 ],'LineStyle','-'); 

xlim([1,3])

% legend([h;a;hlow;hmid],{'Experimental Data','Median of Exp. Data','Synth Low Noise','Synth Mid Noise'});

legend([h;a;hmid;blow],{'Experimental Data','Median of Exp. Data','Synth Mid Noise','Median of Synth. Data'});
% SAVE PLOT 6C
if doSave; saveas(gcf,fullfile(saveTo, strcat(fnameRoot,['.' ext])),ext);  end








% PLOT 6D  % THETA vs ANGLE
% prepare data
tool = 'vista6';
thetas = [];
angles = [];
for nl  = 1:length(useLabels)
    lab = useLabels{nl};
    thetas = [thetas;bylabelsums.(tool).(lab).Th];
    angles = [angles;bylabelsums.(tool).(lab).angle];
end

fnameBegin = 'RealData_AnglevsTheta';
% Create main plot with the ground truth lines
fnameEnd = sprintf('TR-%i_Dur-%is_C.I.-%i',tr,duration,centerPerc);
fnameRoot = strcat(fnameBegin,'-', fnameEnd);
disp(fnameRoot)
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
identityLine(gca)
xlim([-90,90])
ylim([-90,90])
xticks(-90:15:90)
yticks(-90:15:90)
% SAVE PLOT 2b
if doSave; saveas(gcf,fullfile(saveTo, strcat(fnameRoot,['.' ext])),ext);  end
























% REMOVED PLOTS
%{
% PLOT OUT 
fnameBegin = 'RealData_AspectHistogramCumSum';
% Create main plot with the ground truth lines
fnameEnd = sprintf('TR-%i_Dur-%is_C.I.-%i',tr,duration,centerPerc);
fnameRoot = strcat(fnameBegin,'-', fnameEnd);
disp(fnameRoot)
kk = mrvNewGraphWin(fnameRoot);
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  .5  .5]);
for nt=1:length(tools)
    tool = tools{nt};
    subplot(1,1,nt);
    a = [];
    for nl = 1:length(useLabels)
        lab  = useLabels{nl};
        % Remove 10% of extreme aspect cases
        B         = prctile(bylabelsums.(tool).(lab).aspect, [5, 95]);
        inRange   = bylabelsums.(tool).(lab).aspect >= B(1) & bylabelsums.(tool).(lab).aspect <= B(2);
        aspect    = bylabelsums.(tool).(lab).aspect(inRange);
        eccen     = bylabelsums.(tool).(lab).eccen(inRange);
        area      = bylabelsums.(tool).(lab).area(inRange);
        if length(lab)==3;if lab(3)=='d';linestyle='-.';else;linestyle='-';end;end

        [N,EDGES] = histcounts(bylabelsums.(tool).(lab).aspect(inRange));
        a = [a; plot(mean([EDGES(2:end);EDGES(1:end-1)]),100*(cumsum(N)/sum(N)),...
             'Color',Cs(nl,:),'LineStyle',linestyle)];hold on;
        
    end
    
    % Now add synth data
    % TR1 vs TR2
    % aspect1  = A1A2.(tool).aspect(A1A2.TR==1);
    % B1       = prctile(aspect1, [5, 95]);
    % inRange1 = aspect1 >= B1(1) & aspect1 <= B1(2);
    % aspect2  = A1A2.(tool).aspect(A1A2.TR==2);
    % B2       = prctile(aspect2, [5, 95]);
    % inRange2 = aspect2 >= B2(1) & aspect2 <= B2(2);
    
    
    % LOW noise MID noise
    aspect1  = A1A2.(tool).aspect(A1A2.noiseLevel=="low");
    B1       = prctile(aspect1, [5, 95]);
    inRange1 = aspect1 >= B1(1) & aspect1 <= B1(2);
    aspect2  = A1A2.(tool).aspect(A1A2.noiseLevel=="mid");
    B2       = prctile(aspect2, [5, 95]);
    inRange2 = aspect2 >= B2(1) & aspect2 <= B2(2);
    
    
    
    [N1,EDGES1] = histcounts(aspect1(inRange1));
    a = [a; plot(mean([EDGES1(2:end);EDGES1(1:end-1)]),100*(cumsum(N1)/sum(N1)),'k--','LineWidth',3)];
    [N2,EDGES2] = histcounts(aspect2(inRange2));
    a = [a; plot(mean([EDGES2(2:end);EDGES2(1:end-1)]),100*(cumsum(N2)/sum(N2)),'k:','LineWidth',3)];
    
    % Add noise ratio
    % [N,EDGES] = histcounts(limRatio);
    % a = [a; plot(mean([EDGES(2:end);EDGES(1:end-1)]),100*(cumsum(N)/sum(N)),'k-','LineWidth',2)];
    
    % legend(a, [useLabels,{'Synth TR1'},{'Synth TR2'},{'Radius noise'}]);grid on
    legend(a, [useLabels,{'Synth Low Noise'},{'Synth Mid Noise'}],'location','best');grid on
    title(tool)
    xlim([1,4])
    xlabel('Aspect Ratio')
    ylabel('Cumulative Frequency (%)')
    set(gca, 'FontSize', 16) 
end
% SAVE PLOT 4
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,['.' ext])),ext);  



 
%}

%% PLOT 6E: Compare r2 values vista4/vista6
useLabels = {'V1','V2','V3'};
saveTo    = '~/gDrive/STANFORD/PROJECTS/2019_PRF_Validation_methods_(Gari)/__PUBLISH__/ELLIPTICAL/Figures/RAW';
ext       = 'svg';   
fnameRoot = 'RealData_R2diff_histogram_filteredAsTheRest'



% kk = mrvNewGraphWin('R2 vista4/vista6');
% set(kk,'Position',[0.007 0.62  0.6 0.5]);
% compTable.vista6.R2 = abs(sqrt(compTable.vista6.r2));
% compTable.vista4.R2 = abs(sqrt(compTable.vista4.r2));
% % scatter(compTable.vista4.R2,compTable.vista6.R2,3,'b','filled')
% histogram(compTable.vista6.R2 - compTable.vista4.R2)
% xlabel('R2, Vista Elliptical - Circular')
% set(gca,'FontSize',20)


% Use the same restrictions we used in the previous plots
% sMajMIN    = 1;
% sMinMIN    = 1;
% sMajMAX    = 3;
% eccenMIN   = 2;
% eccenMAX   = 6;
% minR2      = 0.25;



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
% title(sprintf('Elliptical median variance explained is %1.2g%% larger than Circular',v64f))
if doSave; saveas(gcf,fullfile(saveTo, strcat(fnameRoot,['.' ext])),ext); end



% Plot scatter of vista4 and vista6 R2 values
%{
kk = mrvNewGraphWin('R2 vista4/vista6 COLORS');
set(kk,'Position',[0.007 0.62  0.6 0.5]);
R2vals        = struct();
R2vals.vista6 = [];
R2vals.vista4 = [];
useLabels  = {'V1' ,   'V2',   'V3'};
C =          [1,0,0;  0,1,0;  0,0,1];
for nt=1:length(tools)
    tool = tools{nt};
    Vcolors= [];
    for nl = 1:length(useLabels)
        lab    = useLabels{nl};
        tmpr2  = abs(sqrt(bylabelsums.(tool).(lab).r2));
        R2vals.(tool) = [R2vals.(tool); tmpr2];
        Vcolors = [Vcolors; repmat(C(nl,:),[length(tmpr2),1])];
    end
end
scatter(R2vals.vista4,R2vals.vista6,5,Vcolors,'filled')
xlabel('R2 (Vista Circular)')
ylabel('R2 (Vista Elliptical)')
title('Circular vs elliptical variance explained in V1-V2-V3 [R-G-B] (R2=abs(sqrt(r2)))')
% legend(C,useLabels)
identityLine(gca)
set(gca,'FontSize',20)
%}


%% PLOT (together): Silson 2018 plot: Real Data 7T

kk = mrvNewGraphWin('R2');
set(kk,'Position',[0.007 0.62  0.6 0.5]);
subplot(1,2,1)
histogram(compTable.afni6.R2);legend({'AFNI6 R2 with calccod'});
subplot(1,2,2)
histogram(abs(sqrt(compTable.vista6.r2)));legend({'VISTA6 sqrt(r2) (rss and rawrss)'});





ext        = 'svg';
centerPerc = 90;
eccenInGT  = true;
xlims      = [0,10];
ylims      = [0,10];
tools      = {'afni6'   , 'vista6'};  
duration   = 300;
tr         = 1;
% Filter results
sMajMIN  = 1; % 1;
sMinMIN  = 0.75; % .75;
sMajMAX  = 4; % 4;
eccenMIN = 2;
eccenMAX = 7;
minR2    = .6;
% How many bins
NeccenBins = 6;
NareaBins  = NeccenBins;
%close all



% AFNI
afni6 = compTable.afni6;
[TH,R]      = cart2pol(afni6.x0, afni6.y0);
afni6.angle = rad2deg(TH);
afni6.eccen = R;
afni6.area  = pmEllipseArea(2*afni6.sMaj, 2*afni6.sMin);
afni6       = afni6(afni6.sMaj > sMajMIN & ...
                    afni6.sMin > sMinMIN & ...
                    afni6.sMaj < sMajMAX & ...
                    afni6.eccen > eccenMIN & ...
                    afni6.eccen < eccenMAX & ...
                    afni6.R2 > minR2,:);
eccen     = afni6.eccen;
% Theta can only be [-90,+90]
% Vista and Afni treat it differently it seems
% I added 90 deg to AFNI, but I still don't know if I need it or not. Remove it
theta     = rad2deg(afni6.Th - deg2rad(90));
theta(theta>180) = theta(theta>180) -180;
theta(theta>90)  = theta(theta>90) -180;
% We can express the theta in the same way, because we only care about the
% radiality, not the exact angle
angle     = afni6.angle;
angle(angle>180) = angle(angle>180) -180;
angle(angle>90)  = angle(angle>90) -180;

aspect    = afni6.sMaj  ./ afni6.sMin;
area      = afni6.area;
% Remove 10% of extreme aspect cases
B         = prctile(aspect, [5, 95]);
inRange   = aspect >= B(1) & aspect <= B(2);
aspect    = aspect(inRange);
eccen     = eccen(inRange);
area      = area(inRange);

% VISTA
vista6         = compTable.vista6;
[TH,R]         = cart2pol(vista6.x0, vista6.y0);
vista6.angle   = rad2deg(TH);
vista6.eccen   = R;
vista6.area    = pmEllipseArea(2*vista6.sMaj, 2*vista6.sMin);
vista6         = vista6(vista6.sMaj > sMajMIN & ...
                                  vista6.sMin > sMinMIN & ...
                                  vista6.sMaj < sMajMAX & ...
                                  vista6.eccen > eccenMIN & ...
                                  vista6.eccen < eccenMAX &  ...
                                  vista6.R2 > minR2,:);
eccenvista     = vista6.eccen;
thetavista     = rad2deg(vista6.Th);
thetavista(thetavista>180) = thetavista(thetavista>180) -180;
thetavista(thetavista>90)  = thetavista(thetavista>90) -180;
anglevista     = vista6.angle;
anglevista(anglevista>180) = anglevista(anglevista>180) -180;
anglevista(anglevista>90)  = anglevista(anglevista>90) -180;
areavista      = vista6.area;
aspectvista    = vista6.sMaj  ./ vista6.sMin;
% Remove 10% of extreme aspect cases
B              = prctile(aspectvista, [5, 95]);
inRange        = aspectvista >= B(1) & aspectvista <= B(2);
aspectvista    = aspectvista(inRange);
eccenvista     = eccenvista(inRange);
areavista      = areavista(inRange);
close all

% Obtain the same eccentricities as in the simulations
% xvalues = linspace(1,9,8);
% xvalues(end) = max([max(eccen);max(eccenvista)]);

% Discretize, to bin the eccentricities for afni6
[Y,E] = discretize(eccen,NeccenBins); 
% Discretize, to bin the eccentricities for Vista6
[Yv,E] = discretize(eccenvista,E);

% Discretize, to bin the eccentricities for afni6
[Z,F] = discretize(area,NareaBins); 
% Discretize, to bin the eccentricities for Vista6
[Zv,F] = discretize(areavista,F);


% Apply percentiles and plot individually
% Create the vectors and then plot all together


areamedsize = zeros(1,length(E)-1);
areaminsize = zeros(1,length(E)-1);
areamaxsize = zeros(1,length(E)-1);
 
arvistmedsize = zeros(1,length(E)-1);
arvistminsize = zeros(1,length(E)-1);
arvistmaxsize = zeros(1,length(E)-1);
 
aspectmedecc = zeros(1,length(E)-1);
aspectminecc = zeros(1,length(E)-1);
aspectmaxecc = zeros(1,length(E)-1);
 
aspectvistmedecc = zeros(1,length(E)-1);
aspectvistminecc = zeros(1,length(E)-1);
aspectvistmaxecc = zeros(1,length(E)-1);
 
aspectmedsize = zeros(1,length(E)-1);
aspectminsize = zeros(1,length(E)-1);
aspectmaxsize = zeros(1,length(E)-1);
 
aspectvistmedsize = zeros(1,length(E)-1);
aspectvistminsize = zeros(1,length(E)-1);
aspectvistmaxsize = zeros(1,length(E)-1);
  
for ne=1:(length(E)-1)
    % ECC - SIZE
    % AFNI6
    arsize          = area(Y==ne);
    % Median and std
    areamedsize(ne) = median(arsize);
    areaminsize(ne) = min(arsize);
    areamaxsize(ne) = max(arsize);
    
    % VISTA6    
    arvistsize      = areavista(Yv==ne);
    % Median and std
    arvistmedsize(ne) = median(arvistsize);
    arvistminsize(ne) = min(arvistsize);
    arvistmaxsize(ne) = max(arvistsize);
    
    
    
    
    % ECC - ASPECT
    % AFNI6
    aspecc       = aspect(Y==ne);
    % Median and std
    aspectmedecc(ne) = median(aspecc);
    aspectminecc(ne) = min(aspecc);
    aspectmaxecc(ne) = max(aspecc);
    
    % VISTA6    
    aspvistecc   = aspectvista(Yv==ne);
    % Median and std
    aspectvistmedecc(ne) = median(aspvistecc);
    aspectvistminecc(ne) = min(aspvistecc);
    aspectvistmaxecc(ne) = max(aspvistecc);
    
    
    
    
    
    % SIZE - ASPECT
    % AFNI6
    aspsize       = aspect(Z==ne);
    % Median and std
    aspectmedsize(ne) = median(aspsize);
    aspectminsize(ne) = min(aspsize);
    aspectmaxsize(ne) = max(aspsize);
    
    % VISTA6    
    aspvistsize   = aspectvista(Zv==ne);
    % Median and std
    aspectvistmedsize(ne) = median(aspvistsize);
    aspectvistminsize(ne) = min(aspvistsize);
    aspectvistmaxsize(ne) = max(aspvistsize);
end






% PLOT 1
fnameBegin = 'RealData_SilsonEcc&Size';
% Create main plot with the ground truth lines
fnameEnd = sprintf('TR-%i_Dur-%is_C.I.-%i',tr,duration,centerPerc);
fnameRoot = strcat(fnameBegin,'-', fnameEnd);
disp(fnameRoot)
kk = mrvNewGraphWin(fnameRoot);
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  1 0.5]);


subplot(1,3,1)  % ECCEN vs SIZE
% Plot it
Emidpoints = mean([E(2:end);E(1:end-1)]);
as = scatter(Emidpoints,areamedsize,80,'k','filled');hold on
a  = plot([Emidpoints;Emidpoints] ,...
          [areaminsize  ; areamaxsize], ...
          'Color','k','LineStyle','-','LineWidth',3);  % 0.75*[0 1 0]

bs = scatter(Emidpoints+.1,arvistmedsize,80,'k^','filled');
b  = plot(([Emidpoints;Emidpoints] +.1), ...
    [arvistminsize ;arvistmaxsize], ...
    'Color','k','LineStyle',':','LineWidth',2);
legend([as,bs],{'AFNI', 'VISTA'})
title(strrep(sprintf('%s_TR-%i_Dur-%is_C.I.-%i',...
    tool,tr,duration,centerPerc),'_','\_'))

xlabel('Eccentricity (deg)')
ylabel('Area (deg^2)')
xlim([2,Emidpoints(end)+.2]);
xticks(Emidpoints)
set(gca, 'FontSize', 16) 


subplot(1,3,2)  % ECCEN vs ASPECT
% Plot it
Emidpoints = mean([E(2:end);E(1:end-1)]);
as = scatter(Emidpoints,aspectmedecc,80,'k','filled');hold on
a  = plot([Emidpoints;Emidpoints] ,...
          [aspectminecc  ; aspectmaxecc], ...
          'Color','k','LineStyle','-','LineWidth',3);  % 0.75*[0 1 0]

bs = scatter(Emidpoints+.1,aspectvistmedecc,80,'k^','filled');
b  = plot(([Emidpoints;Emidpoints] +.1), ...
    [aspectvistminecc ;aspectvistmaxecc], ...
    'Color','k','LineStyle',':','LineWidth',2);
legend([as,bs],{'AFNI', 'VISTA'})
title(strrep(sprintf('%s_TR-%i_Dur-%is_C.I.-%i',...
    tool,tr,duration,centerPerc),'_','\_'))

xlabel('Eccentricity (deg)')
ylabel('pRF aspect ratio')
xlim([2,Emidpoints(end)+.2]);
ylim([0,3])
xticks(Emidpoints)
set(gca, 'FontSize', 16) 

subplot(1,3,3)  % SIZE vs ASPECT
% Plot it
Fmidpoints = mean([F(2:end);F(1:end-1)]);
as = scatter(Fmidpoints,aspectmedsize,80,'k','filled');hold on
a  = plot([Fmidpoints;Fmidpoints],...
     [aspectminsize ;aspectmaxsize], ...
     'Color','k','LineStyle','-','LineWidth',3);  % 0.75*[0 1 0]

bs = scatter(Fmidpoints+5,aspectvistmedsize,80,'k^','filled');
b  = plot(([Fmidpoints;Fmidpoints]+5),...
    [aspectvistminsize  ;aspectvistmaxsize], ...
    'Color','k','LineStyle',':','LineWidth',2);
legend([as,bs],{'AFNI', 'VISTA'})
title(strrep(sprintf('%s_TR-%i_Dur-%is_C.I.-%i',...
    tool,tr,duration,centerPerc),'_','\_'))

xlabel('Area (deg^2)')
ylabel('pRF aspect ratio')
xlim([0,Fmidpoints(end)+8]);
ylim([0,3])
xticks(Fmidpoints)
set(gca, 'FontSize', 16) 



% SAVE PLOT 1
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,['.' ext])),ext);  

 







% PLOT 2  % THETA vs ANGLE
fnameBegin = 'RealData_AnglevsTheta';
% Create main plot with the ground truth lines
fnameEnd = sprintf('TR-%i_Dur-%is_C.I.-%i',tr,duration,centerPerc);
fnameRoot = strcat(fnameBegin,'-', fnameEnd);
disp(fnameRoot)
kk = mrvNewGraphWin(fnameRoot);
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  .6  0.5]);

subplot(1,2,1)
% Plot it
as = plot(angle,theta,'ko');hold on; identityLine(gca)
legend({'AFNI'})
xlabel('Angle (deg)')
ylabel('Theta (deg)')
set(gca, 'FontSize', 16) 

% close all

subplot(1,2,2)
bs = plot(anglevista,thetavista,'k^');hold on; identityLine(gca)
legend({'VISTA'})
title(strrep(sprintf('%s_TR-%i_Dur-%is_C.I.-%i',...
    tool,tr,duration,centerPerc),'_','\_'))

xlabel('Angle (deg)')
ylabel('Theta (deg)')
set(gca, 'FontSize', 16) 

% SAVE PLOT 2
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,['.' ext])),ext);  









% PLOT 3
fnameBegin = 'RealData_AspectHistogram';
% Create main plot with the ground truth lines
fnameEnd = sprintf('TR-%i_Dur-%is_C.I.-%i',tr,duration,centerPerc);
fnameRoot = strcat(fnameBegin,'-', fnameEnd);
disp(fnameRoot)
kk = mrvNewGraphWin(fnameRoot);
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  .5  1]);
subplot(2,1,1);
    histogram(aspect,100);legend({'AFNI'})
subplot(2,1,2);
    histogram(aspectvista,100);legend({'VISTA'})
% SAVE PLOT 3
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,['.' ext])),ext);  


% PLOT 4
fnameBegin = 'RealData_AspectHistogramCumSum';
% Create main plot with the ground truth lines
fnameEnd = sprintf('TR-%i_Dur-%is_C.I.-%i',tr,duration,centerPerc);
fnameRoot = strcat(fnameBegin,'-', fnameEnd);
disp(fnameRoot)
kk = mrvNewGraphWin(fnameRoot);
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  .5  .5]);

    [N,EDGES] = histcounts(aspectvista);
    plot(mean([EDGES(2:end);EDGES(1:end-1)]),(cumsum(N)/sum(N)));
    legend({'VISTA'});grid on
% SAVE PLOT 4
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,['.' ext])),ext);  

