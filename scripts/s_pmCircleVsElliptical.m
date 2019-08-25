%% s_pmCircleVsElliptical.m
%{
AFNI test to check circular versus elliptical
- create hundreds of synthetic signals.
- use the same HRF they use for solving it to generate the signal
- vary the positions randomly in all of them, not controlling for them.
Control this:
- increasing levels of noise: 4 different levels for example
- size of the prf: sigmaMinor = sigmaMajor = [1,2,3,4] for example. 

This is 4 (noise) x 4 (sigma size) x 2 (oldAfni,newAfni) groups. 
In all of them, ideally, the sigma ratio (sigmaMajor/sigmaMinor) should 
be 1+- epsilon. If it is not, or if it is for the old software and not for 
the new, we'll have learned something. 
%}
addpath(genpath('/home/glerma/soft/afni_matlab'));

%% Create dataset for testing circular versus elliptical in Afni
clear all;
COMBINE_PARAMETERS.RF.Centerx0        = [3];  % [-6 5 4 3 2 1 0 1 2 3 4 5 6];
COMBINE_PARAMETERS.RF.Centery0        = [3];  % [-6 5 4 3 2 1 0 1 2 3 4 5 6];
COMBINE_PARAMETERS.RF.Theta           = [0]; %, deg2rad(45)];
COMBINE_PARAMETERS.RF.sigmaMajor      = [0.5,1,2,4,8]  % [1,2,3,4];
COMBINE_PARAMETERS.RF.sigmaMinor      = 'same'; % 'same' for making it the same to Major
COMBINE_PARAMETERS.TR                 = [2];
    HRF(1).Type                       = 'afni_GAM';
    % HRF(2).Type                     = 'friston';
COMBINE_PARAMETERS.HRF                = HRF;
% TODO: implement a more complex noise addition system. 
% Right now only the parameter for white noise can be edited. 
COMBINE_PARAMETERS.Noise.noise2signal = [0, 0.05, 0.1, 0.2];
synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS, 'mult',100);
synthDT = pmForwardModelCalculate(synthDT);
% Visually check that all the combinations we specified are there
% [synthDT.RF(:,{'Centerx0','Centery0','Theta','sigmaMajor','sigmaMinor'}), ...
%  synthDT(:,'TR'), ...
%  synthDT.HRF(:,'Type')...
% ]



% Analyze it
results_AFNI    = pmModelFit(synthDT,'afni_6');

% save('resultsAFNI_multipleCenters.mat', 'results_AFNI')
% save('synthDT_multipleCenters.mat', 'synthDT')
save(fullfile(pmRootPath,'local','resultsAFNI_oneCenter.mat'), 'results_AFNI')
save(fullfile(pmRootPath,'local','synthDT_oneCenter.mat'), 'synthDT')

%% Plot the results
% load('resultsAFNI_multipleCenters.mat');
% load('synthDT_multipleCenters.mat');
load(fullfile(pmRootPath,'local','resultsAFNI_oneCenter.mat'));
load(fullfile(pmRootPath,'local','synthDT_oneCenter.mat'));


paramDefaults = {'Centerx0','Centery0','Theta','sigmaMinor','sigmaMajor'};
[compTable, tSeries] = pmResultsCompare(synthDT, ... % Defines the input params
                            {'afni'}, ... % Analysis names we want to see: 'aPRF','vista',
                            {results_AFNI}, ... % results_analyzePRF,results_vista,
                            'params', paramDefaults, ...
                            'shorten names',true); 

noise2sigs = unique(compTable.noise2sig);
prfsizes   = unique(compTable.synth.sMaj);
mrvNewGraphWin('AFNI comparisons');
plotIndex = 0;
for ns=1:length(noise2sigs)
    noise2sig = noise2sigs(ns);
    for np=1:length(prfsizes)
        prfsize = prfsizes(np);
        % Now we can create the subplots
        plotIndex = plotIndex + 1;
        subplot(length(noise2sigs),length(prfsizes),plotIndex);
        X = compTable.afni.sMin(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
        Y = compTable.afni.sMaj(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
        scatter(X,Y)
        axis equal;
        xlabel('sMin'); ylabel('sMaj')
        title(sprintf('rfsize:%1.1f | noise2sig:%0.2f',prfsize,noise2sig))
        grid on; 
        xlim([-0,8]); ylim([-0,8])
        xticks([0:1:8]);yticks([0:1:8])
        identityLine(gca);
        % h1 = lsline;
        % h1.Color = 'r';
        hold on; 
        % [x0,y0] = fitEllipse(X,Y,'r');
        x0  = median(X); y0 = median(Y);
        sdx = std(X); sdy = std(Y);
        plot([x0-sdx, x0+sdx],[y0 y0],'r');
        plot([x0 x0],[y0-sdy, y0+sdy],'r');
        text(2,0.5,sprintf('Median:[%1.2f(%1.2f),%1.2f(%1.2f)], ratio:%1.2f',x0,sdx,y0,sdy,y0/x0));
    end
end


% Now create the histograms with the thetas
noise2sigs = unique(compTable.noise2sig);
prfsizes   = unique(compTable.synth.sMaj);
mrvNewGraphWin('AFNI comparisons');
plotIndex = 0;
for ns=1:length(noise2sigs)
    noise2sig = noise2sigs(ns);
    for np=1:length(prfsizes)
        prfsize = prfsizes(np);
        % Now we can create the subplots
        plotIndex = plotIndex + 1;
        subplot(length(noise2sigs),length(prfsizes),plotIndex);
        % X = rad2deg(compTable.afni.Th(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize));
        X = (compTable.afni.Th(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize));
        X(X>0) = rad2deg(X(X>0));
        X(X<0) = rad2deg(X(X<0))-360;
        % Y = compTable.afni.sMaj(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
        hist(X)
        axis equal;
        xlabel('Theta in deg'); 
        % ylabel('sMaj')
        title(sprintf('rfsize:%1.1f | noise2sig:%0.2f',prfsize,noise2sig))
        % grid on; 
        % xlim([-0,8]); ylim([-0,8])
        xlim([-90,90]); ylim([0,40]);
        xticks([-80:20:80]);
        % identityLine(gca);
        % h1 = lsline;
        % h1.Color = 'r';s
        hold on; 
        % [x0,y0] = fitEllipse(X,Y,'r');
        x0  = median(X); 
        % y0 = median(Y);
        sdx = std(X); sdy = std(Y);
        % plot([x0-sdx, x0+sdx],[y0 y0],'r');
        % plot([x0 x0],[y0-sdy, y0+sdy],'r');
        % text(2,0.5,sprintf('Median:[%1.2f(%1.2f),%1.2f(%1.2f)], ratio:%1.2f',x0,sdx,y0,sdy,y0/x0));
    end
end



% Now plot the solutions

%{
Given stimulus images over time s(x,y,t), find x0, y0, sigma, R and
theta values that produce a best fit of the model to the data.  Here
x0, y0 are taken to be the center of the population receptive field,
sigma is the minor width of it (sigma_x, below), sigrat R is the ratio
(sigma_y / sigma_x), and theta is the rotation from the y-direction
major axis (so zero is in the positive y-direction).

We assume sigma_y >= sigma_x and refer to sigrat >= 1, since that
sufficiently represents all possibilities.  The reciprocol would
come from the negative complimentary angle, and would therefore be a
redundant solution.

parameter domains:
x,y        : [-1,1], scaled by the mask, itself
sigma      : (0,1], where 1 means the mask radius
R (sigrat) : [1,inf), since sigma defines the smaller size
theta      : [-PI/2, PI/2), since rotation by PI has no effect

%}

noise2sigs = unique(compTable.noise2sig);
prfsizes   = unique(compTable.synth.sMaj);
mrvNewGraphWin('AFNI comparisons');
plotIndex = 0;
for ns=1:length(noise2sigs)
    noise2sig = noise2sigs(ns);
    for np=1:length(prfsizes)
        prfsize = prfsizes(np);
        % Now we can create the subplots
        plotIndex = plotIndex + 1;
        subplot(length(noise2sigs),length(prfsizes),plotIndex);
        
        Th   = compTable.afni.Th(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
        x0   = compTable.afni.x0(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
        y0   = compTable.afni.y0(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
        sMaj = compTable.afni.sMaj(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
        sMin = compTable.afni.sMin(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
        
        for ii=1:1:length(sMin)
            x=x0(ii);
            y=y0(ii);
            a = sMaj(ii);
            b = sMin(ii);
            t = linspace(0,2*pi);
            X = a*cos(t);
            Y = b*sin(t);
            w = Th(ii);
            XX = x + X*cos(w) - Y*sin(w);
            YY = y + X*sin(w) + Y*cos(w);
            plot(XX,YY,'k','LineWidth',1)
            axis equal
            hold on
        end
        
        xlabel('Degrees'); 
        ylabel('Degrees')
        title(sprintf('rfsize:%1.1f | noise2sig:%0.2f',prfsize,noise2sig))
        % grid on; 
%         xlim([-10,10]); xticks(-10:2:10);
%         ylim([-10,10]); yticks(-10:2:10);
        % xlim([-15,15]); xticks(-14:2:14);
        % ylim([-15,15]); yticks(-14:2:14);
        xlim([1,5]); xticks(1:1:5);
        ylim([1,5]); yticks(1:1:5);
    end
end



% Do just one
noise2sig=0.1;
prfsize = 1;


mrvNewGraphWin('AFNI comparisons, one case');

Th   = compTable.afni.Th(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
x0   = compTable.afni.x0(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
y0   = compTable.afni.y0(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
sMaj = compTable.afni.sMaj(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
sMin = compTable.afni.sMin(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);

for ii=1:1:length(sMin)
    x=x0(ii);
    y=y0(ii);
    a = sMaj(ii);
    b = sMin(ii);
    t = linspace(0,2*pi);
    X = a*cos(t);
    Y = b*sin(t);
    w = Th(ii);
    XX = x + X*cos(w) - Y*sin(w);
    YY = y + X*sin(w) + Y*cos(w);
    % Define some colours
    lightblue = [0.75 0.75 1];
    blue      = [0    0    1];
    lightred  = [1    0.75 0.6];
    red       = [1    0    0];
    black     = [0    0    0];
    
    if w < -1.571/2
        colour = red;
    elseif w > -1.571/2 && w < -deg2rad(5)
        colour = lightred;
    elseif w < 1.571/2 && w > deg2rad(5)
        colour = lightblue;
    elseif w > 1.571/2
        colour = blue;
    else
        colour = black;
    end
    
    if a/b >= 1.5
        width  = 2;
        estilo = '-';
    else
        width  = 1;
        estilo = '--';
    end
    
    
    plot(XX,YY,'color',colour,'LineWidth',width,'LineStyle',estilo)
    axis equal
    hold on
end
xlabel('Degrees');
ylabel('Degrees')
title(sprintf('rfsize:%1.1f | noise2sig:%0.2f',prfsize,noise2sig))
% grid on;
%         xlim([-10,10]); xticks(-10:2:10);
%         ylim([-10,10]); yticks(-10:2:10);
% xlim([-10,10]); xticks(-10:2:10);
% ylim([-10,10]); yticks(-10:2:10);
% text(-9.5,9.5,'RED(-90>-45), lightRED(-45>-5), BLACK(+-5deg), lightBLUE(4>45), BLUE(45>90)')
% text(-9.5,9,'THIN DASHED if sigRat < 1.5, THICK SOLID is sigRat >= 1.5')
xlim([1,5]); xticks(1:1:5);
ylim([1,5]); yticks(1:1:5);
text(1.5,4.75,'RED(-90>-45), lightRED(-45>-5), BLACK(+-5deg), lightBLUE(4>45), BLUE(45>90)')
text(1.5,4.5,'THIN DASHED if sigRat < 1.5, THICK SOLID is sigRat >= 1.5')


% Now plot theta versus sigrat
noise2sigs = unique(compTable.noise2sig);
prfsizes   = unique(compTable.synth.sMaj);
mrvNewGraphWin('AFNI theta vs sigrat');
plotIndex = 0;
for ns=1:length(noise2sigs)
    noise2sig = noise2sigs(ns);
    for np=1:length(prfsizes)
        prfsize = prfsizes(np);
        % Now we can create the subplots
        plotIndex = plotIndex + 1;
        subplot(length(noise2sigs),length(prfsizes),plotIndex);
        
        Th   = compTable.afni.Th(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
        Th(Th>0) = rad2deg(Th(Th>0));
        Th(Th<0) = rad2deg(Th(Th<0))-360;
        sMaj = compTable.afni.sMaj(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
        sMin = compTable.afni.sMin(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
        sigRat = sMaj ./ sMin;
        
        scatter(Th,sigRat)
        xlabel('Theta in Degrees'); 
        ylabel('sigRat (unitless)')
        title(sprintf('rfsize:%1.1f | noise2sig:%0.2f',prfsize,noise2sig))
        % grid on; 
%         xlim([-10,10]); xticks(-10:2:10);
%         ylim([-10,10]); yticks(-10:2:10);
        xlim([-90,90]); xticks(-80:20:80);
        ylim([1,6]); yticks(1:1:6);
        hold on;
        ls = lsline;
        ls.Color = 'r';
    end
end





%% Now do the same with mrVista
load(fullfile(pmRootPath,'local','synthDT_oneCenter.mat'));
results_vista      = pmModelFit(synthDT,'vistasoft', ...
                                        'model','one oval gaussian', ...
                                        'grid', false, ... % if true, returns gFit
                                        'wSearch', 'coarse to fine');  %  and hrf

save(fullfile(pmRootPath,'local','resultsVISTA_oneCenter.mat'), 'results_vista')


paramDefaults = {'Centerx0','Centery0','Theta','sigmaMinor','sigmaMajor'};
[compTable, tSeries] = pmResultsCompare(synthDT, ... % Defines the input params
                            {'vista'}, ... % Analysis names we want to see: 'aPRF','vista',
                            {results_vista}, ... % results_analyzePRF,results_vista,
                            'params', paramDefaults, ...
                            'shorten names',true); 

noise2sigs = unique(compTable.noise2sig);
prfsizes   = unique(compTable.synth.sMaj);
mrvNewGraphWin('VISTA comparisons');
plotIndex = 0;
for ns=1:length(noise2sigs)
    noise2sig = noise2sigs(ns);
    for np=1:length(prfsizes)
        prfsize = prfsizes(np);
        % Now we can create the subplots
        plotIndex = plotIndex + 1;
        subplot(length(noise2sigs),length(prfsizes),plotIndex);
        X = compTable.vista.sMin(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
        Y = compTable.vista.sMaj(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
        scatter(X,Y)
        axis equal;
        xlabel('sMin'); ylabel('sMaj')
        title(sprintf('rfsize:%1.1f | noise2sig:%0.2f',prfsize,noise2sig))
        grid on; 
        xlim([-0,8]); ylim([-0,8])
        xticks([0:1:8]);yticks([0:1:8])
        identityLine(gca);
        % h1 = lsline;
        % h1.Color = 'r';
        hold on; 
        % [x0,y0] = fitEllipse(X,Y,'r');
        x0  = median(X); y0 = median(Y);
        sdx = std(X); sdy = std(Y);
        plot([x0-sdx, x0+sdx],[y0 y0],'r');
        plot([x0 x0],[y0-sdy, y0+sdy],'r');
        text(2,0.5,sprintf('Median:[%1.2f(%1.2f),%1.2f(%1.2f)], ratio:%1.2f',x0,sdx,y0,sdy,y0/x0));
    end
end


% Now create the histograms with the thetas
noise2sigs = unique(compTable.noise2sig);
prfsizes   = unique(compTable.synth.sMaj);
mrvNewGraphWin('VISTA comparisons');
plotIndex = 0;
for ns=1:length(noise2sigs)
    noise2sig = noise2sigs(ns);
    for np=1:length(prfsizes)
        prfsize = prfsizes(np);
        % Now we can create the subplots
        plotIndex = plotIndex + 1;
        subplot(length(noise2sigs),length(prfsizes),plotIndex);
        % X = rad2deg(compTable.afni.Th(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize));
        X = (compTable.vista.Th(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize));
        X(X>0) = rad2deg(X(X>0));
        X(X<0) = rad2deg(X(X<0))-360;
        % Y = compTable.vista.sMaj(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
        hist(X)
        axis equal;
        xlabel('Theta in deg'); 
        % ylabel('sMaj')
        title(sprintf('rfsize:%1.1f | noise2sig:%0.2f',prfsize,noise2sig))
        % grid on; 
        % xlim([-0,8]); ylim([-0,8])
        xlim([-90,90]); ylim([0,40]);
        xticks([-80:20:80]);
        % identityLine(gca);
        % h1 = lsline;
        % h1.Color = 'r';s
        hold on; 
        % [x0,y0] = fitEllipse(X,Y,'r');
        x0  = median(X); 
        % y0 = median(Y);
        sdx = std(X); sdy = std(Y);
        % plot([x0-sdx, x0+sdx],[y0 y0],'r');
        % plot([x0 x0],[y0-sdy, y0+sdy],'r');
        % text(2,0.5,sprintf('Median:[%1.2f(%1.2f),%1.2f(%1.2f)], ratio:%1.2f',x0,sdx,y0,sdy,y0/x0));
    end
end



% Now plot the solutions

%{
Given stimulus images over time s(x,y,t), find x0, y0, sigma, R and
theta values that produce a best fit of the model to the data.  Here
x0, y0 are taken to be the center of the population receptive field,
sigma is the minor width of it (sigma_x, below), sigrat R is the ratio
(sigma_y / sigma_x), and theta is the rotation from the y-direction
major axis (so zero is in the positive y-direction).

We assume sigma_y >= sigma_x and refer to sigrat >= 1, since that
sufficiently represents all possibilities.  The reciprocol would
come from the negative complimentary angle, and would therefore be a
redundant solution.

parameter domains:
x,y        : [-1,1], scaled by the mask, itself
sigma      : (0,1], where 1 means the mask radius
R (sigrat) : [1,inf), since sigma defines the smaller size
theta      : [-PI/2, PI/2), since rotation by PI has no effect

%}

noise2sigs = unique(compTable.noise2sig);
prfsizes   = unique(compTable.synth.sMaj);
mrvNewGraphWin('VISTA comparisons');
plotIndex = 0;
for ns=1:length(noise2sigs)
    noise2sig = noise2sigs(ns);
    for np=1:length(prfsizes)
        prfsize = prfsizes(np);
        % Now we can create the subplots
        plotIndex = plotIndex + 1;
        subplot(length(noise2sigs),length(prfsizes),plotIndex);
        
        Th   = compTable.vista.Th(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
        x0   = compTable.vista.x0(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
        y0   = compTable.vista.y0(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
        sMaj = compTable.vista.sMaj(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
        sMin = compTable.vista.sMin(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
        
        for ii=1:1:length(sMin)
            x=x0(ii);
            y=y0(ii);
            a = sMaj(ii);
            b = sMin(ii);
            t = linspace(0,2*pi);
            X = a*cos(t);
            Y = b*sin(t);
            w = Th(ii);
            XX = x + X*cos(w) - Y*sin(w);
            YY = y + X*sin(w) + Y*cos(w);
            plot(XX,YY,'k','LineWidth',1)
            axis equal
            hold on
        end
        
        xlabel('Degrees'); 
        ylabel('Degrees')
        title(sprintf('rfsize:%1.1f | noise2sig:%0.2f',prfsize,noise2sig))
        % grid on; 
%         xlim([-10,10]); xticks(-10:2:10);
%         ylim([-10,10]); yticks(-10:2:10);
        % xlim([-15,15]); xticks(-14:2:14);
        % ylim([-15,15]); yticks(-14:2:14);
        xlim([1,5]); xticks(1:1:5);
        ylim([1,5]); yticks(1:1:5);
    end
end



% Do just one
noise2sig=0.1;
prfsize = 1;


mrvNewGraphWin('VISTA comparisons, one case');

Th   = compTable.vista.Th(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
x0   = compTable.vista.x0(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
y0   = compTable.vista.y0(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
sMaj = compTable.vista.sMaj(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
sMin = compTable.vista.sMin(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);

for ii=1:1:length(sMin)
    x=x0(ii);
    y=y0(ii);
    a = sMaj(ii);
    b = sMin(ii);
    t = linspace(0,2*pi);
    X = a*cos(t);
    Y = b*sin(t);
    w = Th(ii);
    XX = x + X*cos(w) - Y*sin(w);
    YY = y + X*sin(w) + Y*cos(w);
    % Define some colours
    lightblue = [0.75 0.75 1];
    blue      = [0    0    1];
    lightred  = [1    0.75 0.6];
    red       = [1    0    0];
    black     = [0    0    0];
    
    if w < -1.571/2
        colour = red;
    elseif w > -1.571/2 && w < -deg2rad(5)
        colour = lightred;
    elseif w < 1.571/2 && w > deg2rad(5)
        colour = lightblue;
    elseif w > 1.571/2
        colour = blue;
    else
        colour = black;
    end
    
    if a/b >= 1.5
        width  = 2;
        estilo = '-';
    else
        width  = 1;
        estilo = '--';
    end
    
    
    plot(XX,YY,'color',colour,'LineWidth',width,'LineStyle',estilo)
    axis equal
    hold on
end
xlabel('Degrees');
ylabel('Degrees')
title(sprintf('rfsize:%1.1f | noise2sig:%0.2f',prfsize,noise2sig))
% grid on;
%         xlim([-10,10]); xticks(-10:2:10);
%         ylim([-10,10]); yticks(-10:2:10);
% xlim([-10,10]); xticks(-10:2:10);
% ylim([-10,10]); yticks(-10:2:10);
% text(-9.5,9.5,'RED(-90>-45), lightRED(-45>-5), BLACK(+-5deg), lightBLUE(4>45), BLUE(45>90)')
% text(-9.5,9,'THIN DASHED if sigRat < 1.5, THICK SOLID is sigRat >= 1.5')
xlim([1,5]); xticks(1:1:5);
ylim([1,5]); yticks(1:1:5);
text(1.5,4.75,'RED(-90>-45), lightRED(-45>-5), BLACK(+-5deg), lightBLUE(4>45), BLUE(45>90)')
text(1.5,4.5,'THIN DASHED if sigRat < 1.5, THICK SOLID is sigRat >= 1.5')


% Now plot theta versus sigrat
noise2sigs = unique(compTable.noise2sig);
prfsizes   = unique(compTable.synth.sMaj);
mrvNewGraphWin('VISTA theta vs sigrat');
plotIndex = 0;
for ns=1:length(noise2sigs)
    noise2sig = noise2sigs(ns);
    for np=1:length(prfsizes)
        prfsize = prfsizes(np);
        % Now we can create the subplots
        plotIndex = plotIndex + 1;
        subplot(length(noise2sigs),length(prfsizes),plotIndex);
        
        Th   = compTable.vista.Th(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
        Th(Th>0) = rad2deg(Th(Th>0));
        Th(Th<0) = rad2deg(Th(Th<0))-360;
        sMaj = compTable.vista.sMaj(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
        sMin = compTable.vista.sMin(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
        sigRat = sMaj ./ sMin;
        
        scatter(Th,sigRat)
        xlabel('Theta in Degrees'); 
        ylabel('sigRat (unitless)')
        title(sprintf('rfsize:%1.1f | noise2sig:%0.2f',prfsize,noise2sig))
        % grid on; 
%         xlim([-10,10]); xticks(-10:2:10);
%         ylim([-10,10]); yticks(-10:2:10);
        xlim([-90,90]); xticks(-80:20:80);
        ylim([1,6]); yticks(1:1:6);
        hold on;
        ls = lsline;
        ls.Color = 'r';
    end
end



%% COMpare old AFNI with new afni

old_oneCenter   = load(fullfile(pmRootPath,'local','results_old_AFNI_oneCenter.mat'));
new_oneCenter   = load(fullfile(pmRootPath,'local','results_new_AFNI_oneCenter.mat'));
old_multiCenter = load(fullfile(pmRootPath,'local','results_old_AFNI_multipleCenters.mat'));
new_multiCenter = load(fullfile(pmRootPath,'local','results_new_AFNI_multipleCenters.mat'));






% Plot for oneCenter
load('synthDT_oneCenter.mat');
paramDefaults = {'Centerx0','Centery0','Theta','sigmaMinor','sigmaMajor'};
OLDcompTable = pmResultsCompare(synthDT, ... % Defines the input params
                            {'afni'}, ... % Analysis names we want to see: 'aPRF','vista',
                            {old_oneCenter.results_AFNI}, ... % results_analyzePRF,results_vista,
                            'params', paramDefaults, ...
                            'shorten names',true); 
NEWcompTable = pmResultsCompare(synthDT, ... % Defines the input params
                            {'afni'}, ... % Analysis names we want to see: 'aPRF','vista',
                            {new_oneCenter.results_AFNI}, ... % results_analyzePRF,results_vista,
                            'params', paramDefaults, ...
                            'shorten names',true); 





noise2sigs = unique(OLDcompTable.noise2sig);
prfsizes   = unique(OLDcompTable.synth.sMaj);
mrvNewGraphWin('AFNI old vs new oneCenter');
plotIndex = 0;
for ns=1:length(noise2sigs)
    noise2sig = noise2sigs(ns);
    for np=1:length(prfsizes)
        prfsize = prfsizes(np);
        % Now we can create the subplots
        plotIndex = plotIndex + 1;
        subplot(length(noise2sigs),length(prfsizes),plotIndex);
        oldsMin = OLDcompTable.afni.sMin(OLDcompTable.noise2sig==noise2sig & OLDcompTable.synth.sMaj==prfsize);
        oldsMaj = OLDcompTable.afni.sMaj(OLDcompTable.noise2sig==noise2sig & OLDcompTable.synth.sMaj==prfsize);
        newsMin = NEWcompTable.afni.sMin(NEWcompTable.noise2sig==noise2sig & NEWcompTable.synth.sMaj==prfsize);
        newsMaj = NEWcompTable.afni.sMaj(NEWcompTable.noise2sig==noise2sig & NEWcompTable.synth.sMaj==prfsize);
        oldSigRat = oldsMaj ./ oldsMin;
        newSigRat = newsMaj ./ newsMin;
        scatter(oldSigRat,newSigRat)
        axis equal;
        xlabel('oldSigmaRatio'); ylabel('newSigmaRatio')
        title(sprintf('rfsize:%1.1f | noise2sig:%0.2f',prfsize,noise2sig))
        grid on; 
        xlim([1,8]); ylim([1,8])
        xticks([1:1:8]);yticks([1:1:8])
        identityLine(gca);
        % h1 = lsline;
        % h1.Color = 'r';
        hold on; 
        % [x0,y0] = fitEllipse(X,Y,'r');
        % x0  = median(X); y0 = median(Y);
        % sdx = std(X); sdy = std(Y);
        % plot([x0-sdx, x0+sdx],[y0 y0],'r');
        % plot([x0 x0],[y0-sdy, y0+sdy],'r');
        % text(2,0.5,sprintf('Median:[%1.2f(%1.2f),%1.2f(%1.2f)], ratio:%1.2f',x0,sdx,y0,sdy,y0/x0));
    end
end




% Plot for multiCenter
load('synthDT_multipleCenters.mat');
paramDefaults = {'Centerx0','Centery0','Theta','sigmaMinor','sigmaMajor'};
OLDcompTable = pmResultsCompare(synthDT, ... % Defines the input params
                            {'afni'}, ... % Analysis names we want to see: 'aPRF','vista',
                            {old_multiCenter.results_AFNI}, ... % results_analyzePRF,results_vista,
                            'params', paramDefaults, ...
                            'shorten names',true); 
NEWcompTable = pmResultsCompare(synthDT, ... % Defines the input params
                            {'afni'}, ... % Analysis names we want to see: 'aPRF','vista',
                            {new_multiCenter.results_AFNI}, ... % results_analyzePRF,results_vista,
                            'params', paramDefaults, ...
                            'shorten names',true); 





noise2sigs = unique(OLDcompTable.noise2sig);
prfsizes   = unique(OLDcompTable.synth.sMaj);
mrvNewGraphWin('AFNI old vs new multiCenter');
plotIndex = 0;
for ns=1:length(noise2sigs)
    noise2sig = noise2sigs(ns);
    for np=1:length(prfsizes)
        prfsize = prfsizes(np);
        % Now we can create the subplots
        plotIndex = plotIndex + 1;
        subplot(length(noise2sigs),length(prfsizes),plotIndex);
        oldsMin = OLDcompTable.afni.sMin(OLDcompTable.noise2sig==noise2sig & OLDcompTable.synth.sMaj==prfsize);
        oldsMaj = OLDcompTable.afni.sMaj(OLDcompTable.noise2sig==noise2sig & OLDcompTable.synth.sMaj==prfsize);
        newsMin = NEWcompTable.afni.sMin(NEWcompTable.noise2sig==noise2sig & NEWcompTable.synth.sMaj==prfsize);
        newsMaj = NEWcompTable.afni.sMaj(NEWcompTable.noise2sig==noise2sig & NEWcompTable.synth.sMaj==prfsize);
        oldSigRat = oldsMaj ./ oldsMin;
        newSigRat = newsMaj ./ newsMin;
        scatter(oldSigRat,newSigRat)
        axis equal;
        xlabel('oldSigmaRatio'); ylabel('newSigmaRatio')
        title(sprintf('rfsize:%1.1f | noise2sig:%0.2f',prfsize,noise2sig))
        grid on; 
        xlim([1,8]); ylim([1,8])
        xticks([1:1:8]);yticks([1:1:8])
        identityLine(gca);
        % h1 = lsline;
        % h1.Color = 'r';
        hold on; 
        % [x0,y0] = fitEllipse(X,Y,'r');
        % x0  = median(X); y0 = median(Y);
        % sdx = std(X); sdy = std(Y);
        % plot([x0-sdx, x0+sdx],[y0 y0],'r');
        % plot([x0 x0],[y0-sdy, y0+sdy],'r');
        % text(2,0.5,sprintf('Median:[%1.2f(%1.2f),%1.2f(%1.2f)], ratio:%1.2f',x0,sdx,y0,sdy,y0/x0));
    end
end


