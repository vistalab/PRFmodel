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


%% Create dataset for testing circular versus elliptical in Afni

clear all;
COMBINE_PARAMETERS.RF.Centerx0        = [-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6];
COMBINE_PARAMETERS.RF.Centery0        = [-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6];
COMBINE_PARAMETERS.RF.Theta           = [0]; %, deg2rad(45)];
COMBINE_PARAMETERS.RF.sigmaMajor      = [1,2,3,4];
COMBINE_PARAMETERS.RF.sigmaMinor      = 'same'; % 'same' for making it the same to Major
COMBINE_PARAMETERS.TR                 = [2];
    HRF(1).Type                       = 'afni_GAM';
    % HRF(2).Type                     = 'friston';
COMBINE_PARAMETERS.HRF                = HRF;
% TODO: implement a more complex noise addition system. 
% Right now only the parameter for white noise can be edited. 
COMBINE_PARAMETERS.Noise.noise2signal = [0, 0.05, 0.1, 0.2];


synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS);
synthDT = pmForwardModelCalculate(synthDT);
% Visually check that all the combinations we specified are there
% [synthDT.RF(:,{'Centerx0','Centery0','Theta','sigmaMajor','sigmaMinor'}), ...
%  synthDT(:,'TR'), ...
%  synthDT.HRF(:,'Type')...
% ]



% Analyze it
results_AFNI    = pmModelFit(synthDT,'afni_6');

save('results_AFNI.mat', 'results_AFNI')

%% Plot the results
load('results_AFNI.mat');
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
        subplot(length(prfsizes),length(noise2sigs),plotIndex);
        X = compTable.afni.sMin(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
        Y = compTable.afni.sMaj(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
        scatter(X,Y)
        axis equal;
        xlabel('sMin'); ylabel('sMaj')
        title(sprintf('rfsize:%i | noise2sig:%0.2f',prfsize,noise2sig))
        grid on; 
        xlim([-0,8]); ylim([-0,8])
        xticks([0:1:8]);yticks([0:1:8])
        identityLine(gca);
        % h1 = lsline;
        % h1.Color = 'r';
        hold on; [x0,y0] = fitEllipse(X,Y,'r');
        text(2,0.5,sprintf('Center:[%1.2f,%1.2f], ratio:%1.2f',x0,y0,y0/x0));
    end
end


























