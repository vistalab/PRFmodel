%% s_pmAccVsprec.m
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


%% Create dataset for testing noiseless and noise values in all datasets. 
% USE IDEAL HRF FOR EACH OF ONE
clear all;

% DO aPRF
COMBINE_PARAMETERS                    = struct();
COMBINE_PARAMETERS.RF.Centerx0        = [0];  % [-6 5 4 3 2 1 0 1 2 3 4 5 6];
COMBINE_PARAMETERS.RF.Centery0        = [0];  % [-6 5 4 3 2 1 0 1 2 3 4 5 6];
COMBINE_PARAMETERS.RF.Theta           = [0]; %, deg2rad(45)];
COMBINE_PARAMETERS.RF.sigmaMajor      = [0.5,1,2,4,8];  % [1,2,3,4];
COMBINE_PARAMETERS.RF.sigmaMinor      = 'same'; % 'same' for making it the same to Major
COMBINE_PARAMETERS.TR                 = [1.5];
    HRF(1).Type                       = 'canonical';
COMBINE_PARAMETERS.HRF                = HRF;
% Right now only the parameter for white noise can be edited. 
COMBINE_PARAMETERS.Noise.noise2signal = [0:0.05:0.45];
synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS, 'mult',100);
synthDT = pmForwardModelCalculate(synthDT);
% Save it
aprfsynthDTfName = ['synthDT_aprf_oneCenter_' datestr(datetime,'yyyymmddTHHMMSS','local') '.mat'];
save(fullfile(pmRootPath,'local',aprfsynthDTfName), 'synthDT');

% Analyze it, it takes a lot of time
results    = pmModelFit(synthDT, 'aprf', 'useParallel', true);
% Save the result, we don't want to lose it because matlab crashes...
% Use the same datetime to associate synthDT and result
aprfresultfName = ['result_aprf_oneCenter_' datestr(datetime,'yyyymmddTHHMMSS','local') '.mat'];
% if exist(fullfile(pmRootPath,'local',aprfresultfName),'file')
%     error('Change file name or delete existing file.')
    % delete(fullfile(pmRootPath,'local',resultfName))
% else
    save(fullfile(pmRootPath,'local',aprfresultfName), 'results');
% end







% DO POPEYE
COMBINE_PARAMETERS                    = struct();
COMBINE_PARAMETERS.RF.Centerx0        = [0];  % [-6 5 4 3 2 1 0 1 2 3 4 5 6];
COMBINE_PARAMETERS.RF.Centery0        = [0];  % [-6 5 4 3 2 1 0 1 2 3 4 5 6];
COMBINE_PARAMETERS.RF.Theta           = [0]; %, deg2rad(45)];
COMBINE_PARAMETERS.RF.sigmaMajor      = [0.5,1,2,4,8];  % [1,2,3,4];
COMBINE_PARAMETERS.RF.sigmaMinor      = 'same'; % 'same' for making it the same to Major
COMBINE_PARAMETERS.TR                 = [1.5];
    HRF(1).Type                       = 'popeye_twogammas';
COMBINE_PARAMETERS.HRF                = HRF;
% Right now only the parameter for white noise can be edited. 
COMBINE_PARAMETERS.Noise.noise2signal = [0:0.05:0.45];
synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS, 'mult',100);
synthDT = pmForwardModelCalculate(synthDT);
% Save it
popsynthDTfName = ['synthDT_pop_oneCenter_' datestr(datetime,'yyyymmddTHHMMSS','local') '.mat'];
save(fullfile(pmRootPath,'local',popsynthDTfName), 'synthDT');

% Analyze it, it takes a lot of time
results    = pmModelFit(synthDT, 'popeye_onegaussian');
% Save the result, we don't want to lose it because matlab crashes...
% Use the same datetime to associate synthDT and result
popresultfName = ['results_pop_oneCenter_' datestr(datetime,'yyyymmddTHHMMSS','local') '.mat'];
% if exist(fullfile(pmRootPath,'local',popresultfName),'file')
%     error('Change file name or delete existing file.')
    % delete(fullfile(pmRootPath,'local',resultfName))
% else
    save(fullfile(pmRootPath,'local',popresultfName), 'results');
% end



 


% DO mrVISTA
COMBINE_PARAMETERS                    = struct();
COMBINE_PARAMETERS.RF.Centerx0        = [0];  % [-6 5 4 3 2 1 0 1 2 3 4 5 6];
COMBINE_PARAMETERS.RF.Centery0        = [0];  % [-6 5 4 3 2 1 0 1 2 3 4 5 6];
COMBINE_PARAMETERS.RF.Theta           = [0]; %, deg2rad(45)];
COMBINE_PARAMETERS.RF.sigmaMajor      = [0.5,1,2,4,8];  % [1,2,3,4];
COMBINE_PARAMETERS.RF.sigmaMinor      = 'same'; % 'same' for making it the same to Major
COMBINE_PARAMETERS.TR                 = [1.5];
    HRF(1).Type                       = 'vista_twogammas';
COMBINE_PARAMETERS.HRF                = HRF;
% TODO: implement a more complex noise addition system. 
% Right now only the parameter for white noise can be edited. 
COMBINE_PARAMETERS.Noise.noise2signal = [0:0.05:0.45];
synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS, 'mult',100);
synthDT = pmForwardModelCalculate(synthDT);
% Save it
vistasynthDTfName = ['synthDT_vista_oneCenter_' datestr(datetime,'yyyymmddTHHMMSS','local') '.mat'];
save(fullfile(pmRootPath,'local',vistasynthDTfName), 'synthDT');

% Analyze it, it takes a lot of time
results    = pmModelFit(synthDT, 'mrvista','model','one gaussian', ...
                                        'grid', false, ... % if true, returns gFit
                                        'wSearch', 'coarse to fine');
% Save the result, we don't want to lose it because matlab crashes...
% Use the same datetime to associate synthDT and result
vistaresultfName = ['results_vista_oneCenter_' datestr(datetime,'yyyymmddTHHMMSS','local') '.mat'];
% if exist(fullfile(pmRootPath,'local',resultfName),'file')
%     error('Change file name or delete existing file.')
    % delete(fullfile(pmRootPath,'local',resultfName))
% else
    save(fullfile(pmRootPath,'local',vistaresultfName), 'results');
% end










% DO AFNI
COMBINE_PARAMETERS                    = struct();
COMBINE_PARAMETERS.RF.Centerx0        = [0];  % [-6 5 4 3 2 1 0 1 2 3 4 5 6];
COMBINE_PARAMETERS.RF.Centery0        = [0];  % [-6 5 4 3 2 1 0 1 2 3 4 5 6];
COMBINE_PARAMETERS.RF.Theta           = [0]; %, deg2rad(45)];
COMBINE_PARAMETERS.RF.sigmaMajor      = [0.5,1,2,4,8];  % [1,2,3,4];
COMBINE_PARAMETERS.RF.sigmaMinor      = 'same'; % 'same' for making it the same to Major
COMBINE_PARAMETERS.TR                 = [1.5];
    HRF(1).Type                       = 'afni_spm';
COMBINE_PARAMETERS.HRF                = HRF;
% TODO: implement a more complex noise addition system. 
% Right now only the parameter for white noise can be edited. 
COMBINE_PARAMETERS.Noise.noise2signal = [0:0.05:0.45];
synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS, 'mult',100);
synthDT = pmForwardModelCalculate(synthDT);
% Save it
afnisynthDTfName = ['synthDT_afnispm_oneCenter_' datestr(datetime,'yyyymmddTHHMMSS','local') '.mat'];
save(fullfile(pmRootPath,'local',afnisynthDTfName), 'synthDT');

% Analyze it, it takes a lot of time
results    = pmModelFit(synthDT, 'afni_4');
% Save the result, we don't want to lose it because matlab crashes...
% Use the same datetime to associate synthDT and result
afniresultfName = ['results_afni4_oneCenter_' datestr(datetime,'yyyymmddTHHMMSS','local') '.mat'];
% if exist(fullfile(pmRootPath,'local',resultfName),'file')
%    error('Change file name or delete existing file.')
    % delete(fullfile(pmRootPath,'local',resultfName))
% else
    save(fullfile(pmRootPath,'local',afniresultfName), 'results');
% end




% {

%% Identify the results we want to use manually
prfSynth   = load(fullfile(pmRootPath,'local','synthDT_aprf_oneCenter_20190825T212722.mat'));
prfRes     = load(fullfile(pmRootPath,'local','result_aprf_oneCenter_20190825T221545.mat'));

vistaSynth = load(fullfile(pmRootPath,'local','synthDT_vista_oneCenter_20190826T015219.mat'));
vistaRes   = load(fullfile(pmRootPath,'local','results_vista_oneCenter_20190826T020241.mat'));

afni4Synth = load(fullfile(pmRootPath,'local','synthDT_afnispm_oneCenter_20190826T023028.mat'));
afni4Res   = load(fullfile(pmRootPath,'local','results_afni4_oneCenter_20190826T094209.mat'));

popSynth   = load(fullfile(pmRootPath,'local','synthDT_pop_oneCenter_20190825T224055.mat'));
popRes     = load(fullfile(pmRootPath,'local','results_pop_oneCenter_20190826T012732.mat'));


%% (single tool) Plot the results

paramDefaults = {'Centerx0','Centery0','Theta','sigmaMinor','sigmaMajor'};
[compTable, tSeries] = pmResultsCompare(popSynth.synthDT, ... % Defines the input params
                            {'tool'}, ... % Analysis names we want to see: 'aPRF','vista',
                            {popRes.results}, ... % results_analyzePRF,results_vista,
                            'params', paramDefaults, ...
                            'shorten names',true); 
                        
pmTseriesPlot(tSeries, prfSynth.synthDT(:,'TR'), ...
                            'to compare', {'synth','tool'}, ...
                            'voxel',[1], ... 
                            'newWin',true)

% 2 rows, top row accuracy, bottom precision
% each col will be different rf sizes
% each plot will have noise in the x

%% (single tool) ACCURACY AND PRECISION PLOT: by noise levels
prfsizes   = unique(compTable.synth.sMaj);
mrvNewGraphWin([useTool ' accuracy and precision']);

for np=1:length(prfsizes)
    subplot(length(prfsizes), 1, np);
    prfsize      = prfsizes(np);
    % Reduce the table for only the RF size-s we want
    DT           = compTable(compTable.synth.sMaj==prfsize,:);
    % Call the function for the distribution plot
    % subplot(length(prfsizes), 1, np);
    [result,noiseVals,sds] = pmNoise2AccPrec(DT, 'both','plotIt',false);
    % errorbar(noiseVals,result, sds,'Color','b','LineStyle','-','LineWidth',2);hold on;
    plot(noiseVals,result,'Color','b','LineStyle','-','LineWidth',2);hold on;
    jbfill(noiseVals,result+sds,result-sds,'b','b',1,.5);hold on;
    h1 = plot([noiseVals(1),noiseVals(end)],prfsize*[1,1],'Color','k','LineStyle','-.','LineWidth',2);
    xticks([])
    set(gca,'box','off','color','none')
    ylabel('Acc. & prec. (deg)','FontSize',14,'FontWeight','bold','Color','k');
    % title(['ACCURACY AND PRECISION (vs noise and RF size)']);
    % legend(num2str(prfsizes),'Location','northwest');
end
xticks(noiseVals)
xlabel('Noise levels','FontSize',18,'FontWeight','bold','Color','k');

%% (multiple tool) Plot the results
%{
paramDefaults        = {'Centerx0','Centery0','Theta','sigmaMinor','sigmaMajor'};
[compTable, tSeries] = pmResultsCompare(prfSynth.synthDT, ... % Defines the input params
                            {'aprf','vista','afni4','pop'}, ... % Analysis names we want to see: 'aPRF','vista',
                            {prfRes.results, vistaRes.results, afni4Res.results, popRes.results}, ... % results_analyzePRF,results_vista,
                            'params', paramDefaults, ...
                            'shorten names',true); 
                        
pmTseriesPlot(tSeries, prfSynth.synthDT(:,'TR'), ...
                            'to compare', {'synth','aprf','vista','afni4','pop'}, ...
                            'voxel',[1:9], ... 
                            'newWin',true);
%}

% NO POPEYE
paramDefaults        = {'Centerx0','Centery0','Theta','sigmaMinor','sigmaMajor'};
[compTable, tSeries] = pmResultsCompare(prfSynth.synthDT, ... % Defines the input params
                            {'aprf','vista','afni4'}, ... % Analysis names we want to see: 'aPRF','vista',
                            {prfRes.results, vistaRes.results, afni4Res.results, popRes.results}, ... % results_analyzePRF,results_vista,
                            'params', paramDefaults, ...
                            'shorten names',true); 
                        
pmTseriesPlot(tSeries, prfSynth.synthDT(:,'TR'), ...
                            'to compare', {'synth','aprf','vista','afni4'}, ...
                            'voxel',[1:9], ... 
                            'newWin',true);

% 2 rows, top row accuracy, bottom precision
% each col will be different rf sizes
% each plot will have noise in the x

%% (multiple tool) ACCURACY AND PRECISION PLOT: by noise levels
prfsizes = unique(compTable.synth.sMaj);
% tools    = {'aprf','vista','afni4','pop'};
tools    = {'aprf','vista','afni4'};
colors   = distinguishable_colors(length(tools),'w');

mrvNewGraphWin(['Accuracy and precision']);

for np=1:length(prfsizes)
    subplot(length(prfsizes), 1, np);
    prfsize      = prfsizes(np);
   a = [];b = [];
   for nt=1:length(tools)
       tool = tools{nt};
       % Reduce the table for only the RF size-s we want, and the tool we want
       DT           = compTable(compTable.synth.sMaj==prfsize,:);
       % Call the function for the distribution plot
       % subplot(length(prfsizes), 1, np);
       [result,noiseVals,sds] = pmNoise2AccPrec(DT, 'both','plotIt',false,'tool',tool);
       % errorbar(noiseVals,result, sds,'Color','b','LineStyle','-','LineWidth',2);hold on;
       % a = [a;plot(noiseVals,result,'Color',colors(nt,:),'LineStyle',':','LineWidth',.5)];
       a = [a;plot(noiseVals,result+sds,'Color',colors(nt,:),'LineStyle','-','LineWidth',2)];hold on;
       b = [b;plot(noiseVals,result-sds,'Color',colors(nt,:),'LineStyle','-','LineWidth',2)];
       
       jbfill(noiseVals,result+sds,result-sds,colors(nt,:),colors(nt,:),1,.05); hold on;
    end
    h1 = plot([noiseVals(1),noiseVals(end)],prfsize*[1,1],'Color','k','LineStyle','-.','LineWidth',1);
    grid
    xticks([])
    set(gca,'box','off','color','none')
    ylabel('Acc. & prec. (deg)','FontSize',14,'FontWeight','bold','Color','k');
    % title(['ACCURACY AND PRECISION (vs noise and RF size)']);
    legend([h1;a],['synth',tools],'Location','southwest');
    
end
xticks(noiseVals)
% set(gca,'TickLength',[0 0])
ax = gca;
ax.XGrid = 'off';
xlabel('Noise levels','FontSize',18,'FontWeight','bold','Color','k');
x0    = unique(compTable.synth.x0);
y0    = unique(compTable.synth.x0);
Theta = unique(compTable.synth.Th);
TR    = unique(compTable.TR);
dr_suptitle(sprintf('Ideal HRF | Location=[%i,%i] | Theta=%i | TR=%1.2f | Band= 1 SD', x0, y0, Theta, TR));


%


%}

%% RANDOM HRFs
% DO RANDOM HRF
COMBINE_PARAMETERS                    = struct();
COMBINE_PARAMETERS.RF.Centerx0        = [0];  % [-6 5 4 3 2 1 0 1 2 3 4 5 6];
COMBINE_PARAMETERS.RF.Centery0        = [0];  % [-6 5 4 3 2 1 0 1 2 3 4 5 6];
COMBINE_PARAMETERS.RF.Theta           = [0]; %, deg2rad(45)];
COMBINE_PARAMETERS.RF.sigmaMajor      = [0.5,1,2,4,8];  % [1,2,3,4];
COMBINE_PARAMETERS.RF.sigmaMinor      = 'same'; % 'same' for making it the same to Major
COMBINE_PARAMETERS.TR                 = [1.5];
    HRF(1).Type                       = 'random';
COMBINE_PARAMETERS.HRF                = HRF;
% TODO: implement a more complex noise addition system. 
% Right now only the parameter for white noise can be edited. 
COMBINE_PARAMETERS.Noise.noise2signal = [0:0.05:0.45];
synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS, 'mult',100);
synthDT = pmForwardModelCalculate(synthDT);
% Save it
randomSynthDTfName = ['synthDT_RANDOM_HRF_oneCenter_' datestr(datetime,'yyyymmddTHHMMSS','local') '.mat'];
save(fullfile(pmRootPath,'local',randomSynthDTfName), 'synthDT');





% AnalyzePRF
results    = pmModelFit(synthDT, 'aprf', 'useParallel', true);
aprfresultfName = ['result_rnd_aprf_oneCenter_' datestr(datetime,'yyyymmddTHHMMSS','local') '.mat'];
save(fullfile(pmRootPath,'local',aprfresultfName), 'results');
% popeye
results    = pmModelFit(synthDT, 'popeye_onegaussian');
popresultfName = ['result_rnd_pop_oneCenter_' datestr(datetime,'yyyymmddTHHMMSS','local') '.mat'];
save(fullfile(pmRootPath,'local',popresultfName), 'results');
% mrvista
results    = pmModelFit(synthDT, 'mrvista','model','one gaussian', ...
                                        'grid', false, ... % if true, returns gFit
                                        'wSearch', 'coarse to fine');
vistaresultfName = ['results_rnd_vista_oneCenter_' datestr(datetime,'yyyymmddTHHMMSS','local') '.mat'];
save(fullfile(pmRootPath,'local',vistaresultfName), 'results');
% Afni
results    = pmModelFit(synthDT, 'afni_4');
afniresultfName = ['results_rnd_afni_oneCenter_' datestr(datetime,'yyyymmddTHHMMSS','local') '.mat'];
save(fullfile(pmRootPath,'local',afniresultfName), 'results');




%% TODO
% Averiguar que le pasa a popeye
% Conseguir correr todo en black
% Hacer los plots independientemente como antes, o sea, dos.
% Ver la diferencia entre los tools entre usar ideal HRF o random.
% Kendrick:
%    pedir HRFs
%    comentarle lo de los malos resultados para big RFs
% Write
% Noah: OSF and Dockers
% Noise models: more tests
% Location tests?
% HRF: test with params? enough if we test with DDBB and randoms






























%{


%% ACCURACY AND PRECISION PLOT: by prfsize levels
noiseVals    = unique(compTable.noise2sig); 
rednoiseVals = noiseVals(2) : 0.1 : noiseVals(10);
mrvNewGraphWin([useTool ' accuracy and precision']);

for np = 1:length(rednoiseVals)
    subplot(length(rednoiseVals), 1, np);
    noise      = rednoiseVals(np);
    % Reduce the table for only the RF size-s we want
    DT         = compTable(compTable.noise2sig==noise,:);
    % Call the function for the distribution plot
    % subplot(length(prfsizes), 1, np);
    [result,noiseVals,sds] = pmNoise2AccPrec(DT, 'both','plotIt',false);
    % errorbar(noiseVals,result, sds,'Color','b','LineStyle','-','LineWidth',2);hold on;
    plot(noiseVals,result,'Color','b','LineStyle','-','LineWidth',2);hold on;
    jbfill(noiseVals,result+sds,result-sds,'b','b',1,.5);hold on;
    h1 = plot([noiseVals(1),noiseVals(end)],prfsize*[1,1],'Color','k','LineStyle','-.','LineWidth',2);
    xlabel('Noise levels');
    ylabel('Acc. & prec. (deg)');
    % title(['ACCURACY AND PRECISION (vs noise and RF size)']);
    % legend(num2str(prfsizes),'Location','northwest');
end


%% DISTRIBUTION PLOTS
prfsizes   = unique(compTable.synth.sMaj);
mrvNewGraphWin([useTool ' comparisons']);

for np=1:length(prfsizes)
    prfsize      = prfsizes(np);
    % Reduce the table for only the RF size-s we want
    DT           = compTable(compTable.synth.sMaj==prfsize,:);
    % Select the variable to compare the distributions with
    varToCompare = 'noise2sig';
    % Call the function for the distribution plot
    subplot(length(prfsizes), 1, np);
    pmAccPrecDistributions(DT, varToCompare,'normalksdensity','ksdensity');
    % Add text to the plots
    xlabel('Prediction of RF size (deg)');
    ylabel('Count');
    title(['Encoded RF size:' num2str(prfsize) ' deg']);
%     xlim([-3,10])
%     xticks(-3:1:10)
    xlim([-0.5,3])
    xticks(-0.5:0.5:3)
end
 


%% ACCURACY AND PRECISION PLOTS
prfsizes   = unique(compTable.synth.sMaj);
mrvNewGraphWin([useTool ' accuracy and precision']);

% Do precision
subplot(1,2,1)
for np=1:length(prfsizes)
    prfsize      = prfsizes(np);
    % Reduce the table for only the RF size-s we want
    DT           = compTable(compTable.synth.sMaj==prfsize,:);
    % Call the function for the distribution plot
    % subplot(length(prfsizes), 1, np);
    [result,noiseVals] = pmNoise2AccPrec(DT, 'precision','plotIt',false);
    plot(noiseVals,result);hold on;
end
xlabel('Noise levels'); 
ylabel('SD of the results (deg)'); 
title(['PRECISION for different noise and RF size']);
legend(num2str(prfsizes),'Location','northwest');

% Do accuracy
subplot(1,2,2)
for np=1:length(prfsizes)
    prfsize      = prfsizes(np);
    % Reduce the table for only the RF size-s we want
    DT           = compTable(compTable.synth.sMaj==prfsize,:);
    % Call the function for the distribution plot
    % subplot(length(prfsizes), 1, np);
    [result,noiseVals] = pmNoise2AccPrec(DT, 'accuracy','plotIt',false);
    plot(noiseVals,result);hold on;
end
xlabel('Noise levels'); 
ylabel('Abs difference from the mean (deg)'); 
title(['ACCURACY for different noise and RF size']);
legend(num2str(prfsizes),'Location','northwest');



%% TO DELETE
%     for ns=1:length(noise2sigs)
%         noise2sig = noise2sigs(ns);
%         
%         % Now we can create the subplots
%         
%         subplot(length(noise2sigs),length(prfsizes),plotIndex);
%         X = compTable.afni.sMin(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
%         Y = compTable.afni.sMaj(compTable.noise2sig==noise2sig & compTable.synth.sMaj==prfsize);
%         scatter(X,Y)
%         axis equal;
%         xlabel('sMin'); ylabel('sMaj')
%         title(sprintf('rfsize:%1.1f | noise2sig:%0.2f',prfsize,noise2sig))
%         grid on;
%         xlim([-0,8]); ylim([-0,8])
%         xticks([0:1:8]);yticks([0:1:8])
%         identityLine(gca);
%         % h1 = lsline;
%         % h1.Color = 'r';
%         hold on;
%         % [x0,y0] = fitEllipse(X,Y,'r');
%         x0  = median(X); y0 = median(Y);
%         sdx = std(X); sdy = std(Y);
%         plot([x0-sdx, x0+sdx],[y0 y0],'r');
%         plot([x0 x0],[y0-sdy, y0+sdy],'r');
%         text(2,0.5,sprintf('Median:[%1.2f(%1.2f),%1.2f(%1.2f)], ratio:%1.2f',x0,sdx,y0,sdy,y0/x0));
%     end
% end

%}
