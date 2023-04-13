%% s_pmHRFtests.m



% DO aPRF
COMBINE_PARAMETERS                    = struct();
COMBINE_PARAMETERS.RF.Centerx0        = [0,5,-5];
COMBINE_PARAMETERS.RF.Centery0        = [0,5,-5];  % [-6 5 4 3 2 1 0 1 2 3 4 5 6];
COMBINE_PARAMETERS.RF.Theta           = [0]; %, deg2rad(45)];
COMBINE_PARAMETERS.RF.sigmaMajor      = [2]; % [0.5,1,2,4,8];  % [1,2,3,4];
COMBINE_PARAMETERS.RF.sigmaMinor      = 'same'; % 'same' for making it the same to Major
COMBINE_PARAMETERS.TR                 = [1.5];
    HRF(1).Type                       = 'friston';
    % HRF(2).Type                       = 'boynton';
    % HRF(3).Type                       = 'canonical';
    % HRF(4).Type                       = 'vista_twogammas';
    % HRF(5).Type                       = 'popeye_twogammas';
    % HRF(6).Type                       = 'afni_gam';
     %HRF(7).Type                       = 'afni_spm';
COMBINE_PARAMETERS.HRF                = HRF;
% Right now only the parameter for white noise can be edited. 
% COMBINE_PARAMETERS.Noise.noise2signal = [0:0.05:0.45];
synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS, 'repeats',5);
synthDT = pmForwardModelCalculate(synthDT);
% Save it
aprfsynthDTfName = ['synthDT_ALL_' datestr(datetime,'yyyymmddTHHMMSS','local') '.mat'];
save(fullfile(pmRootPath,'local',aprfsynthDTfName), 'synthDT');





% Now we need to solve per each tool
% AnalyzePRF
aprfresults         = pmModelFit(synthDT, 'aprf', 'useParallel', true);
aprfresultfName = ['result_ALL_aprf_oneCenter_' datestr(datetime,'yyyymmddTHHMMSS','local') '.mat'];
save(fullfile(pmRootPath,'local',aprfresultfName), 'aprfresults');

% mrvista
options.vista = struct('model'        ,'one gaussian'   , ...
                       'grid'         , false           , ...
                       'wSearch'      , 'coarse to fine', ...
                       'detrend'      , 1               , ...
                       'keepAllPoints', false           , ...
                       'obtainPreds'  , false           , ...
                       'fixcssexp'    , 0               , ...
                       'numberStimulusGridPoints',   50);
vistaresults    = pmModelFit(synthDT, 'mrvista', options);
vistaresultfName = ['results_ALL_vista_oneCenter_' datestr(datetime,'yyyymmddTHHMMSS','local') '.mat'];
save(fullfile(pmRootPath,'local',vistaresultfName), 'vistaresults');

% Afni
results    = pmModelFit(synthDT, 'afni_4');
afniresultfName = ['results_ALL_afni_oneCenter_' datestr(datetime,'yyyymmddTHHMMSS','local') '.mat'];
save(fullfile(pmRootPath,'local',afniresultfName), 'results');
% popeye
results    = pmModelFit(synthDT, 'popeye_onegaussian');
popresultfName = ['result_ALL_pop_oneCenter_' datestr(datetime,'yyyymmddTHHMMSS','local') '.mat'];
save(fullfile(pmRootPath,'local',popresultfName), 'results');







% 







%% Copy from below if needed, delete all at the end
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



%% IDEAL HRF: Create dataset for testing noiseless and noise values in all datasets. 

%{
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
%}

% Now that they have been generated, and saved locally, load them to FW. 
% Next, we will check if they exist locally, otherwise download from FW.
%  ...and upload it to the collection


% THESE WILL BE 
allFileNames = {
fullfile(pmRootPath,'local','synthDT_aprf_oneCenter_20190825T212722.mat')
fullfile(pmRootPath,'local','result_aprf_oneCenter_20190825T221545.mat')
fullfile(pmRootPath,'local','synthDT_vista_oneCenter_20190826T015219.mat')
fullfile(pmRootPath,'local','results_vista_oneCenter_20190826T020241.mat')
fullfile(pmRootPath,'local','synthDT_afnispm_oneCenter_20190826T023028.mat')
fullfile(pmRootPath,'local','results_afni4_oneCenter_20190826T094209.mat')
fullfile(pmRootPath,'local','synthDT_pop_oneCenter_20190825T224055.mat')
fullfile(pmRootPath,'local','results_pop_oneCenter_20190826T012732.mat')
};

% UPLOAD LOCAL FILES
%{
st   = scitran('stanfordlabs'); st.verify;
cc   = st.search('collection','collection label exact','PRF_StimDependence');

for nf=1:length(allFileNames)
    localfname = allFileNames{nf};
    stts = st.fileUpload(localfname, cc{1}.collection.id, 'collection');
    % Check that the data is there
    [~,fname,ext] = fileparts(localfname);
    data          = load(st.fw.downloadFileFromCollection(cc{1}.collection.id,...
        [fname ext],localfname));
end
%}





% Identify the results we want to use manually
% If the file exists locally, read it, otherwise download from the server locally
data = {};
for nf=1:length(allFileNames)
    localfname = allFileNames{nf};
    if exist(localfname,'file')
        data{nf}          = load(localfname);
    else
        % Read the data
        [~,fname,ext] = fileparts(localfname);
        data{nf}          = load(st.fw.downloadFileFromCollection(cc{1}.collection.id,...
            [fname ext],localfname));
    end
end

prfSynth   = load(fullfile(pmRootPath,'local','synthDT_aprf_oneCenter_20190825T212722.mat'));    % 1
prfRes     = load(fullfile(pmRootPath,'local','result_aprf_oneCenter_20190825T221545.mat'));     % 2 

vistaSynth = load(fullfile(pmRootPath,'local','synthDT_vista_oneCenter_20190826T015219.mat'));   % 3
vistaRes   = load(fullfile(pmRootPath,'local','results_vista_oneCenter_20190826T020241.mat'));   % 4

afni4Synth = load(fullfile(pmRootPath,'local','synthDT_afnispm_oneCenter_20190826T023028.mat')); % 5
afni4Res   = load(fullfile(pmRootPath,'local','results_afni4_oneCenter_20190826T094209.mat'));   % 6

popSynth   = load(fullfile(pmRootPath,'local','synthDT_pop_oneCenter_20190825T224055.mat'));     % 7
popRes     = load(fullfile(pmRootPath,'local','results_pop_oneCenter_20190826T012732.mat'));     % 8


% Plot the results ACCURACY AND PRECISION PLOT: by noise levels
paramDefaults        = {'Centerx0','Centery0','Theta','sigmaMinor','sigmaMajor'};
[compTable, tSeries] = pmResultsCompare(data{1}.synthDT, ... % Defines the input params
                            {'aprf','vista','afni4','pop'}, ... % Analysis names we want to see: 'aPRF','vista',
                            {data{2}.results, data{4}.results, data{6}.results, data{8}.results}, ... % results_analyzePRF,results_vista,
                            'params', paramDefaults, ...
                            'shorten names',true); 
                        
% pmTseriesPlot(tSeries, prfSynth.synthDT(:,'TR'), 'voxel',[1:9], 'newWin',true, ...
%                             'to compare', {'synth','aprf','vista','afni4','pop'});


% Plot the tool comparison plots for noise and rfSize ACCURACY AND PRECISION PLOT: by noise levels
pmNoiseSizePlots(compTable, {'aprf','vista','afni4','pop'})

% CHECK the results, CLEAN, and PLOT again
% Well, the thing is that at least aPRF and popeye have several outliers that should be cleaned. 
% So, this is now a outlier cleaning tool, that will output:
%   - Cleaned result table. 
%   - A table with outlier counts, that will evaluate quality of analysis 
[newDT,outlierReport] = pmResultsCleanOutliers(compTable,'outlierLevel',20);

% Plot again, cleaned dataset
pmNoiseSizePlots(newDT, {'aprf','vista','afni4','pop'},'randomHRF',false)



%% RANDOM HRFs
%{
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




%RUN

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
%}




allFileNames = {
fullfile(pmRootPath,'local','synthDT_RANDOM_HRF_oneCenter_20190826T100808.mat')
fullfile(pmRootPath,'local','result_rnd_aprf_oneCenter_20190826T110119.mat')
fullfile(pmRootPath,'local','results_rnd_vista_oneCenter_20190826T124614.mat')
fullfile(pmRootPath,'local','results_rnd_afni_oneCenter_20190826T185712.mat')
fullfile(pmRootPath,'local','result_rnd_pop_oneCenter_20190826T123508.mat')
};

% UPLOAD LOCAL FILES
%{
st   = scitran('stanfordlabs'); st.verify;
cc   = st.search('collection','collection label exact','PRF_StimDependence');

for nf=1:length(allFileNames)
    localfname = allFileNames{nf};
    stts = st.fileUpload(localfname, cc{1}.collection.id, 'collection');
    % Check that the data is there
    [~,fname,ext] = fileparts(localfname);
    data          = load(st.fw.downloadFileFromCollection(cc{1}.collection.id,...
        [fname ext],localfname));
end
%}


% Identify the results we want to use manually
% If the file exists locally, read it, otherwise download from the server locally
data = {};
for nf=1:length(allFileNames)
    localfname = allFileNames{nf};
    if exist(localfname,'file')
        data{nf}          = load(localfname);
    else
        % Read the data
        [~,fname,ext] = fileparts(localfname);
        data{nf}          = load(st.fw.downloadFileFromCollection(cc{1}.collection.id,...
            [fname ext],localfname));
    end
end



% Identify and read the results we want to use manually
%  TODO: Backup them in FW, run in the server and read them here
rand_Synth      = load(fullfile(pmRootPath,'local','synthDT_RANDOM_HRF_oneCenter_20190826T100808.mat'));

rand_prfRes     = load(fullfile(pmRootPath,'local','result_rnd_aprf_oneCenter_20190826T110119.mat'));
rand_vistaRes   = load(fullfile(pmRootPath,'local','results_rnd_vista_oneCenter_20190826T124614.mat'));
rand_afni4Res   = load(fullfile(pmRootPath,'local','results_rnd_afni_oneCenter_20190826T185712.mat'));
rand_popRes     = load(fullfile(pmRootPath,'local','result_rnd_pop_oneCenter_20190826T123508.mat'));

% Plot the results
paramDefaults                  = {'Centerx0','Centery0','Theta','sigmaMinor','sigmaMajor'};
[rand_compTable, rand_tSeries] = pmResultsCompare(data{1}.synthDT, ... % Defines the input params
                            {'aprf','vista','afni4','pop'}, ... % Analysis names we want to see: 'aPRF','vista',
                            {data{2}.results, data{3}.results, data{4}.results, data{5}.results}, ... % results_analyzePRF,results_vista,
                            'params', paramDefaults, ...
                            'shorten names',true); 
% Plot some example voxels to see time series                        
%pmTseriesPlot(tSeries, rand_Synth.synthDT(:,'TR'), 'voxel',[1:9], 'newWin',true, ...
%                            'to compare', {'synth','aprf','vista','afni4','pop'});
% Plot the tool comparison plots for noise and rfSize ACCURACY AND PRECISION PLOT: by noise levels
pmNoiseSizePlots(rand_compTable, {'aprf','vista','afni4','pop'},'randomHRF',true)


% CHECK the results, CLEAN, and PLOT again
% Well, the thing is that at least aPRF and popeye have several outliers that should be cleaned. 
% So, this is now a outlier cleaning tool, that will output:
%   - Cleaned result table. 
%   - A table with outlier counts, that will evaluate quality of analysis 
[rand_newDT,rand_outlierReport] = pmResultsCleanOutliers(rand_compTable,'outlierLevel',20);

% Plot again, cleaned dataset
pmNoiseSizePlots(rand_newDT, {'aprf','vista','afni4','pop'},'randomHRF',true)


%% Difference between ideal and HRF
% meterlo en pmNoiseSizePlots, DT sera cell, if cell comparalos, si no hazlo uno a uno
%     comparacion 

pmNoiseSizePlots({newDT,rand_newDT}, {'aprf','vista','afni4','pop'})

%% TODO
% (done) Investigate what happens to popeye and aPRF
%    --- Outliers!
% Run everything in black (or GCP)
%     (done) Docker is still failing after the reboot: reinstalled and rebooted
%     (done) Matlab with python crashes, segfault as Ariel said: hardcode the
%            popeye HRFs in Linux
%     ()     Actually run it
% (done) Plot the ideal and random hrf plots
% (ongoing) Plot difference between tools when using ideal HRF or random.
% () Kendrick:
%    () ask for HRF database?
%    () comment bad results big HRFs
% () Write manuscript
% () Noah: OSF and Dockers
% () Noise models: more tests
% () Location tests?
% () HRF: test with params? enough if we test with DDBB and randoms






























%% OLD CODE
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



%}