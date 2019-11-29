%% s_pmPaperResultsFigures.m
% First the nifti with the data was created using prf-synthesize
% The json params files used is in scripts/params_big_test_for_paper_v01.json

% Create the names that will be used in this script

close all;clear all;clc
tbUse prfmodel;

% Control execution of the script file
defineFileNames          = true;
generateInputFiles       = false;
uploadInputFiles         = false;
analysisInputFiles       = false;
downloadLoadFiles        = true;
createComptTables        = true;
plotGroupAnalysis        = false;
plotCloudPoints          = false;
plotHRFwidthtestsREALS   = false;
plotHRFwidthtestsBOYNTON = true;
plotHRFwidthtestsVISTA   = false;

%% defineFileNames
if defineFileNames
    % INPUTS
    inputJSON    = fullfile(pmRootPath, 'local/output/paper', ...
                           'params_big_test_for_paper_3-3_v02BOLD.json');
    outputFolder = fullfile(pmRootPath, 'local/output/paper');
    
    niftiBOLDfile  = fullfile(outputFolder, ...
        'BIDS/sub-paperResults/ses-prfsynthBOLDx3y3V02/func',...
        'sub-paperResults_ses-prfsynthBOLDx3y3V02_task-prf_acq-normal_run-01_bold.nii.gz');
    jsonSynthFile  = fullfile(outputFolder, ...
        'BIDS/derivatives/prfsynth/sub-paperResults/ses-prfsynthBOLDx3y3V02',...
        'sub-paperResults_ses-prfsynthBOLDx3y3V02_task-prf_acq-normal_run-01_bold.json');
    stimNiftiFname = fullfile(outputFolder, ...
        'BIDS','stimuli',...
        'sub-paperResults_ses-prfsynthBOLDx3y3V02_task-prf_apertures.nii.gz');
    inputNiftis = {niftiBOLDfile, jsonSynthFile, stimNiftiFname};
    % OUTPUTS
    aprfcssresultfName     = fullfile(pmRootPath,'local',['paper04_result_aprfcss.mat']);
    aprfresultfName        = fullfile(pmRootPath,'local',['paper04_result_aprf.mat']);
    vistaresultfName       = fullfile(pmRootPath,'local',['paper04_result_vista.mat']); % 2 is BOLD, 1 is contrast
    vistaandhrfresultfName = fullfile(pmRootPath,'local',['paper04_result_vistaandhrf.mat']);
    popresultfName         = fullfile(pmRootPath,'local',['paper04_result_pop.mat']);
    popnohrfresultfName    = fullfile(pmRootPath,'local',['paper04_result_popnohrf.mat']);
    afni4resultfName       = fullfile(pmRootPath,'local',['paper04_result_afni4.mat']);
    afni6resultfName       = fullfile(pmRootPath,'local',['paper04_result_afni6.mat']);
    
    outputNiftis = {aprfcssresultfName,aprfresultfName, ...
        vistaresultfName,vistaandhrfresultfName, ...
        popresultfName, popnohrfresultfName, ...
        afni4resultfName, afni6resultfName};
    % Connect to fw
    st   = scitran('stanfordlabs'); st.verify;
    cc   = st.search('collection','collection label exact','PRF_StimDependence');
    
    % Save files to
    saveTo = '/Users/glerma/gDrive/STANFORD/PROJECTS/2019_PRF_Validation_methods_(Gari)/__PUBLISH__/PAPER_SUBMISSION01/Figures/RAW';
end

%% generateInputFiles
if generateInputFiles
    cmd = [fullfile(pmRootPath,'gear','prfsynth','run_prfsynth.sh ') ...
           inputJSON ' ' outputFolder];
    system(cmd);
end

%% uploadInputFiles
if (uploadInputFiles)
    for ii=1:length(inputNiftis)
        stts = st.fileUpload(inputNiftis{ii}, cc{1}.collection.id, 'collection');
    end
end

%% analysisInputFiles
if (analysisInputFiles)
% SOLVE
%  - aPRF CSS=true
disp('START aPRF CSS=true ... ')
options = struct('seedmode',[0,1,2], 'display','off', 'usecss',true);
results = pmModelFit(inputNiftis, 'aprf', 'options', options);
save(aprfcssresultfName, 'results'); clear results
disp(' ... END aPRF CSS=true')

%  - aPRF CSS=false
disp('START aPRF CSS=false ...')
options = struct('seedmode',[0,1,2], 'display','off', 'usecss',false);
results = pmModelFit(inputNiftis, 'aprf', 'options', options);
save(aprfresultfName, 'results'); clear results
disp(' ... END aPRF CSS=false')
% {
%  - mrVista
disp('START mrVista ... ')
results = pmModelFit(inputNiftis, 'mrvista','model','one gaussian', ...
                  'grid', false, 'wSearch', 'coarse to fine');
save(vistaresultfName, 'results'); clear results
disp(' ... END mrVista')

%  - mrvista and hrf
disp('START mrvista and hrf ... ')
results = pmModelFit(inputNiftis, 'mrvista','model','one gaussian', ...
                    'grid', false, 'wSearch', 'coarse to fine and hrf');
save(vistaandhrfresultfName, 'results'); clear results
disp(' ... END mrvista and hrf')
%}
%  - popeye, default with hrf search
disp('START popeye, default with hrf search ... ')
results = pmModelFit(inputNiftis, 'popeye'); 
save(popresultfName, 'results'); clear results
disp(' ... END popeye, default with hrf search')

%  - popeye no hrf
disp('START popeye no hrf ... ')
results= pmModelFit(inputNiftis, 'popeyenohrf'); 
save(popnohrfresultfName, 'results'); clear results
disp(' ... END popeye no hrf')

%  - Afni 4
disp('START Afni 4 ... ')
results = pmModelFit(inputNiftis, 'afni_4','afni_hrf','SPM'); 
save(afni4resultfName, 'results'); clear results
disp(' ... END Afni 4')

%  - Afni 6
disp('START Afni 6 ... ')
results = pmModelFit(inputNiftis, 'afni_6','afni_hrf','SPM'); 
save(afni6resultfName, 'results'); clear results
disp(' ... END Afni 6')

% UPLOAD
for nu=1:length(outputNiftis)
    resultfName = outputNiftis{nu}
    stts      = st.fileUpload(resultfName, cc{1}.collection.id, 'collection');
end
end

%% downloadLoadFiles
if (downloadLoadFiles)
    % (Download and) read the input nifti-s (only the json for now)
    for ni=1:length(inputNiftis)
        [fpath,fname,ext] = fileparts(inputNiftis{ni});
        if ~exist(inputNiftis{ni},'file')
            mkdir(fpath)
            st.fw.downloadFileFromCollection(cc{1}.collection.id,[fname ext],inputNiftis{ni});
        end
        if strcmp(ext, '.json')
            SynthDT  = struct2table(jsonread(inputNiftis{ni}));
            for na=1:width(SynthDT),if isstruct(SynthDT{:,na})
                    SynthDT.(SynthDT.Properties.VariableNames{na}) = struct2table(SynthDT{:,na});
                end,end
        end
    end
    
    % (Download and) read the output nifti-s
    res = struct();
    for outputNifti=outputNiftis
        rname = outputNifti{:}(1:end-4);
        rname = split(rname,'_');
        rname = rname{end,:};
        if exist(outputNifti{:},'file')
            load(outputNifti{:});
        else
            % TODO Check that the data is there before trying to download
            [~,fname,ext] = fileparts(outputNifti{:});
            load(st.fw.downloadFileFromCollection(cc{1}.collection.id,...
                [fname ext],outputNifti{:}));
        end
        res.(rname) = results;
    end
end

%% createCompTables
if createComptTables
    paramDefaults = {'Centerx0','Centery0','Theta','sigmaMinor','sigmaMajor'};
    % Create the concatenated result table
    
    % Select analysis names we want to see: 'aPRF','vista',... by default, all of them
    %     {'aprf','aprfcss','vista','vistahrf','pop','popno','afni4','afni6'}
    anNames = fieldnames(res);
    % Select the result files
    ress = [{res.(anNames{1})}];
    for an =2:length(anNames)
        ress = [ress, {res.(anNames{an})}];
    end
    compTable  = pmResultsCompare(SynthDT, anNames, ress, ...
        'params', paramDefaults, ...
        'shorten names',true, ...
        'dotSeries', false);
    % Now create the new plots
    sortHRFlike = {'friston','afni_gam','boynton','afni_spm',...
        'popeye_twogammas','canonical'};  % sorted according noise=0, RFsize=2
    
end

%% plotGroupAnalysis
if (plotGroupAnalysis) 
    % PLOTS
    tools = {'vista','afni4','popnohrf','aprf'}';
    % Or select them all
    % tools = anNames;
    
    
    % After seeing the plots we decided to do noisefree, 3 sizes
%     pmNoisePlotsByHRF(compTable, tools, 'x0y0',[3,3],...
%         'sortHRF',sortHRFlike,'usemetric','rfsize', ...
%         'noisevalues',{'none'}, 'userfsize',0.5, ...
%         'ylims',[0,2.5], 'CIrange',50,...
%         'saveTo',saveTo,'saveToType','svg','fontSize',16)
    pmNoisePlotsByHRF(compTable, tools, 'x0y0',[3,3],...
        'sortHRF',sortHRFlike,'usemetric','rfsize', ...
        'noisevalues',{'none'}, 'userfsize',2, ...
        'ylims',[0.5,3.5], 'CIrange',50,'lines',true,...
        'saveTo',saveTo,'saveToType','svg','fontSize',16)
%     pmNoisePlotsByHRF(compTable, tools, 'x0y0',[3,3],...
%         'sortHRF',sortHRFlike, 'usemetric','rfsize', ...
%         'noisevalues',{'none'}, 'userfsize',4, ...
%         'ylims',[2,6], 'CIrange',50,...
%         'saveTo',saveTo,'saveToType','svg','fontSize',16)
    
    
    
    
    %{
    % Noise free
    pmNoisePlotsByHRF(compTable, tools, ... 
        'sortHRF',sortHRFlike,'usemetric','eccentricity', ...
        'noisevalues',{'none'}, 'userfsize',2, ...
        'ylims',[4.2,4.7], 'CIrange',50,...
        'saveTo',saveTo,'saveToType','svg','fontSize',16);
    
    pmNoisePlotsByHRF(compTable, tools, ... 
        'sortHRF',sortHRFlike,'usemetric','polarangle', ...
        'noisevalues',{'none'}, 'userfsize',2, ...
        'ylims',[44,46], 'CIrange',50,...
        'saveTo',saveTo,'saveToType','svg','fontSize',16)
    
    pmNoisePlotsByHRF(compTable, tools, 'x0y0',[3,3],...
        'sortHRF',sortHRFlike,'usemetric','rfsize', ...
        'noisevalues',{'none'}, 'userfsize',2, ...
        'ylims',[0,4], 'CIrange',50,...
        'saveTo',saveTo,'saveToType','svg','fontSize',16)
    
    % Three noise levels
    pmNoisePlotsByHRF(compTable, tools, ... 
        'sortHRF',sortHRFlike,'usemetric','eccentricity', ...
        'noisevalues',{'mid','high','low'}, 'userfsize',2, ...
        'ylims',[3.8,4.8], 'CIrange',50,...
        'saveTo',saveTo,'saveToType','svg','fontSize',16)
    
    pmNoisePlotsByHRF(compTable, tools, ... 
        'sortHRF',sortHRFlike,'usemetric','polarangle', ...
        'noisevalues',{'mid','high','low'}, 'userfsize',2, ...
        'ylims',[38,52], 'CIrange',50,...
        'saveTo',saveTo,'saveToType','svg','fontSize',16)
    
    pmNoisePlotsByHRF(compTable, tools, 'x0y0',[3,3],...
        'sortHRF',sortHRFlike,'usemetric','rfsize', ...
        'noisevalues',{'mid','high','low'}, 'userfsize',2, ...
        'ylims',[0,4], 'CIrange',50,...
        'saveTo',saveTo,'saveToType','svg','fontSize',16)

    %}
end
  
%% plotCloudPoints
if plotCloudPoints
% Optional params to be used as varargin when creating the function
% tools = {'aprf','aprfcss','vista','vistahrf','pop','popnohrf','afni4','afni6'};

kk = mrvNewGraphWin('NoiselessCloudPoints','wide');
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(kk,'Position',[0.007 0.62  0.8  0.3]);
subplot(1,4,1)
tools  = {'vista'};
useHRF = 'friston';
nslvl  = 'none';
pmCloudOfResults(compTable   , tools ,'onlyCenters',false ,'userfsize' , 2, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 2     ,'noiselevel' ,nslvl , ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')

subplot(1,4,2)
tools  = {'afni4'};
useHRF = 'afni_spm';
nslvl  = 'none';
pmCloudOfResults(compTable   , tools ,'onlyCenters',false ,'userfsize' , 2, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 2     ,'noiselevel' ,nslvl , ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')

subplot(1,4,3)
tools  = {'popnohrf'};
useHRF = 'popeye_twogammas';
nslvl  = 'none';
pmCloudOfResults(compTable   , tools ,'onlyCenters',false ,'userfsize' , 2, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 2     ,'noiselevel' ,nslvl , ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')

subplot(1,4,4)
tools  = {'aprf'};
useHRF = 'canonical';
nslvl  = 'none';
pmCloudOfResults(compTable   , tools ,'onlyCenters',false ,'userfsize' , 2, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-', ...
                 'lineWidth' , 1.5   ,'noiselevel' ,nslvl , ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')

fnameRoot = 'Noisefree_accuracy';
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.svg')),'svg');



% LOW NOISE
mm = mrvNewGraphWin('NoiselessCloudPoints');
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(mm,'Position',[0.007 0.62  0.8  0.8]);
tools   = {'vista','afni4','popnohrf','aprf'};
useHRFs = {'friston','afni_spm','popeye_twogammas','canonical'};
nslvl   = 'low';
np      = 0;
for tool = tools; for useHRF = useHRFs
    np=np+1;
    subplot(4,4,np)
    pmCloudOfResults(compTable   , tool ,'onlyCenters',false ,'userfsize' , 2, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF{:},'lineStyle' , '-', ...
                 'lineWidth' , .7     ,'noiselevel' ,nslvl , 'addtext',false, ...
                 'color', [0.5,0.5,0.5], 'xlims',[.5, 5.5],'ylims',[.5, 5.5],...
                 'xtick',[1:5],'ytick',[1:5], 'addcibar', true, ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')
end;end
fnameRoot = ['CloudPlots_4x4_Noise_' nslvl];
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.svg')),'svg');


% MID NOISE
mm = mrvNewGraphWin('NoiselessCloudPoints');
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(mm,'Position',[0.007 0.62  0.8  0.8]);
tools   = {'vista','afni4','popnohrf','aprf'};
useHRFs = {'friston','afni_spm','popeye_twogammas','canonical'};
nslvl   = 'mid';
np      = 0;
for useHRF = useHRFs; for tool = tools
    np=np+1;
    subplot(4,4,np)
    % figure
    pmCloudOfResults(compTable   , tool ,'onlyCenters',false ,'userfsize' , 2, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF{:},'lineStyle' , '-', ...
                 'lineWidth' , .7     ,'noiselevel' ,nslvl , 'addtext',false, ...
                 'color', [0.5,0.5,0.5], 'xlims',[.5, 5.5],'ylims',[.5, 5.5],...
                 'xtick',[1:5],'ytick',[1:5], 'addcibar', false, ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')
end;end
fnameRoot = ['CloudPlots_4x4_Noise_' nslvl];
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.svg')),'svg');


% HIGH NOISE
mm = mrvNewGraphWin('NoiselessCloudPoints');
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(mm,'Position',[0.007 0.62  0.8  0.8]);
tools   = {'vista','afni4','popnohrf','aprf'};
useHRFs = {'friston','afni_spm','popeye_twogammas','canonical'};
nslvl   = 'high';
np      = 0;
for tool = tools; for useHRF = useHRFs
    np=np+1;
    subplot(4,4,np)
    pmCloudOfResults(compTable   , tool ,'onlyCenters',false ,'userfsize' , 2, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF{:},'lineStyle' , '-', ...
                 'lineWidth' , .7     ,'noiselevel' ,nslvl , 'addtext',false, ...
                 'color', [0.5,0.5,0.5], 'xlims',[.5, 5.5],'ylims',[.5, 5.5],...
                 'xtick',[1:5],'ytick',[1:5], 'addcibar', true, ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')
end;end
fnameRoot = ['CloudPlots_4x4_Noise_' nslvl];
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.svg')),'svg');






% MID NOISE, ALL MIXED HRFs
mm = mrvNewGraphWin('MidNoiseMixHRFCloudPoints');
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(mm,'Position',[0.007 0.62  0.8  0.3]);
tools   = {'vista','afni4','popnohrf','aprf'};
useHRF  = 'mix';
nslvl   = 'mid';
np      = 0;
for tool = tools
    np=np+1;
    subplot(1,4,np)
    % figure
    pmCloudOfResults(compTable   , tool ,'onlyCenters',false ,'userfsize' , 2, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF ,'lineStyle' , '-', ...
                 'lineWidth' , .7     ,'noiselevel' ,nslvl , 'addtext',true, ...
                 'color', [0.5,0.5,0.5], 'xlims',[.5, 5.5],'ylims',[.5, 5.5],...
                 'xtick',[1:5],'ytick',[1:5], 'addcibar', true, ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')
end
fnameRoot = ['CloudPlots_MixHRF_Noise_' nslvl];
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.svg')),'svg');






end

%% plotHRFwidthtestsREALS
if plotHRFwidthtests
    COMBINE_PARAMETERS                       = struct();
    COMBINE_PARAMETERS.RF.Centerx0           = [3];
    COMBINE_PARAMETERS.RF.Centery0           = [3];  
    COMBINE_PARAMETERS.RF.sigmaMajor         = [2];  
    COMBINE_PARAMETERS.RF.sigmaMinor         = 'same';
    COMBINE_PARAMETERS.TR                    = 1;
    %{
    HRF                                      = struct();
    HRF(1).Type                              = 'vista_twogammas';  
    HRF(1).normalize                         = 'absarea'; 
        
    HRF(2).Type                              = 'afni_spm';
    HRF(2).normalize                         = 'absarea'; 
    
    HRF(3).Type                              = 'popeye_twogammas';
    HRF(3).normalize                         = 'absarea'; 
    
    HRF(4).Type                              = 'canonical'; 
    HRF(4).normalize                         = 'absarea'; 
    %}
    HRF                                      = struct();
    HRF(1).Type                              = 'boynton';  
    HRF(1).normalize                         = 'absarea'; 
    HRF(1).params.n = 3;
    HRF(1).params.tau = 1.08;
    HRF(1).params.delay = 2.05;
        
    HRF(2).Type                              = 'boynton';
    HRF(2).normalize                         = 'absarea'; 
    HRF(2).params.n = 3;
    HRF(2).params.tau = 1.38;
    HRF(2).params.delay = 2;
    
    HRF(3).Type                              = 'boynton';
    HRF(3).normalize                         = 'absarea'; 
    HRF(3).params.n = 3;
    HRF(3).params.tau = 1.68;
    HRF(3).params.delay = 1.75;
    
    HRF(4).Type                              = 'boynton';
    HRF(4).normalize                         = 'absarea'; 
    HRF(4).params.n = 3;
    HRF(4).params.tau = 1.935;
    HRF(4).params.delay = 1.65;
    
    HRF(5).Type                              = 'canonical'; 
    HRF(5).normalize                         = 'absarea'; 
    
    COMBINE_PARAMETERS.HRF                   = HRF;
        Noise                                = struct();
        Noise(1).seed                        = 'none';
    COMBINE_PARAMETERS.Noise                 = Noise;
    synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS, 'repeats',1);
    synthDT = pmForwardModelCalculate(synthDT);
    sDT = synthDT;
    
    %% Solve it
    results = pmModelFit(sDT, 'aprf')
    
    %% Plot it
    mrvNewGraphWin('HRF comparison','tall')
    subplot(3,1,1)
    a = [];
    for ii = 1:height(sDT)
        thispm = sDT.pm(ii);
        a = [a;thispm.HRF.plot('window',false,'dots',false,'addwidth',true,'xlims',[0,20])];
    end
    legend(a,{'vista\_twogammas','afni\_spm','pop\_twogammas','aprf\_canonical'})
    
    %% Create 5 data points to simulate the two-sigmas
    pm    = prfModel;
    pm.TR = 1;
    pm.HRF.Type = 'vista_twogammas';
    pm.HRF.compute
    
    subplot(1,2,1)
    pm.Stimulus.plot('window',window)
    
    subplot(1,2,2)
    pm.computeBOLD
    pm.plot('what','nonoisetimeseries','window',window,'centerzero',true); hold on
    ylim([-.1,.2]);xlim([0,30])
    pm.HRF.Type = 'popeye_twogammas';
    pm.HRF.compute
    pm.computeBOLD
    pm.plot('what','nonoise','window',window,'color','r','addtext',false,'centerzero',true)
    title('middle bar one frame, two HRF')
    legend({'timeseries','vista\_twogammas','popeye\_twogammas'})
end


%% plotHRFwidthtestsBOYNTON
if plotHRFwidthtests
    COMBINE_PARAMETERS                       = struct();
    COMBINE_PARAMETERS.RF.Centerx0           = [3];
    COMBINE_PARAMETERS.RF.Centery0           = [3];  
    COMBINE_PARAMETERS.RF.sigmaMajor         = [4];  
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
    
    HRF(5).Type                              = 'canonical'; 
    HRF(5).normalize                         = 'height'; 
    
    COMBINE_PARAMETERS.HRF                   = HRF;
        Noise                                = struct();
        Noise(1).seed                        = 'none';
    COMBINE_PARAMETERS.Noise                 = Noise;
    synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS, 'repeats',1);
    synthDT = pmForwardModelCalculate(synthDT);
    sDT = synthDT;
    
    %% Solve it
    boyntonresults = pmModelFit(sDT, 'aprf');
    
    %% Create comptTable
    boyntoncompTable  = pmResultsCompare(sDT, {'aprf'}, {boyntonresults}, ...
        'params', paramDefaults, ...
        'shorten names',true, ...
        'dotSeries', false);
    
    %% Plot it
    hh = mrvNewGraphWin('HRF comparison');
    set(hh,'Position',[0.007 0.62  0.8  0.8]);

    Cs  = 0.65 * distinguishable_colors(6,'w');
    
    % Create the fit plots with the ground truth
    tools  = {'aprf'}; nslvl  = 'none';
    HRFs   = {'boynton','boynton','boynton','boynton','canonical'};
    for ii=1:height(boyntoncompTable)
        subplot(2,5,ii)
        useHRF = HRFs{ii};
        ttable = boyntoncompTable(ii,:);
        pmCloudOfResults(ttable   , tools ,'onlyCenters',false ,'userfsize' , 4, ...
            'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-','color',Cs(ii+1,:), ...
            'lineWidth' , 2     ,'noiselevel' ,nslvl , ...
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
    
    fnameRoot = 'HRF_and_width';
    saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.svg')),'svg');
    
    
    
end


%% plotHRFwidthtestsVISTA
if plotHRFwidthtestsVISTA
    %%
    COMBINE_PARAMETERS                       = struct();
    COMBINE_PARAMETERS.RF.Centerx0           = [3];
    COMBINE_PARAMETERS.RF.Centery0           = [3];  
    COMBINE_PARAMETERS.RF.sigmaMajor         = [2];  
    COMBINE_PARAMETERS.RF.sigmaMinor         = 'same';
    COMBINE_PARAMETERS.TR                    = 1;

    HRF                                      = struct();
    HRF(1).Type                              = 'friston';  
    HRF(1).normalize                         = 'absarea'; 
    HRF(1).params.a = [6 12];
    HRF(1).params.b = [.9 .9];
    HRF(1).params.c = 0.4;
        
    HRF(2).Type                              = 'friston';
    HRF(2).normalize                         = 'absarea'; 
    HRF(2).params.a = [7 12];
    HRF(2).params.b = [.9 .9];
    HRF(2).params.c = 0.3;
    
    HRF(3).Type                              = 'friston';
    HRF(3).normalize                         = 'absarea'; 
    HRF(3).params.a = [8 12];
    HRF(3).params.b = [.9 .9];
    HRF(3).params.c = 0.2;
    
    HRF(4).Type                              = 'friston';
    HRF(4).normalize                         = 'absarea'; 
    HRF(4).params.a = [9 12];
    HRF(4).params.b = [.9 .9];
    HRF(4).params.c = 0.1;
    
    HRF(5).Type                              = 'friston';
    HRF(5).normalize                         = 'absarea'; 
    HRF(5).params.a = [10 12];
    HRF(5).params.b = [.9 .9];
    HRF(5).params.c = 0.0;
    
    COMBINE_PARAMETERS.HRF                   = HRF;
        Noise                                = struct();
        Noise(1).seed                        = 'none';
    COMBINE_PARAMETERS.Noise                 = Noise;
    synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS, 'repeats',1);
    synthDT = pmForwardModelCalculate(synthDT);
    sDT = synthDT;
    
    %% Solve it
    fristonresults = pmModelFit(sDT, 'aprf');
    
    %% Create comptTable
    fristoncompTable  = pmResultsCompare(sDT, {'aprf'}, {fristonresults}, ...
        'params', paramDefaults, ...
        'shorten names',true, ...
        'dotSeries', false);
    
    %% Plot it
    hh = mrvNewGraphWin('HRF comparison');
    set(hh,'Position',[0.007 0.62  0.8  0.8]);

    Cs  = 0.65 * distinguishable_colors(6,'w');
    
    % Create the fit plots with the ground truth
    tools  = {'aprf'}; nslvl  = 'none';
    HRFs   = {'friston','friston','friston','friston','friston'};
    for ii=1:height(fristoncompTable)
        subplot(2,5,ii)
        useHRF = HRFs{ii};
        ttable = fristoncompTable(ii,:);
        pmCloudOfResults(ttable   , tools ,'onlyCenters',false ,'userfsize' , 2, ...
            'centerPerc', 90    ,'useHRF'     ,useHRF,'lineStyle' , '-','color',Cs(ii+1,:), ...
            'lineWidth' , 2     ,'noiselevel' ,nslvl , ...
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
            thisleg = {sprintf('friston %i',ii)};
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
    
    fnameRoot = 'HRF_and_width';
    saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.svg')),'svg');
    
    
    
end


%% 
pm = prfModel;
pm.compute;
% Make the RF independent as a figure, otherwise it breaks Designer
rfrf = mrvNewGraphWin('RF');
set(rfrf,'Position',[0.007 0.62  .7  .7]);
pm.RF.plot('window',false); view([5 5 13])
grid off; set(gca, 'xtick',[],'ytick',[])
colormap('gray')
fnameRoot = ['ForwardModel_justRF'];
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.png')),'png');

nrows = 4; ncols = 3;
xx = mrvNewGraphWin('ForwardModel','wide');
set(xx,'Position',[0.007 0.62  1  1]);
subplot(nrows,ncols,1)
subplot(nrows,ncols,2)
pm.Stimulus.plot('window',false)
grid off; set(gca, 'xtick',[],'ytick',[])
subplot(nrows,ncols,3)
pm.HRF.plot('window',false,'dots',false,'color','k');
grid off; set(gca, 'xtick',[],'ytick',[])
subplot(nrows,ncols,4)
pm.plot('window',false,'what','timeseries','color','k');
grid off; set(gca, 'xtick',[],'ytick',[])
subplot(nrows,ncols,5)
pm.plot('what','noiseless','window',false,'color','k')
grid off; set(gca, 'xtick',[],'ytick',[])
subplot(nrows,ncols,6)
pm.plot('what','withnoise','window',false,'color','k')
grid off; set(gca, 'xtick',[],'ytick',[])
subplot(nrows,ncols,7)
pm.Noise.plot('window',false,'color','k')
grid off; set(gca, 'xtick',[],'ytick',[])

fnameRoot = ['ForwardModel'];
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.svg')),'svg');