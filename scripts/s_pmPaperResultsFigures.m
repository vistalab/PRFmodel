%% s_pmPaperResultsFigures.m
% First the nifti with the data was created using prf-synthesize
% The json params files used is in scripts/params_big_test_for_paper_v01.json

% Create the names that will be used in this script

close all;clear all;clc
tbUse prfmodel;

% Control execution of the script file
defineFileNames    = true;
generateInputFiles = false;
uploadInputFiles   = false;
analysisInputFiles = false;
downloadLoadFiles  = true;
plotGroupAnalysis  = true;
plotCloudPoints    = true;

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

%% plotGroupAnalysis
if (plotGroupAnalysis)
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
  
    
    % PLOTS
    tools = {'vista','afni4','popnohrf','aprf'}';
    % Or select them all
    % tools = anNames;
    
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
                 'lineWidth' , .7     ,'noiselevel' ,nslvl , ...
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
for tool = tools; for useHRF = useHRFs
    np=np+1;
    subplot(4,4,np)
    pmCloudOfResults(compTable   , tool ,'onlyCenters',false ,'userfsize' , 2, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF{:},'lineStyle' , '-', ...
                 'lineWidth' , .7     ,'noiselevel' ,nslvl , ...
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
                 'lineWidth' , .7     ,'noiselevel' ,nslvl , ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')
end;end
fnameRoot = ['CloudPlots_4x4_Noise_' nslvl];
saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.svg')),'svg');


end


