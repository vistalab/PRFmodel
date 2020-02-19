function prfReportWrapper(json, output_dir)
% Takes a json file with file names used to generate reports
% If instead of a jason file an empty string is passed or it cannot find the file, then it will write a default json file into output_dir 
% 
% REQUIRED INPUTS:
%      output_dir is where BIDS folder should exist already, otherwise will
%      fail
% 
% HELP: 
%       If 'help', '-h', '--help', or nothing (nargin==0), is passed in
%       this help will be displayed.
% 
% 
% USAGE:
%       Pass in a JSON file, a JSON text string, or a path to a directory
%       containing a JSON file to the docker container to initiate a
%       dtiInit processing run (see INPUT section for JSON schema):
% 
% 
% Use compile.sh for compiling
% Use this command to launch in matlab
%{
    % Create files
    jsonPath   = fullfile(pmRootPath,'local','test2','prfreport-configuration-defaults.json');
    output_dir = fullfile(pmRootPath,'local','test2');
    jsonPath   = fullfile(pmRootPath,'local','ellipse','prfreport-configuration-ellipse-sess02.json');
    output_dir = fullfile(pmRootPath,'local','ellipse');
    prfReportWrapper(jsonPath, output_dir);

%}
% Use this command to run the docker in the directory
% 
% (C) Vista Lab, Stanford University, 2019
% 

%% Initial checks

% If nothing was passed in, display help and return
if nargin == 0
    help_file = '/opt/help.txt';
    if exist(help_file, 'file')
        system(['cat ', help_file]);
    else
        help(mfilename);
    end
    return
end

% Assume the user wanted to see the help, and show it
if ischar(json) 
    if strcmpi(json, 'help') || strcmpi(json, '-help') || strcmpi(json, '-h') || strcmpi(json, '--help')
        help(mfilename);
    end
end

%% Parse the JSON file or object
% Create a pm instance, we will use it in both cases
% Check if we need to read a json or provide a default one
if exist(json, 'file') == 2
    fprintf('[prfreport] Reading json file...')
    J = loadjson(json);
    if iscell(J)
        J=J{:};
    end
    fprintf('done\n')
else
    DEFAULTS             = struct();
    % We can only compare to one synthetic dataset
    DEFAULTS.subjectName = "editDefaultSubjectName";
    DEFAULTS.sessionName = "editDefaultSessionName";
    % But we can compare several diferent analyses to the same data
        analyze             = struct(); 
        analyze(1).Type     = "valid: aprf popeye vista afni";
        analyze(2).Type     = "there can be more than one";
    DEFAULTS.analyze        = analyze;
    DEFAULTS.createplots    = false;
    % Other params
        resultParams        = struct();
        resultParams(1).name = "Centerx0";
        resultParams(2).name = "Centery0";
        resultParams(3).name = "Theta";
        resultParams(4).name = "sigmaMinor";
        resultParams(5).name = "sigmaMajor";
        resultParams(6).name = "R2";
    DEFAULTS.resultParams      = resultParams;    
        resultParamsShort(1).shortname = "x0";
        resultParamsShort(2).shortname = "y0";
        resultParamsShort(3).shortname = "Th";
        resultParamsShort(4).shortname = "sMin";
        resultParamsShort(5).shortname = "sMaj";
        resultParamsShort(6).shortname = "R2";
    DEFAULTS.resultParamsShort = resultParamsShort;    
    DEFAULTS.shortenParamNames = true;
    DEFAULTS.doTSeries         = false;
    % pmCloudOfResults params
        pmCloudOfResultsParams = struct();
        pmCloudOfResultsParams.onlyCenters = false;
        pmCloudOfResultsParams.userfsize   = 2;
        pmCloudOfResultsParams.centerPerc  = 90;
        pmCloudOfResultsParams.useHRF      = "{}";
        pmCloudOfResultsParams.lineStyle   = "-";
        pmCloudOfResultsParams.lineWidth   = 0.7;
        pmCloudOfResultsParams.fontsize    = 14;
        pmCloudOfResultsParams.noiselevel  = "mid";
        pmCloudOfResultsParams.addtext     = true;
        pmCloudOfResultsParams.useellipse  = false;
        pmCloudOfResultsParams.color       = [0.5,0.5,0.5];
        pmCloudOfResultsParams.xlims       = [0, 5.5];
        pmCloudOfResultsParams.ylims       = [0, 5.5];
        pmCloudOfResultsParams.xtick       = [1:5];
        pmCloudOfResultsParams.ytick       = [1:5];
        pmCloudOfResultsParams.addcihist   = true;
        pmCloudOfResultsParams.addcibar    = false;
        pmCloudOfResultsParams.newWin      = false;
        pmCloudOfResultsParams.saveToType  = 'svg';
    DEFAULTS.pmCloudOfResultsParams        = pmCloudOfResultsParams;  
    % Select filename to be saved
    fname = fullfile(output_dir, 'prfreport-configuration-defaults.json');
    % Encode json
    jsonString = jsonencode(DEFAULTS);
    % Format a little bit
    jsonString = strrep(jsonString, ',', sprintf(',\n'));
    jsonString = strrep(jsonString, '[{', sprintf('[\n{\n'));
    jsonString = strrep(jsonString, '}]', sprintf('\n}\n]'));

    % Write it
    fid = fopen(fname,'w');if fid == -1,error('Cannot create JSON file');end
    fwrite(fid, jsonString,'char');fclose(fid);
    % Permissions
    if isdeployed
        fileattrib(fname,'+w +x', 'o g'); 
    end
    disp('prfreport-configuration-defaults.json written, edit it and pass it to the prfreport container as the configuration json.')
    return
end

%% Find the BIDS structure and files
fprintf('[prfreport] Finding BIDS structure and files ...')
% Noah let me know how you want to handle this, I can remove it from here
BIDSdir = fullfile(output_dir, 'BIDS');
if ~exist(BIDSdir,'dir');error("No BIDS structure in provided output_dir %s",output_dir);end

% Find the synthetic json file
synthDefFile = fullfile(BIDSdir,'derivatives','prfsynth',...
                        ['sub-' J.subjectName],['ses-' J.sessionName], ...
                        ['sub-' J.subjectName '_ses-' J.sessionName '_task-prf_acq-normal_run-01_bold.json']);
if ~exist(synthDefFile,'file');error("Can't find synthetic definition file %s",synthDefFile);end

% Generate the table with the synthetic data
SynthDT  = struct2table(jsonread(synthDefFile));
for na=1:width(SynthDT); if isstruct(SynthDT{:,na});
        SynthDT.(SynthDT.Properties.VariableNames{na}) = struct2table(SynthDT{:,na});
end; end




% Find and load the result file
resultsFile = {};
resultsNames= {};
for nr=1:length(J.analyze)
    antype           = J.analyze{nr}.Type;
    resultsNames{nr} = antype; 
    resultDir = fullfile(BIDSdir,'derivatives',['prfanalyze-' antype],...
                                  ['sub-' J.subjectName],['ses-' J.sessionName]);
    if ~exist(resultDir,'dir');error("Can't find results directory %s",resultDir);end
    cd(resultDir)
	
    %% Read the results back to Matlab
    % Read all the nifti-s for this type of tool 
	niftis = dir('*.nii.gz');
    if length(niftis) < 6
		error('Not enough nifti result files in output, at least there should be x0,y0,sigmamajor,sigmaminor,theta and modelpred')
	end
	filesRead = 0;
    
    pmEstimates = table();
	for nn=1:length(niftis)
		fname = niftis(nn).name;
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
	    % asign it to the table
        % Noah used other names in the filenames, leave this hardcoded here until we decide what to do
        % noahnames = {'modelpred', 'sigmamajor','sigmaminor', 'testdata', 'theta', 'x0' ,'y0'};
        % {'Centerx0','Centery0','Theta','sigmaMinor','sigmaMajor'}
		switch resname
            case 'sigmaminor', resname = 'sigmaMinor';
            case 'sigmamajor', resname = 'sigmaMajor';
            case 'theta'     , resname = 'Theta';
            case 'centerx0'  , resname = 'Centerx0';
            case 'centery0'  , resname = 'Centery0';
            case 'x0'        , resname = 'Centerx0';
            case 'y0'        , resname = 'Centery0';
        end
        pmEstimates.(resname) = data;	        
    end
    
    % We are always going to add by default the R2
    % This will be calculated between the testdata and modelpred
    pmEstimates.R2 = calccod(pmEstimates.testdata,  pmEstimates.modelpred,2);
 
    % If one of the elements is missing, fill it with zeros. 
    % Popeye was not providing Theta.
    if ~contains('Theta', pmEstimates.Properties.VariableNames)
        pmEstimates.Theta = zeros(size(pmEstimates.Centerx0));
    end

   
    % Return the pmEstimates for this tool
    resultsFile{nr} = pmEstimates;
end
fprintf('files loaded and table created.\n')



%% Concatenate synthetic data params and the actual results
% Define how to manage this, as input in the config?
% anNames = fieldnames(res);
% Select the result files
% ress = [{res.(anNames{1})}];
% for an =2:length(anNames)
%     ress = [ress, {res.(anNames{an})}];
% end
fprintf('[prfreport] Concatenating results from tools with synthetic ground-truth...')
if length(J.resultParams) ~= length(J.resultParamsShort)
    error('The number of long and short names needs to be same and in the same order.')
end
for np=1:length(J.resultParams)
    paramDefaults{np} = J.resultParams{np}.name;
    shortnames{np}    = J.resultParamsShort{np}.shortname;
end

compTable  = pmResultsCompare(SynthDT, resultsNames, resultsFile, ...
    'params', paramDefaults, ...
    'shorten names', J.shortenParamNames, ...
    'short names', shortnames, ...
    'addsnrcol',false, ...
    'dotSeries', J.doTSeries);
fprintf('done\n')

%% Generate the output directory
fprintf('[prfreport] Generating output directory and saving files...')
reportDir = fullfile(BIDSdir,'derivatives','prfreport', ...
                        ['sub-' J.subjectName],['ses-' J.sessionName]);
if ~exist(reportDir,'dir');mkdir(reportDir);end

%% Generate the output files
% The .mat files first. Decide final format. I think json. 
fname_trunk = ['sub-' J.subjectName '_ses-' J.sessionName '-prf_acq-normal_run-01_bold'];

reportFile  = fullfile(reportDir,[fname_trunk '.mat']);
save(reportFile, 'compTable')

% Save it as  json file as well
% Select filename to be saved
fname = fullfile(reportDir, [fname_trunk '.json']);
% Encode json
jsonString = jsonencode(compTable);
% Format a little bit
jsonString = strrep(jsonString, ',', sprintf(',\n'));
jsonString = strrep(jsonString, '[{', sprintf('[\n{\n'));
jsonString = strrep(jsonString, '}]', sprintf('\n}\n]'));

% Write it
fid = fopen(fname,'w');if fid == -1,error('Cannot create JSON file');end
fwrite(fid, jsonString,'char');fclose(fid);
% Permissions
if isdeployed
    fileattrib(fname,'+w +x', 'o g');
end
fprintf('done\n')




%% Generate the output figures
if J.createplots
% Change some defaults for the plots
set(0,'defaultAxesFontName', 'Helvetica')
set(0,'defaultTextFontName', 'Helvetica')
% MID NOISE, ALL MIXED HRFs
% {
fprintf('[prfreport] Generating and saving output figures...')
% Generic values coming from the config.json
% TODO: add defaults if variables not set in the json file
onlyCenters = J.pmCloudOfResultsParams.onlyCenters;
userfsize   = J.pmCloudOfResultsParams.userfsize;
useHRF      = J.pmCloudOfResultsParams.useHRF;
centerPerc  = J.pmCloudOfResultsParams.centerPerc;
lineStyle   = J.pmCloudOfResultsParams.lineStyle;
lineWidth   = J.pmCloudOfResultsParams.lineWidth;
fontsize    = J.pmCloudOfResultsParams.fontsize;
noiselevel  = J.pmCloudOfResultsParams.noiselevel;
addtext     = J.pmCloudOfResultsParams.addtext;
useellipse  = J.pmCloudOfResultsParams.useellipse;
color       = J.pmCloudOfResultsParams.color;
xlims       = J.pmCloudOfResultsParams.xlims;
ylims       = J.pmCloudOfResultsParams.ylims;
xtick       = J.pmCloudOfResultsParams.xtick;
ytick       = J.pmCloudOfResultsParams.ytick;
addcihist   = J.pmCloudOfResultsParams.addcihist;
addcibar    = J.pmCloudOfResultsParams.addcibar;
saveToType  = J.pmCloudOfResultsParams.saveToType;
% load(fullfile(pmRootPath,'local','results.mat'));
saveTo = reportDir;
%% FIGURE 5
fnameRoot = 'Noisefree_accuracy';
kk = mrvNewGraphWin(fnameRoot,'wide','off');
numanalysis = length(J.analyze);
set(kk,'Units','centimeters','Position',[0 0 10*numanalysis 10]);
nslvl  = 'none';
for na=1:numanalysis
    subplot(1,numanalysis,na)
    tools  = J.analyze{na}.Type;
    switch tools
        case {'vista','mrvista','vistasoft'}
            useHRF = 'vista_twogammas';
        case {'pop','popeye'}
            useHRF = 'popeye_twogammas';
        case {'afni','afni4','afni6','afnidog'}
            useHRF = 'afni_spm';
        case {'aprf','analyzeprf'}
            useHRF = 'canonical';
        otherwise
            warning('%s not recorded, using vista_twogammas as default',J.analyze{na})
    end    
    pmCloudOfResults(compTable   , {tools} ,'onlyCenters',onlyCenters ,'userfsize' , userfsize, ...
                 'centerPerc', centerPerc ,'useHRF'     ,useHRF,'lineStyle' , lineStyle, ...
                 'lineWidth' , lineWidth, 'noiselevel' ,nslvl , 'fontsize', fontsize, ...
                 'useellipse', useellipse, ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType', saveToType)
end
set(gca,'FontName', 'Arial')
saveas(kk,fullfile(saveTo, strcat(fnameRoot,'.',saveToType)),saveToType);



%% FIGURE 7
tools   = {};
useHRFs = {};
for nj=1:length(J.analyze)
    tool = J.analyze{nj}.Type;
    switch tool
        case {'vista','mrvista','vistasoft'}
            useHRF = 'vista_twogammas';
        case {'pop','popeye'}
            useHRF = 'popeye_twogammas';
        case {'afni','afni4','afni6','afnidog'}
            useHRF = 'afni_spm';
        case {'aprf','analyzeprf'}
            useHRF = 'canonical';
        otherwise
            warning('%s not recorded, using vista_twogammas as default',J.analyze{na})
    end    
    tools{nj}   = tool;
    useHRFs{nj} = useHRF;
end

nslvls   = {'low','mid'};
for nslvl = nslvls
    fnameRoot = ['CloudPlots_4x4_Noise_' nslvl{:}];
    mm = mrvNewGraphWin(fnameRoot,[],'off');
    set(mm,'Units','centimeters','Position',[0 0 10*numanalysis 10*numanalysis]);
    np      = 0;
    for tool = tools; for useHRF = useHRFs
        np=np+1;
        subplot(numanalysis,numanalysis,np)
        pmCloudOfResults(compTable   , tool ,'onlyCenters',onlyCenters ,'userfsize' , userfsize, ...
                     'centerPerc', centerPerc    ,'useHRF'     ,useHRF{:},'lineStyle' , lineStyle, ...
                     'lineWidth' , lineWidth     ,'noiselevel' ,nslvl{:} , 'useellipse', useellipse, ...
                     'addtext',addtext, 'adddice',false,'addsnr',false,...
                     'color', color, 'xlims',xlims,'ylims',ylims,'fontsize', fontsize, ...
                     'xtick',xtick,'ytick',ytick, 'addcibar', addcibar,'addcihist', addcihist,  ...
                     'newWin'    , false ,'saveTo'     ,'','saveToType',saveToType)
    end;end
    set(gca,'FontName', 'Arial')
    saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.',saveToType)),saveToType);
end


%% FIGURE 8
useHRF  = 'mix';
nslvls  = {'low','mid'};
for nslvl=nslvls   
    fnameRoot = ['CloudPlots_MixHRF_Noise_HIST' nslvl{:}];
    mm = mrvNewGraphWin(fnameRoot,[],'off');
    set(mm,'Units','centimeters','Position',[0 0 10*numanalysis 10]);
    np      = 0;
    for tool = tools
        np=np+1;
        subplot(1,numanalysis,np)
        % figure
        pmCloudOfResults(compTable   , tool ,'onlyCenters',onlyCenters ,'userfsize' , userfsize, ...
                     'centerPerc', centerPerc    ,'useHRF'     ,useHRF ,'lineStyle' , lineStyle, ...
                     'lineWidth' , lineWidth     ,'noiselevel' ,nslvl{:} , 'addtext',addtext, ...
                     'color', color, 'xlims',xlims,'ylims',ylims,'fontsize', fontsize, ...
                     'useellipse', useellipse, ...
                     'xtick',xtick,'ytick',ytick, 'addcihist', addcihist, ...
                     'newWin'    , false ,'saveTo'     ,'','saveToType',saveToType)
    end
    set(gca,'FontName', 'Arial')
    saveas(gcf,fullfile(saveTo, strcat(fnameRoot,'.',saveToType)),saveToType);
end


fprintf('done\n')
end  % createplots


%% Change file attributes and close
if isdeployed
    fprintf('[prfreport] Changing file attributes and exiting...')
    fileattrib(fullfile(BIDSdir,'derivatives','prfreport'),'+w +x', 'o'); 
    fileattrib(reportDir,'+w +x', 'o'); 
    fprintf('done.\n')
end
fprintf('[prfreport] Exiting.\n')
return 



