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
    jsonPath   = fullfile(pmRootPath,'local','output','alex', ...
                                           'prfreport_alex5.json');
    output_dir = fullfile(pmRootPath,'local','output','alex');
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
    J = loadjson(json);
    if iscell(J)
        J=J{:};
    end
else
    DEFAULTS             = struct();
    % We can only compare to one synthetic dataset
    DEFAULTS.subjectName = "editDefaultSubjectName";
    DEFAULTS.sessionName = "editDefaultSessionName";
    % But we can compare several diferent analyses to the same data
        analyze             = struct(); 
        analyze(1).Type     = "valid: aprf popeye mrvista afni";
        analyze(2).Type     = "example_afni_or_mrvista";
    DEFAULTS.analyze        = analyze;
    % Other params
        resultParams        = struct();
        resultParams(1).name   = "Centerx0";
        resultParams(2).name   = "Centery0";
        resultParams(3).name   = "Theta";
        resultParams(4).name   = "sigmaMinor";
        resultParams(5).name   = "sigmaMajor";
    DEFAULTS.resultParams   = resultParams;    
    DEFAULTS.shortenParamNames = true;
    DEFAULTS.doTSeries      = false;
    
    % Select filename to be saved
    fname = fullfile(output_dir, 'defaultParams_ToBeEdited.json');
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
    fileattrib(fname,'+w +x', 'o g'); 
    disp('defaultParams_ToBeEdited.json written, edit it and pass it to the container to generate synthetic data.')
    return
end

%% Find the BIDS structure and files
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
    % I asked noah to change the name, revert this back
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
		tmp       = niftiRead(fname);
        data      = squeeze(tmp.data);	
	    % asign it to the table
        % Noah used other names in the filenames, leave this hardcoded here until we decide what to do
        % noahnames = {'modelpred', 'sigmamajor','sigmaminor', 'testdata', 'theta', 'x0' ,'y0'};
        % {'Centerx0','Centery0','Theta','sigmaMinor','sigmaMajor'}
		switch resname
            case 'sigmaminor',resname='sigmaMinor';
            case 'sigmamajor',resname='sigmaMajor';
            case 'theta',resname='Theta';
            case 'x0',resname='Centerx0';
            case 'y0',resname='Centery0';
        end
        pmEstimates.(resname) = data;	
	end
    resultsFile{nr} = pmEstimates;
end

%% Obtain the result files
% Right now is reading the matlab .mat, but when the prfanalyze is
% implemented, here we will do the conversion to the .mat from the result
% niftis
% res = struct();
% for outputNifti=outputNiftis
%     rname = outputNifti{:}(1:end-4);
%     rname = split(rname,'_');
%     rname = rname{end,:};
%     if exist(outputNifti{:},'file')
%         load(outputNifti{:});
%     else
%         % TODO Check that the data is there before trying to download
%         [~,fname,ext] = fileparts(outputNifti{:});
%         load(st.fw.downloadFileFromCollection(cc{1}.collection.id,...
%             [fname ext],outputNifti{:}));
%     end
%     res.(rname) = results;
% end

%% Concatenate synthetic data params and the actual results
% Define how to manage this, as input in the config?
% anNames = fieldnames(res);
% Select the result files
% ress = [{res.(anNames{1})}];
% for an =2:length(anNames)
%     ress = [ress, {res.(anNames{an})}];
% end

paramDefaults = {};
for np=1:length(J.resultParams)
    paramDefaults{np} = J.resultParams{np}.name;
end



compTable  = pmResultsCompare(SynthDT, resultsNames, resultsFile, ...
    'params', paramDefaults, ...
    'shorten names', J.shortenParamNames, ...
    'dotSeries', J.doTSeries);

%% Generate the output directory
reportDir = fullfile(BIDSdir,'derivatives','prfreport', ...
                        ['sub-' J.subjectName],['ses-' J.sessionName]);
if ~exist(reportDir,'dir');mkdir(reportDir);end

%% Generate the output files
% The .mat files first. Decide final format. I think json. 
reportFile = fullfile(reportDir,['results.mat']);
save(reportFile, 'compTable')

% Save it as  json file as well
% Select filename to be saved
fname = fullfile(reportDir, 'results.json');
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
fileattrib(fname,'+w +x', 'o g');




% Can we save .svg plots
% MID NOISE, ALL MIXED HRFs
%{
mm = pmNewGraphWin('MidNoiseMixHRFCloudPoints',[],'off');
% Fig size is relative to the screen used. This is for laptop at 1900x1200
set(mm,'Position',[0.007 0.62  0.8  0.3]);
tools   =  resultsNames;
useHRF  = 'mix';
nslvl   = {'high'};
np      = 0;
for tool = tools
    np=np+1;
    subplot(1,length(tools),np)
    % figure
    pmCloudOfResults(compTable   , tool ,'onlyCenters',false , ...% 'userfsize' , 2, ...
                 'centerPerc', 90    ,'useHRF'     ,useHRF ,'lineStyle' , '-', ...
                 'lineWidth' , .7     ,'noiselevel' ,nslvl , 'addtext',true, ...
                 'color', [0.5,0.5,0.5], 'xlims',[.5, 5.5],'ylims',[.5, 5.5],...
                 'xtick',[1:5],'ytick',[1:5], 'addcibar', true, ...
                 'newWin'    , false ,'saveTo'     ,'','saveToType','svg')
end
fnameRoot = ['CloudPlots_MixHRF_Noise_' nslvl];
saveas(gcf,fullfile(reportDir, strcat(fnameRoot,'.svg')),'svg');
%}


%% Permissions
fileattrib(fullfile(BIDSdir,'derivatives','prfreport'),'+w +x', 'o'); 
fileattrib(reportDir,'+w +x', 'o'); 

return 



