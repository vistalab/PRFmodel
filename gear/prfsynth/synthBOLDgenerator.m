function synthBOLDgenerator(json, output_dir)
% Takes a json file with parameters required to generate:
%     1/ name.nii.gii    : nifti file with the synthetic BOLD signal
%     2/ name.json       : json file with the parameters of the BOLD tSeries
%     3/ name_Stim.nii.gz: nifti file with the stimulus used to generate the BOLD tSeries
% If instead of a jason file an empty string is passed or it cannot find the file, then it will write a default json file into output_dir 
% 
% REQUIRED INPUTS:
%       jsonParams file. Use the synthBOLDgenerator_paramsExample.json as a
%       template always, and see instructions for editing. 
% 
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
%       % Using a JSON file
%        docker run --rm -ti -v `pwd`/input:/input -v `pwd`/output:/output vistalab/dtiinit /input/dtiInit.json
% 
%       % Using a JSON string
%        docker run --rm -ti -v `pwd`/input:/input -v `pwd`/output:/output vistalab/dtiinit '{"input_dir":"/input", "output_dir": "/output"}'
% 
%       % Using a directory (in the container), containing a JSON (.json)
%        docker run --rm -ti -v `pwd`/input:/input -v `pwd`/output:/output vistalab/dtiinit /input/
% 
% Use compile.sh for compiling
% Use this command to launch in matlab
%{
    % Create files
    jsonPath   = fullfile(pmRootPath,'local','output','paper', ...
                          'params_big_test_for_paper_3-3_v02BOLD.json');
    output_dir = fullfile(pmRootPath,'local','output','paper');
    synthBOLDgenerator(jsonPath, output_dir);








    % Run the analysis
    results = pmModelFit({'synthBOLD_20190917T222512/synthBOLD.nii.gz', ...
                          'synthBOLD_20190917T222512/synthBOLD.json', ...
                          'synthBOLD_20190917T222512/synthBOLD_Stim.nii.gz'}, ...
                          'aPRF');
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
pm = prfModel;
% Check if we need to read a json or provide a default one
if exist(json, 'file') == 2
    J = loadjson(json);
    if iscell(J)
        J=J{:};
    end
else
    % return the default json and exit
    DEFAULTS = pm.defaultsTable;
    % Add two params required to generate the file
    DEFAULTS.subjectName = "editDefaultSubjectName";
    DEFAULTS.sessionName = "editDefaultSessionName";
    DEFAULTS.repeats  = 2;
    % Reorder fieldnames
    DEF_colnames = DEFAULTS.Properties.VariableNames;
    DEF_colnames = DEF_colnames([end-2:end,1:end-3]);
    DEFAULTS     = DEFAULTS(:,DEF_colnames);
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

%% Create an output subfolder for the outputs 
% outputSubFolder = [J.subjectName '_' datestr(datetime,'yyyymmddTHHMMSS','local')];
outputSubFolder = [J.subjectName '_' J.sessionName];
output_dir = fullfile(output_dir, outputSubFolder);
mkdir(output_dir);

%% Generate the table with the synthetic data
% This is a working example, to use in the testing process against the json data
%{
COMBINE_PARAMETERS.RF.Centerx0        = [0,5];
COMBINE_PARAMETERS.RF.Centery0        = [5];
COMBINE_PARAMETERS.RF.Theta           = [0];
COMBINE_PARAMETERS.RF.sigmaMajor      = [1,2];
COMBINE_PARAMETERS.RF.sigmaMinor      = 'same';
COMBINE_PARAMETERS.Noise.noise2signal = [0, 0.2];  % By default only white noise added
COMBINE_PARAMETERS.TR                 = [1.5];
    HRF(1).Type = 'friston';
    HRF(2).Type = 'canonical';
COMBINE_PARAMETERS.HRF                = HRF;
TEST = pmForwardModelTableCreate(COMBINE_PARAMETERS, 'repeats', 2);
TEST = pmForwardModelCalculate(TEST);
%}
% Create one pm instance to obtain defaults
PARAMETERS = J;
PARAMETERS = rmfield(PARAMETERS,'subjectName');
PARAMETERS = rmfield(PARAMETERS,'sessionName');
PARAMETERS = rmfield(PARAMETERS,'repeats');
% Add first the simple parameters

% Add now the sub-table params
% Do the required conversions before creating the table
% To concatenate structs they need to have the same num of fields
for nh=1:length(PARAMETERS.HRF)
    thisHRF = PARAMETERS.HRF(nh); 
        PARAMETERS.HRF(nh) = {pmParamsCompletenessCheck(thisHRF{:}, ...
                              table2struct(pm.defaultsTable.HRF))};
end
PARAMETERS.HRF       = [PARAMETERS.HRF{:}];
PARAMETERS.RF        = [PARAMETERS.RF{:}];
PARAMETERS.Stimulus  = [PARAMETERS.Stimulus{:}];
% Noise can have different amount of variables. Complete them and concatenate
for nh=1:length(PARAMETERS.Noise)
    thisNoise = PARAMETERS.Noise(nh); 
    if isfield(thisNoise{:},'voxel')
        completeNoise = pmParamsCompletenessCheck(thisNoise{:}, ...
            table2struct(pm.Noise.defaultsGet('voxel',thisNoise{:}.voxel)));
    else
        completeNoise = pmParamsCompletenessCheck(thisNoise{:}, ...
            table2struct(pm.defaultsTable.Noise));
    end
    PARAMETERS.Noise(nh) = {completeNoise};
end
PARAMETERS.Noise     = [PARAMETERS.Noise{:}];
% Convert some char-s to string-s, char-s are treated as individual elements...
PARAMETERS.Type             = string(PARAMETERS.Type);
PARAMETERS.RF.Type          = string(PARAMETERS.RF.Type);
PARAMETERS.HRF.Type         = string(PARAMETERS.HRF.Type);
PARAMETERS.HRF.normalize    = string(PARAMETERS.HRF.normalize);
PARAMETERS.Noise.seed       = string(PARAMETERS.Noise.seed);
PARAMETERS.Stimulus.expName = string(PARAMETERS.Stimulus.expName);

% Generate the same thing from the json file
synthDT = pmForwardModelTableCreate(PARAMETERS, 'repeats', J.repeats);
synthDT = pmForwardModelCalculate(synthDT);

%% Generate the files

% Save the default niftis with different TR and HRF to be used as tests later on

% BOLD FILE
cd(output_dir)
fname = [J.subjectName '.nii.gz'];
pmForwardModelToNifti(synthDT, 'fname',fname, 'demean',false);

% JSON FILE
jsonSynthFile = [J.subjectName '.json'];
% Encode json
jsonString = jsonencode(synthDT(:,1:(end-1)));
% Format a little bit
jsonString = strrep(jsonString, ',', sprintf(',\n'));
jsonString = strrep(jsonString, '[{', sprintf('[\n{\n'));
jsonString = strrep(jsonString, '}]', sprintf('\n}\n]'));
% Write it
fid = fopen(jsonSynthFile,'w');if fid == -1,error('Cannot create JSON file');end
fwrite(fid, jsonString,'char');fclose(fid);

% STIM FILE
stimNiftiFname = [J.subjectName '_Stim.nii.gz'];
pm1            = synthDT.pm(1);
stimNiftiFname = pm1.Stimulus.toNifti('fname',stimNiftiFname);

% Permissions
fileattrib(output_dir,'+w +x', 'o'); 

return 



