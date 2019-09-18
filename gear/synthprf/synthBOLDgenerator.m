function synthBOLDgenerator(json, output_dir)
% Takes a json file with parameters required to generate:
%     1/ name.nii.gii    : nifti file with the synthetic BOLD signal
%     2/ name.json       : json file with the parameters of the BOLD tSeries
%     3/ name_Stim.nii.gz: nifti file with the stimulus used to generate the BOLD tSeries
% 
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
    jsonPath   = fullfile(pmRootPath,'gear','synthprf','synthBOLDgenerator_paramsExample.json');
    output_dir = fullfile(pmRootPath,'local','output');
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

if exist(json, 'file') == 2
    J = loadjson(json);
elseif exist(json, 'dir') == 7
    jsonFile = dir(fullfile(json, '*.json'));
    jsonFile = fullfile(json, jsonFile.name);
    disp(jsonFile);
    if ~isempty(jsonFile)
        J = loadjson(jsonFile);
    else
        error('No JSON file could be found');
    end
elseif ~isempty(json) && ischar(json)
    try
        J = loadjson(json);
    catch ME
        disp(ME.message); 
        return
    end
else
    error('Could not find/parse the json file/structure');
end

%% Check the json object for required fields
%{
required = {'input_dir', 'output_dir'};
err = false;

for r = 1:numel(required)
    if ~isfield(J, required{r})
        err = true;
        fprintf('%s not found in JSON object!\n', required{r});
    elseif ~exist(J.(required{r}), 'dir')
        fprintf('%s Does not exist!\n', required{r});
        err = true;
    end
end

% If there was a problem, return
if err 
    error('Exiting! There was a problem with the inputs. Please check input_dir and output_dir!');
end
%}


% Create an output subfolder for the outputs 
outputSubFolder = ['synthBOLD_', datestr(datetime,'yyyymmddTHHMMSS','local')];
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
TEST = pmForwardModelTableCreate(COMBINE_PARAMETERS, 'mult', 2);
TEST = pmForwardModelCalculate(TEST);
%}
PARAMETERS = J;
PARAMETERS = rmfield(PARAMETERS,'output_dir');
PARAMETERS = rmfield(PARAMETERS,'fileName');
PARAMETERS = rmfield(PARAMETERS,'mult');
% Do the required conversions before creating the table
PARAMETERS.HRF       = [PARAMETERS.HRF{:}]; 
PARAMETERS.RF        = [PARAMETERS.RF{:}];
PARAMETERS.Stimulus  = [PARAMETERS.Stimulus{:}];
PARAMETERS.Noise     = [PARAMETERS.Noise{:}];
% Convert some char-s to string-s, char-s are treated as individual elements...
PARAMETERS.Type      = string(PARAMETERS.Type);
PARAMETERS.RF.Type   = string(PARAMETERS.RF.Type);
PARAMETERS.Stimulus.expName = string(PARAMETERS.Stimulus.expName);
% Generate the same thing from the json file
synthDT = pmForwardModelTableCreate(PARAMETERS, 'mult', J.mult);
synthDT = pmForwardModelCalculate(synthDT);

%% Generate the files

% Save the default niftis with different TR and HRF to be used as tests later on

% BOLD FILE
cd(output_dir)
fname = [J.fileName '.nii.gz'];
pmForwardModelToNifti(synthDT, 'fname',fname, 'demean',false);

% JSON FILE
jsonSynthFile = [J.fileName '.json'];
% Encode json
jsonString = jsonencode(synthDT(:,1:(end-1)));
% Format a little bit
jsonString = strrep(jsonString, ',', sprintf(',\r'));
jsonString = strrep(jsonString, '[{', sprintf('[\r{\r'));
jsonString = strrep(jsonString, '}]', sprintf('\r}\r]'));
% Write it
fid = fopen(jsonSynthFile,'w');if fid == -1,error('Cannot create JSON file');end
fwrite(fid, jsonString,'char');fclose(fid);

% STIM FILE
stimNiftiFname = [J.fileName '_Stim.nii.gz'];
pm1            = synthDT.pm(1);
stimNiftiFname = pm1.Stimulus.toNifti('fname',stimNiftiFname);

% Permissions
fileattrib(output_dir,'+w +x', 'o'); 

return 



