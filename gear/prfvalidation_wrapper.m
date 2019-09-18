function output_dir = prfvalidation_wrapper(json)
% Main wrapper to the PRFmodel validation software
% Options: 
%     - Generate synthetic data and exit
%     - Run one of the four implemented packages and exit
%     - Generate a report on the results and exit
%    (- Some combinations of above)
% 
%     The options are specified in the config.json file. 
%     The inputs are specified in config.json as well. 



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
J.output_dir = fullfile(J.output_dir, outputSubFolder);
mkdir(J.output_dir);


%% Do the tasks
8




synthBOLDgenerator(json)






end

