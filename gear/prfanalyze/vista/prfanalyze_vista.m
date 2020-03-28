function prfanalyze_vista(opts_file, json_file, bold_file, stim_file, output_dir)
% 
% (C) Vista Lab, Stanford University, 2019
% 
%{
bold_file  = '/data/localhome/glerma/toolboxes/PRFmodel/local/paper02/BIDS/sub-paper/ses-sess02/func/sub-paper_ses-sess02_task-prf_acq-normal_run-01_bold.nii.gz';
json_file  = '/data/localhome/glerma/toolboxes/PRFmodel/local/paper02/BIDS/derivatives/prfsynth/sub-paper/ses-sess02/sub-paper_ses-sess02_task-prf_acq-normal_run-01_bold.json';
stim_file  = '/data/localhome/glerma/toolboxes/PRFmodel/local/paper02/BIDS/stimuli/sub-paper_ses-sess02_task-prf_apertures.nii.gz';
output_dir = '/data/localhome/glerma/toolboxes/PRFmodel/local/paper02';
opts_file  = '/data/localhome/glerma/toolboxes/PRFmodel/local/paper02/prfanalyze_vista_config_paper_sess02.json'
%}
    
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
if ischar(json_file)
    if strcmpi(json_file, 'help') ...
            || strcmpi(json_file, '-help') ...
            || strcmpi(json_file, '-h') ...
            || strcmpi(json_file, '--help')
        help(mfilename);
    end
end

% read in the opts file
%{
opts = {};
if ~isempty(opts_file)
    tmp = loadjson(opts_file);
    if ~isempty(tmp)
        % tmp = tmp.options;
        fs = fields(tmp);
        for ii = 1:numel(fs)
            opts{end+1} = fs{ii};
            opts{end+1} = getfield(tmp, fs{ii});
        end
    end
end
%}

% read in the opts file
if ~isempty(opts_file)
    fprintf('This is the config.json file being read: %s\n',opts_file)
    tmp = loadjson(opts_file);
    disp('These are the contents of the json file:')
    tmp
    if ~isempty(tmp)
        opt.options.vista = tmp;
        opts = {'options', opt};
    else
        opts = {};
    end
else
    opts = {};
end

% Make the output directory
mkdir(output_dir);


%% check that the other relevant files eist
if exist(bold_file, 'file') ~= 2
    disp(sprintf('Given BOLD 4D nifti file does not exist: %s', bold_file))
    return
end
if exist(stim_file, 'file') ~= 2
    disp(sprintf('Given stimulus 3D nifti file does not exist: %s', stim_file))
    return
end

%% Call pmModelFit!
disp('================================================================================');
disp(opts_file);
disp(bold_file);
disp(json_file);
disp(stim_file);
disp('--------------------------------------------------------------------------------');

[pmEstimates, results] = pmModelFit({bold_file, json_file, stim_file}, 'vista', opts{:});
%% Write out the results
estimates_file = fullfile(output_dir, 'estimates.mat');
estimates = struct(pmEstimates);
save(estimates_file, 'estimates', 'pmEstimates');
results_file = fullfile(output_dir, 'results.mat');
save(results_file, 'results');


% Save it as  json file as well
% Select filename to be saved
fname = fullfile(output_dir, ['estimates.json']);
% Encode json
jsonString = jsonencode(pmEstimates);
% Format a little bit
jsonString = strrep(jsonString, ',', sprintf(',\n'));
jsonString = strrep(jsonString, '[{', sprintf('[\n{\n'));
jsonString = strrep(jsonString, '}]', sprintf('\n}\n]'));

% Write it
fid = fopen(fname,'w');if fid == -1,error('Cannot create JSON file');end
fwrite(fid, jsonString,'char');fclose(fid);







% Permissions
fileattrib(output_dir,'+w +x', 'o'); 

return 



