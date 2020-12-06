function pmLaunchDockerCommand(varargin)
%LAUNCHDOCKERCOMMAND Launches Dockers based on config files we already have


% docker = 'prfsynth';
% docker = 'prfanalyze-vista';
% solver = [];
% solver = 'vista6'
% sub    = 'ellipse';
% ses    = 'sizesv2';

% launchDockerCommand(docker,sub,ses,solver)

docker = varargin{1};
sub = varargin{2};
ses = varargin{3};
if length(varargin) < 4
    solver = [];
else
    solver = varargin{4};
end




% Examples of actual calls in black server
% ./run_prfsynth.sh $ELLIPSE /data/localhome/glerma/toolboxes/PRFmodel/local/ellipse/prfsynth-config_sub-ellipse_sess-sizesv2TR1.json
% ./run_prfanalyze.sh vista $WF /data/localhome/glerma/toolboxes/PRFmodel/local/WORDSFOV/prfanalyze-vista4-config_sub-003_ses-1.json
% ./run_prfreport.sh $ELLIPSE /data/localhome/glerma/toolboxes/PRFmodel/local/ellipse/prfreport-configuration_sub-ellipse_ses-eccsv2TR1.json

% What is the osf link for the config files?
if strcmp(sub,'testmode')
    teststr = 'testmode_';
    osfLink = 'https://osf.io/ehws5/download';
else
    teststr = '';
    osfLink = 'https://osf.io/2q9uw/download';
end


% Be sure the config file is there, otherwise download it all
switch sub
    case {'ellipse','testmode'}
        basedir = fullfile(pmRootPath,'local',sub); 
        if ~isfolder(basedir); mkdir(basedir); end
        cd(basedir)
        if ~isfolder(fullfile(pmRootPath,'local',sub,[teststr 'config_files']))
            configszip = websave(fullfile(pmRootPath,'local',[teststr 'config_files']),osfLink);
            unzip(configszip);
        end
    otherwise
        basedir = fullfile(pmRootPath,'local','realdata');
        if ~isfolder(basedir); mkdir(basedir); end
        cd(basedir)
        if ~isfolder(fullfile(pmRootPath,'local','realdata',[teststr 'config_files']))
            configszip = websave(fullfile(pmRootPath,'local',[teststr 'config_files']),osfLink);
            unzip(configszip);
        end
end





 %basedir = fullfile(pmRootPath,'local',sub);
switch docker
    case {'prfsynth'}
        % This is the config file
        config_fname = fullfile(basedir,[teststr 'config_files'],...
                            [docker '-config_sub-' sub '_ses-' ses '.json']);
        % This is the command line
        cmd = [fullfile(pmRootPath,'gear',docker,['run_' docker '.sh']) ' ' ...
                                                     basedir ' ' config_fname];
        resultFile = fullfile(basedir,'BIDS','derivatives',docker,...
                      ['sub-' sub],['ses-' ses],...
                      ['sub-' sub '_ses-' ses '_task-prf_acq-normal_run-01_bold.json']);
        if isfile(resultFile)
            fprintf('\n\nResult file already exists, skipping. \n')
        else
            fprintf('\n\nSending the following command:\n')
            fprintf('---> %s\n', cmd)

            s = system(cmd);
            if s > 0
                error('The command failed')
            else
                fprintf('Running ...\n\n')
            end
        end
    case {'prfanalyze'}
        d = split(docker,'-');
        % We are going to use this config file for the Docker container
       
        config_fname = fullfile(basedir,[teststr 'config_files'],...
            [docker '-' solver(1:end-1) '-config_sub-' sub '_ses-' ses ...
                                                    '_solver-' solver '.json']);
       
        % Generate the command line to launch the docker container
        cmd = [fullfile(pmRootPath,'gear',docker,'run_prfanalyze.sh --version 2.0.0') ...
               ' ' solver(1:end-1) ' ' basedir ' ' config_fname]
        % This is one example file in the output, to check if it has already
        % been run
        resultFile = fullfile(basedir,'BIDS','derivatives',docker,...
                      ['sub-' sub],['ses-' ses],...
                      ['sub-' sub '_ses-' ses '_task-prf_acq-normal_run-01_bold.json']);
        % If we run it already, skip it, otherwise run
        if isfile(resultFile)
            fprintf('\n\nResult file already exists, skipping. \n')
        else
            fprintf('\n\nSending the following command:\n')
            fprintf('---> %s\n', cmd)

            s = system(cmd);
            if s > 0
                error('The command failed')
            else
                fprintf('Running ...\n\n')
            end
        end    
    case {'prfreport'}
        % This is the config file
        config_fname = fullfile(basedir,[teststr 'config_files'],...
            [docker '-config_sub-' sub '_ses-' ses '.json']);
        % This is the command line
        cmd = [fullfile(pmRootPath,'gear',docker,['run_' docker '.sh']) ' ' basedir ' ' config_fname];
        resultFile = fullfile(basedir,'BIDS','derivatives',docker,...
                      ['sub-' sub],['ses-' ses],...
                      ['sub-' sub '_ses-' ses '_task-prf_acq-normal_run-01_bold.json']);
        if isfile(resultFile)
            fprintf('\n\nResult file already exists, skipping. \n')
        else
            fprintf('\n\nSending the following command:\n')
            fprintf('---> %s\n', cmd)

            s = system(cmd);
            if s > 0
                error('The command failed')
            else
                fprintf('Running ...\n\n')
            end
        end        

    otherwise
        error('%s not recognized, valid options are prfsynth, prfreport, prfanalyze-vista, prfanalyze-afni');
end


end

