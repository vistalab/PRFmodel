function checkEnvironment(repeatCalculations)
%CHECKENVIRONMENT Prepares environment to create figures from downloaded or
%                 calculated data

% Check if ToolboxToolbox is installed or not
doYouUseToolboxToolbox = false;
[~,f,~] = fileparts(which('tbUse'));
if strcmp(f,'tbUse')
    doYouUseToolboxToolbox = true;
end


% Load dependencies
if doYouUseToolboxToolbox
    % Check were the configurations are stores
    tbUse('sample-repo');
    [p,~,~] = fileparts(which('master.txt'));
    yourToolboxToolboxRegistryConfigurationsFolder = ...
           fullfile(strrep(p,'sample-repo','ToolboxRegistry'),'configurations');
    % Download and add this config file
    configFile = fullfile(yourToolboxToolboxRegistryConfigurationsFolder,'PRFmodel.json');
    if ~exist(configFile,'file')
        websave(configFile,'https://github.com/garikoitz/ToolboxRegistry/blob/master/configurations/PRFmodel.json')
    end
    tbUse PRFmodel;
    
    fprintf('\n\n------------------------------------------------------------------------\n',...
             fullfile(strrep(p,'sample-repo','PRFmodel'),'local'));
    fprintf('--> NOTE: all required data will be downloaded to the gitignored folder:\n');
    fprintf('-->       %s',fullfile(strrep(p,'sample-repo','PRFmodel'),'local'));
    fprintf('\n------------------------------------------------------------------------\n\n',...
             fullfile(strrep(p,'sample-repo','PRFmodel'),'local'));
else
    % If not a ToolboxToolbox user:
    softDir = '~/softForEllipsePaper/';
    if ~exist(softDir,'dir');mkdir(softDir);end
    cd(softDir);
    % PRFmodel
    system('git clone https://github.com/vistalab/PRFmodel.git')
    addpath(genpath(fullfile(softDir,'PRFmodel')));
    % Vistasoft
    system('git clone https://github.com/vistalab/vistasoft.git')
    addpath(genpath(fullfile(softDir,'vistasoft')));
    % freesurfer_mrtrix_afni_matlab_tools
    system('git clone https://github.com/garikoitz/freesurfer_mrtrix_afni_matlab_tools.git')
    addpath(genpath(fullfile(softDir,'freesurfer_mrtrix_afni_matlab_tools')));
    
    
    
    fprintf('\n\n------------------------------------------------------------------------\n',...
             fullfile(strrep(p,'sample-repo','PRFmodel'),'local'));
    fprintf('--> NOTE 1: all required software will be installed in:\n');
    fprintf('-->         %s\n\n',softDir);
    
    fprintf('--> NOTE 2: all required data will be downloaded to the gitignored folder:\n');
    fprintf('-->         %s',fullfile(softDir,'PRFmodel','local'));
    fprintf('\n------------------------------------------------------------------------\n\n',...
             fullfile(strrep(p,'sample-repo','PRFmodel'),'local'));
      
    
end

% Check for the Docker containers
if repeatCalculations
    % prfsynth
    dockerPullCommand('garikoitz/prfsynth:2.0.0')
    
    % prfanalyze-vista
    dockerPullCommand('garikoitz/prfanalyze-vista:2.0.0')
    
    % prfanalyze-afni
    dockerPullCommand('garikoitz/prfanalyze-afni:2.0.0')
        
    % prfreport
    dockerPullCommand('garikoitz/prfreport:2.0.0')
    
end


end