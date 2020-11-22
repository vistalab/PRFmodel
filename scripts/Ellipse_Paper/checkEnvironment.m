function checkEnvironment(repeatCalculations)
% Prepares the compute environment for reproducing PRF elliptical 
%
% Synopsis
%  checkEnvironment(repeatCalculations)  (maybe pmSetEnvironment)
%
% Brief description
%
%   Installs the Matlab repositories and the Docker containers needed to
%   reproduce the figures in the PRF elliptical model paper.
% 
% Inputs
%   repeatCalculations:  If true, the Docker containers are also downloaded
%                        (Logical) 
%
% Outputs
%
% See also
%

%% Check if ToolboxToolbox is installed or not

doYouUseToolboxToolbox = false;
[~,f,~] = fileparts(which('tbUse'));
if strcmp(f,'tbUse')
    doYouUseToolboxToolbox = true;
end


%% Load dependencies

if doYouUseToolboxToolbox
    % If the user has ToolboxToolbox from the mighty Brainard, then we
    % check were the configurations are stored
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
    % A disappointing person, like Wandell, who does not yet use
    % ToolboxToolbox.  We clone the two repositories into local
    
    localDir = fullfile(pmRootPath,'local');
    cd(localDir);
    % The PRFmodel is already installed.
    
    % Vistasoft
    system('git clone https://github.com/vistalab/vistasoft.git')
    addpath(genpath(fullfile(localDir,'vistasoft')));
    disp('Vistasoft downloaded and added to path')
    
    % freesurfer_mrtrix_afni_matlab_tools
    system('git clone https://github.com/garikoitz/freesurfer_mrtrix_afni_matlab_tools.git')
    addpath(genpath(fullfile(localDir,'freesurfer_mrtrix_afni_matlab_tools')));
    disp('Freesurfer tools downloaded and added to path')
    
end

%% Check for the Docker containers
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