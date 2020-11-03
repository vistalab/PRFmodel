% SCRIPT s00_MainFiguresScript.m
%        Downloads the data and
%        Generates the figures in the publication
%        (add citation here)
% 
% Note: png files are for validations, for the paper the files are saved as svg
%       and then edited in Affinnity Designer into the final form

%% PREPARE ENVIRONMENT
clear all; close all; clc

doYouUseToolboxToolbox = true;

% Load dependencies
if doYouUseToolboxToolbox
    % Download and add this config file
    yourToolboxToolboxRegistryConfigurationsFolder='~/toolboxes/ToolboxRegistry/configurations/';
    configFile = fullfile(yourToolboxToolboxRegistryConfigurationsFolder,'PRFmodel.json');
    if ~exist(configFile,'file')
        websave(configFile,'https://github.com/garikoitz/ToolboxRegistry/blob/master/configurations/PRFmodel.json')
    end
    tbUse PRFmodel
else
    % If not a ToolboxToolbox user:
    softDir = '~/soft/';
    if ~exist(softDir,'dir');mkdir(softDir);end
    cd(softDir);
    % PRFmodel
    system('git clone https://github.com/vistalab/PRFmodel.git')
    addpath(genpath(fullfile(softDir,'PRFmodel')));
    % Vistasoft
    system('git clone https://github.com/vistalab/vistasoft.git')
    addpath(genpath(fullfile(softDir,'vistasoft')));
    % Check more, maybe there are more dependencies
end

%% DOWNLOAD DATA
% Download synthetic data (30Mb)
ellipsezip = websave(fullfile(pmRootPath,'local','ellipse.zip'),'https://osf.io/27axp/download');
unzip(ellipsezip);

% NOTE: although they are not required to generate the figures, 
%       the config files used to generate the data just downloaded 
%       using the the prfsynthesis, prfanalyze, and prfreport docker images
%       can be obtained in the following link: https://osf.io/zh5cs/download
%       To run it will require installing a Docker client and downloading the
%       images. Please check the wiki for indications on how to run install and
%       run the docker containers using these config files. (22Kb)



% Download experimental data (2.3Gb)
realdatazip = websave(fullfile(pmRootPath,'local','realdata.zip'),'https://osf.io/s9fd5/download');
unzip(realdatazip);

% NOTE: same with experimental data. We downloaded the results that we use for
%       analyses purposes. The data was originally downloaded  from the HCP 7T
%       repository, modified, and then analyzed with the docker containers. 
%       We provide the raw, intermediate and processed files in this link:
%           hcp_7T_data_and_analysisWithConfig_00 (5Gb) https://osf.io/az5y6/download
%           hcp_7T_data_and_analysisWithConfig_01 (4Gb) https://osf.io/udzs2/download
%           To link them back together: 
%              cat hcp_7T_data_and_analysisWithConfig_00 hcp_7T_data_and_analysisWithConfig_01 > hcp7T.zip
%       With the provided config files and using the same docker images as above
%       the results can be reproduce in any machine with Docker installed. 

% Only the data for the figures is about 2Gb. The whole thing is around 10Gb


%% FIGURE CALLS
saveFigTo = '~/Downloads/TEST'; if ~isfolder(saveFigTo); mkdir(saveFigTo); end
fileType  = 'png';
ConfInt   = 50;

s01_Ellipse_Fig1(saveFigTo,fileType);
s02_Ellipse_Fig2(saveFigTo,fileType);
s03_Ellipse_Fig3(saveFigTo,fileType); % todo: download files from OSF
s04_5_Ellipse_Fig4_5(saveFigTo,fileType,ConfInt); % todo: download files from OSF
s06_S7_S8_Ellipse_Fig6_S7_S8(saveFigTo,fileType); % todo: download files from OSF
ss06_Ellipse_FigS6(saveFigTo,fileType);
