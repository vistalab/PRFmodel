% SCRIPT s00_MainFiguresScript.m
%
% This script sets up your Matlab environment and downloads data needed to
% reproduce the calculations in
%
%      <PaperName>.  
% 
% There are two ways to use the script.  In both cases we
%
%    * Check the compute environment
%    * Download the relevant data
%
% A: Repeat calculations from pre-processed data.  This method takes 15
% minutes and requires about 3GB of data.
%
%    1) Download the processed data from the Open Science Foundation (OSF)
%       site (2.3GB) 
%    2) Run the Matlab scripts to produce the figures
%
% B: Repeat all calculations from scratch.  This method takes several hours
% and requires a substantial amount of disk space (-15GB-) 
%
%    1) Use prfSynthesize to create the ground-truth data (Docker)
%    2) Use prfAnalyze to estimate the parameters (Docker)
%    3) Use prfReport to generate the summary statistics (Docker)
%    4) Run the Matlab scripts that create the published figures
%
% See also
%

%% Choose your option: 
repeatCalculations = false; % It will download calculated data (2.4Gb) from OSF
%                              and plot figures
% repeatCalculations = true; % (1) It will check you have Docker installed and the
%                                  required containers are downloaded. 
%                              (2) It will start a long process of data synthesis
%                                  and analysis. 
%                              (3) It will download the experimental data and
%                                  analyze it (15Gb aprox). 
% 
%                              We recommend to use a server for this option. 
% 
testMode = true;  % Set to true if you want to test your environment with a small 
                  % dataset before running the long calculations

%% CHECK ENVIRONMENT AND PREPARE DATA

% If calculations need to be repeated we need to check the containers are
% installed %   The function
%  'checkEnvironment' has the necessary commands for installing and setting
%  up the working environment.
%
%  To begin, we suggest you verify that your Matlab paths and Docker
%  containers are properly installed.
%
%
pmCheckEnvironment(repeatCalculations)

% Prepare data for the figures (download it or calculate it all)
pmPrepareData(repeatCalculations, testMode)

%% MAKE THE FIGURES 
%{ 
% Options: 
% -----------------------
% 1/ png files are for validations, for the paper the files are saved as svg
%    and then edited in Affinnity Designer into the final form
    fileType  = 'png'; % or 'svg'
    
% 2/ Add a local folder to write the Matlab figures
    saveFigTo = fullfile(pmRootPath,'local','figures');  % Folder path
    if ~exist(saveFigTo,'dir'), mkdir(saveFigTo); end
    
% 3/ Select a confidence interval for the plots. In the manuscript we used just
%    50% and 90%. Select one and repeat plots. 
    ConfInt   = 50;  % or 90
%}

% Run the figure scripts: 
% -----------------------
s01_Ellipse_Fig1;

s02_Ellipse_Fig2;

s03_Ellipse_Fig3; 

s04_5_Ellipse_Fig4_5;

s06_S7_S8_Ellipse_Fig6_S7_S8(saveFigTo,fileType); 

% Supplementary to Figure 6
ss06_Ellipse_FigS6(saveFigTo,fileType);

%% END