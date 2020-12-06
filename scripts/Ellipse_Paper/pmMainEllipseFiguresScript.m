% SCRIPT pmMainEllipseFiguresScript.m
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
%       site (750MB) 
%    2) Run the Matlab scripts to produce the figures
%
% B: Repeat all calculations from scratch.  This method can take days
% and requires a substantial amount of disk space (-more than 15GB-) 
%
%    1) Use prfSynthesize to create the ground-truth data (Docker)
%    2) Use prfAnalyze to estimate the parameters (Docker)
%    3) Use prfReport to generate the summary statistics (Docker)
%    4) Run the Matlab scripts that create the published figures
%
% See also
%

%% Choose your option: 
repeatCalculations = false; % It will download calculated data (750MB)from OSF
%                              and plot figures
% repeatCalculations = true; % (1) It will check you have Docker installed and the
%                                  required containers are downloaded. 
%                              (2) It will start a long process of data synthesis
%                                  and analysis. 
%                              (3) It will download the experimental data and
%                                  analyze it with vista6 and vista4. 
% 
%                              We recommend to use a server for this option. 
% 
testMode = false;  % Set to true if you want to test your environment with a small 
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

%% RUN THE FIGURE SCRIPTS
% There are some options that can be set inside the scripts
% 
pmEllipse_Fig1;
pmEllipse_Fig2;
pmEllipse_Fig3; 
% pmEllipse_Fig3_TestMode;  % Use it to plot the results of the testMode 
pmEllipse_Fig4_5;
pmEllipse_Fig6_S7_S8; 
pmEllipse_FigS6;

%% END