% SCRIPT s00_MainFiguresScript.m

% To obtain the data and replicate the figures, there are two options: 
%    --- A: Repeat all calculations (it takes many hours and disk space -15Gb-) ---
%    1) SYNTHESIZE
%    2) ANALYZE
%    3) REPORT
%    4) MAKE FIGURES

%    --- B: Download pre-calculated data ---
%    1) DOWNLOAD CALCULATIONS FROM OSF (2.3Gb)
%    2) MAKE FIGURES

% Choose your option: 
repeatCalculations = false; % It will download calculated data (2.4Gb and plot
%                               figures)
% repeatCalculations = true; % It will start a long process of data synthesis
%                              and analysis, for the experimental data it will 
%                              download and analyze it (15Gb aprox total). We
%                              recommend to use a server for this option. 

% Other Notes and options: 
% -----------------------
% 1/ png files are for validations, for the paper the files are saved as svg
%    and then edited in Affinnity Designer into the final form
    fileType  = 'png'; % or 'svg'
    
% 2/ Add a local folder to write the figures
    saveFigTo = '~/Downloads/TEST';  % Folder path
    
% 3/ Select a confidence interval for the plots. In the manuscript we used just
%    50% and 90%. Select one and repeat plots. 
    ConfInt   = 50;  % or 90

%% CHECK ENVIRONMENT AND PREPARE DATA
% If calculations need to be repeated we need to check the containers are installed
checkEnvironment(repeatCalculations)

% Either download the calculations from OSF or run everything from scratch
prepareData(repeatCalculations)

%% MAKE THE FIGURES 
s01_Ellipse_Fig1(saveFigTo,fileType);
s02_Ellipse_Fig2(saveFigTo,fileType);
s03_Ellipse_Fig3(saveFigTo,fileType); 
s04_5_Ellipse_Fig4_5(saveFigTo,fileType,ConfInt); 
s06_S7_S8_Ellipse_Fig6_S7_S8(saveFigTo,fileType); 
ss06_Ellipse_FigS6(saveFigTo,fileType);
