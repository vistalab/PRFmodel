function pmPrepareData(repeatCalculations, testMode)
% Download empirical data or calculate synthetic ground-truth data
%
%  Synopsis
%   pmPrepareData(repeatCalculations, testMode)
%
% Input:
%   repeatCalculations:  
%
% These data are used for replicating the figures in the paper
%
%   Title: 
%
%
% See also
%
% TODO: 
% - Convert to switch statements so that we can choose what to run

%% DOWNLOAD DATA OR CALCULATE IT

if ~repeatCalculations
    cd(fullfile(pmRootPath,'local'));
    % Download calculated data from synthetic data (30Mb)
    fname = fullfile(pmRootPath,'local','ellipse.zip');
    if ~isfile(fname)
        ellipsezip = websave(fname,'https://osf.io/27axp/download');
        unzip(ellipsezip);
    end

    % Download calculated data from experimental data (2.3Gb)
    fname = fullfile(pmRootPath,'local','realdata.zip');
    if ~isfile(fname)
        realdatazip = websave(fname,'https://osf.io/s9fd5/download');
        unzip(realdatazip);
    end
    calculateThis = 'downloadData';
else
    if testMode
        calculateThis = 'testData';
    else
        calculateThis = 'allData';
    end
end
    
   
%% CALCULATE RESULTS
switch calculateThis
    case {'downloadData'}
        disp('Data has been downloaded and it is ready to create figures.')
    case {'testData'}
        % --------------
        % 1/   Synthesize
        % --------------
        pmLaunchDockerCommand('prfsynth','testmode','testv1')
        % --------------
        % 2.1/ Analyze-afni
        % --------------
        pmLaunchDockerCommand('prfanalyze','testmode','testv1','afni6')
        % --------------
        % 2.2/ Analyze-vista
        % --------------
        pmLaunchDockerCommand('prfanalyze','testmode','testv1','vista6')
        % --------------
        % 3/   Report
        % --------------
        pmLaunchDockerCommand('prfreport','testmode','testv1')
        
    case {'allData'}
        %% SYNTHETIC DATA
        % {
            % --------------
            % 1/ Synthesize
            % --------------
            pmLaunchDockerCommand('prfsynth','ellipse','eccsv2')
            pmLaunchDockerCommand('prfsynth','ellipse','eccsv2TR1')
            pmLaunchDockerCommand('prfsynth','ellipse','sizesv2')
            pmLaunchDockerCommand('prfsynth','ellipse','sizesv2TR1')
            pmLaunchDockerCommand('prfsynth','ellipse','tr1dur300v2')
            pmLaunchDockerCommand('prfsynth','ellipse','thetasv2')

            % --------------
            % 2.1/ Analyze-afni
            % --------------
            pmLaunchDockerCommand('prfanalyze','ellipse','eccsv2','afni6')
            pmLaunchDockerCommand('prfanalyze','ellipse','eccsv2TR1','afni6')
            pmLaunchDockerCommand('prfanalyze','ellipse','sizesv2','afni6')
            pmLaunchDockerCommand('prfanalyze','ellipse','sizesv2TR1','afni6')
            pmLaunchDockerCommand('prfanalyze','ellipse','thetasv2','afni6')
            
            % --------------
            % 2.2/ Analyze-vista
            % --------------
            pmLaunchDockerCommand('prfanalyze','ellipse','eccsv2','vista6')
            pmLaunchDockerCommand('prfanalyze','ellipse','eccsv2TR1','vista6')
            pmLaunchDockerCommand('prfanalyze','ellipse','sizesv2','vista6')
            pmLaunchDockerCommand('prfanalyze','ellipse','sizesv2TR1','vista6')
            pmLaunchDockerCommand('prfanalyze','ellipse','thetasv2','vista6')
            pmLaunchDockerCommand('prfanalyze','ellipse','tr1dur300v2','vista6')
            pmLaunchDockerCommand('prfanalyze','ellipse','tr1dur300v2','vista4')

            % --------------
            % 3/ Report
            % --------------
            pmLaunchDockerCommand('prfreport','ellipse','eccsv2')
            pmLaunchDockerCommand('prfreport','ellipse','eccsv2TR1')
            pmLaunchDockerCommand('prfreport','ellipse','sizesv2')
            pmLaunchDockerCommand('prfreport','ellipse','sizesv2TR1')
            pmLaunchDockerCommand('prfreport','ellipse','tr1dur300v2')
            pmLaunchDockerCommand('prfreport','ellipse','thetasv2')
        %}

        %% EXPERIMENTAL DATA

        % --------------
        % 1/ Download the HCP 7T data from OSF
        % --------------
        % {
            fname   = fullfile(pmRootPath,'local','hcp7T.zip');
            if ~isfile(fname)
                hcp1    = websave(fullfile(pmRootPath,'local','hcp1.zip'),...
                                               'https://osf.io/az5y6/download');
                hcp2    = websave(fullfile(pmRootPath,'local','hcp2.zip'),...
                                               'https://osf.io/udzs2/download');
                concat  = system(['cat ' hcp1 ' ' hcp2 ' > ' fname]);
                unzip(fname);
            end
            
        %}    

        % NOTE: same with experimental data. We downloaded the results that we use for
        %       analyses purposes. The data was originally downloaded  from the HCP 7T
        %       repository, modified, and then analyzed with the docker containers. 
        %       We provide the raw, intermediate and processed files in this link:
        %           hcp_7T_data_and_analysisWithConfig_00 (5Gb) https://osf.io/az5y6/download
        %           hcp_7T_data_and_analysisWithConfig_01 (4Gb) https://osf.io/udzs2/download
        %           To link them back together: 
        %              cat hcp_7T_data_and_analysisWithConfig_00 hcp_7T_data_and_analysisWithConfig_01 > hcp7T.zip


        % --------------
        % 2.2/ Analyze-vista
        % --------------
        pmLaunchDockerCommand('prfanalyze','realdata','tr1dur300v3','vista6')            


        % --------------
        % 3/ Report
        % --------------
        % launchDockerCommand('garikoitz/prfreport:2.0.0','eccv2')
        % launchDockerCommand('garikoitz/prfreport:2.0.0','eccv2')
        % launchDockerCommand('garikoitz/prfreport:2.0.0','eccv2')
end  % END switch calculateThis
end  % END function
    

