function pmPrepareData(repeatCalculations, testMode)
%PREPAREDATA Downloads or calculates the data for the ellipse paper



%% DOWNLOAD DATA OR CALCULATE IT

if ~repeatCalculations
    cd(fullfile(pmRootPath,'local'));
    % Download calculated data from synthetic data (30Mb)
    ellipsezip = websave(fullfile(pmRootPath,'local','ellipse.zip'),'https://osf.io/27axp/download');
    unzip(ellipsezip);

    % Download calculated data from experimental data (2.3Gb)
    realdatazip = websave(fullfile(pmRootPath,'local','realdata.zip'),'https://osf.io/s9fd5/download');
    unzip(realdatazip);
    
else
    %% SYNTHETIC DATA
    
    % --------------
    % 1/ Synthesize
    % --------------
    pmLaunchDockerCommand('prfsynth','ellipse','sizesv2')
    % launchDockerCommand('garikoitz/prfsynth:2.0.0','eccv2')
    % launchDockerCommand('garikoitz/prfsynth:2.0.0','eccv2')
    
    % --------------
    % 2.1/ Analyze-afni
    % --------------
    pmLaunchDockerCommand('prfanalyze','ellipse','sizesv2','afni6')
    % launchDockerCommand('garikoitz/prfanalyze-afni:2.0.0','eccv2')
    % launchDockerCommand('garikoitz/prfanalyze-afni:2.0.0','eccv2')
    
    % --------------
    % 2.2/ Analyze-vista
    % --------------
    pmLaunchDockerCommand('prfanalyze','ellipse','sizesv2','vista6')
    % launchDockerCommand('garikoitz/prfanalyze-vista:2.0.0','eccv2')
    % launchDockerCommand('garikoitz/prfanalyze-vista:2.0.0','eccv2')
    
    
    % --------------
    % 3/ Report
    % --------------
    pmLaunchDockerCommand('prfreport','ellipse','sizesv2')
    % launchDockerCommand('garikoitz/prfreport:2.0.0','eccv2')
    % launchDockerCommand('garikoitz/prfreport:2.0.0','eccv2')
    
    
    %% EXPERIMENTAL DATA
    
    % --------------
    % 1/ Download the HCP 7T data from OSF
    % --------------
    hcp1 = websave(fullfile(pmRootPath,'local','hcp1.zip'),'https://osf.io/az5y6/download');
    hcp2 = websave(fullfile(pmRootPath,'local','hcp2.zip'),'https://osf.io/udzs2/download');
    concat = system(['cat ' hcp1 ' ' hcp2 ' > ' fullfile(pmRootPath,'local','hcp7T.zip')]);
    unzip(ellipsezip);
    

    % NOTE: same with experimental data. We downloaded the results that we use for
    %       analyses purposes. The data was originally downloaded  from the HCP 7T
    %       repository, modified, and then analyzed with the docker containers. 
    %       We provide the raw, intermediate and processed files in this link:
    %           hcp_7T_data_and_analysisWithConfig_00 (5Gb) https://osf.io/az5y6/download
    %           hcp_7T_data_and_analysisWithConfig_01 (4Gb) https://osf.io/udzs2/download
    %           To link them back together: 
    %              cat hcp_7T_data_and_analysisWithConfig_00 hcp_7T_data_and_analysisWithConfig_01 > hcp7T.zip


        
    % --------------
    % 2.1/ Analyze-afni
    % --------------
    % launchDockerCommand('garikoitz/prfanalyze-afni:2.0.0','eccv2')
    % launchDockerCommand('garikoitz/prfanalyze-afni:2.0.0','eccv2')
    % launchDockerCommand('garikoitz/prfanalyze-afni:2.0.0','eccv2')
    
    % --------------
    % 2.2/ Analyze-vista
    % --------------
    % launchDockerCommand('garikoitz/prfanalyze-vista:2.0.0','eccv2')
    % launchDockerCommand('garikoitz/prfanalyze-vista:2.0.0','eccv2')
    % launchDockerCommand('garikoitz/prfanalyze-vista:2.0.0','eccv2')
    
    
    % --------------
    % 3/ Report
    % --------------
    % launchDockerCommand('garikoitz/prfreport:2.0.0','eccv2')
    % launchDockerCommand('garikoitz/prfreport:2.0.0','eccv2')
    % launchDockerCommand('garikoitz/prfreport:2.0.0','eccv2')
    
end
    
end

