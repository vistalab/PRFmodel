function regressors = createNewRegressors(inputPathName,outputPathName)
% Reads the regressors .tsv from fmriprep and writes back the selected columns
% TODO: make colums an optional variable

%{
inputPathName  = fullfile(prfPath,'local','sub-14magno7806','fmriprep','sub-14MAGNO7806','ses-20190303','func','sub-14MAGNO7806_ses-20190303_task-ret_run-01_desc-confounds_regressors.txt');
outputPathName = fullfile(prfPath,'local','sub-14magno7806','fmriprep','sub-14MAGNO7806','ses-20190303','func','sub-14MAGNO7806_ses-20190303_task-ret_run-01_desc-confounds_regressors_SELECTED.txt');
    regressors = createNewRegressors(inputPathName,outputPathName)



Garikoitz Lerma-Usabiaga 04.2019 garikoitz@gmail.com Stanford Vista Lab
%}



%% PArse inputs
    p = inputParser;

    addRequired(p, 'datamat');
    addRequired(p, 'confoundsmat');
    addOptional(p, 'writeNifti', false, @islogical);
    addOptional(p, 'add2name'  , '_REGRESSED', @ischar);

    parse(p,datamat,confoundsmat,varargin{:});

    writeNifti = p.Results.writeNifti;
    add2name   = p.Results.add2name;
    
%% Do the thing
[FILEPATH,NAME,EXT] = fileparts(inputPathName);
if EXT=='.tsv'
    EXT='.txt';
    inputFile = fullfile(FILEPATH,[NAME EXT]);
    copyfile(inputPathName, inputFile);
end



% Read the whole file
regressors = readtable(inputFile);
% Select the columns we are interested in

motRegressors   = {'trans_x','trans_y','trans_z','rot_x','rot_y','rot_z'};
acompRegressors = {'a_comp_cor_00','a_comp_cor_01','a_comp_cor_02',...
                    'a_comp_cor_03', 'a_comp_cor_04', 'a_comp_cor_05'};
csfwm           = {'csf','white_matter'};
% aComp, mot, (mot-1), mot^2,(mot-1)^2                
% regressors = [regressors{:,acompRegressors}, ...
%               regressors{:,motRegressors}, ... 
%               [0 0 0 0 0 0; regressors{1:end-1,motRegressors}], ...
%               regressors{:,motRegressors}.^2, ...
%               [0 0 0 0 0 0; regressors{1:end-1,motRegressors}] .^2];

regressors = [regressors{:,csfwm}, ...
              regressors{:,motRegressors}, ... 
              [0 0 0 0 0 0;diff([regressors{:,motRegressors}],1,1)] ...
              ];
%           , ...
%               diff([regressors{:,motRegressors}],2,1), ... % Second derivative
%               regressors{:,motRegressors}.^2, ...
%               [0 0 0 0 0 0; regressors{1:end-1,motRegressors}] .^2];          
dlmwrite(outputPathName, regressors);
end