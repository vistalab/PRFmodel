function newts=regressCfdsfromTS(datamat, confoundsmat, varargin)
% Regresses Friston 24 confounds out from fmriprep output. 
% Based on: https://www.mail-archive.com/hcp-users@humanconnectome.org/msg07291.html
% It accepts matrices directly or nifti files, and it can write back the
% corrected file if writeNifti is set to 1
%{
datamat: can be a matrix or a path to a nifti/mgh file, if it is nifti file it will
rewrite with same name _REGRESSED
confoundsmat: can be a matrix or a path name to the .tsv output file
Optional parameters:

datamat      = '/Users/glerma/soft/PRF/local/sub-23MAGNO6134T2/fmriprep/sub-23MAGNO6134/ses-T2/func/sub-23MAGNO6134_ses-T2_task-ret_run-03_space-fsnative_hemi-L.func.gii';
confoundsmat = '/Users/glerma/soft/PRF/local/sub-23MAGNO6134T2/fmriprep/sub-23MAGNO6134/ses-T2/func/sub-23MAGNO6134_ses-T2_task-ret_run-03_desc-confounds_regressors.tsv';
writeNifti   = true;
add2name     = '_REGRESSED2';

    newts = regressCfdsfromTS(datamat, confoundsmat, 'writeNifti', true, 'add2name', add2name);

Garikoitz Lerma-Usabiaga 04.2019 garikoitz@gmail.com Stanford Vista Lab
%}
    %% Parse inputs
    p = inputParser;

    addRequired(p, 'datamat');
    addRequired(p, 'confoundsmat');
    addOptional(p, 'writeNifti', false, @islogical);
    addOptional(p, 'add2name'  , '_REGRESSED', @ischar);

    parse(p,datamat,confoundsmat,varargin{:});

    writeNifti = p.Results.writeNifti;
    add2name   = p.Results.add2name;
    
    %% Do the thing

    % Obtain the matrices
    if ischar(datamat)
        [FILEPATH,NAME,EXT] = fileparts(datamat);
        if strcmp(EXT,'.gii')
            mghFileName = [FILEPATH filesep NAME '.mgh'];
            system(['mri_convert ' datamat ' ' mghFileName])
            datamat = mghFileName;
        end
        mriFile = MRIread(datamat);
        volFile = mriFile.vol;
        mriSize = size(volFile);
        if mriSize(1) == 1
            ts = squeeze(volFile);
        else      
            ts = reshape(volFile, [mriSize(1)*mriSize(2)*mriSize(3),mriSize(4)]);
        end
    end
    
    % If we pass the .tsv from fmriprep, obtain the friston 24, and write them
    % just in case
    % Confoundsmat will be size(ts,2) in size, and with 24 regressor columns
    if ischar(confoundsmat)
        [FILEPATH,NAME] = fileparts(confoundsmat);
        outputPathName = [FILEPATH filesep NAME '_friston24.txt'];
        confoundsmat = createNewRegressors(confoundsmat,outputPathName);
    end
    
    
    % Check it just in case, if it doesn't work throw an error and ask to check
    if ~isequal(size(ts,2), size(confoundsmat,1))
        error('Data and regressors do not have same size, revise')
    end
    
    
    % Do the regression
    demeanmat = confoundsmat - mean(confoundsmat);
    newts     = ts - (demeanmat * (pinv(demeanmat) * ts'))';
    
    % Write the file back if asked to
    if writeNifti && ischar(datamat)
        mriFile.vol = reshape(newts, mriSize);
        [FILEPATH,NAME,EXT] = fileparts(datamat);
        if strcmp(EXT,'.gz')
            NAME = NAME(1:end-4);
            EXT  = '.nii.gz';
        end
        writeNiftiName = [FILEPATH filesep NAME add2name EXT];
        MRIwrite(mriFile, writeNiftiName);
    end
end
