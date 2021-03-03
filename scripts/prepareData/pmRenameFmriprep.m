function [newNames,newStimNames,configNames] = pmRenameFmriprep(basedir,subname,ses,tr,stimDiam,solver)
%
% Example:
%{

basedir  = '/black/localhome/glerma/TESTDATA/RETPIPE';
subname  = 'ex22163';
ses      = '22163';
tr       = 2;
stimDiam = 24;
solver   = 'vista4';

[newNames,newStimNames,configNames]=pmRenameFmriprep(basedir,subname,ses,tr,stimDiam,solver)

docker = 'prfanalyze';
dockerVer = '2.0.0';
% Launch de docker containers here
for nc=1:length(configNames)
    config_fname = fullfile(basedir,'BIDS','config_files',configNames{nc});
    cmd = [fullfile(pmRootPath,'gear',docker,'run_prfanalyze.sh') ...
               '  --version ' dockerVer ' ' solver(1:end-1) ...
               ' ' basedir ' ' config_fname]
    
    fprintf('\n\nSending the following command:\n')
    fprintf('---> %s\n', cmd)

    s = system(cmd);
    if s > 0
        error('The command failed')
    else
        fprintf('Running ...\n\n')
    end
end





% Stimulus pattern: sub-ellipse_ses-eccsv2_task-prf_apertures.nii.gz
% nifti2 pattern: sub-ellipse_ses-eccsv2_task-prf_acq-normal_run-01_bold.nii.gz

% Match stim files with acqu runs

sub-ex22163_ses-22163_task-HRFPRF_run-1_space-fsnative_hemi-L.func.gii
sub-ex22163_ses-22163_task-HRFPRF_run-2_space-fsnative_hemi-L.func.gii
sub-ex22163_ses-22163_task-shuffledPRF_run-1_space-fsnative_hemi-L.func.gii
sub-ex22163_ses-22163_task-shuffledScotomaPRF_run-1_space-fsnative_hemi-L.func.gii
sub-ex22163_ses-22163_task-shuffledScotomaPRF_run-2_space-fsnative_hemi-L.func.gii
sub-ex22163_ses-22163_task-sweepPRF_run-1_space-fsnative_hemi-L.func.gii
sub-ex22163_ses-22163_task-sweepScotomaPRF_run-1_space-fsnative_hemi-L.func.gii
sub-ex22163_ses-22163_task-sweepScotomaPRF_run-2_space-fsnative_hemi-L.func.gii

% The order of the stimulus should be:
HRFestimate_20200717T102902.mat
HRFestimate_20200717T103458.mat
8barsShuffled_20200717T100058.mat
8barsShuffledUpperLeftScotoma_20200717T101646.mat
8barsShuffledUpperLeftScotoma_20200717T104003.mat
8bars_20200717T102313.mat
8barsUpperLeftScotoma_20200717T101023.mat
8barsUpperLeftScotoma_20200717T104548.mat

stimMats = {'HRFestimate_20200717T102902.mat'
            '8barsShuffled_20200717T100058.mat'
            '8barsShuffledUpperLeftScotoma_20200717T101646.mat'
            '8bars_20200717T102313.mat'
            '8barsUpperLeftScotoma_20200717T101023.mat'
            };

%}
%
 
%  See also:   

% Create folder names
BIDSdir      = fullfile(basedir,'BIDS');
fmriprepDir  = fullfile(BIDSdir,'derivatives','fmriprep',['sub-' subname],['ses-' ses]);
subDir       = fullfile(BIDSdir,['sub-' subname]); if ~isfolder(subDir);mkdir(subDir);end
stimMatsPath = fullfile(BIDSdir,'stimMats');
stimulusDir  = fullfile(BIDSdir,'stimuli'); if ~isfolder(stimulusDir);mkdir(stimulusDir);end
configsDir   = fullfile(BIDSdir,'config_files'); if ~isfolder(configsDir);mkdir(configsDir);end
hemi         = {'L','R'};

% STIMULUS
% stimMats     = dir(fullfile(stimMatsPath,'*mat'));
%{
stimMats = {'HRFestimate_20200717T102902.mat'
            '8barsShuffled_20200717T100058.mat'
            '8barsShuffledUpperLeftScotoma_20200717T101646.mat'
            '8bars_20200717T102313.mat'
            '8barsUpperLeftScotoma_20200717T101023.mat'
            };
imMats   = {...
            'HRFestimate.mat'
            '8barsShuffled.mat'
            '8barsShuffledUpperLeftScotoma.mat'
            '8bars.mat'
            '8barsUpperLeftScotoma.mat'
            };
%}
stimMats = {...
            '8barsShuffled_20200717T100058.mat'
            '8barsShuffledUpperLeftScotoma_20200717T101646.mat'
            '8bars_20200717T102313.mat'
            '8barsUpperLeftScotoma_20200717T101023.mat'
            };
imMats   = {...
            '8barsShuffled.mat'
            '8barsShuffledUpperLeftScotoma.mat'
            '8bars.mat'
            '8barsUpperLeftScotoma.mat'
            };
stimMats = {'HRFestimate_20200717T102902.mat'};
imMats   = {'HRFestimate.mat'};

funcFiles    = {};
for nh=[1:2]
    h=hemi{nh};
    funcFiles{nh}  = dir(fullfile(fmriprepDir,'func',['*' h '.func.gii']));
end
% todo: use logic to create this
% toAverage    = [1,2; 3,0; 4,5; 6,0; 7,8];
% Don't do HRF
toAverage    = [3,0; 4,5; 6,0; 7,8];
toAverage    = [1,2];
% Gifti 2 nifti2
% Remove vistasoft's gifti from the path and use the latest one
rmpath(genpath('/data/localhome/glerma/toolboxes/vistasoft/fileFilters/gifti-1.6'));
% avgdata = cell(2,size(toAverage,1));
newNames     = [];
newStimNames = [];
configNames  = [];
for nh=[1:2]
    h=hemi{nh};
    for nr=1:size(toAverage,1)
        whichOnes = toAverage(nr,:);
        run1 = fullfile(fmriprepDir,'func',funcFiles{nh}(whichOnes(1)).name);
        if isfile(run1) 
            r1  = gifti(run1);
            r1d = r1.cdata;
        end
        if whichOnes(2) > 0
            run2 = fullfile(fmriprepDir,'func',funcFiles{nh}(whichOnes(2)).name);
            if isfile(run2)
                r2  = gifti(run2);
                r2d = r2.cdata; 
            end
            % Check lenght
            lr1 = size(r1.cdata,2);
            lr2 = size(r2.cdata,2);
            if lr1 > lr2; r1d = r1d(:,1:lr2); end
            if lr2 > lr1; r2d = r2d(:,1:lr1); end
            % Obtain the mean of both runs
            avgdata    = (r1d + r2d)/2;
        else
            avgdata    = r1d;
        end
    
        % Filter with a label (or not, for now)

        % Write the nifti2 file in the correct location and name
        [oldLoc,oldName] = fileparts(run1);
        oldName          = replace(oldName,'.func','');
        % Obtain new session name
        components  = split(oldName,'_');
        task = split(components{3},'-');task = task{2};
        components  = components([3,5,6]);
        for nc=1:length(components)
            components{nc} = replace(components{nc},'-','');
        end
        % HARDCODE T01 for ses, right now they have the same as sub in ses's place
        ses    = 'T01';
        rest   = replace(join(components),' ','');
        newSes = [ses rest{1}];
        % Generate the new name of then nifti file, follow the pattern
        newName     = ['sub-' subname '_ses-' newSes '_task-prf_acq-normal_run-01_bold.nii.gz'];
        newNameDir  = fullfile(subDir,['ses-' newSes],'func');
        if ~isfolder(newNameDir); mkdir(newNameDir); end
        newNamePath = fullfile(newNameDir, newName);

        % Do the stimulus thing now
        % We need to do it only once, not separated in hemis
        % Well, the sessio names includes the L and R, so let's create them
        stimFile    = fullfile(stimMatsPath,stimMats{nr});
        imFile      = fullfile(stimMatsPath,imMats{nr});
        stimNii     = ['sub-' subname '_ses-' newSes '_task-prf_apertures.nii.gz'];
        stimNiiPath = fullfile(stimulusDir,stimNii);
        % Read and modify the stim file
        S = load(stimFile);
        % Assert tr is ok
        trparams = S.params.tr;
        assert(tr==trparams)
        % Assert length timepoints is ok
        timepointparams = S.params.scanDuration/tr;
        if ~isequal(size(avgdata,2),timepointparams)
            avgdata = avgdata(:,1:timepointparams);
        end
        
        % Now we can save the nifti with the BOLD
        avgdata = reshape(avgdata,size(avgdata,1),1,1,size(avgdata,2));
        ni = niftiCreate('data', avgdata, 'fname',newNamePath, 'tr',tr);
        ni.version=2;
        % ni.dim    = size(avgdata);
        % ni.pixdim=[1, 1, 1,  tr];
        niftiWrite(ni);
        newNames = [newNames; {newName}];
        
        
        % prescanDuration: 12
        % period: 192
        % 12+192=204secs, == 102 time points. 
        % So, for 12secs = 6 volumes, there are no images
        % Then we have the 8 bars inside the rest 192/2=96 timepoints, = 96
        % images
        % We get what image from the params file
        % Check the shown image every 32 images
        steps = size(S.stimulus.seq)/timepointparams;
        weShowedTheseImages = S.stimulus.seq(1:steps(1):end);
        % Read the image matrix
        im = load(imFile);
        videoStim=im.stimulus.images(:,:,weShowedTheseImages);
        % convert it to apertures, so: binarize
        videoStim(videoStim==255) = 1;
        videoStim(videoStim==0)   = 1;
        videoStim(videoStim==128) = 0;
        % Resize, normalize and binarize
        temp = zeros(101, 101, size(videoStim,3));
        for p=1:size(videoStim, 3)
            temp(:,:,p) = imresize(videoStim(:,:,p),[101 101],'cubic');
        end
        videoStim = temp;
        % Save the file
        nis = niftiCreate('data', videoStim, 'fname',stimNiiPath, 'tr',tr);
        nis.version=2;
        nis.pixdim(end)=tr;
        niftiWrite(nis);
        newStimNames = [newStimNames; {stimNii}];

        % Create and save the events file
        eventsName        = ['sub-' subname '_ses-' newSes '_task-prf_events.tsv'];
        eventsNamePath    = fullfile(newNameDir,eventsName);
        eventsTable       = table();
        eventsTable.onset = [0:timepointparams-1]';
        eventsTable.duration  = tr * ones(timepointparams,1);
        eventsTable.stim_file = categorical(repmat({stimNii},[timepointparams,1]));
        eventsTable.stim_file_index = [1:timepointparams]';
        writetable(eventsTable,eventsNamePath,'Delimiter','\t','FileType','text')
        
        
        
        
        
        
        % Create the config file
        s.subjectName=subname;
        s.sessionName=newSes;
        s.solver=solver;
        s.isPRFSynthData=false;
        s.options.model='css';
        s.options.grid=false,
        s.options.wsearch='coarse to fine';
        s.options.detrend=1;
        s.options.keepAllPoints=false;
        s.options.numberStimulusGridPoints=50;
        s.stimulus.stimulus_diameter=stimDiam;
        % Encode json
        jsonString=jsonencode(s);
        % Format a little bit
        jsonString = strrep(jsonString, ',', sprintf(',\n'));
        jsonString = strrep(jsonString, '[{', sprintf('[\n{\n'));
        jsonString = strrep(jsonString, '}]', sprintf('\n}\n]'));
        % Write it
        jsonsname=['prfanalyze-' s.solver(1:end-1) '-config_sub-' subname ...
                   '_ses-' newSes '_solver-' s.solver '.json'];
        jsonsnamePath=fullfile(configsDir, jsonsname);
        fid = fopen(jsonsnamePath,'w');if fid == -1,error('Cannot create JSON file');end
        fwrite(fid, jsonString,'char');fclose(fid);
        if ~exist(jsonsnamePath,'file'),error('Could not create output json file %s', jsonsnamePath);end  
        configNames = [configNames; {jsonsname}];
    end
end

% From aparc, extract
% lingual, pericalcarine, fusiform, cuneus, lateraloccipital,
% superiorparietal, inferiorparietal

% Reduce files by applying ROIs

% Create the final file in the Nifti dir

%  Create the config file
end

