function [compTable,bylabelsums] = pmEllipse_loadExpData(proj,tools,subs,ses,run)

% pmEllipse_loadExpData Loads Experimental Data
%   Loads experimental data from the HCP 7T 3 subjects for the Ellipse paper


% Filter data:
% Size: both radius at least 1 deg
% Eccen: bigger than 2deg, less than 6

% TASKS:
% Download new data
% 3 subject mixed
% left and right mixed
% - plot size vs aspect as well
% - make range a little bit biiger to find errors predicted by synthetic
% - spearate V1v-d from V2v-d
% plot histograms of aspect ratios
% think about the limitations we add to the data?
% V1-V2 should be around 1-2 deg, check kendrick kay and noah
% filter by R2, bigger than 45%
% check radiality: correlation plot between angle and Theta, if it is organized
% around zero, there is a relationship (the difference is zero)

fprintf('\n\nLoading experimental HCP 7T data')

%{
proj   = 'realdata';
% tools  = {'afni6','afni4','vista6','vista4'};
tools  = {'vista6','vista4'};
subs   = {'115017','164131','536647'}; 
ses    = '01';
run    = '01';
%}


% Read labels from FS (it will read them from fsaverage)
% setenv('SUBJECTS_DIR','/Applications/freesurfer/subjects')
% subjectsDir = fullfile(pmRootPath,'local','realdata','BIDS','derivatives','freesurfer');
% fsavgDir    = fullfile(subjectsDir,'fsaverage');
% if isfolder(fsavgDir)
%     setenv('SUBJECTS_DIR', subjectsDir)
% else
%     error('fsaverage missing in SUBJECTS_DIR: %s', subjectsDir)
% end

hemis  = {'lh','rh'};
% The mask labels where copied to the fsaverage folder, but they are in
% prfmodel/scripts/realDataAnalysis/
labels = {'V1_exvivo','V2_exvivo', ...
          'V1_exvivo.thresh','V2_exvivo.thresh', 'dorsalV1V2V3Mask'};
kklab = struct();
path2label = fullfile(pmRootPath,'scripts','readDataAnalysis','fsaverageLabels');
for nh=[1,2]; for nl=1:length(labels)
    labelNameWithPath = fullfile(path2label,[hemis{nh} '.' labels{nl} '.label']);
    tmp = pmFSread_label('fsaverage',labelNameWithPath,true);
    kklab.(hemis{nh}).(strrep(labels{nl},'.','_')) = tmp(:,1) + 1;  % fs is a zero based index
    
end;end
lhben   = MRIread(fullfile(pmRootPath,'scripts/readDataAnalysis/fsaverageLabels','lh.benson14_varea.v4_0.mgz'));
rhben   = MRIread(fullfile(pmRootPath,'scripts/readDataAnalysis/fsaverageLabels','rh.benson14_varea.v4_0.mgz'));
ltmpmgh = lhben;
rtmpmgh = rhben;
kklab.lh.bensonV1V2V3 = lhben.vol;
kklab.rh.bensonV1V2V3 = rhben.vol;
% Now we need to obtain the final labels, separating dorsal and ventral as well
kklab.lh.V1  = find(kklab.lh.bensonV1V2V3==1);
kklab.lh.V2  = find(kklab.lh.bensonV1V2V3==2);
kklab.lh.V3  = find(kklab.lh.bensonV1V2V3==3);
kklab.rh.V1  = find(kklab.rh.bensonV1V2V3==1);
kklab.rh.V2  = find(kklab.rh.bensonV1V2V3==2);
kklab.rh.V3  = find(kklab.rh.bensonV1V2V3==3);
% Separated dorsal
il1=ismember(kklab.lh.V1,kklab.lh.dorsalV1V2V3Mask); kklab.lh.V1d = kklab.lh.V1(il1);
il2=ismember(kklab.lh.V2,kklab.lh.dorsalV1V2V3Mask); kklab.lh.V2d = kklab.lh.V2(il2);
il3=ismember(kklab.lh.V3,kklab.lh.dorsalV1V2V3Mask); kklab.lh.V3d = kklab.lh.V3(il3);
ir1=ismember(kklab.rh.V1,kklab.lh.dorsalV1V2V3Mask); kklab.rh.V1d = kklab.rh.V1(ir1);
ir2=ismember(kklab.rh.V2,kklab.lh.dorsalV1V2V3Mask); kklab.rh.V2d = kklab.rh.V2(ir2);
ir3=ismember(kklab.rh.V3,kklab.lh.dorsalV1V2V3Mask); kklab.rh.V3d = kklab.rh.V3(ir3);
% Separated ventral
kklab.lh.V1v = kklab.lh.V1(~il1);
kklab.lh.V2v = kklab.lh.V2(~il2);
kklab.lh.V3v = kklab.lh.V3(~il3);
kklab.rh.V1v = kklab.rh.V1(~ir1);
kklab.rh.V2v = kklab.rh.V2(~ir2);
kklab.rh.V3v = kklab.rh.V3(~ir3);

% Read the data, but we need to get one with all the data together, and another
% one with regions separated. 

% TODO: WORST CODE EVER BELOW
% JUST ADD TEH WHOLE fsaverage and then filter it...    

compTable   = struct();
bylabel     = struct();
bylabelsums = struct();
for nt=1:length(tools)
    tool  = tools{nt};
    % disp(tool)   
    subestimates   = table();
    bylabel.(tool) = struct();
    % Do the links here to have life easier afterwards
    bylabelsums.(tool).V1    =  table();
    bylabelsums.(tool).V2    =  table();
    bylabelsums.(tool).V3    =  table();
    % lh-rh-Dorsal
    bylabelsums.(tool).V1d   =  table();
    bylabelsums.(tool).V2d   =  table();
    bylabelsums.(tool).V3d   =  table();
    % lh-rh-Ventral
    bylabelsums.(tool).V1v   =  table();
    bylabelsums.(tool).V2v   =  table();
    bylabelsums.(tool).V3v   =  table();
    for ns=1:length(subs)
        sub   = subs{ns};
        % disp(sub)
        p = fullfile(pmRootPath,'local',proj,'BIDS','derivatives',['prfanalyze-' tool],['sub-' sub],['ses-' ses]);
        bylabel.(tool).(['s' sub]) = struct();
        % Copy the code from prfReportWrapper to create the tables
        %% Read the results back to Matlab
        % Read all the nifti-s for this type of tool 
        cd(p)
        niftis = dir(['*run-' run '*.nii.gz']);
        if length(niftis) < 6
            error('Not enough nifti result files in output, at least there should be x0,y0,sigmamajor,sigmaminor,theta and modelpred')
        end
        filesRead = 0;

        pmEstimates = table();
        bylabel.(tool).(['s' sub]).lhV1_exvivo =  table();
        bylabel.(tool).(['s' sub]).lhV2_exvivo =  table();
        bylabel.(tool).(['s' sub]).rhV1_exvivo =  table();
        bylabel.(tool).(['s' sub]).rhV2_exvivo =  table();
        
        % Benson V1-V2-V3
        bylabel.(tool).(['s' sub]).lhV1  =  table();
        bylabel.(tool).(['s' sub]).lhV2  =  table();
        bylabel.(tool).(['s' sub]).lhV3  =  table();
        bylabel.(tool).(['s' sub]).rhV1  =  table();
        bylabel.(tool).(['s' sub]).rhV2  =  table();
        bylabel.(tool).(['s' sub]).rhV3  =  table();
        % Dorsal
        bylabel.(tool).(['s' sub]).lhV1d =  table();
        bylabel.(tool).(['s' sub]).lhV2d =  table();
        bylabel.(tool).(['s' sub]).lhV3d =  table();
        bylabel.(tool).(['s' sub]).rhV1d =  table();
        bylabel.(tool).(['s' sub]).rhV2d =  table();
        bylabel.(tool).(['s' sub]).rhV3d =  table();
        % Ventral
        bylabel.(tool).(['s' sub]).lhV1v =  table();
        bylabel.(tool).(['s' sub]).lhV2v =  table();
        bylabel.(tool).(['s' sub]).lhV3v =  table();
        bylabel.(tool).(['s' sub]).rhV1v =  table();
        bylabel.(tool).(['s' sub]).rhV2v =  table();
        bylabel.(tool).(['s' sub]).rhV3v =  table();
        % Do the links here to have life easier afterwards
        bylabel.(tool).(['s' sub]).V1    =  table();
        bylabel.(tool).(['s' sub]).V2    =  table();
        bylabel.(tool).(['s' sub]).V3    =  table();
        % lh-rh-Dorsal
        bylabel.(tool).(['s' sub]).V1d   =  table();
        bylabel.(tool).(['s' sub]).V2d   =  table();
        bylabel.(tool).(['s' sub]).V3d   =  table();
        % lh-rh-Ventral
        bylabel.(tool).(['s' sub]).V1v   =  table();
        bylabel.(tool).(['s' sub]).V2v   =  table();
        bylabel.(tool).(['s' sub]).V3v   =  table();
        
        for nn=1:length(niftis)
            fname = niftis(nn).name;
            % Assume always it will be .nii.gz and that the . has not been used in the filename
            trunkname = split(fname,'.');
            if length(trunkname) ~= 3; error('%s should be a filename with no dots and .nii.gz',fname);end	
            trunkname = trunkname{1};
            resname   = split(trunkname, '_');
            resname   = resname{end};
            % read the nifti and squeeze the result matrix
            % fprintf('Attempting to read %s\n',fullfile(pwd,fname))
            tmp       = niftiRead(fname);
            data      = squeeze(tmp.data);	
            % Data can be separated in different labels
            lhtmpfs   = zeros(1,163842);
            rhtmpfs   = zeros(1,163842);
            % Assign our results to a real brain
            lV1l = length(kklab.lh.V1_exvivo);
            lV2l = length(kklab.lh.V2_exvivo);
            rV1l = length(kklab.rh.V1_exvivo);
            rV2l = length(kklab.rh.V2_exvivo);
            lhtmpfs(kklab.lh.V1_exvivo) = data(1:lV1l,1)';
            lhtmpfs(kklab.lh.V2_exvivo) = data(lV1l+1 : lV1l+lV2l,1)';
            rhtmpfs(kklab.rh.V1_exvivo) = data(lV1l+lV2l+1:lV1l+lV2l+rV1l,1)';
            rhtmpfs(kklab.rh.V2_exvivo) = data(lV1l+lV2l+rV1l+1:lV1l+lV2l+rV1l+rV2l,1)';
            % Now is better if we use Noah's separation, this is in kklab.lh.bensonV1V2V3
            % Assign it to the table
            switch resname
                case 'sigmaminor', resname = 'sMin';
                case 'sigmamajor', resname = 'sMaj';
                case 'theta'     , resname = 'Th';
                case 'centerx0'  , resname = 'x0';
                case 'centery0'  , resname = 'y0';
            end
            pmEstimates.(resname) = data;	        
            
            % Now assign it to a structute with all values separated
            bylabel.(tool).(['s' sub]).lhV1_exvivo.(resname) = data(1:lV1l,1);
            bylabel.(tool).(['s' sub]).lhV2_exvivo.(resname) = data(lV1l+1 : lV1l+lV2l,1);
            bylabel.(tool).(['s' sub]).rhV1_exvivo.(resname) = data(lV1l+lV2l+1:lV1l+lV2l+rV1l,1);
            bylabel.(tool).(['s' sub]).rhV2_exvivo.(resname) = data(lV1l+lV2l+rV1l+1:lV1l+lV2l+rV1l+rV2l,1);
            % Benson V1-V2-V3
            bylabel.(tool).(['s' sub]).lhV1.(resname) = lhtmpfs(kklab.lh.V1)';
            bylabel.(tool).(['s' sub]).lhV2.(resname) = lhtmpfs(kklab.lh.V2)';
            bylabel.(tool).(['s' sub]).lhV3.(resname) = lhtmpfs(kklab.lh.V3)';
            bylabel.(tool).(['s' sub]).rhV1.(resname) = rhtmpfs(kklab.rh.V1)';
            bylabel.(tool).(['s' sub]).rhV2.(resname) = rhtmpfs(kklab.rh.V2)';
            bylabel.(tool).(['s' sub]).rhV3.(resname) = rhtmpfs(kklab.rh.V3)';
            % Dorsal
            bylabel.(tool).(['s' sub]).lhV1d.(resname) = lhtmpfs(kklab.lh.V1d)';
            bylabel.(tool).(['s' sub]).lhV2d.(resname) = lhtmpfs(kklab.lh.V2d)';
            bylabel.(tool).(['s' sub]).lhV3d.(resname) = lhtmpfs(kklab.lh.V3d)';
            bylabel.(tool).(['s' sub]).rhV1d.(resname) = rhtmpfs(kklab.rh.V1d)';
            bylabel.(tool).(['s' sub]).rhV2d.(resname) = rhtmpfs(kklab.rh.V2d)';
            bylabel.(tool).(['s' sub]).rhV3d.(resname) = rhtmpfs(kklab.rh.V3d)';
            % Ventral
            bylabel.(tool).(['s' sub]).lhV1v.(resname) = lhtmpfs(kklab.lh.V1v)';
            bylabel.(tool).(['s' sub]).lhV2v.(resname) = lhtmpfs(kklab.lh.V2v)';
            bylabel.(tool).(['s' sub]).lhV3v.(resname) = lhtmpfs(kklab.lh.V3v)';
            bylabel.(tool).(['s' sub]).rhV1v.(resname) = rhtmpfs(kklab.rh.V1v)';
            bylabel.(tool).(['s' sub]).rhV2v.(resname) = rhtmpfs(kklab.rh.V2v)';
            bylabel.(tool).(['s' sub]).rhV3v.(resname) = rhtmpfs(kklab.rh.V3v)';
            % Do the links here to have life easier afterwards
            bylabel.(tool).(['s' sub]).V1.(resname) = [lhtmpfs(kklab.lh.V1),rhtmpfs(kklab.rh.V1)]';
            bylabel.(tool).(['s' sub]).V2.(resname) = [lhtmpfs(kklab.lh.V2),rhtmpfs(kklab.rh.V2)]';
            bylabel.(tool).(['s' sub]).V3.(resname) = [lhtmpfs(kklab.lh.V3),rhtmpfs(kklab.rh.V3)]';
            % lh-rh-Dorsal
            bylabel.(tool).(['s' sub]).V1d.(resname) = [lhtmpfs(kklab.lh.V1d),rhtmpfs(kklab.rh.V1d)]';
            bylabel.(tool).(['s' sub]).V2d.(resname) = [lhtmpfs(kklab.lh.V2d),rhtmpfs(kklab.rh.V2d)]';
            bylabel.(tool).(['s' sub]).V3d.(resname) = [lhtmpfs(kklab.lh.V3d),rhtmpfs(kklab.rh.V3d)]';
            % lh-rh-Ventral
            bylabel.(tool).(['s' sub]).V1v.(resname) = [lhtmpfs(kklab.lh.V1v),rhtmpfs(kklab.rh.V1v)]';
            bylabel.(tool).(['s' sub]).V2v.(resname) = [lhtmpfs(kklab.lh.V2v),rhtmpfs(kklab.rh.V2v)]';
            bylabel.(tool).(['s' sub]).V3v.(resname) = [lhtmpfs(kklab.lh.V3v),rhtmpfs(kklab.rh.V3v)]';
        end
        % Concatenate different subjects
        subestimates = [subestimates;pmEstimates];
        
        bylabelsums.(tool).V1    =  [bylabelsums.(tool).V1; bylabel.(tool).(['s' sub]).V1];
        bylabelsums.(tool).V2    =  [bylabelsums.(tool).V2; bylabel.(tool).(['s' sub]).V2];
        bylabelsums.(tool).V3    =  [bylabelsums.(tool).V3; bylabel.(tool).(['s' sub]).V3];
        % lh-rh-Dorsal
        bylabelsums.(tool).V1d   =  [bylabelsums.(tool).V1d; bylabel.(tool).(['s' sub]).V1d];
        bylabelsums.(tool).V2d   =  [bylabelsums.(tool).V2d; bylabel.(tool).(['s' sub]).V2d];
        bylabelsums.(tool).V3d   =  [bylabelsums.(tool).V3d; bylabel.(tool).(['s' sub]).V3d];
        % lh-rh-Ventral
        bylabelsums.(tool).V1v   =  [bylabelsums.(tool).V1v; bylabel.(tool).(['s' sub]).V1v];
        bylabelsums.(tool).V2v   =  [bylabelsums.(tool).V2v; bylabel.(tool).(['s' sub]).V2v];
        bylabelsums.(tool).V3v   =  [bylabelsums.(tool).V3v; bylabel.(tool).(['s' sub]).V3v];
    end
    compTable.(tool) = subestimates;
    compTable.(tool).R2 = calccod(compTable.(tool).modelpred', compTable.(tool).testdata')';
end

disp ('... done with load')

end

