function DTcalc = pmForwardModelCalculate(DTDT,varargin)
% Creates a table with all parameters (defaults) required to perform a forward
% calculation
% To work with big tables, it can use parfor
% 
%  Inputs: table with one or several rows of parameters to calculate bold series
%          in fwd model
% 
%  Outputs: returns the same table with the last column of pm-s calculated
% 
%  See also: forwardModelTableCreate
% 
%  GLU Vistalab 2019.05


%% Read the inputs
varargin = mrvParamFormat(varargin);
p = inputParser;
p.addRequired('DTDT' ,    @(x)(isa(x,'table')));
% By default do not create multiple copies
p.addParameter('useparallel', false, @islogical);
p.addParameter('writefiles' , false, @islogical);
p.addParameter('outputdir'  , '', @ischar);
p.addParameter('subjectname', '', @ischar);
p.parse(DTDT, varargin{:});
% Assign it
useparallel = p.Results.useparallel;
writefiles  = p.Results.writefiles;
outputdir   = p.Results.outputdir;
subjectName = p.Results.subjectname;
if writefiles && (isempty(outputdir) || isempty(subjectName))
    error('If writefiles is set to true, then outputdir and subjectName must be provided')
end
% Asign output as well, to an empty table, depending on options we won't return the table (too large)
DTcalc = table();

if useparallel
	%%%%%   CLUSTER PARPOOL    %%%%%%
	myclusterLocal = parcluster('local');
	NumWorkers = myclusterLocal.NumWorkers;
	% [st, re] = system('qstat -g c | grep matlab.q');
	% [Tok, Rem] = strtok(re);
	% [Tok, Rem] = strtok(Rem);
	% [Tok, Rem] = strtok(Rem);
	% [Tok, Rem] = strtok(Rem);
	% [available] = strtok(Rem)
	% parpool('ips_base', str2num(available))
	%%%%% END CLUSTER PARPOOL  %%%%%%
else
	NumWorkers = 0;
end
% chksize = ceil(height(DTDT) / (NumWorkers));
% Optimize for Matlab memory problems
chksize = 3000;
if height(DTDT) < chksize
    DTcc{1}       = DTDT;  
    nchcks        = 1;
    startindex{1} = 1;
    endindex{1}   = height(DTDT);
    DTheight{1}   = height(DTDT);
else
    nchcks = ceil(height(DTDT) / chksize);
    for nn=1:nchcks
        startindex{nn} = (nn*chksize) + 1 - chksize;
        if nn == nchcks
            endindex{nn}   = height(DTDT);
        else
            endindex{nn}   = nn*chksize;
        end
        DTcc{nn}     = DTDT(startindex{nn}:endindex{nn},:);
        DTheight{nn} = length(startindex{nn}:endindex{nn});
    end
end
fprintf('There are %d voxels \n', height(DTDT))

%% Create temp folders
tmpName = tempname(fullfile(pmRootPath,'local'));
mkdir(tmpName);
        
tic
parfor (nn=1:nchcks, NumWorkers)
    DT = DTcc{nn};
    % Initialize prev variables, for parallel toolbox
    dtprev = [];
    pmprev = [];
    for ii=1:DTheight{nn}
        
        %%  Do it row to row and parameter to parameter first, for debugging
        dispat = 100; if height(DT)<dispat;dispat=2;end
        if mod(ii,dispat)==0
            disp([num2str(nn) ' -- ' num2str(ii) ' -- ' num2str(height(DT))])
        end
        %% Initialize the basic model with defaults
        dt     = DT(ii,:);
        % we need a fresh copy of the prfModel class here, otherwise it references
        % the same one and changes are not persistent
        pm     = prfModel;
        
        %% High Level Variables
        isprfmodel = @(x)(isa(x,'prfModel'));
        for vn = dt.Properties.VariableNames
            if ~istable(dt.(vn{:})) && ~isprfmodel(dt.(vn{:}))
                if iscell(dt.(vn{:}))
                    pm.(vn{:}) = dt.(vn{:}){:};
                else
                    pm.(vn{:}) = dt.(vn{:});
                end
            end
        end
        
        %% Stimulus
        for jj=1:width(dt.Stimulus)
                paramName               = dt.Stimulus.Properties.VariableNames{jj};
                pm.Stimulus.(paramName) = dt.Stimulus.(paramName);
        end
        
        % Check if did not change from the previous line, calculate again or copy the old one
        if ii > 1
            if isequal(dtprev.Stimulus, dt.Stimulus)
                pm.Stimulus.userVals = pmprev.Stimulus.userVals;
            else
                pm.Stimulus.compute;
                pm.Stimulus.userVals = pm.Stimulus.getStimValues;
            end
        else
            pm.Stimulus.compute;
            pm.Stimulus.userVals = pm.Stimulus.getStimValues;
        end
        
        %% RF
        for jj=1:width(dt.RF)
            paramName         = dt.RF.Properties.VariableNames{jj};
            pm.RF.(paramName) = dt.RF.(paramName);
        end
        pm.RF.compute;
        
        %% HRF
        for jj=1:width(dt.HRF)
            paramName          = dt.HRF.Properties.VariableNames{jj};
            val                = dt.HRF.(paramName);
            if iscell(val)
                pm.HRF.(paramName) = val{:};
            else
                pm.HRF.(paramName) = val;
            end
        end
        pm.HRF.compute;
        
        %% Noise
        for jj=1:width(dt.Noise)
            paramName            = dt.Noise.Properties.VariableNames{jj};
            % voxel: we did not create it as a variable, as it is just to specify a
            % series of noise defaults. I don't like to have specific variable
            % names here, think about making it a noise param
            if strcmp(paramName,'seed')
                [val, status] = str2num(dt.Noise.(paramName){:});
                if status
                    pm.Noise.(paramName) = val;
                else
                    pm.Noise.(paramName) = dt.Noise.(paramName);
                end
            elseif ~strcmp(paramName,'voxel')
                pm.Noise.(paramName) = dt.Noise.(paramName);
            end
            
        end
        pm.Noise.compute;
        
        %% Compute the synthetic signal
        % The compute at the top level computes all the lovel level ones.
        % Just do it once here.
        % By default, pm.compute computes all subclasses. 
        % Here, we want to to calculate or copy from the previous one to save
        % time. Therefore we will set pm.computeSubclasses=false and it will
        % only calculate the last step.
        pm.computeSubclasses = false;
        pm.compute;
        
        %% Assign it to the cell array (or Write back the updated pm model)
        DT.pm(ii) = pm;
        
        %% Save as the previous one
        dtprev = dt;
        pmprev = copy(pm);  % If we don't copy it it is just a pointer
        
        
    end
    %% Write the result to file
    if writefiles
    	% BOLD FILE
    	fname = fullfile(tmpName, sprintf('%s_%04i.nii.gz', subjectName,nn));
    	pmForwardModelToNifti(DT, 'fname',fname, 'demean',false);
    
    	% JSON FILE
    	jsonSynthFile = fullfile(tmpName, sprintf('%s_%04i.json', subjectName,nn));
    	% Encode json
    	jsonString = jsonencode(DT(:,1:(end-1)));
    	% Format a little bit
    	jsonString = strrep(jsonString, ',', sprintf(',\n'));
    	jsonString = strrep(jsonString, '[{', sprintf('[\n{\n'));
    	jsonString = strrep(jsonString, '}]', sprintf('\n}\n]'));
    	% Write it
    	fid = fopen(jsonSynthFile,'w');if fid == -1,error('Cannot create JSON file');end
    	fwrite(fid, jsonString,'char');fclose(fid);
    	
    	% STIM FILE
        % There should be only one, just write it once
    	stimNiftiFname = fullfile(tmpName, sprintf('%s_Stim.nii.gz', subjectName));
    	pm1            = DT.pm(1);
        if nn==1
    	    stimNiftiFname = pm1.Stimulus.toNifti('fname',stimNiftiFname);
        end
    else
        % we dont want to write .mat files, they can be too large, write the definite thing and later concatenate
        fName = fullfile(tmpName, sprintf('tmpDT_%04i.mat',nn));
        m     = matfile(fName,'writable',true);
        m.DT  = DT;
    end
end

toc
disp('The parfor par has been finished')

%% Concatenate back, either files or final table
if writefiles
    disp('Concatenating the nifti and json files back')
    % BOLD files
    fname1 = fullfile(tmpName, sprintf('%s_%04i.nii.gz', subjectName,1));
    BOLDnifti = niftiRead(fname1);
    BOLDdata  = BOLDnifti.data;
    % json file
    jsonSynthFile1 = fullfile(tmpName, sprintf('%s_%04i.json', subjectName,1));
    jsonfile       = jsonread(jsonSynthFile1);
    % Stim file: there should be just one
    stimNiftiFname = fullfile(tmpName, sprintf('%s_Stim.nii.gz', subjectName));
    for nn=2:nchcks
        % BOLD files
    	tmpboldname  = fullfile(tmpName, sprintf('%s_%04i.nii.gz', subjectName,nn));
        tmpboldnifti = niftiRead(tmpboldname);
        tmpbolddata  = tmpboldnifti.data;
        BOLDdata     = [BOLDdata;tmpbolddata];
        % json files    
        tmpjsonname  = fullfile(tmpName, sprintf('%s_%04i.json',subjectName,nn));
        tmpjson      = jsonread(tmpjsonname);
        jsonfile     = [jsonfile;tmpjson]; 
    end
    % BOLD
    BOLDnifti.data   = BOLDdata;
    BOLDnifti.dim    = size(BOLDdata);
    BOLDnifti.fname  = fullfile(outputdir, sprintf('%s.nii.gz', subjectName));
    niftiWrite(BOLDnifti);
    if ~exist(BOLDnifti.fname,'file'),error('Could not create output BOLD %s in %s', BOLDnifti.fname,outputdir);end

    % json
    jsonfname = fullfile(outputdir, sprintf('%s.json', subjectName));
    % Encode json
    jsonString = jsonencode(jsonfile);
    % Format a little bit
    jsonString = strrep(jsonString, ',', sprintf(',\n'));
    jsonString = strrep(jsonString, '[{', sprintf('[\n{\n'));
    jsonString = strrep(jsonString, '}]', sprintf('\n}\n]'));
    % Write it
    fid = fopen(jsonfname,'w');if fid == -1,error('Cannot create JSON file');end
    fwrite(fid, jsonString,'char');fclose(fid);
    if ~exist(jsonfname,'file'),error('Could not create output json file %s in %s', jsonfname, outputdir);end
        
    % Stim: is always de same
    succ = movefile(stimNiftiFname,outputdir);
    if ~succ;error('Could not move %s to %s', stimNiftiFname,outputdir);end

else
    disp('Concatenating the .mat files back')
    for nn=1:nchcks
        fName = fullfile(tmpName, sprintf('tmpDT_%04i.mat',nn));
        tmp   = load(fName,'DT');
        DTcalc = [DTcalc; tmp.DT];
    end
end


%% Remove tmp dir
rmdir(tmpName,'s')

end
