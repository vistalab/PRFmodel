function [pmEstimates, results] = pmModelFit(input, prfimplementation, varargin)
% Select and apply a PRF model to estimate model parameters
% 
% Syntax:
%    [pmEstimates, results] = pmModelFit(inputTable, prfImplementation)
%  
% Brief description:
%    Takes a table whose rows describe the BOLD time series in a
%    single voxel (along with its metadata) and estimates the pRF
%    model parameters for some implementation.
%
% Inputs:
%   input       - It can be a Data table including all the required info 
%                 (stimulus, all parameters and BOLD time series). Each row of 
%                 the data table is a voxel.
%               - It can be a cell array with three paths to physical file
%                 locations: 4D-nifti with synthetic bold series, 3D-nifti with
%                 stimuli info, json file with all the params for all the
%                 voxels. 
%                 -- NOTE: if there is only one entry in the json file, the
%                          same params apply to all voxels in the nifti file. 
%   prfImplementation - String defining the model
%
% Outputs: 
%   pmEstimates: Table format of the pRF model parameters in results
%   results:     The struct from analyzePRF
%
% Key/val parameters (Optional)
%   N/A
%
% GLU Vistalab 05.2019
%
% See also:
%     pmXXX
%

% 

% Examples:
%{
  pmCompute...
%}

%% Read the inputs
% Make varargin lower case, remove white spaces...
varargin = mrvParamFormat(varargin);
% Parse
p = inputParser;
p.addRequired('input');
p.addRequired('prfimplementation',@ischar);
p.addParameter('useparallel'    ,  true        , @islogical);
    options       = struct();
    options.aprf  = struct('seedmode', [0 1 2], ...
                           'display' , 'off'  , ...
                           'usecss'  , true  );
    options.vista = struct('model'        ,'one gaussian'   , ...
                           'grid'         , false           , ...
                           'wsearch'      , 'coarse to fine', ...
                           'detrend'      , 1               , ...
                           'keepAllPoints', false           , ...
                           'obtainPreds'  , false           , ...
                           'numberStimulusGridPoints',   50);
    options.afni  = struct('model','afni4', ...
                           'hrf'  , 'SPM');
    options.mlr   = struct('quickFit', 0 , ...
                           'doParallel'  , 0, ...
                           'rfType','gaussian');
p.addParameter('options', options, @isstruct);
% Parse. Assign result inside each case
p.parse(input,prfimplementation,varargin{:});
% Read here only the generic ones
useParallel = p.Results.useparallel;
allOptions  = p.Results.options;
% We need to be sure that if only some of the params are passed, the rest will
% be taken from the defaults 
allOptions  = pmParamsCompletenessCheck(allOptions, options);

disp('This is the options file comming in:')
allOptions
allOptions.vista

%% Convert the input
% Check if the pm-s come in a table or alone
% Read, preprocess, and check the values

if iscell(input)
    niftiInputs = true;
    % Extract filenames
    BOLDname = input{1};
    JSONname = input{2};
    STIMname = input{3};
    % Check the files
    if ~exist(BOLDname,'file'),error('Cannot find %s',BOLDname),end
    if ~exist(JSONname,'file'),error('Cannot find %s',JSONname),end
    if ~exist(STIMname,'file'),error('Cannot find %s',STIMname),end
    
    % All three files come from the same synthDT table. 
    % Table is separated into a nifti with the bold series and TR info,
    % Into a nifti with the stimuli used, same for every bold tSeries, 
    % and a json file with all the parameters. Some parameters have to be the
    % same for all tSeries (TR, for example), but other are specific (x0y0
    % locations or rfSizes, for example).
    
    % 1: nifti with the synthetic BOLD series
    data     = niftiRead(BOLDname);
    TRdata   = data.pixdim(end);
    % Synthetic datasets larger than 32000 rows come in more dimensions
    % Convert it to 2D again concatenating
    % [s1,s2,s3,s4] = size(data.data);
    % if s1>1 && s2 == 1 && s3==1 && s4>1
        data = squeeze(data.data);
    % elseif s1>1 && s2 > 1 && s3==1 && s4>1
    %     % A reshape might do, but just wanted control...
    %     tmp = zeros(s1*s2,s4);
    %     for jj=1:s2
    %         tmp(((jj*s1) - s1 + 1):(jj*s1),:) = squeeze(data.data(:,jj,s3,:));
    %     end
    %     % We will have many 0 time series, delete them
    %     nonzeroindx = ~(sum(tmp,2) == 0);
    %     data = tmp(nonzeroindx,:);
    % else
    %     error('The dimension of nifti is not recognized')
    % end
    % 2: json with the metadata
    synthDT  = struct2table(jsonread(JSONname));
    if ismember('isPRFSynthData',synthDT.Properties.VariableNames)
        TRsynth           = TRdata;
        signalPercentage  = 'bold';
        % Check that the stim diameter was passed
		stimulus_diameter = synthDT.stimulus_diameter; 

        % 3: nifti with the stimuli
        stimulus = niftiRead(STIMname);
		TRstim   = stimulus.pixdim(end);
        stimulus = squeeze(stimulus.data);
        % Check dimensions
        if ~isequal(size(data,2), size(stimulus,3))
            error('Data and stimulus have different time points')
        end

        % Obtain the main parameters required for analysis
        % TR
        if isclose(TRdata,TRstim,'tolerance',0.001) && isclose(TRstim,TRsynth,'tolerance',0.001)
            TR = TRdata;
        else
            error('Data and stimulus have different TR')
        end
        % Stimulus related
        Stimulus.ResizedHorz = size(stimulus,2);
        Stimulus.ResizedVert = size(stimulus,1);
        Stimulus.fieldofviewHorz = stimulus_diameter;
        Stimulus.fieldofviewVert = stimulus_diameter;
        Stimulus.spatialSampleHorz = Stimulus.fieldofviewHorz / Stimulus.ResizedHorz;
        Stimulus.spatialSampleVert = Stimulus.fieldofviewVert / Stimulus.ResizedVert;
        stimradius    = Stimulus.fieldofviewHorz / 2;
        % BOLD related
        BOLDmeanValue = 10000; % Not used by mrVista, fix for the other tools 
        if (size(stimulus,3) == size(data,2)),timePointsN = size(data,2), end

    else
        for na=1:width(synthDT)
            if isstruct(synthDT{:,na})
                synthDT.(synthDT.Properties.VariableNames{na}) = struct2table(synthDT{:,na});
            end
        end
        TRsynth = unique(synthDT.TR); if (length(TRsynth) ~= 1),error('More than 1 TR in synth json'),end
        signalPercentage = unique(synthDT.signalPercentage); if (length(signalPercentage) ~= 1),error('More than 1 signalPercentage in synth json'),end
        % 3: nifti with the stimuli
        stimulus = niftiRead(STIMname);
        TRstim   = stimulus.pixdim(4);
        stimulus = squeeze(stimulus.data);
        % Check dimensions
        if ~isequal(size(data,2), size(stimulus,3))
            error('Data and stimulus have different time points')
        end
        
        % Obtain the main parameters required for analysis
        % TR
        if isclose(TRdata,TRstim,'tolerance',0.001) && isclose(TRstim,TRsynth,'tolerance',0.001)
            TR = TRdata;
        else
            error('Data and stimulus have different TR')
        end
        % Stimulus related
        Stimulus.ResizedHorz = unique(synthDT.Stimulus.ResizedHorz); if (length(Stimulus.ResizedHorz) ~= 1),error('More than 1 ResizedHorz in synth json'),end
        Stimulus.ResizedVert = unique(synthDT.Stimulus.ResizedVert); if (length(Stimulus.ResizedVert) ~= 1),error('More than 1 ResizedVert in synth json'),end
        Stimulus.fieldofviewHorz = unique(synthDT.Stimulus.fieldofviewHorz);if (length(Stimulus.fieldofviewHorz) ~= 1),error('More than 1 fieldofviewHorz in synth json'),end
        Stimulus.fieldofviewVert = unique(synthDT.Stimulus.fieldofviewVert);if (length(Stimulus.fieldofviewVert) ~= 1),error('More than 1 fieldofviewVert in synth json'),end
        Stimulus.spatialSampleHorz = Stimulus.fieldofviewHorz/Stimulus.ResizedHorz;
        Stimulus.spatialSampleVert = Stimulus.fieldofviewVert / Stimulus.ResizedVert;
        stimradius    = Stimulus.fieldofviewHorz / 2;
        % BOLD related
        BOLDmeanValue = unique(synthDT.BOLDmeanValue);
        if (length(BOLDmeanValue) ~= 1),error('More than 1 BOLDmeanValue in synth json'),end
        if (size(stimulus,3) == size(data,2)),timePointsN = size(data,2), end
    end
else
    niftiInputs = false;
    if ~istable(input)
        temp    = table();
        temp.pm = input;
        input   = temp;
    end
    % Obtain the same parameters as above
    pm1      = input.pm(1);
    stimulus = double(pm1.Stimulus.getStimValues);
    % Check if the TR is the same
    if length(unique(input.TR)) > 1, error('TR has to be the same for all synthetic BOLD series'),end
    TR       = pm1.TR;
    % Check if the signalPercentage is the same
    if length(unique(input.signalPercentage)) > 1, error('signalPercentage has to be the same for all synthetic BOLD series'),end
    signalPercentage = pm1.signalPercentage;
    % Stimulus related
    Stimulus.ResizedHorz = pm1.Stimulus.ResizedHorz;
    Stimulus.ResizedVert = pm1.Stimulus.ResizedVert;
    Stimulus.fieldofviewHorz = pm1.Stimulus.fieldofviewHorz;
    Stimulus.fieldofviewVert = pm1.Stimulus.fieldofviewVert;
    Stimulus.spatialSampleHorz = pm1.Stimulus.spatialSampleHorz;
    Stimulus.spatialSampleVert = pm1.Stimulus.spatialSampleVert;
    stimradius    = pm1.Stimulus.fieldofviewHorz/2;
    % BOLD related
    BOLDmeanValue = pm1.BOLDmeanValue;
    timePointsN   = pm1.timePointsN;
    % Add the time series
    PMs        = input.pm;
    data       = repmat(nan([1,pm1.timePointsN]), [height(input),1]);
    for ii=1:height(input); data(ii,:)=PMs(ii).BOLDnoise;end
end

%% Choose the analysis case
prfimplementation = mrvParamFormat(prfimplementation);

switch prfimplementation
    case {'aprf','analyzeprf'}
        % Read the options
        options = allOptions.aprf;
        
        %% Create and check tables
        % TODO: Let's define the format for the estimates so that this is
        %       the same for all the methods.
        pmEstimates = table();           
        if useParallel
            pmEstimates.testdata = data;
            
            % Calculate PRF
            results  = analyzePRF({stimulus}, {data}, TR, options);
            
            % Prepare output
            Centerx0         = squeeze(results.params(1,2,:));
            Centery0         = squeeze(results.params(1,1,:));
            % Change to the center of the fov
            if all(iseven(Stimulus.ResizedHorz))
                Centerx0 = Centerx0 - Stimulus.ResizedHorz/2;
            else
                Centerx0 = Centerx0 - ((Stimulus.ResizedHorz-1)/2 + 1);
            end
            if iseven(Stimulus.ResizedVert)
                Centery0 = Centery0 - Stimulus.ResizedVert/2;
            else
                Centery0 = Centery0 - ((Stimulus.ResizedVert-1)/2 + 1);
            end
            pmEstimates.Centerx0   = (Stimulus.spatialSampleHorz * Centerx0);
            pmEstimates.Centery0   = (Stimulus.spatialSampleVert * Centery0);
            pmEstimates.Theta      = zeros(size(Centery0));
            pmEstimates.sigmaMinor = Stimulus.spatialSampleHorz * results.rfsize;
            pmEstimates.sigmaMajor = pmEstimates.sigmaMinor;
            pmEstimates.modelpred  = results.modelpred;
            % Scale it back
            % With the new method, the first point will be in BOLDmeanValue; 
            pmEstimates.modelpred  = pmEstimates.modelpred + BOLDmeanValue;
            if pmEstimates.modelpred(1) ~= data(1)
                pmEstimates.modelpred = pmEstimates.modelpred + ...
                                        (data(1) - pmEstimates.modelpred(1));
            end
            
        else
            if niftiInputs
                error('aPRF when using not parallel need a table or a pm, not niftis')
            end
            %% Calculate values voxel to voxel
            for ii=1:height(input)
                % TODO: use parfor if the number of rows if the table is larger than XX
                pm       = input.pm(ii);
                stimulus = double(pm.Stimulus.getStimValues);
                data     = pm.BOLDnoise;
                TR       = pm.TR;
                
                % Calculate PRF
                results  = analyzePRF({stimulus}, {data}, TR, options);
                % This is the explanation of the results
                % Select only those parameters we need to use to compare to the
                % input RF variables we manipulated.
                %{
            % The results structure contains the following fields:
            % <ang> contains pRF angle estimates.  Values range between 0 and 360 degrees.
            %   0 corresponds to the right horizontal meridian, 90 corresponds to the upper vertical
            %   meridian, and so on.
            % <ecc> contains pRF eccentricity estimates.  Values are in pixel units with a lower
            %   bound of 0 pixels.
            % <rfsize> contains pRF size estimates.  pRF size is defined as sigma/sqrt(n) where
            %   sigma is the standard of the 2D Gaussian and n is the exponent of the power-law
            %   function.  Values are in pixel units with a lower bound of 0 pixels.
            % <expt> contains pRF exponent estimates.
            % <gain> contains pRF gain estimates.  Values are in the same units of the data
            %   and are constrained to be non-negative.
            % <R2> contains R^2 values that indicate the goodness-of-fit of the model to the data.
            %   Values are in percentages and generally range between 0% and 100%.  The R^2 values
            %   are computed after projecting out polynomials from both the data and the model fit.
            %   (Because of this projection, R^2 values can sometimes drop below 0%.)  Note that
            %   if cross-validation is used (see <xvalmode>), the interpretation of <R2> changes
            %   accordingly.
            % <resnorms> and <numiters> contain optimization details (residual norms and
            %   number of iterations, respectively).
            % <meanvol> contains the mean volume, that is, the mean of each voxel's time-series.
            % <noisereg> contains a record of the noise regressors used in the model.
            % <params> contains a record of the raw parameter estimates that are obtained internally
            %   in the code.  These raw parameters are transformed to a more palatable format for
            %   the user (as described above).
            % <options> contains a record of the options used in the call to analyzePRF.m.
                %}
                
                % Add a new row of results
                % (Here we need to select what results go to the table and adecuate them to be comparable to the vakues we put into)
                % We will use this to interpret and modify the parameters.
                %{
               The first seed is a generic large pRF that is centered with
               respect to the stimulus, has a pRF size equal to 1/4th of the
               stimulus extent (thus, +/- 2 pRF sizes matches
               the stimulus extent), and has an exponent of 0.5.
               resmx is the res in pixels of the biggest side
              (1+res(1))/2 (1+res(2)()/2 resmx/4*sqrt(0.5)   options.typicalgain 0.5
                %}
                % Check if the results make sense considering the inputs we used
                % Note: for him x and y are row and cols, so y is cols and x is rows
                % Flip it here
                Centerx0         = results.params(2);
                Centery0         = results.params(1);
                %             if (Centerx0 < 0 || Centerx0 > size(pm.RF.values,1) || ...
                %                     Centery0 < 0 || Centery0 > size(pm.RF.values,2))
                %                 error('The parameter estimate cannot be outside RF size limits')
                %             end
                
                % Change to the center of the fov
                if iseven(pm.Stimulus.ResizedHorz)
                    Centerx0 = Centerx0 - pm.Stimulus.ResizedHorz/2;
                else
                    Centerx0 = Centerx0 - ((pm.Stimulus.ResizedHorz-1)/2 + 1);
                end
                if iseven(pm.Stimulus.ResizedVert)
                    Centery0 = Centery0 - pm.Stimulus.ResizedVert/2;
                else
                    Centery0 = Centery0 - ((pm.Stimulus.ResizedVert-1)/2 + 1);
                end
                
                %% Convert the inputs to the same units we used
                Centerx0 = (pm.Stimulus.spatialSampleHorz * Centerx0);
                Centery0 = (pm.Stimulus.spatialSampleVert * Centery0);
                
                % Calculate RMSE
                RMSE     = sqrt(mse(results.testdata, results.modelpred));
                
                %{
                % This is the formula he uses to create the gaussian:
                % makegaussian2d(resmx,p1,p2,p3,p3) /(2*pi*p3^2)
                % p3 is the standard deviation in the vertical/horizontal direction
                % if you want an L1-normalized image, divide the image by 2*pi*p3^2
                % note that this is in reference to the ideal case where the Gaussian has
                % enough room to extend out.  so, if you are constructing a Gaussian that
                % does not fit very well within the image, the actual L1 length of the image
                % that is constructed will not be exactly 1.
                %
                % note that it doesn't matter if <sr> or <sc> are negative, since they
                % are always squared in function evaluation.

                % In his RF function he is normalizing it, then the sigma has not
                % the same meaning...
                % He does this:
                % r = (-1/res) * r + (.5 + .5/res);  % this is faster
                % c = (1/res) * c + (-.5 - .5/res);  % this is faster
                % sr = sr/res;
                % sc = sc/res;
                % In any case, he is providing the rfsize in pixel units, like this:
                % <rfsize> contains pRF size estimates.  pRF size is defined as sigma/sqrt(n) where
                %   sigma is the standard of the 2D Gaussian and n is the exponent of the power-law
                %   function.  Values are in pixel units with a lower bound of 0 pixels.
                % It is calculated in analyzePRF in this line:
                % results.rfsize(options.vxs,:) = permute(abs(paramsA(:,3,:)) ./ sqrt(posrect(paramsA(:,5,:))),[3 1 2]);
                % This means that we can go back to the sigma in pixels, and then
                % convert it to sigma in deg, as we inputed it.
                %
                % After conversation with Jon, adding the /sqrt(n)
                % As the result is always circular, we want to check that rfsize is
                % always the same as sigmaMinor and sigmaMajor.
                % It is, so we are removing the rfsize column from the table and
                % assigning it to sigmaMinor and sigmaMajor withouth the calculation
                % Leave the code here comented to undertand it better in the future
                % sigmaMinor = (pm.Stimulus.spatialSampleVert * abs(results.params(3))/sqrt(posrect(results.params(5))));
                % sigmaMajor = (pm.Stimulus.spatialSampleHorz * abs(results.params(3))/sqrt(posrect(results.params(5))));
                % tmpTable.rfsize = pm.Stimulus.spatialSampleHorz * results.rfsize;
                %}
                
                tmpTable            = struct2table(results,'AsArray',true);
                tmpTable.Centerx0   = Centerx0;
                tmpTable.Centery0   = Centery0;
                tmpTable.Theta      = 0; % Is circular, we can't model it
                % sigmaMajor and sigmaMinor will have same value, is circular RF
                % Convert it to degrees
                tmpTable.sigmaMajor = pm.Stimulus.spatialSampleHorz * results.rfsize;  %  * sqrt(posrect(results.params(5)));
                tmpTable.sigmaMinor = pm.Stimulus.spatialSampleHorz * results.rfsize;  %  * sqrt(posrect(results.params(5)));
                tmpTable.RMSE       = RMSE;
                % Make results table smaller (this was done for Brian The Substractor :) )
                
                % Leave the testdata and modelpred so that we have the fit.
                % Testdata is not exactly the same to the data (=pm.BOLDnoise), he
                % removes I think part of the low frew noise, check later. If we
                % cannot get testdata in other tools, we will just obtain modelpred.
                % We'll need to add pm.BOLDmeanValue to modelpred to have it in the
                % same range as the data (=pm.BOLDnoise)
                tmpTable            = tmpTable(:,{'Centerx0','Centery0', 'Theta' ,...
                                                  'sigmaMinor', 'sigmaMajor' ,...
                                                  'testdata','modelpred','R2','RMSE'});
       
                % tmpTable.testdata  = tmpTable.testdata  + pm.BOLDmeanValue;
                tmpTable.modelpred = tmpTable.modelpred + pm.BOLDmeanValue;

                pmEstimates = [pmEstimates; tmpTable];
            end
        end
        
    case {'vista','mrvista','vistasoft'}
        % Read the options
        options = allOptions.vista;
        %% Create temp folders
        tmpName = tempname(fullfile(pmRootPath,'local'));
        mkdir(tmpName);
        %% Change the inputs
        if niftiInputs
            niftiBOLDfile  = fullfile(tmpName, 'tmp.nii.gz')
            stimNiftiFname = fullfile(tmpName, 'tmpstim.nii.gz');
            
            copyfile(BOLDname,niftiBOLDfile)
            copyfile(STIMname,stimNiftiFname)
        else
            %% Write the stimuli as a nifti
            pm1            = input.pm(1);
            stimNiftiFname = fullfile(tmpName, 'tmpstim.nii.gz');
            stimNiftiFname = pm1.Stimulus.toNifti('fname',stimNiftiFname);

            %% Create a tmp nifti file
            niftiBOLDfile  = pmForwardModelToNifti(input, 'fname', ...
                fullfile(tmpName,'tmp.nii.gz'));
            
        end
        
        %% Prepare the function call
        homedir       = tmpName;
        stimfile      = stimNiftiFname;
        datafile      = niftiBOLDfile;
        warning('mrvista is assuming all stimuli with same radius. Fix this')
        model         = options.model;
        grid          = options.grid;
        obtainPreds   = options.obtainPreds;
        % TODO: fix this
        if isfield(options, 'wSearch'); wSearch = options.wSearch;end
        if isfield(options, 'wsearch'); wSearch = options.wsearch;end
        detrend       = options.detrend;
        keepAllPoints = options.keepAllPoints;
        numberStimulusGridPoints = options.numberStimulusGridPoints;
        
        % Make the call to the function based on Jon's script
        
        % TODO: save some steps if we already start with niftis
        mrvCleanWorkspace
        results = pmVistasoft(homedir, stimfile, datafile, stimradius,...
            'model'  , model, ...
            'grid'   , grid, ...
            'wSearch', wSearch, ...
            'detrend', detrend, ...
            'keepAllPoints', keepAllPoints, ...
            'numberStimulusGridPoints', numberStimulusGridPoints);
        
        %% Prepare the outputs in a table format
        pmEstimates = table();
        pmEstimates.Centerx0   = results.model{1}.x0';
        pmEstimates.Centery0   = results.model{1}.y0';
        pmEstimates.Theta      = results.model{1}.sigma.theta';
        pmEstimates.sigmaMajor = results.model{1}.sigma.major';
        pmEstimates.sigmaMinor = results.model{1}.sigma.minor';
        % Add the time series as well
        if niftiInputs
            data                   = niftiRead(BOLDname);
            pmEstimates.testdata   = squeeze(data.data);
        else
            pmEstimates.testdata   = repmat(ones([1,pm1.timePointsN]), ...
                [height(pmEstimates),1]);
            PMs                    = input.pm;
            for ii=1:height(pmEstimates); pmEstimates{ii,'testdata'}=PMs(ii).BOLDnoise;end
        end
        % Obtain the variance explained
        varexp         = 1-results.model{1}.rss./results.model{1}.rawrss;
        pmEstimates.R2 = varexp';
                
        % Obtain the modelfit if asked
        % Obtaining the prediction is too costly for big datasets
        if obtainPreds
            pmEstimates = pmVistaObtainPrediction(pmEstimates, results);
        else
            % Make the output same as the input
            pmEstimates.modelpred = pmEstimates.testdata;
        end
                
        % pmEstimates.R2         = calccod(pmEstimates.testdata,  pmEstimates.modelpred,2);
        % pmEstimates.R2        = results.model{1}.varExp';
        % pmEstimates.RMSE       = sqrt(mean((pmEstimates.testdata - pmEstimates.modelpred).^2,2));
        % errperf(T,P,'mae')
    
    case {'mrtools','mlrtools','mlr'}
        % Read the options
        options = allOptions.mlr;
        %% Create temp folders
        tmpName = tempname(fullfile(pmRootPath,'local'));
        mkdir(tmpName);
        %% Change the inputs
        if niftiInputs
            niftiBOLDfile  = fullfile(tmpName, 'tmp.nii.gz')
            stimNiftiFname = fullfile(tmpName, 'tmpstim.nii.gz');
            
            copyfile(BOLDname,niftiBOLDfile)
            copyfile(STIMname,stimNiftiFname)
        else
            %% Write the stimuli as a nifti
            pm1            = input.pm(1);
            stimNiftiFname = fullfile(tmpName, 'tmpstim.nii.gz');
            stimNiftiFname = pm1.Stimulus.toNifti('fname',stimNiftiFname);

            %% Create a tmp nifti file
            niftiBOLDfile  = pmForwardModelToNifti(input, 'fname', ...
                fullfile(tmpName,'tmp.nii.gz'));
            
        end
        
        %% Prepare the function call
        homedir       = tmpName;
        stimfile      = stimNiftiFname;
        datafile      = niftiBOLDfile;
        % Make the call to the function based on Jon's script
        
        % TODO: save some steps if we already start with niftis
        stimsize = [stimradius*2, stimradius*2];
        results  = mlrRunPRF(homedir, datafile, stimfile, stimsize,...
                             'quickFit', options.quickFit ,'doParallel',options.doParallel, ...
                             'rfType',options.rfType);
        if results.status==-1;error('mlr failed, check results.errorstring for more info');end
        
        %% Prepare the outputs in a table format
        pmEstimates = table();
        pmEstimates.Centerx0   = results.x;
        pmEstimates.Centery0   = results.y;
        pmEstimates.Theta      = zeros(size(results.y));
        pmEstimates.sigmaMajor = results.rfHalfWidth;
        pmEstimates.sigmaMinor = results.rfHalfWidth;
        
        % Add the time series as well
        if niftiInputs
            data                   = niftiRead(BOLDname);
            pmEstimates.testdata   = squeeze(data.data);
        else
            pmEstimates.testdata   = repmat(ones([1,pm1.timePointsN]), ...
                [height(pmEstimates),1]);
            PMs                    = input.pm;
            for ii=1:height(pmEstimates); pmEstimates{ii,'testdata'}=PMs(ii).BOLDnoise;end
        end
                
                
        % Obtain the modelfit,
        % pmEstimates = pmVistaObtainPrediction(pmEstimates, results);
        pmEstimates.modelpred = pmEstimates.testdata;
        pmEstimates.R2         = results.r2; % calccod(pmEstimates.testdata,  pmEstimates.modelpred,2);
        % pmEstimates.rss        = results.model{1}.rss';
        % pmEstimates.RMSE       = sqrt(mean((pmEstimates.testdata - pmEstimates.modelpred).^2,2));
        % errperf(T,P,'mae')
        
        % Delete the defaults file after using this solver, at least we have
        % observed that it messes up vistasoft's results
        if isfile('~/.mrDefaults.mat')
            delete '~/.mrDefaults.mat'
        else
            warning('mrtool (mlr) was run and it could find ~/.mrDefaults.mat to be deleted, check if it was writtme somewhere else')
        end
               
    case {'afni', 'simpleafni', 'basicafni', 'afni_4', 'afni4', ...
            'afni6', 'afni_6','withsigmaratio', ...
            'afni_dog', 'afnidog', 'dog'}
        % Read the options
        options = allOptions.afni;
        %% AFNI doc: algorithm options (TODO: add other implementations)
        % Afni help commented below
        %{
          Conv_PRF         : 4-param Population Receptive Field Model
                             (A, X, Y, sigma)
                             see model_conv_PRF.c
                  for help : setenv AFNI_MODEL_HELP_CONV_PRF YES
                             3dNLfim -signal bunnies
                  The model is made from parameters A, x0, y0, sigma, and from stimulus
                  time series input (visual field masks over time) by:

                  1. compute a Gaussian curve centered at x0, y0 of with spread sigma
                         g(x,y) = e^-( [(x-x0)^2+(y-y0)^2] / (2*sigma^2) )
                  2. multiply this 2-D image by each 2-D stimulus mask image
                  3. convolve the result with an ideal HRF
                  4. scale by the amplitude A

                  Currently, x0, y0, and sigma are limited to [-1,1], which the stimulus
                  images are evaluated on

          Conv_PRF_6       : 6-param Population Receptive Field Model
                             (A, X, Y, sigma, sigrat, theta)
                             see model_conv_PRF_6.c
                  for help : setenv AFNI_MODEL_HELP_CONV_PRF_6 YES
                             3dNLfim -signal bunnies
                  Given stimulus images over time s(x,y,t), find x0, y0, sigma, R and
                  theta values that produce a best fit of the model to the data.  Here
                  x0, y0 are taken to be the center of the population receptive field,
                  sigma is the minor width of it (sigma_x, below), sigrat R is the ratio
                  (sigma_y / sigma_x), and theta is the rotation from the y-direction
                  major axis (so zero is in the positive y-direction).

                  We assume sigma_y >= sigma_x and refer to sigrat >= 1, since that
                  sufficiently represents all possibilities.  The reciprocol would
                  come from the negative complimentary angle, and would therefore be a
                  redundant solution.

                  parameter domains:
                     x,y        : [-1,1], scaled by the mask, itself
                     sigma      : (0,1], where 1 means the mask radius
                     R (sigrat) : [1,inf), since sigma defines the smaller size
                     theta      : [-PI/2, PI/2), since rotation by PI has no effect

          Conv_PRF_DOG     : 6-param 'Difference of Gaussians' PRF Model
                             (as Conv_PRF, but with second A and sigma)
                             (A, X, Y, sig, A2, sig2)
                             see model_conv_PRF_DOG.c
                  for help : setenv AFNI_MODEL_HELP_CONV_PRF_DOG YES
                             3dNLfim -signal bunnies
        %}
                
        %% Change the inputs
        tmpName = tempname(fullfile(pmRootPath,'local'));
        mkdir(tmpName);
        if niftiInputs
            niftiBOLDfile  = fullfile(tmpName, 'tmp.nii.gz');         
            copyfile(BOLDname,niftiBOLDfile)
            
            % Demean it if pm.signalPercent == 'bold'
            if iscell(signalPercentage);signalPercentage=signalPercentage{:};end
            switch signalPercentage
                case {'bold'}
                    demean = true;
                otherwise
                    demean = false;
            end
            if demean
                demeanData = data;
                for ii=1:size(data,1)
                    tmp = data(ii,:);
                    demeanData(ii,:) = (tmp-mean(tmp)) / mean(tmp);
                end
                demeanData = reshape(demeanData,[size(demeanData,1),1,1,size(demeanData,2)]);
                TMP = niftiRead(niftiBOLDfile);
                TMP.data = demeanData;
                niftiWrite(TMP);
            end
            
        else
            warning('For AFNI analysis, be sure that all options have the same TR and the same stimulus')
            % Create a tmp nifti file and convert it to a tmp AFNI format
            pm1            = input.pm(1);
            if iscell(signalPercentage);signalPercentage=signalPercentage{:};end
            switch signalPercentage
                case {'bold'}
                    demean = true;
                otherwise
                    demean = false;
            end
            niftiBOLDfile  = pmForwardModelToNifti(input, ...
                'fname',fullfile(tmpName,'tmp.nii.gz'), ...
                'demean',demean);
        end
        
        % Create an AFNI file
        setenv('AFNI_NIFTI_TYPE_WARN','YES');
        if exist(fullfile(tmpName,'tmp.nii.gz'), 'file')
            system(['3dcopy ' niftiBOLDfile ' ' fullfile(tmpName,'tmp')]);
        end
        % We need to add this info fields
        % Now that it has been converted to 1D, it doesn't need or accept the
        % -space orig instruction anymore
        system(['3drefit -TR ' num2str(TR) ' ' fullfile(tmpName,'tmp.1D')])
        % system(['3drefit -space orig ' fullfile(tmpName, 'tmp+orig')])
        
        %% STIMULI
        % Prepare stimuli for AFNI
        mkdir(fullfile(tmpName,'images'));
        for ii=1:size(stimulus,3)
            thisImage = stimulus(:,:,ii);
            imwrite(thisImage, ...
                    fullfile(tmpName,'images', ...
                             sprintf('DSIMG_%03d.jpg', ii-1 )));
        end
        % Create a movie now
        % Careful, NIMAGES is one less than what we have
        setenv('NIMAGES', num2str(size(stimulus,3)-1)) % NUMBER OF IMAGES
        system(['3dcopy ' fullfile(tmpName,'images','DSIMG_000.jpg') ...
            ' ' fullfile(tmpName,'images','DSIMAGE_MOVIE')]);
        system(['for NUM in $(count -digits 3 1 $NIMAGES);do 3dTcat -glueto ' ...
            fullfile(tmpName, 'images','DSIMAGE_MOVIE+orig.') ' ' ...
            fullfile(tmpName, 'images','DSIMG_') '${NUM}.jpg;done']);
        % Add the stimuli we just created to env variables
        setenv('AFNI_MODEL_PRF_STIM_DSET', fullfile(tmpName, 'images','DSIMAGE_MOVIE+orig'));
        
        % We need to add this info fields
        system(['3drefit -TR ' num2str(TR) ' ' fullfile(tmpName,'images','DSIMAGE_MOVIE+orig')])
        system(['3drefit -space orig ' fullfile(tmpName, 'images','DSIMAGE_MOVIE+orig')])
        
        %% HRF
        % Obtain the HRF: follow the steps on the 3dNLfim.help file
        % TODO: add to options what HRF to use
        afni_hrf = options.hrf;
        switch afni_hrf
            % TODO: synch it with the HRF creation process in pmHRF.m 
            case {'GAM','gam'}
                % Default GAM normalized to 1
                system(['3dDeconvolve ' ...
                          '-nodata 50 ' num2str(TR) ' ' ...
                          '-polort -1 ' ... % Do not calculate detrending polinomials
                          '-num_stimts 1 ' ...
                          '-stim_times 1 "1D:0" GAM ' ... % k, tname, Rmodel
                          '-x1D ' fullfile(tmpName, 'conv.ref.GAM.1D')]);
                % Set the enviroment variable
                setenv('AFNI_CONVMODEL_REF', fullfile(tmpName, 'conv.ref.GAM.1D'));
            case {'SPM','spm'}
                system(['3dDeconvolve ' ...
                          '-nodata 50 ' num2str(TR) ' ' ...
                          '-polort -1 ' ...
                          '-num_stimts 1 ' ...
                          '-stim_times 1 "1D:0" SPMG1\(0\) ' ...
                          '-x1D ' fullfile(tmpName, 'sisar.conv.ref.SPMG1.1D')]);
                % Set the enviroment variable
                setenv('AFNI_CONVMODEL_REF', fullfile(tmpName, 'sisar.conv.ref.SPMG1.1D'));
            otherwise
                error('%s afni hrf not recognized',afni_hrf)
        end
        
        %% SET OTHER CONTROL ENVIROMENTAL VARIABLES
        % Not sure about the options here, Reynolds sent them to us
        setenv('AFNI_MODEL_DEBUG', '3');
        setenv('AFNI_MODEL_PRF_ON_GRID', 'NO');
        setenv('AFNI_MODEL_PRF_RAM_STATS', 'N');
        % If it takes too much, change this from 100 to 10 (if images too big, for
        % 100x100, 100 is ok)
        setenv('AFNI_MODEL_PRF_SIGMA_NSTEPS','100');
        
        %% NOISE MODELS
        %{
        % Noise Models (see the appropriate model_*.c file for exact details) :
        %
        %   Zero                     : Zero Noise Model
        %                              (no parameters)
        %                              see model_zero.c
        %
        %   Constant                 : Constant Noise Model
        %                              (constant)
        %                              see model_constant.c
        %
        %   Linear                   : Linear Noise Model
        %                              (constant, linear)
        %                              see model_linear.c
        %
        %   Linear+Ort               : Linear+Ort Noise Model
        %                              (constant, linear, Ort)
        %                              see model_linplusort.c
        %
        %   Quadratic                : Quadratic Noise Model
        %                              (constant, linear, quadratic)
        %                              see model_quadratic.c
        %}
        
        %% RUN THE ANALYSIS
        % Create the command line for the fit
        % Apply it to every voxel
        % c = parcluster('local');
        c.NumWorkers = 0;
        % Adding the absolute path is throwing an error. Not Matlab related, same in CLI
        % Use cd() and then launch the command
        cd(fullfile(tmpName))
        
        switch options.model
            case {'afni', 'simpleafni', 'basicafni', 'afni_4', 'afni4','Conv_PRF'}
                % 4-param pRF Field Model 
                % (A, X, Y, sigma)
                modelName   = 'Conv_PRF';
                modelConstr =   [...
                    '-sconstr 0 -10.0 10.0 ' ...
                    '-sconstr 1  -1.0  1.0 ' ...
                    '-sconstr 2  -1.0  1.0 ' ...
                    '-sconstr 3   0.0  1.0 ' ...
                    ];
            case {'afni6', 'afni_6','withsigmaratio','Conv_PRF_6'}
                % 6-param pRF Field Model 
                % (A, X, Y, sigma, sigrat, theta)
                % Theta in radiansm, constraints from -90 to +90
                modelName   = 'Conv_PRF_6';
                modelConstr =   [ ...
                    '-sconstr 0 -10.0 10.0 ' ...
                    '-sconstr 1  -1.0 1.0 ' ...
                    '-sconstr 2  -1.0 1.0 ' ...
                    '-sconstr 3   0.0 1.0 ' ...
                    '-sconstr 4   1.0 5.0 ' ...
                    '-sconstr 5  -1.571 1.570 ' ...
                    ];
            case {'afni_dog', 'afnidog', 'dog','Conv_PRF_DOG'}
                % 6-param 'Difference of Gaussians' PRF Model
                % As Conv_PRF, but with second A and sigma
                % (A, X, Y, sig, A2, sig2)
                modelName   = 'Conv_PRF_DOG';
                modelConstr =   [...
                    '-sconstr 0 -10.0 10.0 ' ...
                    '-sconstr 1  -1.0  1.0 ' ...
                    '-sconstr 2  -1.0  1.0 ' ...
                    '-sconstr 3   0.0  1.0 ' ...
                    '-sconstr 4 -10.0 10.0 ' ...
                    '-sconstr 5   0.0  1.0 ' ...
                    ];
            otherwise
                error('%s not implemented yet',options.model)
        end
        
        % Launch the command
        tic
        system([...
            '3dNLfim -input tmp.1D ' ...
            '-noise Zero ' ... 
            '-signal ' modelName ' ' ...
            modelConstr ....
            '-BOTH ' ...
            '-nrand 10000 ' ...
            '-nbest 5 ' ...
            '-bucket 0 buck_tmp.PRF ' ...
            '-snfit snfit_tmp.PRF ' ...
            '-TR ' num2str(TR) ' ' ...
            '-jobs 36' ... % num2str(c.NumWorkers)...
            ]);
        toc
        
        %% Read the results back to Matlab
        kk = which('BrikLoad');
        if isempty(kk)
            addpath(genpath('~/soft/afni_matlab'));
        end
        [err, V, Info, ErrMessage] = BrikLoad('snfit_tmp.PRF.1D');
        % Plot the fit and the signal of voxel 2,2,2
        % plot(squeeze(myCube(5,5,2,:)),'k-');hold on;plot(squeeze(V(5,5,2,:)),'r-');
        [berr, bV, bInfo, bErrMessage] = BrikLoad('buck_tmp.PRF.1D');
        %% Delete the tmp folder
        cd(fullfile(pmRootPath,'local'));
        rmdir(fullfile(tmpName),'s');
        %% Prepare data, return as results
        % fitSeries = reshape(V,  [height(input), timePointsN]);
        % results   = reshape(bV, [height(input), size(bV,4)]); % Two more params
        results   = bV;
        fits      = V;
        
        % Calculate correct estimates, according to Reynolds:
        %   You can provide constraints on the command line.
        %   But note that [-1,1] is actually the definition of the
        %   mask area.  It is evaluated as a unit field of view,
        %   since the input width is really arbitrary anyway.
        %   It sounds like the results are consistent between
        %   the packages.  If your mask is of width 60, from
        %   -30 to 30, then the X and Y results should be scaled
        %   by that factor of 30 when comparing.
        %   The '30' does not really add information here. A unit FOV made sense to me.

        % Write results common to all of them
        % Params
        pmEstimates = table();
        pmEstimates.A          = results(:,1);
        pmEstimates.Centerx0   = results(:,2) * Stimulus.fieldofviewHorz/2;
        pmEstimates.Centery0   = results(:,3) * Stimulus.fieldofviewVert/2;
        
          
        
        % Fits
        %   - Testdata: same as input in the nifti
        pmEstimates.testdata = data;
        
        %   - modelpred: the fits comming out from AFNI
        % restore back the demeaning
        fitsScaledBack = (BOLDmeanValue * fits) + BOLDmeanValue;
        pmEstimates.modelpred  = fitsScaledBack;
        
        
        switch options.model
            case {'afni', 'simpleafni', 'basicafni', 'afni_4', 'afni4', 'Conv_PRF'}
                % 4-param pRF Field Model
                % (A, X, Y, sigma)
                pmEstimates.Theta      = zeros(size(results(:,3)));
                pmEstimates.sigmaMajor = results(:,4) * Stimulus.fieldofviewHorz/2;
                pmEstimates.sigmaMinor = results(:,4) * Stimulus.fieldofviewVert/2;
            case {'afni6', 'afni_6','withsigmaratio', 'Conv_PRF_6'}
                % 6-param pRF Field Model
                % (A, X, Y, sigma, sigrat, theta)
                % Theta in radiansm, constraints from -90 to +90
                pmEstimates.sigmaMajor = (results(:,4) .* results(:,5)) * Stimulus.fieldofviewHorz/2;
                pmEstimates.sigmaMinor = results(:,4) * Stimulus.fieldofviewHorz/2;
                pmEstimates.Theta      = results(:,6) + deg2rad(90);
            case {'afni_dog', 'afnidog', 'dog', 'Conv_PRF_DOG'}
                % 6-param 'Difference of Gaussians' PRF Model
                % As Conv_PRF, but with second A and sigma
                % (A, X, Y, sig, A2, sig2)
                pmEstimates.A2         = results(:,5);
                pmEstimates.sigmaMajor = results(:,4) * Stimulus.fieldofviewHorz/2;
                pmEstimates.sigmaMinor = results(:,6) * Stimulus.fieldofviewVert/2;
            otherwise
                error('%s not implemented yet',options.model)
        end
        
    case {'popeye','popeye_onegaussian','popeye_CSS','popnohrf','popeyenohrf'}
        %% Create temp folders
        tmpName = tempname(fullfile(pmRootPath,'local'));
        mkdir(tmpName);
        %% Change the inputs
        if niftiInputs
            niftiBOLDfile  = fullfile(tmpName, 'tmp.nii.gz')
            stimNiftiFname = fullfile(tmpName, 'tmpstim.nii.gz');
            copyfile(BOLDname,niftiBOLDfile)
            copyfile(STIMname,stimNiftiFname)   
        else
            warning('For popeye analysis, be sure that all options have the same TR and the same stimulus')
            % Create a tmp nifti file
            pm1            = input.pm(1);
            niftiBOLDfile  = pmForwardModelToNifti(input, ...
                'fname',fullfile(tmpName,'tmp.nii.gz'), ...
                'demean',false);
            %% Create stimulus file in nifti
            pm1            = input.pm(1);
            stimNiftiFname = fullfile(tmpName, 'tmpstim.nii.gz');
            stimNiftiFname = pm1.Stimulus.toNifti('fname',stimNiftiFname);
            
        end
        
        %Params
        params.pixels_per_degree = Stimulus.ResizedHorz / Stimulus.fieldofviewHorz;
        
        %% Prepare the json and function call
        % Create struct  with variables for params.json
        params.stimulus_file     = 'tmpstim.nii.gz';
        params.data_file         = 'tmp.nii.gz';
        params.screen_distance   = 50;
        params.screen_width      = 28.67;
        params.invert_y          = true;
        % Encode into json file 
        jsonStr = jsonencode(params);
        fid = fopen(fullfile(tmpName,'params.json'), 'w');
        if fid == -1, error('Cannot create JSON file'); end
        fwrite(fid, jsonStr, 'char');
        fclose(fid);
        % Run the Docker container
        switch prfimplementation
            case {'popnohrf','popeyenohrf'}
                system(['docker run -it --rm '...
                        '-v "' tmpName ':/input" garikoitz/popeye-mini-nohrf' ...
                        ]);
            otherwise
                system(['docker run -it --rm '...
                        '-v "' tmpName ':/input" garikoitz/popeye-mini' ...
                        ]);
        end
        %% Read the results back to Matlab
        x0    = niftiRead(fullfile(tmpName,'out_x.nii.gz'));
        x0    = x0.data(:,1,1,1);
        y0    = niftiRead(fullfile(tmpName,'out_y.nii.gz'));
        y0    = y0.data(:,1,1,1);
        gain  = niftiRead(fullfile(tmpName,'out_gain.nii.gz'));
        gain  = gain.data(:,1,1,1);
        sigma = niftiRead(fullfile(tmpName,'out_sigma.nii.gz'));
        sigma = sigma.data(:,1,1,1);
        modelpred = niftiRead(fullfile(tmpName,'out_modelpred.nii.gz'));
        modelpred = modelpred.data;
        
        % Write results common to all of them
        % Params
        pmEstimates = table();
        pmEstimates.A          = gain;
        pmEstimates.Centerx0   = x0;
        pmEstimates.Centery0   = y0;
        pmEstimates.sigmaMajor = sigma;
        pmEstimates.sigmaMinor = sigma;
        pmEstimates.Theta      = zeros(size(sigma));
        % Fits
        %   - Testdata: same as input in the nifti
        pmEstimates.testdata   = data;
        pmEstimates.modelpred  = squeeze(modelpred);
        % restore back the demeaning
        
        % for ii=1:size(modelpred,1)
        %     pmEstimates.modelpred(ii,:) = BOLDmeanValue * (pmEstimates.modelpred(ii,:) + 1);
        % end
              
    otherwise
        error('Method %s not implemented yet.', prfimplementation)
end

% Make the order of pmEstimates is always the same, this is required for the unpacking python function in prfanalyze
% This was fixed in python already
% pmEstimates = pmEstimates(:,{'Centerx0','Centery0','Theta','sigmaMinor','sigmaMajor','R2','testdata','modelpred'});


end

