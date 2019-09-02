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
%   inputTable  - Data table including stimulus parameters and BOLD
%                 time series. Each row of the data table is another
%                 voxel
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
p.addParameter('useParallel'    ,  false        , @islogical);

% This options structs are defaults for analyzePRF
options  = struct('seedmode', [0 1], 'display' , 'off');
% Implementation specifics
% AnalyzePRF
p.addParameter('options'    ,  options        , @isstruct);
% Vistasoft
p.addParameter('model'        , 'one gaussian'  , @ischar);
p.addParameter('grid'         , false           , @islogical);
p.addParameter('wsearch'      , 'coarse to fine', @ischar);
p.addParameter('detrend'      , true            , @islogical);
p.addParameter('keepAllPoints', false           , @islogical);
p.addParameter('numberStimulusGridPoints', 50   , @isnumeric);

% AFNI
p.addParameter('afni_hrf'   , 'SPM'           , @ischar);
% Popeye
% p.addParameter('afni_hrf'   , 'SPM'           , @ischar);



% Parse. Assign result inside each case
p.parse(input,prfimplementation,varargin{:});
% Read here only the generic ones
useParallel = p.Results.useParallel;



%% Choose the analysis case
prfimplementation = mrvParamFormat(prfimplementation);

switch prfimplementation
    case {'aprf','analyzeprf'}
        %% Create and check tables
        % TODO: Let's define the format for the estimates so that this is
        %       the same for all the methods.
        pmEstimates = table();
        
        
        % Check if the pm-s come in a table or alone
        if ~istable(input)
            temp = table();
            temp.pm = input;
            input = temp;
        end
        
        if useParallel
            % Check if the TR is the same
            if length(unique(input.TR)) > 1
                error('TR has to be the same for all synthetic BOLD series')
            end
            % Add the time series as well
            pm1        = input.pm(1);
            BOLDdata   = repmat(ones([1,pm1.timePointsN]), [height(input),1]);
            PMs        = input.pm;
            for ii=1:height(input); pmEstimates{ii,'testdata'}=PMs(ii).BOLDnoise;end
            
            % Create variables and function call
            stimulus = double(pm1.Stimulus.getStimValues);
            data     = pmEstimates.testdata;
            TR       = pm1.TR;
            options  = p.Results.options;
            
            % Calculate PRF
            results  = analyzePRF({stimulus}, {data}, TR, options);
            
            % Prepare output
            Centerx0         = squeeze(results.params(1,2,:));
            Centery0         = squeeze(results.params(1,1,:));
            % Change to the center of the fov
            if all(iseven(pm1.Stimulus.ResizedHorz))
                Centerx0 = Centerx0 - pm1.Stimulus.ResizedHorz/2;
            else
                Centerx0 = Centerx0 - ((pm1.Stimulus.ResizedHorz-1)/2 + 1);
            end
            if iseven(pm1.Stimulus.ResizedVert)
                Centery0 = Centery0 - pm1.Stimulus.ResizedVert/2;
            else
                Centery0 = Centery0 - ((pm1.Stimulus.ResizedVert-1)/2 + 1);
            end
            pmEstimates.Centerx0   = (pm1.Stimulus.spatialSampleHorz * Centerx0);
            pmEstimates.Centery0   = (pm1.Stimulus.spatialSampleVert * Centery0);
            pmEstimates.Theta      = zeros(size(Centery0));
            pmEstimates.sigmaMinor = pm1.Stimulus.spatialSampleHorz * results.rfsize;
            pmEstimates.sigmaMajor = pmEstimates.sigmaMinor;
            pmEstimates.modelpred  = results.modelpred;
            % Scale it back
            pmEstimates.modelpred  = pmEstimates.modelpred + pm1.BOLDmeanValue;
            
        else
            %% Calculate values voxel to voxel
            for ii=1:height(input)
                % TODO: use parfor if the number of rows if the table is larger than XX
                pm       = input.pm(ii);
                stimulus = double(pm.Stimulus.getStimValues);
                data     = pm.BOLDnoise;
                TR       = pm.TR;
                options  = p.Results.options;
                
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
        %% Create temp folders
        tmpName = tempname(fullfile(pmRootPath,'local'));
        mkdir(tmpName);
        
        %% Write the stimuli as a nifti
        pm1            = input.pm(1);
        stimNiftiFname = fullfile(tmpName, 'tmpstim.nii.gz');
        stimNiftiFname = pm1.Stimulus.toNifti('fname',stimNiftiFname);
        
        %% Create a tmp nifti file
        niftiBOLDfile  = pmForwardModelToNifti(input, 'fname', ...
            fullfile(tmpName,'tmp.nii.gz'));
        
        %% Prepare the function call
        homedir       = tmpName;
        stimfile      = stimNiftiFname;
        datafile      = niftiBOLDfile;
        warning('mrvista is assuming all stimuli with same radius. Fix this')
        stimradius    = pm1.Stimulus.fieldofviewHorz/2;
        model         = p.Results.model;
        grid          = p.Results.grid;
        wSearch       = p.Results.wsearch;
        detrend       = p.Results.detrend;
        keepAllPoints = p.Results.keepAllPoints;
        numberStimulusGridPoints = p.Results.numberStimulusGridPoints;
        
        % Make the call to the function based on Jon's script
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
        pmEstimates.testdata   = repmat(ones([1,pm1.timePointsN]), ...
                                                       [height(pmEstimates),1]);
        PMs                    = input.pm;
        for ii=1:height(pmEstimates); pmEstimates{ii,'testdata'}=PMs(ii).BOLDnoise;end
        
        % Obtain the modelfit,
        pmEstimates = pmVistaObtainPrediction(pmEstimates, input, results);
        
        pmEstimates.R2         = calccod(pmEstimates.testdata,  pmEstimates.modelpred,2);
        pmEstimates.rss        = results.model{1}.rss';
        pmEstimates.RMSE       = sqrt(mean((pmEstimates.testdata - pmEstimates.modelpred).^2,2));
        % errperf(T,P,'mae')
        
    case {'afni', 'simpleafni', 'basicafni', 'afni_4', 'afni4', ...
          'afni6', 'afni_6','withsigmaratio', ...
          'afni_dog', 'afnidog', 'dog'}
        %% AFNI doc: algorithm options (TODO: add other implementations)
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
        
        warning('For AFNI analysis, be sure that all options have the same TR and the same stimulus')
        tmpName = tempname(fullfile(pmRootPath,'local'));
        mkdir(tmpName);
        % Create a tmp nifti file and convert it to a tmp AFNI format
        pm1            = input.pm(1);
        niftiBOLDfile  = pmForwardModelToNifti(input, ...
                                               'fname',fullfile(tmpName,'tmp.nii.gz'), ...
                                               'demean',true);
        % Create an AFNI file
        setenv('AFNI_NIFTI_TYPE_WARN','YES');
        if exist(fullfile(tmpName,'tmp.nii.gz'), 'file')
            system(['3dcopy ' fullfile(tmpName,'tmp.nii.gz') ' ' fullfile(tmpName,'tmp')]);
        end
        % We need to add this info fields
        % Now that it has been converted to 1D, it doesn't need or accept the
        % -space orig instruction anymore
        system(['3drefit -TR ' num2str(pm1.TR) ' ' fullfile(tmpName,'tmp.1D')])
        % system(['3drefit -space orig ' fullfile(tmpName, 'tmp+orig')])
        
        %% STIMULI
        % Prepare stimuli for AFNI
        pm1      = input.pm(1);
        stimulus = pm1.Stimulus.getStimValues;
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
        system(['3drefit -TR ' num2str(pm1.TR) ' ' fullfile(tmpName,'images','DSIMAGE_MOVIE+orig')])
        system(['3drefit -space orig ' fullfile(tmpName, 'images','DSIMAGE_MOVIE+orig')])
        
        %% HRF
        % Obtain the HRF: follow the steps on the 3dNLfim.help file
        % TODO: add to options what HRF to use
        afni_hrf = p.Results.afni_hrf;
        switch afni_hrf
            % TODO: synch it with the HRF creation process in pmHRF.m 
            case {'GAM'}
                % Default GAM normalized to 1
                system(['3dDeconvolve ' ...
                          '-nodata 50 ' num2str(pm1.TR) ' ' ...
                          '-polort -1 ' ... % Do not calculate detrending polinomials
                          '-num_stimts 1 ' ...
                          '-stim_times 1 "1D:0" GAM ' ... % k, tname, Rmodel
                          '-x1D ' fullfile(tmpName, 'conv.ref.GAM.1D')]);
                % Set the enviroment variable
                setenv('AFNI_CONVMODEL_REF', fullfile(tmpName, 'conv.ref.GAM.1D'));
            case {'SPM'}
                system(['3dDeconvolve ' ...
                          '-nodata 50 ' num2str(pm1.TR) ' ' ...
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
        % Not sure about the options here
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
        
        switch prfimplementation
            case {'afni', 'simpleafni', 'basicafni', 'afni_4', 'afni4'}
                % 4-param pRF Field Model 
                % (A, X, Y, sigma)
                modelName   = 'Conv_PRF';
                modelConstr =   [...
                    '-sconstr 0 -10.0 10.0 ' ...
                    '-sconstr 1  -1.0  1.0 ' ...
                    '-sconstr 2  -1.0  1.0 ' ...
                    '-sconstr 3   0.0  1.0 ' ...
                    ];
            case {'afni6', 'afni_6','withsigmaratio'}
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
            case {'afni_dog', 'afnidog', 'dog'}
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
                error('%s not implemented yet',prfimplementation)
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
            '-TR ' num2str(pm1.TR) ' ' ...
            '-jobs ' num2str(c.NumWorkers)...
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
        fitSeries = reshape(V,  [height(input), pm1.timePointsN]);
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
        pmEstimates.Centerx0   = results(:,2) * pm1.Stimulus.fieldofviewHorz/2;
        pmEstimates.Centery0   = results(:,3) * pm1.Stimulus.fieldofviewVert/2;
        % Fits
        %   - Testdata: same as input in the nifti
        pmEstimates.testdata   = repmat(ones([1,pm1.timePointsN]), ...
            [height(pmEstimates),1]);
        PMs                    = input.pm;
        for ii=1:height(pmEstimates); pmEstimates{ii,'testdata'}=PMs(ii).BOLDnoise;end
        
        %   - modelpred: the fits comming out from AFNI
        % restore back the demeaning
        fitsScaledBack = (pm1.BOLDmeanValue * fits) + pm1.BOLDmeanValue;
        pmEstimates.modelpred  = fitsScaledBack;
        
        
        switch prfimplementation
            case {'afni', 'simpleafni', 'basicafni', 'afni_4', 'afni4'}
                % 4-param pRF Field Model
                % (A, X, Y, sigma)
                pmEstimates.Theta      = zeros(size(results(:,3)));
                pmEstimates.sigmaMajor = results(:,4) * pm1.Stimulus.fieldofviewHorz/2;
                pmEstimates.sigmaMinor = results(:,4) * pm1.Stimulus.fieldofviewVert/2;
            case {'afni6', 'afni_6','withsigmaratio'}
                % 6-param pRF Field Model
                % (A, X, Y, sigma, sigrat, theta)
                % Theta in radiansm, constraints from -90 to +90
                pmEstimates.sigmaMajor = (results(:,4) .* results(:,5)) * pm1.Stimulus.fieldofviewHorz/2;
                pmEstimates.sigmaMinor = results(:,4) * pm1.Stimulus.fieldofviewHorz/2;
                pmEstimates.Theta      = results(:,6);
            case {'afni_dog', 'afnidog', 'dog'}
                % 6-param 'Difference of Gaussians' PRF Model
                % As Conv_PRF, but with second A and sigma
                % (A, X, Y, sig, A2, sig2)
                pmEstimates.A2         = results(:,5);
                pmEstimates.sigmaMajor = results(:,4) * pm1.Stimulus.fieldofviewHorz/2;
                pmEstimates.sigmaMinor = results(:,6) * pm1.Stimulus.fieldofviewVert/2;
            otherwise
                error('%s not implemented yet',prfimplementation)
        end
        
    case {'popeye','popeye_onegaussian','popeye_CSS','popeye_dog'}
        %% Create time series file
        warning('For popeye analysis, be sure that all options have the same TR and the same stimulus')
        tmpName = tempname(fullfile(pmRootPath,'local'));
        mkdir(tmpName);
        % Create a tmp nifti file and convert it to a tmp AFNI format
        pm1            = input.pm(1);
        niftiBOLDfile  = pmForwardModelToNifti(input, ...
            'fname',fullfile(tmpName,'tmp.nii.gz'), ...
            'demean',true);
        %% Create stimulus file in nifti
        pm1            = input.pm(1);
        stimNiftiFname = fullfile(tmpName, 'tmpstim.nii.gz');
        stimNiftiFname = pm1.Stimulus.toNifti('fname',stimNiftiFname);
        
        %% Prepare the json and function call
        % Create struct  with variables for params.json
        params.stimulus_file     = 'tmpstim.nii.gz';
        params.data_file         = 'tmp.nii.gz';
        params.screen_distance   = 50;
        params.screen_width      = 28.67;
        params.pixels_per_degree = pm1.Stimulus.ResizedHorz / pm1.Stimulus.fieldofviewHorz;
        params.invert_y          = true;
        % Encode into json file 
        jsonStr = jsonencode(params);
        fid = fopen(fullfile(tmpName,'params.json'), 'w');
        if fid == -1, error('Cannot create JSON file'); end
        fwrite(fid, jsonStr, 'char');
        fclose(fid);
        % Run the docerk container
        system(['docker run -it --rm '...
                '-v "' tmpName ':/input" nben/popeye-mini' ...
                ]);
        
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
        pmEstimates.testdata   = repmat(ones([1,pm1.timePointsN]), ...
            [height(pmEstimates),1]);
        PMs                    = input.pm;
        for ii=1:height(pmEstimates); pmEstimates{ii,'testdata'}=PMs(ii).BOLDnoise;end
        
        pmEstimates.modelpred  = squeeze(modelpred);
        % restore back the demeaning
        for ii=1:size(modelpred,1)
            pmEstimates.modelpred(ii,:) = pm1.BOLDmeanValue * (pmEstimates.modelpred(ii,:) + 1);
        end
              
    otherwise
        error('Method %s not implemented yet.', prfimplementation)
end


end




function pmEstimates = pmVistaObtainPrediction(pmEstimates, input, results)
    % Initialize the results
    pmEstimates.modelpred = pmEstimates.testdata;



    % start from the function rmPlotGUI_makePrediction
    recompFit = true;  % We want to recompute fit
        
    % paramsOrder is as follow : sigmaMajor,sigmaMinor,theta, x0,y0
    % rfGaussian2d(X, Y,         rfParams(3), rfParams(5), rfParams(6), rfParams(1), rfParams(2));
    % rfParams(4) is not used, at least in this function
    % I need to create rfParams = []; per every voxel it seems
    
    for voxel=1:height(input) 
        rfParams = zeros([1,6]);
        rfParams(1) = results.model{1}.x0(voxel);
        rfParams(2) = results.model{1}.y0(voxel);
        rfParams(3) = results.model{1}.sigma.major(voxel);
        % rfParams(4) =  % Esto se asigna luego, suele ser el beta(1)
        rfParams(5) = results.model{1}.sigma.minor(voxel);
        rfParams(6) = results.model{1}.sigma.theta(voxel);
        
        % %% make RFs
        % RFs = rmPlotGUI_makeRFs(modelName, rfParams, M.params.analysis.X, M.params.analysis.Y);
        RFs = rmPlotGUI_makeRFs(results.model{1}.description, ...
            rfParams, ...
            results.params.analysis.X,results.params.analysis.Y);
        
        %% make predictions for each RF
        % pred = M.params.analysis.allstimimages * RFs;
        pred = results.params.analysis.allstimimages * RFs;
        
        
        % Determine which frames have no stimulus. We may want to use this
        % information to highlight the blanks in the time series plots. We need to
        % determine blanks from the original images, not the images that have been
        % convolved with the hRF.
        %{
        stim = [];
        for ii = 1:length(M.params.stim)
            endframe = size(M.params.stim(ii).images_org, 2);
            frames =  endframe - M.params.stim(ii).nFrames+1:endframe;
            stim = [stim M.params.stim(ii).images_org(:, frames)];
        end
        blanks = sum(stim, 1) < .001;
        %}
        stim = [];
        for ii = 1:length(results.params.stim)
            endframe = size(results.params.stim(ii).images_org, 2);
            frames =  endframe - results.params.stim(ii).nFrames+1:endframe;
            stim = [stim results.params.stim(ii).images_org(:, frames)];
        end
        blanks = sum(stim, 1) < .001;
        
        
        
        %% get/make trends
        [trends, ntrends, dcid] = rmMakeTrends(results.params, 0);
        if isfield(results.params.analysis,'allnuisance')
            trends = [trends results.params.analysis.allnuisance];
        end
        
        %% Compute final predicted time series (and get beta values)
        % we also add this to the rfParams, to report later
        M = results;
        
        % In resutls we don't have M.tSeries(:,voxel)
        % I am going to paste the testdata in vertical form. 
        M.tSeries = pmEstimates.testdata';
        
        
        switch M.model{1}.description,
            case {'2D pRF fit (x,y,sigma, positive only)',...
                    '2D RF (x,y,sigma) fit (positive only)',...
                    '1D pRF fit (x,sigma, positive only)'};
                if recompFit==0,
                    beta = rmCoordsGet(M.viewType, model, 'b', coords);
                    beta = beta([1 dcid+1])';
                    
                else
                    beta = pinv([pred trends(:,dcid)])*M.tSeries(:,voxel);
                    beta(1) = max(beta(1),0);
                    
                end
                
                RFs        = RFs .* (beta(1) .* M.params.analysis.HrfMaxResponse);
                
                rfParams(4) = beta(1);
                
                
            case {'2D pRF fit (x,y,sigma_major,sigma_minor)' ...
                    'oval 2D pRF fit (x,y,sigma_major,sigma_minor,theta)'};
                if recompFit==0,
                    beta = rmCoordsGet(M.viewType, model, 'b', coords);
                    beta = beta([1 dcid+1]);
                else
                    beta = pinv([pred trends(:,dcid)])*M.tSeries(:,voxel);
                    beta(1) = max(beta(1),0);
                    
                end
                
                RFs        = RFs .* (beta(1) .* M.params.analysis.HrfMaxResponse);
                
                rfParams(4) = beta(1);
                
            case 'unsigned 2D pRF fit (x,y,sigma)';
                if recompFit==0,
                    beta = rmCoordsGet(M.viewType, model, 'b', coords);
                    beta = beta([1 dcid+1]);
                else
                    beta = pinv([pred trends(:,dcid)])*M.tSeries(:,voxel);
                end
                
                RFs        = RFs .* (beta(1) .* M.params.analysis.HrfMaxResponse);
                
                rfParams(4) = beta(1);
                
            case {'Double 2D pRF fit (x,y,sigma,sigma2, center=positive)',...
                    'Difference 2D pRF fit (x,y,sigma,sigma2, center=positive)',...
                    'Difference 1D pRF fit (x,sigma, sigma2, center=positive)'},
                if recompFit==0,
                    beta = rmCoordsGet(M.viewType, model, 'b', coords);
                    beta = beta([1 2 dcid+2]);
                    beta = beta';
                else
                    beta = pinv([pred trends(:,dcid)])*M.tSeries(:,voxel);
                    beta(1) = max(beta(1),0);
                    beta(2) = max(beta(2),-abs(beta(1)));
                end
                
                RFs        = RFs * (beta(1:2).*M.params.analysis.HrfMaxResponse);
                
                rfParams(:,4) = beta(1);
                
            case {'Two independent 2D pRF fit (2*(x,y,sigma, positive only))'},
                if recompFit==0,
                    beta = rmCoordsGet(M.viewType, model, 'b', coords);
                    beta = beta([1 2 dcid+2]);
                else
                    beta = pinv([pred trends(:,dcid)])*M.tSeries(:,voxel);
                    beta(1:2) = max(beta(1:2),0);
                end
                
                RFs        = RFs * (beta(1:2) .* M.params.analysis.HrfMaxResponse);
                
                rfParams(:,4) = beta(1:2);
                rfParams = rfParams(1,:);
                
            case {'Mirrored 2D pRF fit (2*(x,y,sigma, positive only))'},
                if recompFit==0,
                    beta = rmCoordsGet(M.viewType, model, 'b', coords);
                    beta = beta([1 dcid+1]);
                else
                    beta = pinv([pred trends(:,dcid)])*M.tSeries(:,voxel);
                    beta(1) = max(beta(1),0);
                end
                
                RFs        = RFs * (beta(1) .* M.params.analysis.HrfMaxResponse);
                
                rfParams(:,4) = beta(1);
                rfParams = rfParams(1,:);
                
            case {'Sequential 2D pRF fit (2*(x,y,sigma, positive only))'},
                if recompFit==0,
                    beta = rmCoordsGet(M.viewType, model, 'b', coords);
                    beta = beta([1 2 dcid+2]);
                else
                    beta = pinv([pred trends(:,dcid)])*M.tSeries(:,voxel);
                    beta(1:2) = max(beta(1:2),0);
                end
                
                RFs        = RFs * (beta(1:2) .* M.params.analysis.HrfMaxResponse);
                
                rfParams(:,4) = beta(1:2);
                rfParams = rfParams(1,:);
            case {'css' '2D nonlinear pRF fit (x,y,sigma,exponent, positive only)'}
                % we-do the prediction with stimulus that has not been convolved
                % with the hrf, and then add in the exponent, and then convolve
                
                % make neural predictions for each RF
                pred = (M.params.analysis.allstimimages_unconvolved * RFs).^rfParams(7);
                % reconvolve with hRF
                for scan = 1:length(M.params.stim)
                    these_time_points = M.params.analysis.scan_number == scan;
                    hrf = M.params.analysis.Hrf{scan};
                    pred(these_time_points,:) = filter(hrf, 1, pred(these_time_points,:));
                end
                
                if recompFit
                    beta = pinv([pred trends(:,dcid)])*M.tSeries(:,voxel);
                    beta(1) = max(beta(1),0);
                    
                else
                    beta = rmCoordsGet(M.viewType, model, 'b', coords);
                    beta = beta([1 dcid+1])';   % scale factor (gain) and DC (mean)
                end
                
                RFs        = RFs .* (beta(1) .* M.params.analysis.HrfMaxResponse);
                
                rfParams(4) = beta(1);
                
                
                
            otherwise,
                error('Unknown modelName: %s', modelName);
        end;
        
        % Calculate the prediction
        prediction = [pred trends(:,dcid)] * beta;
        
        % Add it to the return table
        pmEstimates.modelpred(voxel,:) = prediction';
    end


    
    
    
    
    
end





