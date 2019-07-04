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

%  TODO: Make all the outputs the same so that we can use the same function to
%        compare them to the synthetic data. 
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
% This options structs are defaults for analyzePRF
options  = struct('seedmode', [0 1], 'display' , 'off');
p.addParameter('options'    ,  options , @isstruct);
p.parse(input,prfimplementation,varargin{:});




%% Choose the analysis case
prfimplementation = mrvParamFormat(prfimplementation);


switch prfimplementation
    case {'analyzeprf'}
        % TODO: Let's define the format for the estimates so that this is
        %       the same for all the methods.
        pmEstimates = table();
        
        % Make sure we are using winawerlab analyzePRF (I am using branch gari)
        if ~contains(which('analyzePRF'), 'winawerlab')
            addpath(genpath('~/soft/winawerlab/analyzePRF'));
        end
        
        % Check if the pm-s come in a table or alone
        if ~istable(input)
            temp = table();
            temp.pm = input;
            input = temp;
        end
        
        % Go line by line and compute required values for each pRF model
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
            if (Centerx0 < 0 || Centerx0 > size(pm.RF.values,1) || ...
                Centery0 < 0 || Centery0 > size(pm.RF.values,2))
                error('The parameter estimate cannot be outside RF size limits')
            end
            
            % Convert the inputs to the same units we used
            Centerx0 = (pm.Stimulus.spatialSampleHorz * Centerx0) - (pm.Stimulus.fieldofviewHorz/2);
            Centery0 = (pm.Stimulus.spatialSampleVert * Centery0) - (pm.Stimulus.fieldofviewVert/2);
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
            sigmaMinor = (pm.Stimulus.spatialSampleVert * abs(results.params(3)));
            sigmaMajor = (pm.Stimulus.spatialSampleHorz * abs(results.params(3)));
            
            
            tmpTable            = struct2table(results,'AsArray',true);
            tmpTable.Centerx0   = Centerx0;
            tmpTable.Centery0   = Centery0;
            tmpTable.Theta      = 0; % Is circular, we can't model it
            tmpTable.sigmaMajor = sigmaMinor; % Circular, same value 
            tmpTable.sigmaMinor = sigmaMinor; % Circular, same value
            % Make results table smaller (this is for Brian The Substractor :) )
            % Leave the testdata and modelpred so that we have the fit. 
            % Testdata is not exactly the same to the data (=pm.BOLDnoise), he
            % removes I think part of the low frew noise, check later. If we
            % cannot get testdata in other tools, we will just obtain modelpred.
            % We'll need to add pm.BOLDmeanValue to modelpred to have it in the
            % same range as the data (=pm.BOLDnoise)
            tmpTable            = tmpTable(:,{'R2','Centerx0','Centery0', 'Theta' ,...
                                             'rfsize', 'sigmaMinor', 'sigmaMajor' ,...
                                               'testdata','modelpred'});
            pmEstimates = [pmEstimates; tmpTable];
        end
        
    case {'afni', 'afni6', 'afnidog', 'afni_6', 'afni_dog'}
        % AFNI provides 
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
        niftiBOLDfile  = pmForwardModelToNifti(input, 'fname', ...
                                                fullfile(tmpName,'tmp.nii.gz'));
        % Create an AFNI file
        setenv('AFNI_NIFTI_TYPE_WARN','YES');
        if exist(fullfile(tmpName,'tmp.nii.gz'), 'file')
            system(['3dcopy ' fullfile(tmpName,'tmp.nii.gz') ' ' fullfile(tmpName,'tmp')]);
        end
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
        setenv('NIMAGES', '143') % NUMBER OF IMAGES
        system(['3dcopy ' fullfile(tmpName,'images','DSIMG_000.jpg') ...
            ' ' fullfile(tmpName,'images','DSIMAGE_MOVIE')]);
        system(['for NUM in $(count -digits 3 1 $NIMAGES);do 3dTcat -glueto ' ...
                fullfile(tmpName, 'images','DSIMAGE_MOVIE+orig.') ' ' ...
                fullfile(tmpName, 'images','DSIMG_') '${NUM}.jpg;done']);

        % We need to add this info fields
        system(['3drefit -TR ' num2str(pm1.TR) ' ' fullfile(tmpName,'tmp.1D')])
        % system(['3drefit -space orig ' fullfile(tmpName, 'tmp+orig')])
        system(['3drefit -TR ' num2str(pm1.TR) ' ' fullfile(tmpName,'images','DSIMAGE_MOVIE+orig')])
        % system(['3drefit -space orig ' fullfile(tmpName, 'images','DSIMAGE_MOVIE+orig')])
        
        % Obtain the HRF: follow the steps on the 3dNLfim.help file
        system(['3dDeconvolve -nodata 10 ' num2str(pm1.TR) ' -polort -1 ' ...
            '-num_stimts 1 -stim_times 1 "1D:0" GAM ' ...
            '-x1D ' fullfile(tmpName, 'conv.ref.GAM.1D')]);
        
        % Set the enviroment variables
        % Change HRF to other models
        setenv('AFNI_CONVMODEL_REF', fullfile(tmpName, 'conv.ref.GAM.1D'));
        % Implement these other hrf-s for testing
        %{
            TR=2

            echo "Create GAM"
            3dDeconvolve -nodata 50 $TR -polort -1                \
                               -num_stimts 1 -stim_times 1 '1D:0' GAM \
                               -x1D sisar.conv.ref.GAM.1D

            # La GAM por defecto esta normalizada a pico=1

            echo "Create SPMG1"
            3dDeconvolve -nodata 50 $TR -polort -1                \
                               -num_stimts 1 -stim_times 1 '1D:0' 'SPMG1(0)' \
                               -x1D sisar.conv.ref.SPMG1.1D
        %}
        % Add the stimuli we just created
        setenv('AFNI_MODEL_PRF_STIM_DSET', fullfile(tmpName, 'images','DSIMAGE_MOVIE+orig'));
        % Not sure about the options here
        setenv('AFNI_MODEL_DEBUG', '3');
        setenv('AFNI_MODEL_PRF_ON_GRID', 'YES');
        setenv('AFNI_MODEL_PRF_RAM_STATS', 'Y');
        % If it takes too much, change this from 100 to 10 (if images too big, for
        % 100x100, 100 is ok)
        setenv('AFNI_MODEL_PRF_SIGMA_NSTEPS','100');
        
        % Create the command line for the fit
        % Apply it to every voxel
        c = parcluster('local');
        % Adding the absolute path is throwing an error. Not Matlab related, same in CLI
        % Use cd() and then launch the command
        cd(fullfile(tmpName))
        modelName   = 'Conv_PRF';
        modelConstr =   [...
            '-sconstr 0 -10.0 10.0 ' ...
            '-sconstr 1 -1.0 1.0 ' ...
            '-sconstr 2 -1.0 1.0 ' ...
            '-sconstr 3 0.0 1.0 ' ...
            ];
        if strcmp(prfimplementation,'afni6') || strcmp(prfimplementation,'afni_6')
            modelName   = 'Conv_PRF';
            modelConstr =   [...
                '-sconstr 0 -10.0 10.0 ' ...
                '-sconstr 1 -1.0 1.0 ' ...
                '-sconstr 2 -1.0 1.0 ' ...
                '-sconstr 3 0.0 1.0 '
                '-sconstr 4 1.0 1.0 ' ...
                '-sconstr 5 -1.571 1.570 ' ...
                ];
        end
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

        % Read the results back to matlab
        [err, V, Info, ErrMessage] = BrikLoad('snfit_tmp.PRF.1D');
        % Plot the fit and the signal of voxel 2,2,2
        % plot(squeeze(myCube(5,5,2,:)),'k-');hold on;plot(squeeze(V(5,5,2,:)),'r-');
        [berr, bV, bInfo, bErrMessage] = BrikLoad('buck_tmp.PRF.1D');
        
        % Delete the tmp folder
        cd(fullfile(pmRootPath,'local'));
        rmdir(fullfile(tmpName),'s');
        
        % Reshape the data
        fitSeries = reshape(V,  [height(input), pm1.timePointsN]);
        % results   = reshape(bV, [height(input), size(bV,4)]); % Two more params
        results   = bV;
        
        % Results Table
        plot(fitSeries)
%{
        
        % save('fitSeries_Wang_6nogrid.mat', 'fitSeries');
        % save('results_Wang_6nogrid.mat', 'results');
        
        
        
        % Now we want to compare the measured the time series with the predicted
        % time series in the voxels with the largest ratios.
        measured  = reshape(myCube, [500, 144]);
        predicted = reshape(V , [500, 144]);
        results   = reshape(bV, [500, 14]);
        sigRat    = results(:,5);
        
        % decimate for visualization
        % decFact   = 10;
        % measured  = measured(1:decFact:500,:);
        % predicted = predicted(1:decFact:500,:);
        % sigRat    = sigRat(1:decFact:500,:);
        
        % Create bins for ploting
        [Y,E] = discretize(sigRat,5);
        measuredMinusPredicted = measured - predicted;
        
        % Plot the differences in different colours
        for e=1:(length(E)-1)
            eSet = measuredMinusPredicted(Y==e, :);
            % plot(repmat([1:144], [size(eSet,1),1])', eSet'); hold on;
            plot([1:144], mean(eSet)); hold on;
        end
        legend('1<>1.6','1.6<>2.2','2.2<>2.8','2.8<>3.4','3.4<>4')
        title('Categories come from different sigRat values (sigRat goes from 1 to 4 per constraint in the model.)')
        xlabel('Time Points (144, TR=2)')
        ylabel('Mean(Measured - predicted) values (500 voxels)')
%}
        % TODO: Let's define the format for the estimates so that this is
        %       the same for all the methods.
        pmEstimates = table();
        
        % Go line by line and compute required values for each pRF model
        for ii=1:height(input)
            % Assert that every pm has the same TR and the same stimuli, and hrf
            % TR
            assert(DT.pm(index).TR == pm1.TR, 'All BOLDnoise signals TRs should be equal')
            % Stimuli
            assert(DT.pm(index).Stimuli.getStimValues == pm1.Stimuli.getStimValues, ...
                'All Stimuli should be equal')
            % Hrf
            assert(DT.pm(index).HRF.values == pm1.HRF.values, 'All HRF should be equal')
            
            % TODO: use parfor if the number of rows if the table is larger than XX
            pm       = input.pm(ii);
            stimulus = double(pm.Stimulus.getStimValues);
            data     = pm.BOLDnoise;
            TR       = pm.TR;
            options  = p.Results.options;
            
            % Calculate PRF
            results  = analyzePRF({stimulus}, {data}, TR, options);
            
            % Add a new row of results
            pmEstimates = [pmEstimates; struct2table(results,'AsArray',true)];




        end
        
        




case {'popeye'}
        disp('NYI');
    otherwise
        error('Method %s not implemented yet.', prfimplementation)
end


end

