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
            
            % Add a new row of results
            pmEstimates = [pmEstimates; struct2table(results,'AsArray',true)];
        end
        
    case {'afni'}
        warning('For AFNI analysis, be sure that all options have the same TR and the same stimulus')
        tmpName = tempname(fullfile(pmRootPath,'local'));
        mkdir(tmpName);
        % Create a tmp nifti file and convert it to a tmp AFNI format
        niftiBOLDfile  = pmForwardModelToNifti(input, 'fname', ...
                                                fullfile(tmpName,'tmp.nii.gz'));
        % Create an AFNI file
        if exist(fullfile(tmpName,'tmp.nii.gz'), 'file')
            system(['3dcopy ' fullfile(tmpName,'tmp.nii.gz') ' ' fullfile(tmpName,'tmp+orig')]);
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
        system(['3drefit -TR ' num2str(pm1.TR) ' ' fullfile(tmpName,'tmp+orig')])
        system(['3drefit -space orig ' fullfile(tmpName, 'images','DSIMAGE_MOVIE+orig')])
        
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
        setenv('AFNI_MODEL_PRF_SIGMA_NSTEPS','10');
        
        % Create the command line for the fit
        % Apply it to every voxel
        c = parcluster('local');
        system([...
            '3dNLfim -input ' fullfile(tmpName,'tmp+orig') ' ' ...
            '-noise Zero ' ... 
            '-signal Conv_PRF ' ... % Conv_PRF_6 for the 6 param model, it needs '-sconstr 4 1.0 1.0 ' ...'-sconstr 5 -1.571 1.570 ' ...
            '-sconstr 0 -10.0 10.0 ' ...
            '-sconstr 1 -1.0 1.0 ' ...
            '-sconstr 2 -1.0 1.0 ' ...
            '-sconstr 3 0.0 1.0 ' ...
            '-BOTH ' ...
            '-nrand 10000 ' ...
            '-nbest 5 ' ...
            '-bucket 0 '  fullfile(tmpName,'buck_tmp.PRF ') ...
            '-snfit ' fullfile(tmpName,'snfit_tmp.PRF ') ...
            '-TR ' num2str(pm1.TR) ' ' ...
            '-jobs ' num2str(c.NumWorkers)...
            ]);

        % Read the results back to matlab
        [err, V, Info, ErrMessage] = BrikLoad(fullfile(tmpName,'snfit_tmp.PRF'));
        % Plot the fit and the signal of voxel 2,2,2
        % plot(squeeze(myCube(5,5,2,:)),'k-');hold on;plot(squeeze(V(5,5,2,:)),'r-');
        [berr, bV, bInfo, bErrMessage] = BrikLoad(fullfile(tmpName,'buck_tmp.PRF '));
        



%{
        fitSeries = reshape(V,  [500, 144]);
        results   = reshape(bV, [500, 14]); % Two more params
        save('fitSeries_Wang_6nogrid.mat', 'fitSeries');
        save('results_Wang_6nogrid.mat', 'results');
        
        
        
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

