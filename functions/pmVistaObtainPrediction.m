function pmEstimates = pmVistaObtainPrediction(pmEstimates, results)
    % Initialize the results
    pmEstimates.modelpred = pmEstimates.testdata;
    testdata              = pmEstimates.testdata;
    modelprediction       = zeros(size(testdata));

    numVoxels = size(pmEstimates.modelpred,1);

    % start from the function rmPlotGUI_makePrediction
    recompFit = true;  % We want to recompute fit
        
    % paramsOrder is as follow : sigmaMajor,sigmaMinor,theta, x0,y0
    % rfGaussian2d(X, Y,         rfParams(3), rfParams(5), rfParams(6), rfParams(1), rfParams(2));
    % rfParams(4) is not used, at least in this function
    % I need to create rfParams = []; per every voxel it seems
    
    parfor voxel=1:numVoxels 
        thisTestData = testdata(voxel,:);
        if isfield(results.model{1},'exponent')
            rfParams    = zeros([1,7]);
        else
            rfParams    = zeros([1,6]);
        end
        rfParams(1) = results.model{1}.x0(voxel);
        rfParams(2) = results.model{1}.y0(voxel);
        rfParams(3) = results.model{1}.sigma.major(voxel);
        % rfParams(4) =  % Esto se asigna luego, suele ser el beta(1)
        rfParams(5) = results.model{1}.sigma.minor(voxel);
        rfParams(6) = results.model{1}.sigma.theta(voxel);
        if isfield(results.model{1},'exponent')
            rfParams(7) = results.model{1}.exponent(voxel);
        end
        
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
        
        % In resutls we don't have M.tSeries
        % I am going to paste the testdata in vertical form. 
        % M.tSeries = pmEstimates.testdata';
        M.tSeries = thisTestData';
        
        switch M.model{1}.description
            case {'2D pRF fit (x,y,sigma, positive only)',...
                    '2D RF (x,y,sigma) fit (positive only)',...
                    '1D pRF fit (x,sigma, positive only)'}
                if recompFit==0
                    beta = rmCoordsGet(M.viewType, model, 'b', coords);
                    beta = beta([1 dcid+1])';
                    
                else
                    beta = pinv([pred trends(:,dcid)])*M.tSeries;
                    beta(1) = max(beta(1),0);
                    
                end
                
                RFs        = RFs .* (beta(1) .* M.params.analysis.HrfMaxResponse);
                
                rfParams(4) = beta(1);
                
                
            case {'2D pRF fit (x,y,sigma_major,sigma_minor)' ...
                    'oval 2D pRF fit (x,y,sigma_major,sigma_minor,theta)'}
                if recompFit==0
                    beta = rmCoordsGet(M.viewType, model, 'b', coords);
                    beta = beta([1 dcid+1]);
                else
                    beta = pinv([pred trends(:,dcid)])*M.tSeries;
                    beta(1) = max(beta(1),0);
                    
                end
                
                RFs        = RFs .* (beta(1) .* M.params.analysis.HrfMaxResponse);
                
                rfParams(4) = beta(1);
                
            case 'unsigned 2D pRF fit (x,y,sigma)'
                if recompFit==0,
                    beta = rmCoordsGet(M.viewType, model, 'b', coords);
                    beta = beta([1 dcid+1]);
                else
                    beta = pinv([pred trends(:,dcid)])*M.tSeries;
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
                    beta = pinv([pred trends(:,dcid)])*M.tSeries;
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
                    beta = pinv([pred trends(:,dcid)])*M.tSeries;
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
                    beta = pinv([pred trends(:,dcid)])*M.tSeries;
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
                    beta = pinv([pred trends(:,dcid)])*M.tSeries;
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
                    beta = pinv([pred trends(:,dcid)])*M.tSeries;
                    beta(1) = max(beta(1),0);
                    
                else
                    beta = rmCoordsGet(M.viewType, model, 'b', coords);
                    beta = beta([1 dcid+1])';   % scale factor (gain) and DC (mean)
                end
                
                RFs        = RFs .* (beta(1) .* M.params.analysis.HrfMaxResponse);
                
                rfParams(4) = beta(1);
                
                
                
            otherwise
                error('Unknown modelName: %s', modelName);
        end
        
        % Calculate the prediction
        prediction = [pred trends(:,dcid)] * beta;
        
        % Add it to the return table
        modelprediction(voxel,:) = prediction';
    end
    pmEstimates.modelpred = modelprediction;
end


