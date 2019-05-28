function pm = rfCompute(pm)
% Compute the RF from the current parameters
%
%
% See also

    RF = pm.RF;
    stimulus = pm.stimulus;

    % Check if X and Y have been calculated
    if (isempty(stimulus.binary) || isempty(stimulus.fieldofviewHorz) ...
                                 || isempty(stimulus.fieldofviewVert))
        error('Add stimulus and/or field of view to compute the RF.')
    else
        if isempty(stimulus.X) | isempty(stimulus.Y)
            pm = pm.spatialSampleCompute;
            stimulus = pm.stimulus;
        end
    end

    RF.values = rfGaussian2d(stimulus.X, stimulus.Y,...
                             RF.sigmaMajor,RF.sigmaMinor,RF.theta, ...
                             RF.center(1),RF.center(2));

    pm.RF = RF;

end