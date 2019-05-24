function pm = rfCompute(pm)
% Compute the RF from the current parameters
%
%
% See also

RF = pm.RF;
stimulus = pm.stimulus;

RF.values = ...
    rfGaussian2d(stimulus.X, stimulus.Y,...
      RF.sigmaMajor,RF.sigmaMinor,RF.theta, ...
      RF.center(1),RF.center(2));

pm.RF = RF;

end