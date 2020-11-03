%% KAI TESTS: Create noisefree time series

stimshuffle = true;


COMBINE_PARAMETERS                           = struct();

COMBINE_PARAMETERS.TR                    = [1];
COMBINE_PARAMETERS.RF                    = struct();
COMBINE_PARAMETERS.RF.Centerx0           = [0,3];
COMBINE_PARAMETERS.RF.Centery0           = [0,3];
COMBINE_PARAMETERS.RF.sigmaMajor         = [1,2]; % this is sigma of x deg radius
COMBINE_PARAMETERS.RF.sigmaMinor         = "same"; 
COMBINE_PARAMETERS.Stimulus.durationSecs = 300;
COMBINE_PARAMETERS.signalPercentage = 'frac';

if stimshuffle
    COMBINE_PARAMETERS.Stimulus.Shuffle   = true;
end

HRF(1).Type                      = 'friston';
HRF(1).params.a                  = [6 12];
HRF(1).params.b                  = [0.9 0.9];
HRF(1).params.c                  = 0.35;

HRF(2).Type                      = 'friston';
HRF(2).params.a                  = [7 12];
HRF(2).params.b                  = [0.8 0.8];
HRF(2).params.c                  = 0.4;



COMBINE_PARAMETERS.HRF           = HRF;


Noise(1).seed                    = 'none';
COMBINE_PARAMETERS.Noise         = Noise;

synthDT = pmForwardModelTableCreate(COMBINE_PARAMETERS, 'repeats', 1);
synthDT = pmForwardModelCalculate(synthDT,'useparallel',false);



% To visualize time series
synthDT(1,:).pm.plot('color','b','window',true,'what','nonoise');hold on;
synthDT(2,:).pm.plot('color','r','window',false,'what','nonoise')

% To get the time series
synthDT(1,:).pm.BOLD


% Where are the classes to obtain the formulas
% HRF
edit pmHRF


% RF
edit pmRF
% We use this vistasoft function to generate the gaussian RF
pmGaussian2d(XY{1}, XY{2}, rf.sigmaMajor,rf.sigmaMinor,rf.Theta, rf.Centerx0,rf.Centery0);