function pm = spatialSampleCompute(pm)
% Compute the X,Y and spatialSample from Stimulus values
%
% Syntax
%
% Brief description
%
% Inputs
%
% Outputs
%
% Key/value pairs
%
% GL/BW
%
% See also
%

% Obtain the values if it is a path
stimValues     = pmStimulusRead(pm.stimulus.values);

% Derive variables
spatialSampleHorz = pm.stimulus.fieldofviewHorz/size(stimValues,2);
spatialSampleVert = pm.stimulus.fieldofviewVert/size(stimValues,1);

x = (spatialSampleVert:spatialSampleVert:pm.stimulus.fieldofviewVert);
x = x - mean(x);
y = (spatialSampleHorz:spatialSampleHorz:pm.stimulus.fieldofviewHorz);
y = y - mean(y);

% Set the spatial sampling parameters
[X,Y] = meshgrid(x,y);
pm.stimulus.X = X;
pm.stimulus.Y = Y;

end