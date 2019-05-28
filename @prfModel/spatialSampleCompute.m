function pm = spatialSampleCompute(pm)
% Compute the X,Y and spatialSample from Stimulus values
%
%
% See also

    spatialSampleHorz = pm.stimulus.fieldofviewHorz/size(pm.stimulus.values,2);
    spatialSampleVert = pm.stimulus.fieldofviewVert/size(pm.stimulus.values,1);

    x = (spatialSampleVert:spatialSampleVert:pm.stimulus.fieldofviewVert);
    x = x - mean(x);
    y = (spatialSampleHorz:spatialSampleHorz:pm.stimulus.fieldofviewHorz);
    y = y - mean(y);
    
    % Set the spatial sampling parameters
    [X,Y] = meshgrid(x,y);
    pm.stimulus.X = X;
    pm.stimulus.Y = Y;
end