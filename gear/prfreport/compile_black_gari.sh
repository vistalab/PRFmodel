#!/bin/bash
# module load matlab/2018b

cat > build.m <<END
% We do not want ToolboxToolbox to mess up the compiling
restoredefaultpath();
addpath(genpath('/data/localhome/glerma/toolboxes/jsonlab_v1.2'));
addpath(genpath('/data/localhome/glerma/toolboxes/vistasoft'));
addpath(genpath('/data/localhome/glerma/toolboxes/JSONio'));
addpath(genpath('/data/localhome/glerma/toolboxes/PRFmodel'));

mcc -m -R -nodisplay -d compiled prfReportWrapper.m 
exit
END
/software/matlab/r2018b/bin/matlab -nodisplay -nosplash -r build && rm build.m
