#!/bin/bash
# module load matlab/2018b

cat > build.m <<END
% We do not want ToolboxToolbox to mess up the compiling
restoredefaultpath();
addpath(genpath('/Users/glerma/toolboxes/vistasoft'));
addpath(genpath('/Users/glerma/toolboxes/jsonlab_v1.2'));
addpath(genpath('/Users/glerma/toolboxes/JSONio'));
addpath(genpath('/Users/glerma/toolboxes/garikoitzanalyzePRF'));
addpath(genpath('/Users/glerma/toolboxes/freesurfer_mrtrix_afni_matlab_tools'));
addpath(genpath('/Users/glerma/toolboxes/PRFmodel'));

mcc -m -R -nodisplay -a /Users/glerma/toolboxes/PRFmodel/data -d compiled synthBOLDgenerator.m
exit
END
 # /software/matlab/r2018b/bin/matlab -nodisplay -nosplash -r build && rm build.m
/Applications/MATLAB_R2018b.app/bin/matlab -nodisplay -nosplash -r build && rm build.m
