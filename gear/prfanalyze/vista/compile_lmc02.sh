#!/bin/bash
module load matlab/glerma

cat > build.m <<END
% We do not want ToolboxToolbox to mess up the compiling
restoredefaultpath();
addpath(genpath('/export/home/glerma/glerma/toolboxes/vistasoft'));
 rmpath(genpath('/export/home/glerma/glerma/toolboxes/vistasoft/local'));
addpath(genpath('/export/home/glerma/glerma/toolboxes/jsonlab_v1.2'));
addpath(genpath('/export/home/glerma/glerma/toolboxes/JSONio'));
addpath(genpath('/export/home/glerma/glerma/toolboxes/garikoitzanalyzePRF'));
addpath(genpath('/export/home/glerma/glerma/toolboxes/freesurfer_mrtrix_afni_matlab_tools'));
addpath(genpath('/export/home/glerma/glerma/toolboxes/PRFmodel'));
 rmpath(genpath('/export/home/glerma/glerma/toolboxes/PRFmodel/local'));

mcc -m -R -nodisplay -a /export/home/glerma/glerma/toolboxes/PRFmodel/data -d compiled prfanalyze_vista.m 
exit
END
matlab -nodisplay -nosplash -r build && rm build.m
