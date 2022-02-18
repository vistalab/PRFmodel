#!/bin/bash
module load matlab/glerma

cat > build.m <<END
% We do not want ToolboxToolbox to mess up the compiling
restoredefaultpath();
addpath(genpath('/bcbl/home/home_g-m/glerma/toolboxes/vistasoft'));
addpath(genpath('/bcbl/home/home_g-m/glerma/toolboxes/jsonlab_v1.2'));
addpath(genpath('/bcbl/home/home_g-m/glerma/toolboxes/JSONio'));
addpath(genpath('/bcbl/home/home_g-m/glerma/toolboxes/garikoitzanalyzePRF'));
addpath(genpath('/bcbl/home/home_g-m/glerma/toolboxes/freesurfer_mrtrix_afni_matlab_tools'));
addpath(genpath('/bcbl/home/home_g-m/glerma/toolboxes/PRFmodel'));
rmpath(genpath('/bcbl/home/home_g-m/glerma/toolboxes/PRFmodel/local'));

mcc -m -R -nodisplay -a /bcbl/home/home_g-m/glerma/toolboxes/PRFmodel/data -d compiled prfanalyze_vista.m 
exit
END
matlab -nodisplay -nosplash -r build && rm build.m
