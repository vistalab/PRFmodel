#!/bin/bash
# module load matlab/2017a

cat > build.m <<END

addpath(genpath('/data/localhome/glerma/toolboxes/vistasoft'));
addpath(genpath('/data/localhome/glerma/toolboxes/jsonlab'));
addpath(genpath('/data/localhome/glerma/toolboxes/JSONio'));
addpath(genpath('/data/localhome/glerma/toolboxes/garikoitzanalyzePRF'));
addpath(genpath('/data/localhome/glerma/toolboxes/freesurfer_mrtrix_afni_matlab_tools'));
addpath(genpath('/data/localhome/glerma/toolboxes/PRFmodel'));

% mcc -m -R -nodisplay -a /data/localhome/glerma/soft/encode/mexfiles -a /data/localhome/glerma/soft/vistasoft/mrDiffusion -d compiled AFQ_StandAlone_QMR.m
mcc -m -R -nodisplay -d compiled prfvalidation_wrapper.m
exit
END
Matlabr2017a -nodisplay -nosplash -r build && rm build.m
