#!/bin/bash
# module load matlab/2018b

# Try loading the module if we know about the module command
if   [ -d "$MODULESHOME" ] && [ -z "`which matlab`" ]
then module load matlab/2018b
fi

# This is the directory in which the toolboxes are stored
if   [ -n "$1" ]
then TOOLBOXTOOLBOX_PATH="$1"
elif [ -n "$TOOLBOXTOOLBOX_PATH" ]
then : # noop
elif [ -d /data/localhome/glerma/toolboxes ]
then TOOLBOXTOOLBOX_PATH="/data/localhome/glerma/soft/ToolboxToolbox"
else echo "Cannot find toolboxes path!" >&2
     exit 1
fi

# Look for a version of matlab we can use:
if   [ -x "`which Matlabr2018b`" ]
then MLAB="`which Matlabr2018b`"
elif [ -x "/Applications/MATLAB9.5.app/bin/matlab" ]
then MLAB="/Applications/MATLAB9.5.app/bin/matlab"
elif [ -x "/software/matlab/r2018b/bin/matlab" ]
then MLAB="/software/matlab/r2018b/bin/matlab"
elif [ -x "`which matlab`" ]
then echo "Warning: Unsure if matlab is version 2018b" >&2
     MLAB="`which matlab`"
elif [ -x "`which MATLAB`" ]
then echo "Warning: Unsure if matlab is version 2018b" >&2
     MLAB="`which MATLAB`"
else echo "Cannot find a matlab executable!"
     exit 1
fi

cat > build.m <<EOF

%% Make sure the correct toolbox-toolbox directories are being used:
originalDir = pwd();
toolboxToolboxDir = '${TOOLBOXTOOLBOX_PATH}';
% Set up the path.
apiDir = fullfile(toolboxToolboxDir, 'api');
cd(apiDir);
tbResetMatlabPath('reset', 'full');
cd(originalDir);
% Also make sure the toolboxes dir is right:
setpref('ToolboxToolbox', 'toolboxRoot', fullfile('$TOOLBOXTOOLBOX_PATH', 'toolboxes'));

% add use the appropriate toolboxes so that they get added
tbUse('vistasoft');
tbUse('jsonlab');
tbUse('JSONio');
tbUse('garikoitzanalyzePRF');
tbUse('freesurfer_mrtrix_afni_matlab_tools');
tbUse('PRFmodel');

fprintf('Adding paths...\n');
addpath(genpath('${TOOLBOXTOOLBOX_PATH}/toolboxes/vistasoft'));
addpath(genpath('${TOOLBOXTOOLBOX_PATH}/toolboxes/jsonlab_v1.2'));
addpath(genpath('${TOOLBOXTOOLBOX_PATH}/toolboxes/JSONio'));
addpath(genpath('${TOOLBOXTOOLBOX_PATH}/toolboxes/garikoitzanalyzePRF'));
addpath(genpath('${TOOLBOXTOOLBOX_PATH}/toolboxes/freesurfer_mrtrix_afni_matlab_tools'));
addpath(genpath('${TOOLBOXTOOLBOX_PATH}/toolboxes/PRFmodel'));

fprintf("Running mcc...\n");
% mcc -m -R -nodisplay -a /data/localhome/glerma/soft/encode/mexfiles -a /data/localhome/glerma/soft/vistasoft/mrDiffusion -d compiled AFQ_StandAlone_QMR.m
% mcc -m -R -nodisplay -d compiled prfvalidation_wrapper.m
mcc -mv -R -nodisplay -a ${TOOLBOXTOOLBOX_PATH}/toolboxes/PRFmodel/data -d compiled synthBOLDgenerator.m 
exit
EOF

mkdir -p compiled && "$MLAB" -nodesktop -nodisplay -nosplash -r build && rm build.m

