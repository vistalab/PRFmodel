#! /bin/bash

# The run script for the prfanalyze base docker.
################################################################################
set +o verbose   # Command echo off

# If run in debug mode, just exec bash:
if [ "$1" = "DEBUG" ]
then exec /bin/bash
elif [ "$1" = "-h" ] || [ "$1" = "--help" ] || [ "$1" = "-help" ] || [ "$1" = "help" ]
then cat /opt/help.txt
     exit 0
else . /opt/conda/etc/profile.d/conda.sh
     conda activate base
fi


# Some variables and functions #################################################

CONTAINER="[garikoitz/prfanalyze]"
# Built to flywheel-v0 spec.
FLYWHEEL_BASE=/flywheel/v0
OUTPUT_DIR="$FLYWHEEL_BASE"/output
INPUT_DIR="$FLYWHEEL_BASE"/input
CONFIG_FILE="$INPUT_DIR"/config.json
# How we print to stdout:
function note {
    echo "$CONTAINER" "   " "$*"
}
function err {
    echo "<ERROR>" "$CONTAINER" "   " "$*"
}
function die {
    echo "<ERROR>" "$CONTAINER" "   " "$*"
    exit 1
}


# Main Script ##################################################################

# If there's no solve script, this is the bae and we shouldn't have been called
[ -x /solve.sh ] || {
    err "PRF Solver script not found: /solver.sh does not exist or is not executable."
    die "Is this the prfanalyze-base docker image?"
}

# If no input is given, we dump the default config and exit
[ -r "$CONFIG_FILE" ] || {
    note "No config file found. Writing default config and exiting."
    cp /opt/default_config.json "$OUTPUT_DIR"/config.json
    exit 0
}

# otherwise, we run the following python code to parse the json and run the
# /solve.sh script!
mkdir -p /running
python /scripts/run.py || die "Python startup script failed!"
# At this point, the files should have been exported to the appropriate directory,
# which should be linked to /running/output_bids
[ -d /running/output_bids ] || die "Pythin startup script failed to make output link!"

# go to the output_bids path and extract subject and session...
cd -P /running/output_bids
sesdir=$PWD
subdir=$(dirname $sesdir)
ses=$(basename $sesdir)
ses=${ses:4}
sub=$(basename $subdir)
sub=${sub:4}

# For any nifti, mat, or JSON file in the output directory, we want to BIDSify it:
if compgen -G "/running/output_bids/*.nii" > /dev/null
then for fl in /running/output_bids/*.nii
     do bnm="`basename $fl .nii`"
        dnm="`dirname $fl`"
        mv "$fl" "${dnm}/sub-${sub}_ses-${ses}_task-prf_${bnm}.nii"
     done
fi
if compgen -G "/running/output_bids/*.nii.gz" > /dev/null
then for fl in /running/output_bids/*.nii.gz
     do bnm="`basename $fl .nii.gz`"
        dnm="`dirname $fl`"
        mv "$fl" "${dnm}/sub-${sub}_ses-${ses}_task-prf_${bnm}.nii.gz"
     done
fi
if compgen -G "/running/output_bids/*.mat" > /dev/null
then for fl in /running/output_bids/*.mat
     do bnm="`basename $fl .mat`"
        dnm="`dirname $fl`"
        mv "$fl" "${dnm}/sub-${sub}_ses-${ses}_task-prf_${bnm}.mat"
     done
fi
if compgen -G "/running/output_bids/*.json" > /dev/null
then for fl in /running/output_bids/*.json
     do bnm="`basename $fl .json`"
        dnm="`dirname $fl`"
        mv "$fl" "${dnm}/sub-${sub}_ses-${ses}_task-prf_${bnm}.json"
     done
fi

# we don't have any post-processing to do at this point (but later we might)
exit 0



