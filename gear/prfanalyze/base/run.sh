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
VERBOSE=0
FORCE=0 # 1 for force overwrite
FIELDS="task-prf_acq-normal" # default BIDS fields in the filenames
# Built to flywheel-v0 spec.
FLYWHEEL_BASE=/flywheel/v0
OUTPUT_DIR="$FLYWHEEL_BASE"/output
INPUT_DIR="$FLYWHEEL_BASE"/input
DEFAULT_CONFIG_FILE="$INPUT_DIR"/config.json
CONFIG_FILE=""

# How we print to stdout:
function note {
    [ "$VERBOSE" = 1 ] && echo "$CONTAINER" "   " "$*"
}
function err {
    echo "<ERROR>" "$CONTAINER" "   " "$*" >2
}
function die {
    echo "<ERROR>" "$CONTAINER" "   " "$*" >2
    exit 1
}


# Process Arguments ############################################################

while [ "$#" -gt 0 ]
do   case "$1"
     in "--help"|"-h")
            cat /opt/help.txt
            exit 0
            ;;
        "--force"|"-f")
            FORCE=1
            ;;
        "--verbose"|"-v")
            VERBOSE=1
            ;;
        "--bids-fields"|"-b")
            [ "$#" -gt 2 ] || die "Option $1 requires an argument."
            shift
            FIELDS="$1"
            ;;
        *)
            if [ -z "$CONFIG_FILE" ]
            then CONFIG_FILE="$1"
            else die "Too many arguments given to docker"
            fi
            ;;
     esac
     shift
done
[ -z "$CONFIG_FILE" ] && CONFIG_FILE="$DEFAULT_CONFIG_FILE"
([ "$FIELDS" = "-" ] || [ "$FIELDS" = "_" ] || [ "$FIELDS" = " " ]) && FIELDS=""


# Main Script ##################################################################

# If there's no solve script, this is the bae and we shouldn't have been called
[ -x /solve.sh ] || {
    err "PRF Solver script not found: /solver.sh does not exist or is not executable."
    die "Is this the prfanalyze-base docker image?"
}

# If no input is given, we dump the default config and exit
[ -r "$CONFIG_FILE" ] || {
    note "No config file found. Writing default JSON file and exiting."
    cp /opt/default_config.json "$CONFIG_FILE"
    chmod 777 "$CONFIG_FILE"
    exit 0
}

# otherwise, we run the following python code to parse the json and run the
# /solve.sh script!
mkdir -p /running
export FORCE
export VERBOSE
export FIELDS
python /scripts/run.py "$CONFIG_FILE" || die "Python startup script failed!"
# At this point, the files should have been exported to the appropriate directory,
# which should be linked to /running/out
[ -d /running/out ] || die "Python startup script failed to make output link!"

# go to the output_bids path and extract subject and session...
cd -P /running/out
sesdir=$PWD
subdir=$(dirname $sesdir)
ses=$(basename $sesdir)
ses=${ses:4}
sub=$(basename $subdir)
sub=${sub:4}

# For any nifti, mat, or JSON file in the output directory, we want to BIDSify it:
prefix="sub-${sub}_ses-${ses}"
[ -n "$FIELDS" ] && prefix="${prefix}_${FIELDS}"

nn=${#prefix}
if compgen -G "/running/out/*.nii" > /dev/null
then for fl in /running/out/*.nii
     do bnm="`basename $fl .nii`"
        dnm="`dirname $fl`"
        [ "${bnm:0:$nn}" = "$prefix" ] || mv "$fl" "${dnm}/${prefix}_${bnm}.nii"
     done
fi
if compgen -G "/running/out/*.nii.gz" > /dev/null
then for fl in /running/out/*.nii.gz
     do bnm="`basename $fl .nii.gz`"
        dnm="`dirname $fl`"
        [ "${bnm:0:$nn}" = "$prefix" ] || mv "$fl" "${dnm}/${prefix}_${bnm}.nii.gz"
     done
fi
if compgen -G "/running/out/*.mat" > /dev/null
then for fl in /running/out/*.mat
     do bnm="`basename $fl .mat`"
        dnm="`dirname $fl`"
        [ "${bnm:0:$nn}" = "$prefix" ] || mv "$fl" "${dnm}/${prefix}_${bnm}.mat"
     done
fi
if compgen -G "/running/out/*.json" > /dev/null
then for fl in /running/out/*.json
     do bnm="`basename $fl .json`"
        dnm="`dirname $fl`"
        [ "${bnm:0:$nn}" = "$prefix" ] || {
            cat "$fl" | python -m json.tool > "${dnm}/${prefix}_${bnm}.json"
            rm "$fl"
        }
     done
fi

# Handle permissions of the outputs
cd /flywheel/v0/output
find "$OUTPUT_DIR" -type d -exec chmod 777 '{}' ';'
find "$OUTPUT_DIR" -type f -exec chmod 666 '{}' ';'

# we don't have any post-processing to do at this point (but later we might)
exit 0



