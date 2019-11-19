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

# we don't have any post-processing to do at this point (but later we might)
exit 0



