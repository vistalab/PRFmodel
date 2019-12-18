#! /bin/bash

set +o verbose   # Command echo off

# If run in debug mode, just exec bash:
if [ "$1" = "DEBUG" ]
then exec /bin/bash
else . /etc/profile.d/conda.sh
     conda activate base
fi

# How we print to stdout:
function note {
    echo "$CONTAINER" "   " "$*"
}
function die {
    echo "<ERROR>" "$CONTAINER" "   " "$*"
    exit 1
}


# all we have to do is exec python...
export PRF_SOLVER="aprf"
MCR_ROOT=/opt/mcr/v95/

time /compiled/run_prfanalyze_afni.sh "$MCR_ROOT" "$4" "$2" "$3" "$5"
# Check exit status
[ $? = 0 ] || die "An error occurred during execution of the Matlab executable. Exiting!"

exit 0
