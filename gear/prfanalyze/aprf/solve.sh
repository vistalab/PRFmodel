#! /bin/bash

set +o verbose   # Command echo off

# If run in debug mode, just exec bash:
if [ "$1" = "DEBUG" ]
then exec /bin/bash
else
    source /opt/conda/etc/profile.d/conda.sh
    conda activate scientific
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
MCR_ROOT=/opt/mcr/v99/

time /compiled/run_prfanalyze_aprf.sh "$MCR_ROOT" "$1" "$4" "$2" "$3" "$5"
# Check exit status
[ $? = 0 ] || die "An error occurred during execution of the Matlab executable. Exiting!"

exit 0
