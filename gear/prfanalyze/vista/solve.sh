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
export PRF_SOLVER="vista"
MCR_ROOT=/opt/mcr/v99/




export LD_LIBRARY_PATH="/opt/mcr/v99/runtime/glnxa64:/opt/mcr/v99/bin/glnxa64:/opt/mcr/v99/sys/os/glnxa64:/opt/mcr/v99/extern/bin/glnxa64"



time /compiled/run_prfanalyze_vista.sh "$MCR_ROOT" "$1" "$4" "$2" "$3" "$5"
# Check exit status
[ $? = 0 ] || die "An error occurred during execution of the Matlab executable. Exiting!"

exit 0
