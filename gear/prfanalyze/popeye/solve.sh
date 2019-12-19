#! /bin/bash

# all we have to do is exec python...
export PRF_SOLVER="popeye"
exec python /scripts/run_popeye.py "$1" "$2" "$3" "$4" "$5"
