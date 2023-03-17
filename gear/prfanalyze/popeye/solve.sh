#! /bin/bash

# all we have to do is exec python...
export PRF_SOLVER="popeye"

source /opt/conda/etc/profile.d/conda.sh
conda activate scientific

exec python /scripts/run_popeye.py "$1" "$2" "$3" "$4" "$5"
