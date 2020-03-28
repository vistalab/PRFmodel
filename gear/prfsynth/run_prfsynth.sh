#!/bin/bash
#
# run_prfsynth.sh $basedir $basedir/config.json
#
# For example:
#
#
docker run --rm -it \
        -v $2:/flywheel/v0/input/config.json \
        -v $1:/flywheel/v0/output \
        garikoitz/prfsynth:latest $3
