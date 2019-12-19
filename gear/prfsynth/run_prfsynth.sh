#!/bin/bash
docker run --rm -it \
        -v $2:/flywheel/v0/input/config.json \
        -v $3:/flywheel/v0/output \
        garikoitz/prfsynth:$1
