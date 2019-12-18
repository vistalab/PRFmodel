#!/bin/bash
docker run --rm -it \
        -v $1:/flywheel/v0/input/config.json \
        -v $2:/flywheel/v0/output \
        garikoitz/prfreport:1.0.0
