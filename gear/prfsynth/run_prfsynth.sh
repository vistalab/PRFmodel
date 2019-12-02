#!/bin/bash
# docker run --rm -it \
#        -v $(pwd)/synthBOLDgenerator_paramsExample.json:/flywheel/v0/input/config.json \
#        -v $(pwd)/../../local/output:/flywheel/v0/output \

#        --entrypoint=bash \

docker run --rm -it \
        -v $1:/flywheel/v0/input/config.json \
        -v $2:/flywheel/v0/output \
        garikoitz/prfsynth:1.0.3
