#!/bin/bash
docker run --rm -it \
        -v $(pwd)/synthBOLDgenerator_paramsExample.json:/flywheel/v0/input/config.json \
        -v $(pwd)/../../local/output:/flywheel/v0/output \
        garikoitz/synthprf:1.0.0
