#!/bin/bash
docker run --rm -it \
        --entrypoint=bash \
		-v synthBOLDgenerator_paramsExample.json:/flywheel/v0/input/config.json \
        -v ~/toolboxes/PRFmodel/local/output:/flywheel/v0/output \
        vistalab/synthprf:1.0.0
