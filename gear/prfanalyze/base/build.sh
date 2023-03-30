#! /bin/bash


# latest=v95(2018b)
# 3.0.0 = v99(2020b)


USERTAG=garikoitz
SOLVER=base
VERSION=3.1.1
# Versions
# 3.0.0 first one using matlab 2020b, but maintaining old Python
# 3.0.1 using the install method with the yml file and conda 3
# 3.0.2 after merging david changes to use experimental data and prfprepare
# 3.0.3 allow for list of subjects and sessions '[001,002]' or 'all' in config.json
# 3.1.0 move from conda to mamba
# 3.1.1 edit the chhmod in the end of run.py; do not dump config.json

SCRIPTPATH="$( cd "$(dirname "$0")" && pwd -P )"

echo "Building $USERTAG/prfanalyze-$SOLVER:$VERSION"
echo "  Path: $SCRIPTPATH"
echo "------------------------------------------------------------"
echo ""

# docker build --no-cache "$@" --tag "$USERTAG/prfanalyze-$SOLVER:$VERSION" "$SCRIPTPATH"
docker build  "$@" --tag "$USERTAG/prfanalyze-$SOLVER:$VERSION" "$SCRIPTPATH"
