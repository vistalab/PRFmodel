#! /bin/bash

USERTAG=garikoitz
SOLVER=afni
VERSION=2.1.0_3.1.1

# 2.1.0_3.1.1: added the possibility to use custom HRFs




SCRIPTPATH="$( cd "$(dirname "$0")" && pwd -P )"

echo "Building $USERTAG/prfanalyze-$SOLVER:$VERSION"
echo "  Path: $SCRIPTPATH"
echo "------------------------------------------------------------"
echo ""

docker build "$@" --tag "$USERTAG/prfanalyze-$SOLVER:$VERSION" "$SCRIPTPATH"
