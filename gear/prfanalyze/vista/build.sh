#! /bin/bash

[ -z "$USERTAG" ] && USERTAG=garikoitz
[ -z "$SOLVER"  ] && SOLVER=vista
[ -z "$VERSION" ] && VERSION=2.1.0

# 2.0.5: edited in master, fix in base/run.py for singularity

SCRIPTPATH="$( cd "$(dirname "$0")" && pwd -P )"

echo "Building $USERTAG/prfanalyze-$SOLVER:$VERSION"
echo "  Path: $SCRIPTPATH"
echo "------------------------------------------------------------"
echo ""

docker build "$@" --tag "$USERTAG/prfanalyze-$SOLVER:$VERSION" "$SCRIPTPATH"
