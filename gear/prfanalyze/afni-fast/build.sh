#! /bin/bash

[ -z "$USERTAG" ] && USERTAG=garikoitz
[ -z "$SOLVER"  ] && SOLVER=afni
[ -z "$VERSION" ] && VERSION=2.0.0

# 2.0.0 base 3.0.3

SCRIPTPATH="$( cd "$(dirname "$0")" && pwd -P )"

echo "Building $USERTAG/prfanalyze-$SOLVER:$VERSION"
echo "  Path: $SCRIPTPATH"
echo "------------------------------------------------------------"
echo ""

docker build "$@" --tag "$USERTAG/prfanalyze-$SOLVER:$VERSION" "$SCRIPTPATH"
