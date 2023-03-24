#! /bin/bash

USERTAG=garikoitz
SOLVER=vista
VERSION=2.2.1_3.1.0

# 2.0.5: edited in master, fix in base/run.py for singularity
# 2.0.8: with the fixed exp
# 2.0.9: after elines change in css exp range boundaries
# 2.1.0: same as 2.0.9 but with Matlab 2020b (using base 3.0.0)
# 2.1.1: same as 2.0.9 but with Matlab 2020b (using base 3.0.1)
# 2.2.0: includes david's changes to read experimental data from new prfprepare
# 2.2.1: updated base to 3.0.3
# 2.2.1_3.1.0: just changed the naming convention to equal it with the other containers. Removed the checks at the top of this file

SCRIPTPATH="$( cd "$(dirname "$0")" && pwd -P )"

echo "Building $USERTAG/prfanalyze-$SOLVER:$VERSION"
echo "  Path: $SCRIPTPATH"
echo "------------------------------------------------------------"
echo ""

docker build "$@" --tag "$USERTAG/prfanalyze-$SOLVER:$VERSION" "$SCRIPTPATH"
