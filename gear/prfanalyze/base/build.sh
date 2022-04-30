#! /bin/bash


# latest=v95(2018b)
# 3.0.0 = v99(2020b)


[ -z "$USERTAG" ] && USERTAG=garikoitz
[ -z "$SOLVER"  ] && SOLVER=base
[ -z "$VERSION" ] && VERSION=3.0.2
# Versions
# 3.0.0 first one using matlab 2020b, but maintaining old Python
# 3.0.1 using the install method with the yml file and conda 3
# 3.0.2 after merging david changes to use experimental data and prfprepare

SCRIPTPATH="$( cd "$(dirname "$0")" && pwd -P )"

echo "Building $USERTAG/prfanalyze-$SOLVER:$VERSION"
echo "  Path: $SCRIPTPATH"
echo "------------------------------------------------------------"
echo ""

# docker build --no-cache "$@" --tag "$USERTAG/prfanalyze-$SOLVER:$VERSION" "$SCRIPTPATH"
docker build  "$@" --tag "$USERTAG/prfanalyze-$SOLVER:$VERSION" "$SCRIPTPATH"
