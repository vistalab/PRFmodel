#! /bin/bash

USERTAG=garikoitz
SOLVER=aprf
VERSION=2.1.2_3.1.1

# 2.1.0: using Kendrick's latest version instead of the garikoitzanalyzePRF
#        that version had an additionnal option usecss to make it fixed or not
#        the new version has a different way of doing it, we will test 2.1.0 for now
#        and next probably edit the config so that we can control the seeds, kendrick
#        edited the software so that now there is more control over them and we want to exploit that

# 2.1.1: added more option to the config file so that we can use th new stuff, had to fix a couple other
#        things in PRFmodel" now we make sure that aPRF only gets binary stims and we are sure that by default
#        it uses the getcanonicalhrf hrf-s.


# 2.1.2: had to go back to garikoitzanalyzeprf, branch extract_modelpred so that we can have more values in the output
#        now is printing the HRF as well to check that it is the taking the right one 


SCRIPTPATH="$( cd "$(dirname "$0")" && pwd -P )"

echo "Building $USERTAG/prfanalyze-$SOLVER:$VERSION"
echo "  Path: $SCRIPTPATH"
echo "------------------------------------------------------------"
echo ""

docker build "$@" --tag "$USERTAG/prfanalyze-$SOLVER:$VERSION" "$SCRIPTPATH"
