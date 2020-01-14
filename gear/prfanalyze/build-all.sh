#! /bin/bash

# Script for building all of the prfanalyze docker images.
################################################################################

# First, find the directory that this script is in:
path="$(cd "$(dirname "$0")" && pwd -P)"
# each directory that contains a build script should be run, starting with base
cd "$path"
echo "Building base"
echo "================================================================================"
base/build.sh --no-cache || {
    echo ""
    echo "Failed to build base!"
    exit 1
}
for ss in */build.sh
do name=$(dirname $ss)
   [ $name = base ] && continue
   echo ""
   echo ""
   echo "Building $name"
   echo "================================================================================"
   $ss --no-cache || {
       echo ""
       echo "Failed to build $name!"
       exit 1
   }
done

echo ""
echo "All builds complete."

exit 0
