#! /bin/bash

# Script for retagging all of the latest prfanalyze dockers with a version id.
################################################################################

function help {
    cat <<EOF
SYNTAX: tag-all.sh <version>
For example, `tag-all.sh 1.0.4`
Creates tagged images with the given version number out of all the latest docker
container images. The latest images are not deleted or changed. If images with
the tag latest have not been built, this script will fail.
EOF
}

# First, find the directory that this script is in:
path="$(cd "$(dirname "$0")" && pwd -P)"
# each directory that contains a build script should be run, starting with base
cd "$path"
for ss in */build.sh
do name=$(dirname $ss)
   echo "garikoitz/prfanalyze-${name}:latest --> garikoitz/prfanalyze-${name}:${ver}"
   docker tag garikoitz/prfanalyze-${name}:latest garikoitz/prfanalyze-${name}:${ver} || \
       echo " * ERROR: Failed to retag $name!"
done

exit 0
