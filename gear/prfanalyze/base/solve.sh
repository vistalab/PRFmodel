#! /bin/bash

echo ""
echo "/solve.sh called with arguments: "
for arg in "$@"
do echo "  *  $arg"
done

echo ""
echo ""
echo "This is a dummy script in prfanalyze-base. You should replace this script with your own"
echo "script to solve the pRFs specified in the above arguments. These arguments are: "
echo "  * options JSON file"
echo "  * input BOLD NIFTI image"
echo "  * input stimulus NIFTI image"
echo "  * stimulus metadata JSON file"
echo "  * output directory"
echo ""
echo ""

exit 1
