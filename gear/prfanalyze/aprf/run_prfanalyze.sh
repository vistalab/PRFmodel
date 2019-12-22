docker run --rm -it -v "$2":/flywheel/v0/input -v "$3":/flywheel/v0/output garikoitz/prfanalyze-aprf:$1 $4
