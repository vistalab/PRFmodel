docker run --rm -it $4 \
      -v $2:/flywheel/v0/input \
      -v $2:/flywheel/v0/output \
      garikoitz/prfanalyze-afni:$1 $3
