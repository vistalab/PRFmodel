#! /bin/bash
################################################################################

# I/O ##########################################################################

function help {
    cat <<EOF
prfanalyze.sh is a script for running the various prfanalyze-* Docker images
that were produced in conjunction with the paper "A validation framework for
neuroimaging software: the case of population receptive fields" by
Lerma-Usabiaga, Benson, Winawer, and Wandell (DOI: 10.1101/2020.01.07.897991).

In order to run, this script needs three things:
 (1) the name of the prfanalyze docker that is to be run
 (2) the path of the config.json file for running 
 (3) the INPUT directory into which a BIDS/ directory with the results
     should be placed.
Additionally, an input directory (if separate from the INPUT directory) may be
provided.
    
The syntax for running prfanalyze.sh is as follows:
 > prfanalyze.sh <vista|afni|aprf|popeye> <INPUT directory> <config.json path>
 - The input directory may be provided as an optional fourth argument.
 - The following optional arguments may also be provided:
   * --help or -h prints this message then exits.
   * --verbose or -v indicates that this script should print verbose information
     while running.
   * --force or -f may be provided to indicate that files should be overwritten
     on INPUT; otherwise they are placed in a temporary directory in the INPUT
     directory.
   * --version <arg> or -V <arg> specifies that the specific version, <arg> of
     the docker image should be run; by default this is "latest".
   * --base <arg> or -b <arg> indicates that the string <arg> should be used as
     the docker base prefix; by default this is "garikoitz/prfanalyze-".
   * --pull or -p indicates that the docker image should be pulled prior to
     running the docker command.
   * --debug opens the docker in DEBUG mode.
 - If the config.json path is not provided, then a default config.json file is
   usually placed in the INPUT directory by the docker.
EOF
}
function error {
    echo "$@" >&2
}
function die {
    [ "$#" -gt 0 ] && error "$@"
    exit 1
}
function note {
    [ "$VERBOSE" = 1 ] && echo "$@"
}
# This function ensures that a command is run according to the arguments
# provided for redirecting stdout/stderr; if a file has been opened for
# capturing INPUT, commands run later are run in append mode.
PIPE_OPEN=0
function pipecmd {
    if [ "$PIPE_OPEN" = 0 ]
    then if   [ -z "$STDOUT" ] && [ -z "$STDERR" ]
         then "$@"
         elif [ -z "$STDOUT" ] && [ -n "$STDERR" ]
         then "$@" 2>"$STDERR"
         elif [ -n "$STDOUT" ] && [ -z "$STDERR" ]
         then "$@" 1>"$STDOUT"
         elif [ "$STDOUT" = "$STDERR" ]
         then "$@" &>"$STDOUT"
         else "$@" 1>"$STDOUT" 2>"$STDERR"
         fi
         PIPE_OPEN=1
    else if   [ -z "$STDOUT" ] && [ -z "$STDERR" ]
         then "$@"
         elif [ -z "$STDOUT" ] && [ -n "$STDERR" ]
         then "$@" 2>>"$STDERR"
         elif [ -n "$STDOUT" ] && [ -z "$STDERR" ]
         then "$@" 1>>"$STDOUT"
         else "$@" 1>>"$STDOUT" 2>>"$STDERR"
         fi
    fi
}


# Utilities ####################################################################

# `abspath x` puts the resulting filename in ABSPATH (but don't call it with the
# ticks because it will put it in ABSPATH in the child script)
function abspath {
    ABSPATH="$(cd "$(dirname "$1")" && pwd -P)/$(basename "$1")"
}


# Default Options ##############################################################

DOCKERBASE="" # the docker string prefix
VERBOSE=0   # quiet by default
FORCE=0     # don't force overwrite
DEBUG=0     # don't do debug mode by default
PULL=0      # don't pull images by default
STDOUT=""   # these are empty until we parse them from the arguments
STDERR=""
VERSION=""  
SOLVER="" 
INPUT=""
OUTPUT=""
CONFIG=""
STIMDIR=""
SINGIMG=""


OPTS=()



# Arguments ####################################################################

[ "$#" = 0 ] && {
    help
    exit 1
}

while [ "$#" -gt 0 ]
do   case "$1"
     in "--help"|"-h")
            help
            exit 0
            ;;
        "--verbose"|"-v")
            VERBOSE=1
            ;;
        "--out"|"-o")
            [ "$#" -lt 2 ] && die "Option $1 requies an argument."
            shift
            STDOUT="$1"
            ;;
        "--err"|"-e")
            [ "$#" -lt 2 ] && die "Option $1 requies an argument."
            shift
            STDERR="$1"
            ;;
        "--debug"|"--DEBUG")
            DEBUG=1
            ;;
        *)  # this is one of config, INPUT, input...
            if   [ -z "$SOLVER" ]
            then SOLVER="$1"
            elif [ -z "$INPUT" ]
            then INPUT="$1"
            elif [ -z "$CONFIG" ]
            then CONFIG="$1"
            elif [ -z "$STIMDIR" ]
            then STIMDIR="$1"
            elif [ -z "$OUTPUT" ]
            then OUTPUT="$1"
            elif [ -z "$SINGIMG" ]
            then SINGIMG="$1"
            else # we've reached the point where everything that follows is
                 # just extra arguments for the solver...
                 OPTS=( "${OPTS[@]}" "$1" )
            fi
            ;;
     esac
     shift
done
# If not all required inputs were provided, error out.
[ -z "$SOLVER" ] && {
    error 'The pRF solver is a required argument!'
    die 'Use --help flag to see usage documentation'
}
[ -z "$INPUT" ] && {
    error 'The INPUT directory is a required argument!'
    die 'Use --help flag to see usage documentation'
}


# we can construct the docker name now:
# get the abspath for the config file and INPUT directories
abspath "$INPUT"
INPUT="$ABSPATH"
[ -z "$CONFIG" ] || {
    abspath "$CONFIG"
    CONFIG="$ABSPATH"
}
[ -z "$STIMDIR" ] || {
    abspath "$STIMDIR"
    STIMDIR="$ABSPATH"
}

[ -z "$OUTPUT" ] || {
    abspath "$OUTPUT"
    OUTPUT="$ABSPATH"
}

[ -z "$SINGIMG" ] || {
    abspath "$SINGIMG"
    SINGIMG="$ABSPATH"
}


tmp_dir=$(mktemp -d -t ci-XXXXXXXXXX)
#abspath "$tmp_dir"
#TMP="$ABSPATH"


# Okay, let's print a diagnostic message (if verbose is on)
note "Running prfanalyze.sh with the following options:"
note "   solver:   $SOLVER"
note "   singimg:  $SINGIMG"
note "   config:   $CONFIG"
note "   input:    $INPUT"
note "   stimdir:  $STIMDIR"
note "   output:   $OUTPUT"

if   [ ${#OPTS[@]} -gt 0 ]
then note "  options: ${OPTS[@]}"
else note "  options: <none>"
fi
if   [ $DEBUG = 1 ]
then note "  using DEBUG mode"
fi


# Run the Singularity image ###############################################################

# Build up the command line arguments, piece by piece



if [ -d "$CONFIG" ]
then ARGS_IN=("-B" "${CONFIG}:/flywheel/v0/input:rw"
              "-B" "${INPUT}:/flywheel/v0/input/BIDS:rw"
              "-B" "${STIMDIR}:/flywheel/v0/input/BIDS/stimuli:rw"
              "-B" "${tmp_dir}:/running/:rw")

fi
ARGS_OUT=("-B" "${OUTPUT}:/flywheel/v0/output")



if [ $VERBOSE = 1 ]
then OPTS=("${OPTS[@]}" "--verbose")  
fi

# Put them all together:

if  [ $DEBUG = 1 ]
then
ARGS=("${OPTS[@]}" "shell" "${ARGS_IN[@]}" "${ARGS_OUT[@]}" "$SINGIMG")
else
ARGS=("${OPTS[@]}" "run" "${ARGS_IN[@]}" "${ARGS_OUT[@]}" "$SINGIMG")
fi 


note "Running Singularity command:"
note "  > singularity run \\"
[ "${#ARGS_IN[@]}" = 0 ] || note "  |  " "${ARGS_IN[@]}" " \\"
note "  |   " "${ARGS_OUT[@]}" " \\"
note "  |   \"$SINGIMG\" ${OPTS[@]}"
note ""
note "------------------------------------------------------------"
note ""

pipecmd singularity "${ARGS[@]}"
SINGEXIT="$?"

note "Singularity has finished with exit code ${SINGEXIT}."

exit $SINGEXIT
