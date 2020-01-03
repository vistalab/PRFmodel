FROM garikoitz/prfanalyze-base:1.0.0

MAINTAINER Garikoitz Lerma-Usabiaga <glerma@stanford.edu>


# To avoid this runtime error:
# Error using signal.internal.sigcheckfloattype
#   Can't reload '/opt/mcr/v95/bin/glnxa64/libmwcoder_types.so'
RUN apt install -y libxt6 libxcomposite1 libfontconfig1 libasound2


ADD /compiled /compiled

# the solve script and related files
COPY solve.sh /solve.sh
RUN chmod 755 /solve.sh
COPY default_config.json /opt/default_config.json
ENV PRF_SOLVER vista

