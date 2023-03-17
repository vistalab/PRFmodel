FROM garikoitz/prfanalyze-base:3.0.3

LABEL MAINTAINER Garikoitz Lerma-Usabiaga <glerma@stanford.edu>

ADD /compiled /compiled

# the solve script and related files
COPY solve.sh /solve.sh
RUN chmod 755 /solve.sh
COPY default_config.json /opt/default_config.json
ENV PRF_SOLVER aprf

