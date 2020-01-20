FROM nben/prfanalyze-afni:latest

MAINTAINER Garikoitz Lerma-Usabiaga <glerma@stanford.edu>

ADD /compiled /compiled
RUN mkdir -p /compiled/data

# the solve script and related files
COPY solve.sh /solve.sh
RUN chmod 755 /solve.sh
COPY default_config.json /opt/default_config.json
ENV PRF_SOLVER afni

