FROM garikoitz/prfanalyze-base:latest

MAINTAINER Garikoitz Lerma-Usabiaga <glerma@stanford.edu>

ADD /compiled /compiled
RUN mkdir -p /compiled/data


# There were no models in the docker container... add them from my local machine
# models where copied from afni/afni docker container 2019-12-30
COPY ./afnimodels/* /opt/afni/


# the solve script and related files
COPY solve.sh /solve.sh
RUN chmod 755 /solve.sh
COPY default_config.json /opt/default_config.json
ENV PRF_SOLVER afni

