FROM garikoitz/prfanalyze-base:3.1.1

LABEL MAINTAINER Garikoitz Lerma-Usabiaga <garikoitz@gmail.com>

ADD /compiled /compiled
RUN chmod -R 777 /compiled

# the solve script and related files
COPY solve.sh /solve.sh
RUN chmod 755 /solve.sh
COPY default_config.json /opt/default_config.json
ENV PRF_SOLVER aprf

