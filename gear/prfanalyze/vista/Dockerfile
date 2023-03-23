FROM garikoitz/prfanalyze-base:3.0.3

MAINTAINER Garikoitz Lerma-Usabiaga <garikoitz@gmail.com>

ADD /compiled /compiled

# the solve script and related files
COPY solve.sh /solve.sh
RUN chmod 755 /solve.sh
RUN chmod 755 /compiled/run_prfanalyze_vista.sh
RUN chmod 755 /compiled/prfanalyze_vista
COPY default_config.json /opt/default_config.json
ENV PRF_SOLVER vista

