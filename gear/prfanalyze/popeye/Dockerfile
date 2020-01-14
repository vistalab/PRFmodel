FROM garikoitz/prfanalyze-base:latest

MAINTAINER Garikoitz Lerma-Usabiaga <glerma@stanford.edu>

# Install popeye...
RUN . /opt/conda/etc/profile.d/conda.sh \
 && conda activate base \
 && pip install popeye

COPY default_config.json /opt/default_config.json
COPY solve.sh /solve.sh
RUN chmod 755 /solve.sh
COPY run_popeye.py /scripts/run_popeye.py
ENV PRF_SOLVER popeye

