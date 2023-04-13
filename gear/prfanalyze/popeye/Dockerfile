FROM garikoitz/prfanalyze-base:3.1.1

LABEL MAINTAINER Garikoitz Lerma-Usabiaga <glerma@stanford.edu>

# Install popeye...
RUN . /opt/conda/etc/profile.d/conda.sh \
 && conda activate scientific \
 && mamba install cython \
 && pip install git+https://github.com/dlinhardt/popeye

COPY default_config.json /opt/default_config.json
COPY solve.sh /solve.sh
RUN chmod 755 /solve.sh
COPY run_popeye.py /scripts/run_popeye.py
ENV PRF_SOLVER popeye

