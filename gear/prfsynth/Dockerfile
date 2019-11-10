FROM flywheel/matlab-mcr:v95

MAINTAINER Garikoitz Lerma-Usabiaga <glerma@stanford.edu>

# Install jq to parse the JSON config file
RUN apt-get update -qq && apt-get install -y jq

# Install Conda:
RUN apt-get update --fix-missing \
 && apt-get install -y wget bzip2 ca-certificates \
      libglib2.0-0 libxext6 libsm6 libxrender1 \
      git mercurial subversion curl grep sed dpkg gcc
RUN curl -L -o ~/miniconda.sh 'https://repo.anaconda.com/miniconda/Miniconda2-4.5.11-Linux-x86_64.sh'
RUN /bin/bash ~/miniconda.sh -b -p /opt/conda \
 && rm ~/miniconda.sh \
 && ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh \
 && echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc \
 && echo "conda activate base" >> ~/.bashrc
RUN { . /opt/conda/etc/profile.d/conda.sh; conda update -n base -y conda; }
RUN { . /opt/conda/etc/profile.d/conda.sh; conda activate base; conda install -y numpy scipy cython; pip install neuropythy popeye; }

# Make directory for flywheel spec (v0)
ENV FLYWHEEL /flywheel/v0
RUN mkdir -p ${FLYWHEEL}
# COPY manifest.json ${FLYWHEEL}/manifest.json
ADD /compiled /compiled

# Run script
COPY run ${FLYWHEEL}/run
RUN chmod +x ${FLYWHEEL}/run

# Configure entrypoints-
ENTRYPOINT ["/flywheel/v0/run"]
