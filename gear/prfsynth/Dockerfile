# FROM flywheel/matlab-mcr:v95
# We want to use v99(2020b) now

FROM ubuntu:xenial


MAINTAINER Garikoitz Lerma-Usabiaga <garikoitz@gmail.com>


# Install jq to parse the JSON config file
RUN apt-get -qq update && apt-get -qq install -y \
    apt-utils \
    unzip \
    xorg \
    wget \
    curl \
    jq

# Install the MCR dependencies and some things we'll need and download the MCR
# from Mathworks -silently install it
RUN mkdir /mcr-install && \
    mkdir /opt/mcr && \
    cd /mcr-install && \
    wget https://ssd.mathworks.com/supportfiles/downloads/R2020b/Release/6/deployment_files/installer/complete/glnxa64/MATLAB_Runtime_R2020b_Update_6_glnxa64.zip && \
    cd /mcr-install && \
    unzip -q MATLAB_Runtime_R2020b_Update_6_glnxa64.zip && \
    ./install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent && \
    cd / && rm -rf mcr-install

# Configure environment variables for MCR
ENV LD_LIBRARY_PATH /opt/mcr/v99/runtime/glnxa64:/opt/mcr/v99/bin/glnxa64:/opt/mcr/v99/sys/os/glnxa64:/opt/mcr/v99/extern/bin/glnxa64

# # Install Conda:
# RUN apt-get -qq update --fix-missing \
#  && apt-get -qq install -y bzip2 ca-certificates \
#       libglib2.0-0 libxext6 libsm6 libxrender1 \
#       git mercurial subversion grep sed dpkg gcc pip
# # # RUN curl -L -o ~/miniconda.sh 'https://repo.anaconda.com/miniconda/Miniconda2-4.5.11-Linux-x86_64.sh'
# RUN curl -L -o ~/miniconda.sh 'https://repo.anaconda.com/miniconda/Miniconda3-py39_4.10.3-Linux-x86_64.sh'
# RUN /bin/bash ~/miniconda.sh -b -p /opt/conda \
#  && rm ~/miniconda.sh \
#  && ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh \
#  && echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc \
#  && echo "conda activate base" >> ~/.bashrc
# # fix profile to contain the conda startup...
# RUN  . /etc/profile.d/conda.sh
# RUN { ./ opt/conda/etc/profile.d/conda/sh; conda activate base; conda install -y numpy scipy cython; pip install neuropythy; }

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
