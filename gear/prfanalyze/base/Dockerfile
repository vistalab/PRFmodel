# This Dockerfile constructs a docker image, based on the vistalab/freesurfer
# docker image to execute recon-all as a Flywheel Gear.
#
# Example build:
#   docker build --no-cache --tag scitran/freesurfer-recon-all `pwd`
#
# Example usage:
#   docker run -v /path/to/your/subject:/input scitran/freesurfer-recon-all
#
FROM ubuntu:xenial

RUN apt-get update --fix-missing \
 && apt-get install -y wget bzip2 ca-certificates apt-utils \
      libglib2.0-0 libxext6 libsm6 libxrender1 \
      git mercurial subversion curl grep sed dpkg gcc g++ libeigen3-dev zlib1g-dev libqt4-opengl-dev libgl1-mesa-dev libfftw3-dev libtiff5-dev
RUN apt-get install -y libxt6 libxcomposite1 libfontconfig1 libasound2


############################
# Install dependencies
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y \
    apt-utils \
    xvfb \
    xfonts-100dpi \
    xfonts-75dpi \
    xfonts-cyrillic \
    zip \
    unzip \
    python \
    imagemagick \
    wget \
    subversion\
    vim \
    bsdtar

RUN rm /bin/sh && ln -s /bin/bash /bin/sh

############################
# install mamba and env
ENV CONDA_DIR /opt/conda
ENV MAMBA_ROOT_PREFIX="/opt/conda"
RUN wget --quiet https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh -O ~/mamba.sh && \
      /bin/bash ~/mamba.sh -b -p /opt/conda

# Put conda in path so we can use conda activate
ENV PATH=$CONDA_DIR/bin:$PATH

RUN mamba update -n base --all

# install conda env
COPY scientific.yml ./
RUN mamba env create -f scientific.yml
RUN apt-get update && apt-get install -y jq


############################
# Install the MCR dependencies and some things we'll need and download the MCR
# from Mathworks -silently install it
RUN mkdir /mcr-install && mkdir /opt/mcr && cd /mcr-install && \
    wget https://ssd.mathworks.com/supportfiles/downloads/R2020b/Release/6/deployment_files/installer/complete/glnxa64/MATLAB_Runtime_R2020b_Update_6_glnxa64.zip && \
    cd /mcr-install && \
    unzip -q MATLAB_Runtime_R2020b_Update_6_glnxa64.zip && \
    ./install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent && \
    cd / && \
    rm -rf mcr-install

# Configure environment variables for MCR
ENV LD_LIBRARY_PATH /opt/mcr/v99/runtime/glnxa64:/opt/mcr/v99/bin/glnxa64:/opt/mcr/v99/sys/os/glnxa64:/opt/mcr/v99/extern/bin/glnxa64


# Make directory for flywheel spec (v0)
ENV FLYWHEEL /flywheel/v0
RUN mkdir -p ${FLYWHEEL}
WORKDIR ${FLYWHEEL}


# Manage files
COPY run.sh ${FLYWHEEL}/run.sh
RUN chmod 755 ${FLYWHEEL}/run.sh
COPY run.py /scripts/run.py
COPY default_config.json /opt/default_config.json
COPY help.txt /opt/help.txt
COPY solve.sh /solve.sh
RUN chmod 755 /solve.sh

RUN mkdir -p /scripts && chmod 755 /scripts
RUN mkdir -p /running && chmod 755 /running


# Configure entrypoints-
ENTRYPOINT ["/flywheel/v0/run.sh"]

#make it work under singularity
# RUN ldconfig: it fails in BCBL, check Stanford
#https://wiki.ubuntu.com/DashAsBinSh
# RUN rm /bin/sh && ln -s /bin/bash /bin/sh
