FROM garikoitz/prfanalyze-base:3.1.0

LABEL MAINTAINER Garikoitz Lerma-Usabiaga <glerma@stanford.edu>

# Install AFNI:
RUN apt-get install -y \
            make cmake \
            zlib1g-dev libxt-dev libxtst-dev \
            build-essential bzip2 ca-certificates \
            curl freeglut3-dev g++ gcc git \
            libglew-dev libglib2.0-dev libglu1-mesa-dev \
            libglw1-mesa-dev libgsl-dev libjpeg62 \
            libmotif-dev libnetcdf-dev libtool libxi-dev \
            libxmhtml-dev libxmu-dev libxpm-dev libxt-dev \
            m4 ncurses-dev ninja-build pkg-config r-base \
            rsync tcsh vim wget xvfb
RUN mkdir -p /opt/code && cd /opt/code && git clone https://github.com/afni/afni
RUN cd /opt/code/afni/src \
 && cat Makefile.linux_ubuntu_16_64 | sed -E 's|INSTALLDIR = .+|INSTALLDIR = /opt/afni|g' > Makefile \
 && make all \
 && make install
RUN echo 'export PATH="$PATH:/opt/afni"' >> ~/.bashrc

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

