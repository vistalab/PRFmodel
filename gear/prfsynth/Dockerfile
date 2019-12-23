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

# Install AFNI:
# RUN apt-get install -y \
#             make cmake \
#             zlib1g-dev libxt-dev libxtst-dev \
#             build-essential bzip2 ca-certificates \
#             curl freeglut3-dev g++ gcc git \
#             libglew-dev libglib2.0-dev libglu1-mesa-dev \
# libglw1-mesa-dev libgsl-dev libjpeg62 \
#             libmotif-dev libnetcdf-dev libtool libxi-dev \
#             libxmhtml-dev libxmu-dev libxpm-dev libxt-dev \
#             m4 ncurses-dev ninja-build pkg-config r-base \
#             rsync tcsh vim wget xvfb
# RUN mkdir -p /opt/code && cd /opt/code && git clone https://github.com/afni/afni
# RUN cd /opt/code/afni && chmod 755 `find . -name 'configure*'`
# RUN cd /opt/code/afni/src \
#  && cat Makefile.linux_ubuntu_16_64 | sed -E 's|INSTALLDIR = .+|INSTALLDIR = /opt/afni|g' > Makefile \
#  && make all \
#  && make install
# RUN echo 'export PATH="/opt/afni:$PATH"' >> ~/.bashrc

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
