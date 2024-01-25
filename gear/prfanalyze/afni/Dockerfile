FROM garikoitz/prfanalyze-base:3.1.1

LABEL MAINTAINER Garikoitz Lerma-Usabiaga <garikoitz@gmail.com>

ENV PATH="/opt/afni:$PATH" \
    AFNI_PLUGINPATH="/opt/afni"

RUN apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
           ca-certificates \
           cmake \
           curl \
           ed \
           gsl-bin \
           libcurl4-openssl-dev \
           libgl1-mesa-dri \
           libglib2.0-0 \
           libglu1-mesa-dev \
           libglw1-mesa \
           libgomp1 \
           libjpeg-turbo8-dev \
           libjpeg62 \
           libssl-dev \
           libudunits2-dev \
           libxm4 \
           multiarch-support \
           netpbm \
           tcsh \
           xfonts-base \
           xvfb \
    && rm -rf /var/lib/apt/lists/* \
    && _reproenv_tmppath="$(mktemp -t tmp.XXXXXXXXXX.deb)" \
    && curl -fsSL --retry 5 -o "${_reproenv_tmppath}" http://archive.ubuntu.com/ubuntu/pool/main/g/glibc/multiarch-support_2.27-3ubuntu1.6_amd64.deb \
    && apt-get install --yes -q "${_reproenv_tmppath}" \
    && rm "${_reproenv_tmppath}" \
    && _reproenv_tmppath="$(mktemp -t tmp.XXXXXXXXXX.deb)" \
    && curl -fsSL --retry 5 -o "${_reproenv_tmppath}" https://apt.ligo-wa.caltech.edu:8443/debian/pool/stretch/libxp6/libxp6_1.0.2-2_amd64.deb \
    && apt-get install --yes -q "${_reproenv_tmppath}" \
    && rm "${_reproenv_tmppath}" \
    && _reproenv_tmppath="$(mktemp -t tmp.XXXXXXXXXX.deb)" \
    && curl -fsSL --retry 5 -o "${_reproenv_tmppath}" http://debian.snt.utwente.nl/debian-archive/debian/pool/main/libp/libpng/libpng12-0_1.2.49-1+deb7u2_amd64.deb \
    && apt-get install --yes --allow-downgrades -q "${_reproenv_tmppath}" \
    && rm "${_reproenv_tmppath}" \
    %% apt-get update -qq \
    && apt-get install --yes --quiet --fix-missing \
    && rm -rf /var/lib/apt/lists/* \
    && gsl_path="$(find / -name 'libgsl.so.??' || printf '')" \
    && if [ -n "$gsl_path" ]; then \
         ln -sfv "$gsl_path" "$(dirname $gsl_path)/libgsl.so.0"; \
    fi \
    && ldconfig \
    && mkdir -p /opt/afni \
    && echo "Downloading AFNI ..." \
    && curl -fL https://afni.nimh.nih.gov/pub/dist/tgz/linux_openmp_64.tgz \
    | tar -xz -C /opt/afni --strip-components 1 \
    &&   apt-get update -qq \
    && apt-get install -y -q --no-install-recommends \
           libnlopt-dev \
           r-base \
           r-base-dev \
    && rm -rf /var/lib/apt/lists/*

RUN echo 'export PATH="/usr/lib/afni/bin:$PATH:/opt/afni"' >> ~/.bashrc
RUN echo 'export LD_LIBRARY_PATH="/usr/lib/afni/lib:$LD_LIBRARY_PATH"' >> ~/.bashrc

ADD /compiled /compiled
RUN chmod -R 755 /compiled

RUN mkdir -p /compiled/data

# There were no models in the docker container... add them from my local machine
# models where copied from afni/afni docker container 2019-12-30
# COPY ./afnimodels/* /opt/afni/
# RUN chmod -R 755 /compiled

# the solve script and related files
COPY solve.sh /solve.sh
RUN chmod 755 /solve.sh
COPY default_config.json /opt/default_config.json
ENV PRF_SOLVER afni

