FROM flywheel/matlab-mcr:v95

MAINTAINER Garikoitz Lerma-Usabiaga <glerma@stanford.edu>

# Install jq to parse the JSON config file
RUN apt-get update -qq && apt-get install -y jq

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
