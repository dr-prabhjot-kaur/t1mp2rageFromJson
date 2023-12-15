#FROM ubuntu:jammy #Prabhjot
FROM ubuntu:22.04

ARG  DEBIAN_FRONTEND=noninteractive
ENV  TZ="America/New_York"

#RUN echo 'Acquire::http::Proxy "http://proxy.tch.harvard.edu:3128"; Acquire::https::Proxy "https://proxy.tch.harvard.edu:3128";' > /etc/apt/proxy.conf
#RUN cat /etc/apt/proxy.conf

RUN apt-get update 
RUN apt-get install libxt-dev -y
RUN apt-get install --no-install-recommends --no-install-suggests --yes \
  git vim unzip curl

RUN apt-get update && \ 
    apt-get install --no-install-recommends --no-install-suggests --yes apt-transport-https ca-certificates curl software-properties-common

# Let's pull down Update 5 of the Runtime for Matlab R2022a

RUN mkdir /mcr-install && \
    mkdir /opt/mcr && \
    cd /mcr-install && \
    curl -LJO \
        https://ssd.mathworks.com/supportfiles/downloads/R2023a/Release/5/deployment_files/installer/complete/glnxa64/MATLAB_Runtime_R2023a_Update_5_glnxa64.zip && \
    unzip -q MATLAB_Runtime_R2023a_Update_5_glnxa64.zip && \
    rm -f MATLAB_Runtime_R2023a_Update_5_glnxa64.zip && \
    ./install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent && \
    cd / && \
    rm -rf mcr-install

# Cleanup files not required after installation
RUN apt-get remove git -y \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Build for jammy jellyfish
# Use network=host option to provide proxy for apt and pip3
# DOCKER_BUILDKIT=1 docker build --network=host -t crl:t1mp2rage -f Dockerfile .    

# Prabhjot added below- Yet to know
# ALso changed line 13 here for r2022b. Later this commented line 1 & 2 and inserted line 3
COPY . /app
WORKDIR /app
ENV LD_LIBRARY_PATH /opt/mcr/R2022b/runtime/glnxa64:/opt/mcr/R2022b/bin/glnxa64:/opt/mcr/R2022b/sys/os/glnxa64:/opt/mcr/R2022b/extern/bin/glnxa64

# Run the application
#CMD sh run_EstimateT1fromMP2RAGE.sh
