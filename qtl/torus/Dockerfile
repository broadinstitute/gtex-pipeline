# Dockerfile for Torus (https://github.com/xqwen/torus)
FROM ubuntu:16.04
MAINTAINER Francois Aguet

RUN apt-get -qq update && apt-get install -qqy \
        build-essential \
        curl \
        git-all \
        lbzip2 \
        libboost-all-dev \
        libcurl3-dev \
        libgsl-dev \
        python3 \
        python3-pip \
        unzip \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/


# Torus
RUN cd /opt && \
    git clone https://github.com/xqwen/torus.git && \
    cd torus/src && make && mkdir ../bin && cp torus ../bin/ && make clean
ENV PATH /opt/torus/bin:$PATH

# copy scripts
COPY src src/
