# Dockerfile for GTEx RNA-seq pipeline dependencies
FROM ubuntu:18.04
MAINTAINER Aaron Graubert

RUN apt-get update && apt-get install -y \
        software-properties-common \
        bcftools \
        bedtools \
        build-essential \
        cmake \
        curl \
        git \
        libboost-all-dev \
        libbz2-dev \
        libcurl3-dev \
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        python2.7 \
        python-pip \
        python3-pip \
        unzip \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*


#-----------------------------
# Pipeline components
#-----------------------------

# Samtools (phaser)
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.5/samtools-1.5.tar.bz2 && \
    tar -xjf samtools-1.5.tar.bz2 && rm samtools-1.5.tar.bz2 && cd samtools-1.5 && \
    ./configure && make && make install && cd htslib-1.5 && make && make install

# Get Cython (phaser)
RUN python2.7 -m pip install cython pandas scipy pysam intervaltree

# exon map dependencies
# phASER
# f15e83a : Latest commit at time of writing
#RUN cd /opt && \
    #git clone https://github.com/secastel/phaser.git && cd phaser && \
    #git checkout f15e83a && cd phaser && python2.7 setup.py build_ext --inplace
RUN cd /opt && \
    git clone https://github.com/secastel/phaser.git && cd phaser && \
    git fetch origin pull/36/head:subprocess && git checkout subprocess && \
    cd phaser && python2.7 setup.py build_ext --inplace

COPY wrapper.py /opt/phaser/

# clean up
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/
