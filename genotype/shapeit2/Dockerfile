FROM ubuntu:20.04
MAINTAINER Francois Aguet

RUN apt-get update && apt-get install -y software-properties-common && \
    apt-get update && apt-get install -y \
        build-essential \
        cmake \
        curl \
        libboost-all-dev \
        libbz2-dev \
        libcurl3-dev \
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        python3 \
        python3-pip \
        sudo \
        unzip \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/

# htslib
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.11/htslib-1.11.tar.bz2 && \
    tar xf htslib-1.11.tar.bz2 && rm htslib-1.11.tar.bz2 && cd htslib-1.11 && \
    ./configure --enable-libcurl --enable-s3 --enable-plugins --enable-gcs && \
    make && make install && make clean

# samtools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 && \
    tar -xf samtools-1.11.tar.bz2 && rm samtools-1.11.tar.bz2 && cd samtools-1.11 && \
    ./configure --with-htslib=/opt/htslib-1.11 && make && make install && make clean

# bcftools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/bcftools/releases/download/1.11/bcftools-1.11.tar.bz2 && \
    tar -xf bcftools-1.11.tar.bz2 && rm bcftools-1.11.tar.bz2 && cd bcftools-1.11 && \
    ./configure --with-htslib=/opt/htslib-1.11 && make && make install && make clean

# python3
RUN pip3 install --upgrade pip && pip3 install pandas numpy

# SHAPEIT2
RUN cd /opt && \
    wget https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r904.glibcv2.17.linux.tar.gz && \
    tar -xf shapeit.v2.r904.glibcv2.17.linux.tar.gz && rm shapeit.v2.r904.glibcv2.17.linux.tar.gz
ENV PATH /opt/shapeit.v2.904.3.10.0-693.11.6.el7.x86_64/bin:$PATH

# extractPIRs
RUN cd /opt && \
   wget https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/files/extractPIRs.v1.r68.x86_64.tgz && \
   tar xf extractPIRs.v1.r68.x86_64.tgz && rm extractPIRs.v1.r68.x86_64.tgz
ENV PATH /opt/extractPIRs.v1.r68.x86_64:$PATH

# scripts
COPY src src/
ENV PATH /src:$PATH
