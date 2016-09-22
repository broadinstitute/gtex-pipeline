# Dockerfile for GTEx RNA-seq pipeline. Version 0.3
FROM ubuntu:16.04
MAINTAINER Francois Aguet

RUN apt-get update && apt-get install -y software-properties-common && add-apt-repository -y ppa:openjdk-r/ppa && \
    apt-get update && apt-get install -y \
        build-essential \
        cmake \
        libboost-dev \
        libboost-program-options-dev \
        libncurses5-dev \
        openjdk-7-jdk \
        openjdk-8-jdk \
        python3 \
        python3-pip \
        unzip \
        vim-common \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*


#-----------------------------
# Pipeline components
#-----------------------------

# samtools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 && \
    tar -xf samtools-1.3.1.tar.bz2 && rm samtools-1.3.1.tar.bz2 && cd samtools-1.3.1 && make && make install && make clean

# htslib
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.3.1/htslib-1.3.1.tar.bz2 && \
    tar -xf htslib-1.3.1.tar.bz2 && rm htslib-1.3.1.tar.bz2 && cd htslib-1.3.1 && make && make install && make clean

# bamtools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/pezmaster31/bamtools/archive/v2.4.0.tar.gz && \
    tar -xf v2.4.0.tar.gz && rm v2.4.0.tar.gz && cd bamtools-2.4.0 && mkdir build && cd build && cmake .. && make && make install && make clean
ENV LD_LIBRARY_PATH /usr/local/lib/bamtools:$LD_LIBRARY_PATH

# Picard tools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/broadinstitute/picard/releases/download/2.5.0/picard-tools-2.5.0.zip && \
    unzip picard-tools-2.5.0.zip && mv picard-tools-2.5.0 picard-tools && rm picard-tools-2.5.0.zip

# RSEM v1.2.22
RUN cd /opt && \
    wget --no-check-certificate https://github.com/deweylab/RSEM/archive/v1.2.22.tar.gz && \
    tar -xvf v1.2.22.tar.gz && rm v1.2.22.tar.gz && cd RSEM-1.2.22 && make
ENV PATH /opt/RSEM-1.2.22:$PATH

# STAR v2.4.2a
RUN cd /opt && \
    wget --no-check-certificate https://github.com/alexdobin/STAR/archive/STAR_2.4.2a.tar.gz && \
    tar -xf STAR_2.4.2a.tar.gz && rm STAR_2.4.2a.tar.gz && mv STAR-STAR_2.4.2a STAR_2.4.2a && \
    make STAR -C STAR_2.4.2a/source && make STARlong -C STAR_2.4.2a/source && \
    mv STAR_2.4.2a/source/STAR* STAR_2.4.2a/bin/Linux_x86_64/
ENV PATH /opt/STAR_2.4.2a/bin/Linux_x86_64:$PATH

# RNA-SeQC
RUN cd /opt && \
    wget --no-check-certificate https://github.com/francois-a/rnaseqc/releases/download/v1.1.9/RNA-SeQC_1.1.9.zip && \
    unzip RNA-SeQC_1.1.9.zip -d RNA-SeQC_1.1.9 && rm RNA-SeQC_1.1.9.zip

# bamsync
COPY bamsync /opt/bamsync
RUN cd /opt/bamsync && make
ENV PATH /opt/bamsync:$PATH

# python modules
RUN pip3 install pandas

# scripts
COPY src src/

# clean up
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/
