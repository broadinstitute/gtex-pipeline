# Dockerfile for GTEx RNA-seq pipeline
FROM ubuntu:18.04
MAINTAINER Francois Aguet

RUN apt-get update && apt-get install -y software-properties-common && \
    apt-get update && apt-get install -y \
        build-essential \
        cmake \
        curl \
        git \
        libboost-all-dev \
        libbz2-dev \
        libcurl3-dev \
        libhdf5-serial-dev \
        liblzma-dev \
        libncurses5-dev \
        libssl-dev \
        openjdk-8-jdk \
        python3 \
        python3-pip \
        unzip \
        vim-common \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/


#-----------------------------
# Pipeline components
#-----------------------------

# htslib
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    tar -xf htslib-1.9.tar.bz2 && rm htslib-1.9.tar.bz2 && cd htslib-1.9 && \
    ./configure --enable-libcurl --enable-s3 --enable-plugins --enable-gcs && \
    make && make install && make clean

# samtools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    tar -xf samtools-1.9.tar.bz2 && rm samtools-1.9.tar.bz2 && cd samtools-1.9 && \
    ./configure --with-htslib=/opt/htslib-1.9 && make && make install && make clean

# STAR v2.6.1d
RUN cd /opt && \
    wget --no-check-certificate https://github.com/alexdobin/STAR/archive/2.6.1d.tar.gz && \
    tar -xf 2.6.1d.tar.gz && rm 2.6.1d.tar.gz && \
    make STAR -C STAR-2.6.1d/source && make STARlong -C STAR-2.6.1d/source && \
    mv STAR-2.6.1d/source/STAR* STAR-2.6.1d/bin/Linux_x86_64/
ENV PATH /opt/STAR-2.6.1d/bin/Linux_x86_64:$PATH

# RSEM v1.3.1
RUN cd /opt && \
    wget --no-check-certificate https://github.com/deweylab/RSEM/archive/v1.3.1.tar.gz && \
    tar -xvf v1.3.1.tar.gz && rm v1.3.1.tar.gz && cd RSEM-1.3.1 && make
ENV PATH /opt/RSEM-1.3.1:$PATH

# bamtools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/pezmaster31/bamtools/archive/v2.4.1.tar.gz && \
    tar -xf v2.4.1.tar.gz && rm v2.4.1.tar.gz && cd bamtools-2.4.1 && mkdir build && cd build && cmake .. && make && make install && make clean
ENV LD_LIBRARY_PATH /usr/local/lib/bamtools:$LD_LIBRARY_PATH

# bamsync
COPY bamsync /opt/bamsync
RUN cd /opt/bamsync && make
ENV PATH /opt/bamsync:$PATH

# Picard tools
RUN mkdir /opt/picard-tools && \
    wget --no-check-certificate -P /opt/picard-tools/ https://github.com/broadinstitute/picard/releases/download/2.21.4/picard.jar

# kallisto
RUN cd /opt && \
    wget https://github.com/pachterlab/kallisto/releases/download/v0.46.1/kallisto_linux-v0.46.1.tar.gz && \
    tar -xf kallisto_linux-v0.46.1.tar.gz && rm kallisto_linux-v0.46.1.tar.gz
ENV PATH $PATH:/opt/kallisto_linux-v0.46.1

# bedtools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/arq5x/bedtools2/releases/download/v2.29.0/bedtools-2.29.0.tar.gz && \
    tar -xf bedtools-2.29.0.tar.gz && rm bedtools-2.29.0.tar.gz && \
    cd bedtools2 && make && make install && make clean

# UCSC tools
RUN mkdir /opt/ucsc && \
    wget --no-check-certificate -P /opt/ucsc/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph && \
    wget --no-check-certificate -P /opt/ucsc/ http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig && \
    chmod 755 /opt/ucsc/*
ENV PATH /opt/ucsc:$PATH

# python modules
RUN pip3 install --upgrade pip setuptools
RUN pip3 install tables numpy pandas
# numpy dependencies:
RUN pip3 install pyBigWig

# RNA-SeQC
RUN mkdir /opt/rnaseqc && cd /opt/rnaseqc && \
    wget https://github.com/getzlab/rnaseqc/releases/download/v2.3.6/rnaseqc.v2.3.6.linux.gz && \
    gunzip rnaseqc.v2.3.6.linux.gz && mv rnaseqc.v2.3.6.linux rnaseqc && chmod 775 rnaseqc
ENV PATH /opt/rnaseqc:$PATH

# gcloud
RUN export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update -y && apt-get install google-cloud-sdk -y

# scripts
COPY src src/
