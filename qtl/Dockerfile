# Dockerfile for GTEx QTL pipeline
FROM ubuntu:18.04
MAINTAINER Francois Aguet
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get install -y software-properties-common && \
    apt-get update && apt-get install -y \
        build-essential \
        curl \
        lbzip2 \
        libboost-all-dev \
        libcurl3-dev \
        libgsl-dev \
        libhdf5-serial-dev \
        openjdk-8-jdk \
        python3 \
        python3-pip \
        r-base-core \
        unzip \
        vim-common \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/


# workaround for PEER, see https://github.com/mz2/peer/issues/4
RUN apt-get update \
    && apt-get install -y \
        gcc-5 \
        g++-5 \
        gfortran-5 \
        cmake \
    && rm -rf /var/lib/apt/lists/* \
    && update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-5 50 --slave /usr/bin/g++ g++ /usr/bin/g++-5

# R
RUN wget https://github.com/downloads/PMBio/peer/R_peer_source_1.3.tgz && \
    R CMD INSTALL R_peer_source_1.3.tgz && \
    rm R_peer_source_1.3.tgz && \
    echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile && \
    Rscript -e "install.packages(c('argparser'), dependencies=TRUE)" && \
    Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("qvalue"); biocLite("sva"); biocLite("edgeR");'

# htslib
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.15/htslib-1.15.tar.bz2 && \
    tar -xf htslib-1.15.tar.bz2 && rm htslib-1.15.tar.bz2 && cd htslib-1.15 && make && make install && make clean

# bcftools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/bcftools/releases/download/1.15/bcftools-1.15.tar.bz2 && \
    tar -xf bcftools-1.15.tar.bz2 && rm bcftools-1.15.tar.bz2 && cd bcftools-1.15 && \
    ./configure --with-htslib=/opt/htslib-1.15 && make && make install && make clean

# samtools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.15/samtools-1.15.tar.bz2 && \
    tar -xf samtools-1.15.tar.bz2 && rm samtools-1.15.tar.bz2 && cd samtools-1.15 && \
    ./configure --with-htslib=/opt/htslib-1.15 && make && make install && make clean

# PLINK 1.9
RUN mkdir /opt/plink && cd /opt/plink && \
    wget --no-check-certificate https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20220305.zip && \
    unzip plink_linux_x86_64_20220305.zip && rm plink_linux_x86_64_20220305.zip
ENV PATH $PATH:/opt/plink

# METASOFT
RUN mkdir /opt/metasoft && cd /opt/metasoft && \
    wget http://genetics.cs.ucla.edu/meta/repository/2.0.1/Metasoft.zip && \
    unzip Metasoft.zip && rm Metasoft.zip

# Python
RUN pip3 install --upgrade pip setuptools
RUN pip3 install numpy tables pandas scipy matplotlib h5py pysam statsmodels scikits.bootstrap qtl
# numpy dependencies:
RUN pip3 install pyBigWig

# aFC
RUN cd /opt && \
    wget https://github.com/francois-a/aFC/archive/2189fbf403b3d1ced54da21421d00b2d4bf44310.tar.gz && \
    tar -xf 2189fbf403b3d1ced54da21421d00b2d4bf44310.tar.gz && mv aFC-2189fbf403b3d1ced54da21421d00b2d4bf44310 aFC && \
    rm 2189fbf403b3d1ced54da21421d00b2d4bf44310.tar.gz

# copy scripts
COPY src src/
