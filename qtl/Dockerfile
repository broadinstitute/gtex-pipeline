# Dockerfile for GTEx xQTL pipeline
FROM ubuntu:16.04
MAINTAINER Francois Aguet

RUN apt-get -qq update && apt-get install -qqy \
    build-essential \
    curl \
    lbzip2 \
    libboost-all-dev \
    libcurl3-dev \
    libgsl-dev \
    openjdk-8-jdk \
    python3 \
    python3-pip \
    r-base-core=3.2.3-4 \
    unzip \
    vim-common \
    wget \
    zlib1g-dev

# R
RUN wget https://github.com/downloads/PMBio/peer/R_peer_source_1.3.tgz && \
    R CMD INSTALL R_peer_source_1.3.tgz && \
    rm R_peer_source_1.3.tgz && \
    echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile && \
    Rscript -e "install.packages(c('argparser'), dependencies=TRUE)" && \
    Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("qvalue"); biocLite("sva"); biocLite("edgeR");'

# standalone R math library
RUN wget https://cran.r-project.org/src/base/R-3/R-3.2.5.tar.gz && tar -xf R-3.2.5.tar.gz && rm R-3.2.5.tar.gz \
    && cd R-3.2.5 && ./configure --with-x=no && cd src/nmath/standalone && make
ENV RMATH /R-3.2.5/src

# htslib
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.4.1/htslib-1.4.1.tar.bz2 && \
    tar -xf htslib-1.4.1.tar.bz2 && rm htslib-1.4.1.tar.bz2 && cd htslib-1.4.1 && make && make install && make clean

# bcftools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/bcftools/releases/download/1.4.1/bcftools-1.4.1.tar.bz2 && \
    tar -xf bcftools-1.4.1.tar.bz2 && rm bcftools-1.4.1.tar.bz2 && cd bcftools-1.4.1 && \
    make && make install && make clean

# FastQTL
RUN cd /opt && \
    wget https://github.com/francois-a/fastqtl/archive/v7.tar.gz && tar -xf v7.tar.gz && rm v7.tar.gz && mv fastqtl-7 fastqtl && \
    cd fastqtl && mkdir obj && sed -i 's/RMATH=/#RMATH/' Makefile && make
ENV PATH /opt/fastqtl/bin:$PATH

# Python
RUN pip3 install --upgrade pip && pip3 install numpy tables pandas scipy matplotlib pyBigWig h5py feather-format rpy2

# PLINK
RUN mkdir /opt/plink && cd /opt/plink && \
    wget https://www.cog-genomics.org/static/bin/plink170815/plink_linux_x86_64.zip && \
    unzip plink_linux_x86_64.zip && rm plink_linux_x86_64.zip
ENV PATH $PATH:/opt/plink

# METASOFT
RUN mkdir /opt/metasoft && cd /opt/metasoft && \
    wget http://genetics.cs.ucla.edu/meta/repository/2.0.1/Metasoft.zip && \
    unzip Metasoft.zip && rm Metasoft.zip

# copy scripts
COPY src src/
