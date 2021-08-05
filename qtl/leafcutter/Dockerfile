# Dockerfile for LeafCutter
FROM ubuntu:18.04
MAINTAINER Francois Aguet
ENV DEBIAN_FRONTEND noninteractive

RUN apt-get update && apt-get install -y software-properties-common && \
    apt-get update && apt-get install -y \
        build-essential \
        cmake \
        curl \
        lbzip2 \
        libboost-all-dev \
        libcurl3-dev \
        libgsl-dev \
        libssl-dev \
        libxml2-dev \
        openjdk-8-jdk \
        python \
        python-pip \
        r-base-core \
        unzip \
        vim-common \
        wget \
        zlib1g-dev \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get clean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/

# R
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile && \
    Rscript -e "install.packages(c('argparser', 'devtools'), dependencies=TRUE)" && \
    Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("qvalue"); biocLite("sva"); biocLite("edgeR");'

# htslib
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2 && \
    tar -xf htslib-1.13.tar.bz2 && rm htslib-1.13.tar.bz2 && cd htslib-1.13 && \
    ./configure --enable-libcurl --enable-s3 --enable-plugins --enable-gcs && \
    make && make install && make clean

# samtools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.13/samtools-1.13.tar.bz2 && \
    tar -xf samtools-1.13.tar.bz2 && rm samtools-1.13.tar.bz2 && cd samtools-1.13 && \
    ./configure --with-htslib=/opt/htslib-1.13 && make && make install && make clean

# bcftools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/bcftools/releases/download/1.13/bcftools-1.13.tar.bz2 && \
    tar -xf bcftools-1.13.tar.bz2 && rm bcftools-1.13.tar.bz2 && cd bcftools-1.13 && \
    ./configure --with-htslib=/opt/htslib-1.13 && make && make install && make clean

# Python
RUN pip install --upgrade pip && pip install numpy tables pandas scipy scikit-learn matplotlib

# LeafCutter
RUN cd /opt && \
    wget https://github.com/davidaknowles/leafcutter/archive/3f875f1ac05a0f0413884c0d016b3c7ae17713c9.tar.gz && \
    tar -xf 3f875f1ac05a0f0413884c0d016b3c7ae17713c9.tar.gz && \
    rm 3f875f1ac05a0f0413884c0d016b3c7ae17713c9.tar.gz && \
    ln -s leafcutter-3f875f1ac05a0f0413884c0d016b3c7ae17713c9 leafcutter

# regtools
RUN cd /opt && \
    wget https://github.com/griffithlab/regtools/archive/refs/tags/0.5.2.zip && \
    unzip 0.5.2.zip && rm 0.5.2.zip && cd regtools-0.5.2 && mkdir build && cd build && \
    cmake .. && make && mv regtools .. && make clean
ENV PATH /opt/regtools-0.5.2:$PATH

# copy scripts
COPY src src/
