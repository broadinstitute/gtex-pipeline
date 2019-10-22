# Dockerfile for LeafCutter
FROM ubuntu:17.04
MAINTAINER Francois Aguet

RUN apt-get -qq update && apt-get install -qqy \
    build-essential \
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
    r-base-core=3.3.2-1 \
    unzip \
    vim-common \
    wget \
    zlib1g-dev

# R
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile && \
    Rscript -e "install.packages(c('argparser', 'devtools'), dependencies=TRUE)" && \
    Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("qvalue"); biocLite("sva"); biocLite("edgeR");'

# htslib
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.6/htslib-1.6.tar.bz2 && \
    tar -xf htslib-1.6.tar.bz2 && rm htslib-1.6.tar.bz2 && cd htslib-1.6 && make && make install && make clean

# bcftools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/bcftools/releases/download/1.6/bcftools-1.6.tar.bz2 && \
    tar -xf bcftools-1.6.tar.bz2 && rm bcftools-1.6.tar.bz2 && cd bcftools-1.6 && \
    make && make install && make clean

# samtools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2 && \
    tar -xf samtools-1.6.tar.bz2 && rm samtools-1.6.tar.bz2 && cd samtools-1.6 && \
    ./configure --with-htslib=/opt/htslib-1.6 && make && make install && make clean

# Python
RUN pip install --upgrade pip && pip install numpy tables pandas scipy scikit-learn matplotlib feather-format

# LeafCutter
RUN cd /opt && \
    wget https://github.com/davidaknowles/leafcutter/archive/v0.2.7.tar.gz && \
    tar -xf v0.2.7.tar.gz && rm v0.2.7.tar.gz && ln -s leafcutter-0.2.7 leafcutter

# clean up
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/

# copy scripts
COPY src src/
