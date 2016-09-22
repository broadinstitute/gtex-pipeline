# Dockerfile for GTEx xQTL pipeline. Version 0.1
FROM ubuntu:16.04
MAINTAINER Francois Aguet

RUN apt-get -qq update && apt-get install -qqy \
    build-essential \
    lbzip2 \
    libboost-dev \
    libboost-iostreams-dev \
    libboost-program-options-dev \
    libboost-thread-dev \
    libgsl-dev \
    python3 \
    python3-pip \
    r-base-core=3.2.3-4 \
    unzip \
    vim-common \
    wget \
    zlib1g-dev

# Python
RUN pip3 install --upgrade pip numpy pandas scipy

# R
RUN wget https://github.com/downloads/PMBio/peer/R_peer_source_1.3.tgz && \
    R CMD INSTALL R_peer_source_1.3.tgz && \
    rm R_peer_source_1.3.tgz && \
    echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile && \
    Rscript -e "install.packages('argparser')" && \
    Rscript -e 'source("http://bioconductor.org/biocLite.R");biocLite("qvalue");'

# standalone R math library
RUN wget https://cran.r-project.org/src/base/R-3/R-3.2.5.tar.gz && tar -xf R-3.2.5.tar.gz && rm R-3.2.5.tar.gz \
    && cd R-3.2.5 && ./configure --with-x=no && cd src/nmath/standalone && make
ENV RMATH /R-3.2.5/src

# htslib
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.3.1/htslib-1.3.1.tar.bz2 && \
    tar -xf htslib-1.3.1.tar.bz2 && rm htslib-1.3.1.tar.bz2 && cd htslib-1.3.1 && make && make install && make clean

# Install FastQTL
RUN cd /opt && \
    wget https://github.com/francois-a/fastqtl/archive/v6p.zip && unzip v6p.zip && rm v6p.zip && mv fastqtl-6p fastqtl && \
    cd fastqtl && mkdir obj && sed -i 's/RMATH=/#RMATH/' Makefile && make
ENV PATH /opt/fastqtl/bin:$PATH

# copy scripts
COPY src src/
