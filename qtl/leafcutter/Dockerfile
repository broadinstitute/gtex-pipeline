# Dockerfile for LeafCutter
FROM ubuntu:22.04
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

# R
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile && \
    Rscript -e "install.packages(c('argparser', 'devtools', 'foreach'), dependencies=TRUE)" && \
    Rscript -e 'BiocManager::install(c("qvalue", "sva", "edgeR", "dplyr"));'

# htslib
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.19.1/htslib-1.19.1.tar.bz2 && \
    tar -xf htslib-1.19.1.tar.bz2 && rm htslib-1.19.1.tar.bz2 && cd htslib-1.19.1 && \
    ./configure --enable-libcurl --enable-s3 --enable-plugins --enable-gcs && \
    make && make install && make clean

# samtools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.19.2/samtools-1.19.2.tar.bz2 && \
    tar -xf samtools-1.19.2.tar.bz2 && rm samtools-1.19.2.tar.bz2 && cd samtools-1.19.2 && \
    ./configure --with-htslib=/opt/htslib-1.19.1 && make && make install && make clean

# bcftools
RUN cd /opt && \
    wget --no-check-certificate https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2 && \
    tar -xf bcftools-1.19.tar.bz2 && rm bcftools-1.19.tar.bz2 && cd bcftools-1.19 && \
    ./configure --with-htslib=/opt/htslib-1.19.1 && make && make install && make clean

# Python 3
RUN pip3 install --upgrade pip setuptools
RUN pip3 install numpy pyarrow pandas scipy scikit-learn matplotlib qtl

# LeafCutter
RUN cd /opt && \
    wget https://github.com/francois-a/leafcutter/archive/2488118d377baff3354dab85de1f31b03a813c92.tar.gz && \
    tar -xf 2488118d377baff3354dab85de1f31b03a813c92.tar.gz && \
    rm 2488118d377baff3354dab85de1f31b03a813c92.tar.gz && \
    ln -s leafcutter-2488118d377baff3354dab85de1f31b03a813c92 leafcutter

# regtools
RUN cd /opt && \
    wget https://github.com/griffithlab/regtools/archive/refs/tags/0.5.2.zip && \
    unzip 0.5.2.zip && rm 0.5.2.zip && cd regtools-0.5.2 && mkdir build && cd build && \
    cmake .. && make && mv regtools .. && make clean
ENV PATH /opt/regtools-0.5.2:$PATH

# copy scripts
COPY src src/
