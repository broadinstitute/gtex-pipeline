FROM ubuntu:16.04
MAINTAINER Xiao Li

RUN apt-get -qq update && apt-get install -qqy \
    apt-utils \
    build-essential \
    libboost-dev \
    libboost-iostreams-dev \
    libboost-program-options-dev \
    libboost-thread-dev \
    lbzip2 \
    libgsl-dev \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    unzip \
    vim-common \
    python3 \
    python3-pip \
    wget \
    less \
    make

# HTSLIB
RUN cd /opt && \
    wget https://github.com/samtools/htslib/releases/download/1.6/htslib-1.6.tar.bz2 && \
    tar xf htslib-1.6.tar.bz2 && rm htslib-1.6.tar.bz2 && cd htslib-1.6 && \
    ./configure && make && make install

# SAMTOOLS
RUN cd /opt && \
    wget https://github.com/samtools/samtools/releases/download/1.6/samtools-1.6.tar.bz2 && \
    tar xjf samtools-1.6.tar.bz2 && rm samtools-1.6.tar.bz2 && cd samtools-1.6 && ./configure --without-curses && make && make install

# BCFTOOLS
RUN cd /opt && \
    wget https://github.com/samtools/bcftools/releases/download/1.6/bcftools-1.6.tar.bz2 && \
    tar xjf bcftools-1.6.tar.bz2 && rm bcftools-1.6.tar.bz2 && cd bcftools-1.6 && \
    ./configure && make && make install

# python3
RUN pip3 install --upgrade pip && pip3 install argparse pandas numpy

# SHAPEIT2
RUN cd /opt && \
    wget https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz && \
    mkdir shapeit_v2_r837 && tar -xf shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz -C shapeit_v2_r837 && rm shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz
ENV PATH /opt/shapeit_v2_r837/bin:$PATH

# extractPIRs
RUN cd /opt && \
   wget https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/files/extractPIRs.v1.r68.x86_64.tgz && \
   tar zxvf extractPIRs.v1.r68.x86_64.tgz && rm extractPIRs.v1.r68.x86_64.tgz
ENV PATH /opt/extractPIRs.v1.r68.x86_64:$PATH

# clean up
RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get autoclean && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}/

# scripts
COPY src src/
ENV PATH /src:$PATH
