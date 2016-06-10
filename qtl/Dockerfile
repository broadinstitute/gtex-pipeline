# Dockerfile for GTEx xQTL pipeline. Version 0.1
FROM ubuntu:16.04
MAINTAINER Francois Aguet
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get -qq update && apt-get install -qqy wget unzip build-essential vim-common

# Python
RUN apt-get install -qqy python3 python3-pip
RUN pip3 install --upgrade pip numpy pandas scipy

# R
RUN apt-get install -qqy r-base-core=3.2.3-4
RUN wget https://github.com/downloads/PMBio/peer/R_peer_source_1.3.tgz && R CMD INSTALL R_peer_source_1.3.tgz && rm R_peer_source_1.3.tgz
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('argparser')"
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R");biocLite("qvalue");'

# standalone R math library
RUN wget https://cran.r-project.org/src/base/R-3/R-3.2.5.tar.gz && tar -xf R-3.2.5.tar.gz && rm R-3.2.5.tar.gz \
    && cd R-3.2.5 && ./configure --with-x=no && cd src/nmath/standalone && make
ENV RMATH /R-3.2.5/src

# Install htslib
RUN wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.3.1/htslib-1.3.1.tar.bz2 \
	&& tar -xf htslib-1.3.1.tar.bz2 && cd htslib-1.3.1 && make && make install && cd .. && rm -r htslib-1.3.1*

# Install FastQTL
RUN apt-get install -qqy lbzip2 zlib1g-dev libboost-dev libboost-iostreams-dev libboost-thread-dev libboost-program-options-dev libgsl-dev
RUN wget https://github.com/francois-a/fastqtl/archive/master.zip && unzip master.zip && rm master.zip && mv fastqtl-master fastqtl && \
    cd fastqtl && mkdir obj && sed -i 's/RMATH=/#RMATH/' Makefile && make
ENV PATH /fastqtl/bin:$PATH

# copy scripts
COPY src src/
