# Dockerfile for GTEx RNA-seq pipeline. Version 0.1
FROM ubuntu:15.10
MAINTAINER Francois Aguet
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get -qq update && apt-get install -qqy wget unzip python3 python3-pip

# Install Java 8 (required for Picard tools)
RUN apt-get install -qqy openjdk-8-jdk
ENV JAVA_HOME /usr/lib/jvm/java-8-openjdk-amd64

RUN apt-get install -qqy cmake zlib1g-dev libncurses5-dev libboost-dev libboost-program-options-dev

#-----------------------------
# Pipeline components
#-----------------------------

# Install Samtools
RUN wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.3/samtools-1.3.tar.bz2 \
    && tar -xf samtools-1.3.tar.bz2 && cd samtools-1.3 && make && make install && cd .. && rm -r samtools-1.3*

# Install htslib
RUN wget --no-check-certificate https://github.com/samtools/htslib/releases/download/1.3/htslib-1.3.tar.bz2 \
	&& tar -xf htslib-1.3.tar.bz2 && cd htslib-1.3 && make && make install && cd .. && rm -r htslib-1.3*

# Install Bamtools
RUN wget --no-check-certificate https://github.com/pezmaster31/bamtools/archive/v2.4.0.tar.gz \
	&& tar -xf v2.4.0.tar.gz && rm v2.4.0.tar.gz && cd bamtools-2.4.0 && mkdir build && cd build && cmake .. && make && make install
ENV LD_LIBRARY_PATH /usr/local/lib/bamtools:$LD_LIBRARY_PATH

# Install Picard tools
RUN wget --no-check-certificate https://github.com/broadinstitute/picard/releases/download/2.1.1/picard-tools-2.1.1.zip \
    && unzip picard-tools-2.1.1.zip && mv picard-tools-2.1.1 picard-tools && rm picard-tools-2.1.1.zip

# Install STAR v2.4.2a
RUN wget --no-check-certificate https://github.com/alexdobin/STAR/archive/STAR_2.4.2a.tar.gz && tar -xf STAR_2.4.2a.tar.gz \
    && mkdir STAR_2.4.2a && cp /STAR-STAR_2.4.2a/bin/Linux_x86_64_static/* STAR_2.4.2a && rm -r STAR-STAR_2.4.2a && rm STAR_2.4.2a.tar.gz
ENV PATH /STAR_2.4.2a:$PATH

# Install RSEM v1.2.22
RUN wget --no-check-certificate https://github.com/deweylab/RSEM/archive/v1.2.22.tar.gz && tar -xvf v1.2.22.tar.gz && rm v1.2.22.tar.gz && cd RSEM-1.2.22 && make
ENV PATH /RSEM-1.2.22:$PATH

# Install bamsync
COPY bamsync bamsync/
RUN cd bamsync && make && cd ..
ENV PATH /bamsync:$PATH

COPY src src/
