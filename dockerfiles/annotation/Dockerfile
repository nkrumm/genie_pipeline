FROM ubuntu:latest
MAINTAINER Nik Krumm <nkrumm@uw.edu>

RUN rm /bin/sh && ln -s /bin/bash /bin/sh

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=US/Pacific

RUN apt-get update

RUN apt-get install --yes \
 build-essential \
 apt-utils \
 default-jre \
 default-jdk \
 wget \
 unzip

# snpEff and SnpSift.jar in /snpEff
# Get latest source from releases
RUN wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
RUN unzip snpEff_latest_core.zip

# Get and install vcfanno
# PATH to vcfanno executable set

RUN wget https://github.com/brentp/vcfanno/releases/download/v0.3.2/vcfanno_linux64
RUN chmod +x vcfanno_linux64

ENV PATH /snpEff:/:$PATH


## Install python dependencies
RUN apt-get install --yes \
	python3 \
	python3-pip

ADD requirements.txt requirements.txt
RUN pip3 install -r requirements.txt
