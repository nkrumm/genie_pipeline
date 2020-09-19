FROM ubuntu:latest
MAINTAINER Nik Krumm <nkrumm@uw.edu>

RUN rm /bin/sh && ln -s /bin/bash /bin/sh

ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=US/Pacific

## Install python dependencies
RUN apt-get update && \
    apt-get install --yes \
	python3 \
	python3-pip \
	r-base

RUN Rscript -e 'install.packages("BiocManager"); BiocManager::install("CGHcall",  suppressUpdates=TRUE)'

ADD requirements.txt requirements.txt
RUN pip3 install -r requirements.txt
