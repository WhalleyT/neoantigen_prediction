FROM ubuntu:16.04

LABEL authors="whalleyt@cardiff.ac.uk"
LABEL description="Docker image containing all requirements for the neoantigen pipeline"

#create bash
RUN rm /bin/sh && ln -s /bin/bash /bin/sh

#basis stuff
RUN apt-get update && apt-get -y upgrade && \
    apt-get install -y autoconf curl default-jre apt-utils \
    git gsl-bin libcurl4-openssl-dev libgsl0-dev \ 
    libssl-dev python python-pip python3 python3-pip \
    automake make libconfig-inifiles-perl pp-popularity-contest \ 
    wget tar g++ 

#install anaconda
RUN wget http://repo.continuum.io/archive/Anaconda3-5.0.1-Linux-x86_64.sh
RUN bash Anaconda3-5.0.1-Linux-x86_64.sh -b
RUN rm Anaconda3-5.0.1-Linux-x86_64.sh
ENV PATH=/root/anaconda3/bin:$PATH


#updating Anaconda packages
RUN ~/anaconda3/bin/conda update --all


COPY environment.yml /
RUN ~/anaconda3/bin/conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/neoantigen/bin:$PATH
