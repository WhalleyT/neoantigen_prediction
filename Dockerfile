FROM continuumio/miniconda3

LABEL authors="whalleyt@cardiff.ac.uk"
LABEL description="Docker image containing all requirements for the neoantigen pipeline"

ADD environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml

# Pull the environment name out of the environment.yml
RUN echo "source activate neoantigen" > ~/.bashrc
ENV PATH /opt/conda/envs/$(head -1 /tmp/environment.yml | cut -d' ' -f2)/bin:$PATH
RUN . activate neoantigen