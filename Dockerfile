
FROM nfcore/base
LABEL authors="Tom Whalley: whalleyt@cardiff.ac.uk" \
      description="Docker image containing all requirements for the neoantigen prediction pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/neoantigen/bin:$PATH
RUN snpEff download GRCh38.86
