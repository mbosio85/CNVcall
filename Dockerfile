FROM nfcore/base:1.7
LABEL authors="Mattia Bosio" \
      description="Docker image containing all requirements for nf-core/cnvcall pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-cnvcall-1.0dev/bin:$PATH
