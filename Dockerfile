FROM continuumio/miniconda3:latest

LABEL maintainer="João Vitor"
LABEL description="Pipeline p53-regulatory-network"

RUN apt-get update && apt-get install -y procps

WORKDIR /workspace

COPY environment.yml .

RUN conda env create -f environment.yml && conda clean -a

ENV PATH=/opt/conda/envs/p53-regulatory-network/bin:$PATH

CMD ["/bin/bash"]