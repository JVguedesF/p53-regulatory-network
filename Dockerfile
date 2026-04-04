FROM mambaorg/micromamba:1.5.8-debian12-slim

LABEL maintainer="João Vitor Guedes"
LABEL description="p53-regulatory-network pipeline"

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    awscli \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*

WORKDIR /workspace

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml .

RUN micromamba install -y -n base -f environment.yml \
 && micromamba clean -afy

ENV PATH=/opt/conda/bin:$PATH

COPY src/ ./src/

USER $MAMBA_USER

CMD ["/bin/bash"]