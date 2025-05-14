FROM mambaorg/micromamba:2.1.1

USER root
WORKDIR /app

ENV PATH="/opt/conda/bin:$PATH"
ENV CONDA_ENV=base

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    libeigen3-dev \
    wget \
    unzip \
    git

RUN micromamba install -y -n base -c conda-forge \
    python=3.12.9 \
    pip \
    pandas=2.2.3 \
    biopython=1.85 \
    numpy=2.2.4

RUN pip install \
    pyfastx==2.2.0 \
    pyfamsa==0.5.3.post1 \
    pyhmmer==0.11.1 \
    pytrimal==0.7.0

SHELL ["/bin/bash", "-c"]