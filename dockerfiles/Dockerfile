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
    python=3.13.5 \
    pip

RUN pip install \
    biopython==1.85 \
    pyfamsa==0.5.3.post1 \
    pyhmmer==0.11.1 \
    pytrimal==0.8.1 \
    numpy==2.3.0

SHELL ["/bin/bash", "-c"]