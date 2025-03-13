# Use an official Miniconda3 base image
FROM continuumio/miniconda3:latest

# Set the working directory
WORKDIR /app

RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda \
    && conda config --set channel_priority strict

# Install build tools
RUN apt-get update && apt-get install -y build-essential

# Install Python dependencies using pip
RUN pip install pandas pyfastx pyfamsa biopython \
    && pip install git+https://github.com/althonos/pyhmmer.git

# Set environment variable to ensure Conda is properly configured
ENV PATH /opt/conda/bin:$PATH

# Set default shell
SHELL ["/bin/bash", "-c"]