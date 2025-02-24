# Use an official Miniconda3 base image
FROM continuumio/miniconda3:latest

# Set the working directory
WORKDIR /app

# Configure Conda to only use conda-forge and bioconda channels
RUN conda config --add channels conda-forge \
    && conda config --add channels bioconda \
    && conda config --set channel_priority strict
    
RUN conda install -y \
    pandas \
    numpy \
    bioconda::mafft=7.520 \
    bioconda::hmmer=3.4.0 \
    bioconda::easel=0.49

# Create a virtual environment (optional)
RUN conda create -n myenv pandas numpy bioconda::mafft=7.520 bioconda::hmmer=3.4.0 bioconda::easel=0.49

# Activate the environment and set it as the default for the Docker container
SHELL ["conda", "run", "-n", "myenv", "/bin/bash", "-c"]

# Install Biopython using pip (since it's not available in the conda channels)
RUN pip install biopython

# Set environment variable to use conda environment by default
ENV PATH /opt/conda/envs/myenv/bin:$PATH
