#!/bin/bash

## This setup is for RNA-seq analysis

# Create conda environment
conda create -y -n mUC_ICB_RNA python=3.12 r-base

# Activate environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mUC_ICB_RNA

# Install Python packages
conda install -y -c conda-forge \
    pandas \
    numpy \
    matplotlib \
    seaborn \
    scipy \
    statsmodels \
    scanpy \
    jupyter \
    ipykernel

# Install R packages from CRAN and Bioconductor
conda install -y -c conda-forge -c bioconda \
    r-data.table \
    r-ggplot2 \
    r-rcolorbrewer \
    r-ggrepel \
    bioconductor-deseq2 \
    bioconductor-edger \
    bioconductor-biocparallel

# Register Jupyter kernel
python -m ipykernel install --user --name=mUC_ICB_RNA

echo "Environment mUC_ICB_RNA setup complete."