#!/bin/bash

set -euo pipefail

## This setup is for RNA-seq analysis

# Create conda environment
conda create -y -n mUC_ICB_RNA -c conda-forge \
    python=3.12.13 \
    r-base=4.4.3

# Activate environment
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate mUC_ICB_RNA

# Install Python packages
conda install -y -c conda-forge \
    pandas=2.3.3 \
    numpy=2.4.6 \
    matplotlib=3.10.9 \
    seaborn=0.13.2 \
    scipy=1.17.1 \
    statsmodels=0.14.6 \
    scanpy=1.12.1 \
    jupyter=1.1.1 \
    ipykernel=7.2.0 \
    adjustText

# Install R packages from CRAN and Bioconductor
conda install -y -c conda-forge -c bioconda \
    r-data.table=1.17.8 \
    r-ggplot2=4.0.3 \
    r-rcolorbrewer=1.1_3 \
    r-ggrepel=0.9.8 \
    bioconductor-deseq2=1.46.0 \
    bioconductor-edger=4.4.0 \
    bioconductor-biocparallel=1.40.0

# Register Jupyter kernel
python -m ipykernel install --user --name=mUC_ICB_RNA

echo "Environment mUC_ICB_RNA setup complete."
