#!/bin/bash

# This shell script runs the UMAP analysis on the aggregated feature set profiles.
# Steps:
# 1. Concatenate all profiles into a single file using the `concat-profiles.py` script.
# 2. Generate shuffled and non-shuffled versions of the concatenated profile.
# 3. Apply UMAP in the 'umap-analysis.r' notebook to the concatenated profiles.
# 4. Save the results in the 'results/umap' directory.

# Activate the conda environment for CFReT map analysis
conda activate cfret-map

# Convert all Jupyter notebooks in the current directory to Python scripts
# and save them in the 'nbconverted/' directory
jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb

# Execute the profile concatenation script
python nbconverted/concat-profiles.py

# Switch to the R environment to run the UMAP analysis
conda activate cfret-r-analysis

# Execute the UMAP analysis script
Rscript nbconverted/umap-analysis.r
