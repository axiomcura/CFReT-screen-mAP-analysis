#!/bin/bash

# This shell script runs the pairwise-compare analysis on the aggregated feature set profiles.
# Steps:
# 0. converts all notebooks in to python and R scripts
# 1. Executes pairwise-compare generating scores when comparing controls with each other and other where we us controls as areference and compare it to all treated failing cells
# 2. r_plot-pairwise-compare.r generates the plots for the pairwise comparison

# Activate the conda environment for CFReT map analysis
conda activate cfret-map

# Convert all Jupyter notebooks in the current directory to Python scripts
# and save them in the 'nbconverted/' directory
jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb

# Execute the profile concatenation script
python nbconverted/pairwise-compare.py

# Switch to the R environment to run the UMAP analysis
conda activate cfret-r-analysis

# Execute the UMAP analysis script
Rscript nbconverted/r_plot-pairwise-scores.r
