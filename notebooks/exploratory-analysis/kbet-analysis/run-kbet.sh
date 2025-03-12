#!/bin/bash

# This script runs the kbet analysis on the CFReT screen data to see if there is any
# batch effects in the data. The kbet analysis is a non-parametric method that compares
# the distribution of the data in the different batches. The analysis is performed on the
# feature selected aggregated profiles.

# Activate R environment to run the kBET analysis
conda activate cfret-r-analysis

# Convert all Jupyter notebooks in the current directory to scripts
# and save them in the 'nbconverted/' directory
jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb

# execute the kbet analysis script
Rscript nbconverted/kbet-analysis.r
