#!/bin/bash

# This shell script runs the correlation analysis on the aggregated feature set profiles.

# Activate the conda environment for CFReT map analysis
conda activate cfret-map

# Convert all Jupyter notebooks in the current directory to Python scripts
# and save them in the 'nbconverted/' directory
jupyter nbconvert --to script --output-dir=nbconverted/ *.ipynb

# Execute the profile concatenation script
python nbconverted/correlation-analysis.py
