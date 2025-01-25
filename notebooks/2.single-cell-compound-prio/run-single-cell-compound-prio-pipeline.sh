#!/bin/bash

# This shell scripts runs the single-cell compound prioritization pipeline

# Activate the conda environment
conda activate cfret-map

# Convert all Jupyter notebooks in the current directory to Python scripts
# and save them in the 'nbconverted/' directory
jupyter nbconvert --to python --output-dir=nbconverted/ *.ipynb

# Execute the data download script
python nbconverted/1.on_off_morphology_signatures.py
