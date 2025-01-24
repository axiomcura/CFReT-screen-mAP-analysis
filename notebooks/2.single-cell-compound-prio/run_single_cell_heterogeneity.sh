#!/bin/bash

# Shell script executes the pipeline to identify

# Activate the conda environment
conda activate cfret-map

# Convert all Jupyter notebooks in the current directory to Python scripts
# and save them in the 'nbconverted/' directory
jupyter nbconvert --to python --output-dir=nbconverted/ *.ipynb

# Execute the data download script
python nbconverted/2.single_cell_heterogeneity.py
