#!/bin/bash

# Shell script to download data required for this analysis notebook

# Activate the conda environment
conda activate cfret-map

# Convert all Jupyter notebooks in the current directory to Python scripts
# and save them in the 'nbconverted/' directory
jupyter nbconvert --to python --output-dir=nbconverted/ *.ipynb

# Execute the data download script
python nbconverted/1.run-map.py
python nbconverted/2.map-analysis.py
