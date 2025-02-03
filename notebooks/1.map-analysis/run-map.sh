#!/bin/bash

# This script converts Jupyter notebooks to Python scripts and activates the analytical
# pipeline for calculating and analyzing mAP scores for the CFReT dataset.

# Activate the conda environment
conda activate cfret-map

# Convert all Jupyter notebooks in the current directory to Python scripts
# and save them in the 'nbconverted/' directory
jupyter nbconvert --to python --output-dir=nbconverted/ *.ipynb

# Execute the data mAP analysis scripts
python nbconverted/1.run-map.py
python nbconverted/2.map-analysis.py
