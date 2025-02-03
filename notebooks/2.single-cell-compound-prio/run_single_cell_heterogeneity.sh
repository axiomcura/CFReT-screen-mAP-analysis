#!/bin/bash

# This shell script converts all individual steps required in the compound prioritization
# pipeline to Python scripts and executes them in sequence. This is done to ensure that
# the pipeline can be executed in a single command.

# Activate the conda environment
conda activate cfret-map

# Convert all Jupyter notebooks in the current directory to Python scripts
# and save them in the 'nbconverted/' directory
jupyter nbconvert --to python --output-dir=nbconverted/ *.ipynb

# Execute the data download script
python nbconverted/1.on_off_morphology_signatures.py
python nbconverted/2.single_cell_heterogeneity.py
python nbconverted/4.selecting_hits.py
