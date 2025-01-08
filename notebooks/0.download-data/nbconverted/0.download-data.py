#!/usr/bin/env python

# # Module 0: Downloading Metadata
#
# In this module, we download the metadata required for the CFReT screen project.
# All downloaded files and processed data will be placed in the `./data` folder.
#
# The specific data being downloaded in this script is the metadata associated
# with the CFReT screen, hosted on GitHub. The metadata folder's repository link is:
# https://github.com/WayScience/targeted_fibrosis_drug_screen/tree/main/metadata
#
#

# In[1]:


import pathlib
import shutil
import sys
import tempfile
import zipfile

import requests

# importing analysis modules/imports
sys.path.append("../../")
from src.io_utils import load_config

# Setting input and output paths

# In[2]:


# creating a data folder
data_dir_path = pathlib.Path("../data").resolve()
data_dir_path.mkdir(exist_ok=True)

# setting config path
config_path = pathlib.Path("../config.yaml").resolve(strict=True)

# setting compressed file path
zip_path = (data_dir_path / "metadata.zip").resolve()

# create a temp file:
temp_dir = tempfile.mkdtemp()
temp_dir_path = pathlib.Path(temp_dir).resolve()


# Loading in the project configuration file

# In[3]:


# loading config file
loaded_configs = load_config(config_path)
download_configs = loaded_configs["download_configs"]


# Downloading the metadata data and storing it into the `./data` directory

# In[4]:


# Download the repository as a ZIP file
response = requests.get(download_configs["github_metadata_url"], stream=True)
response.raise_for_status()

# once we get a  to github, we start download the metadata directory as a zip file
with open(zip_path, "wb") as zip_file:
    for chunk in response.iter_content(chunk_size=download_configs["download_chunk_size"]):
        zip_file.write(chunk)
print(f"Contents from the GitHub folder have been downloaded: {str(zip_path)}")

# Extract the ZIP file into the temporary directory
with zipfile.ZipFile(zip_path, "r") as zip_ref:
    zip_ref.extractall(temp_dir_path)
print(f"Repository contents extracted into a temporary folder: {str(temp_dir_path)}")

# Locate metadata folder within the extracted repository
extracted_repo_dir = next(temp_dir_path.glob("*"))  # Assuming there's only one top-level folder
source_folder = extracted_repo_dir / "metadata"
target_folder = data_dir_path / "metadata"

if source_folder.exists():
    # Copy the `metadata` folder to the final data directory
    if target_folder.exists():
        shutil.rmtree(target_folder)  # Clear any existing folder
    shutil.move(str(source_folder), str(target_folder))
    print(f"Metadata folder moved to {str(target_folder)}")
else:
    print("Metadata folder not found in the repository.")

# Clean up temporary directory and ZIP file
shutil.rmtree(temp_dir_path)
zip_path.unlink()
print("Temporary files and directories cleaned up.")
