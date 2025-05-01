#!/usr/bin/env python

# # Preprocessing data
#
# In this notebook, we will preprocess the dataset by loading all profiles in batches, adding additional labels, and concatenating them into a single dataset for further exploration.

# In[1]:


import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.append("../../../")
from utils import data_utils, io_utils

# ## Helper functions
# These are helper functions that will be used only in this notebook

# In[2]:


def add_control_type(profile: pd.DataFrame) -> pd.DataFrame:
    """Add control type metadata to the dataframe based on cell type.

    Parameters
    ----------
    profile : pandas.DataFrame
        DataFrame containing the profiles with 'Metadata_cell_type' column.

    Returns
    -------
    pandas.DataFrame
        DataFrame with an additional 'Metadata_control_type' column.
    """
    # add a new column to the dataframe
    profile.insert(2, "Metadata_control_type", np.nan)

    # this adds the label "positive" to wells that contains healthy cells and treated with DMSO
    profile.loc[(profile["Metadata_cell_type"] == "healthy") & (profile["Metadata_treatment"] == "DMSO"), "Metadata_control_type"] = (
        "positive"
    )

    # this adds the label "negative" to wells that contains failing CF cells and treated with DMSO
    profile.loc[(profile["Metadata_cell_type"] == "failing") & (profile["Metadata_treatment"] == "DMSO"), "Metadata_control_type"] = (
        "negative"
    )

    # this adds the label "trt" to wells that contains failing CF cells and treated with a compound
    profile.loc[(profile["Metadata_cell_type"] == "failing") & (profile["Metadata_treatment"] != "DMSO"), "Metadata_control_type"] = (
        "trt"
    )
    return profile

def update_control_treatment(profile: pd.DataFrame) -> pd.DataFrame:
    """Update the Metadata_treatment column based on the Metadata_control_type column.

    If Metadata_control_type is positive, Metadata_treatment is updated to DMSO-positive.
    If Metadata_control_type is negative, Metadata_treatment is updated to DMSO-negative.

    Parameters
    ----------
    profile : pd.DataFrame
        profiles dataframe with Metadata_control_type column.

    Returns
    -------
    pd.DataFrame
        profiles dataframe with updated Metadata_treatment column based on Metadata_control_type.
    """
    profile.loc[profile["Metadata_control_type"] == "positive", "Metadata_treatment"] = "DMSO-positive"
    profile.loc[profile["Metadata_control_type"] == "negative", "Metadata_treatment"] = "DMSO-negative"
    return profile


# Setting up paths on what to load and output directories

# In[3]:


# setting in input paths
data_dir_path = pathlib.Path("../../data")

# selecting aggregated feature selected files
list_of_paths = list(
    (data_dir_path / "aggregated_profiles/").resolve(strict=True).glob("*.parquet")
)

# shared features columns
shared_features_path = pathlib.Path(
    "../../1.map-analysis/results/shared_features.json"
).resolve(strict=True)

# set configs path
config_path = pathlib.Path("../../config.yaml").resolve(strict=True)

# creating a results output directory
results_dir = pathlib.Path("./results/concat_data").resolve()
results_dir.mkdir(exist_ok=True, parents=True)


# Next, we load the configuration file that contains the shared features across all plates within the batch. Then, we load each plate within the batch and add new metadata columns for downstream analysis.

# In[4]:


# loading config
config = io_utils.load_config(config_path)
plate_name_lookup = config["general_configs"]["plate_name_lookup"]["batch_1"]

# loading shared features
shared_features = io_utils.load_config(shared_features_path)["shared_features"]

# loading all feature selected aggregated profiles and updating it with the shared features
loaded_aggregated_profiles = []
loaded_shuffled_profiles = []
for plate_idx, profile_path in enumerate(list_of_paths):
    # getting the plate name
    plate_name = profile_path.stem.split("_")[0]

    # loading aggregated profiles
    aggregated_profiles = pd.read_parquet(profile_path)

    # updating the profile with the shared features
    aggregated_profiles = aggregated_profiles[shared_features]

    # inserting the plate name at the first column
    aggregated_profiles.insert(0, "Metadata_plate_barcode", plate_name)
    aggregated_profiles.insert(1, "Metadata_plate_name", plate_name_lookup[plate_name])


    # next is to shuffled the data
    shuffled_aggregated_profiles = data_utils.shuffle_features(aggregated_profiles)

    # append it to the list
    loaded_aggregated_profiles.append(aggregated_profiles)
    loaded_shuffled_profiles.append(shuffled_aggregated_profiles)

# concatenating all the profiles
loaded_aggregated_profiles = pd.concat(loaded_aggregated_profiles).reset_index(
    drop=True
)
shuffled_aggregated_profiles = pd.concat(loaded_shuffled_profiles).reset_index(
    drop=True
)

# add metadata into the dmso profile where if Metadata_cell_type == "healthy" then Metadata_control_type == "positive"
# add if Metadata_cell_type == "failing" then Metadata_control_type == "negative"
# Apply to both aggregated and shuffled profiles
loaded_aggregated_profiles = add_control_type(loaded_aggregated_profiles)
shuffled_aggregated_profiles = add_control_type(shuffled_aggregated_profiles)

# update Metadata_treatment column based on Metadata_control_type column
loaded_aggregated_profiles = update_control_treatment(loaded_aggregated_profiles)
shuffled_aggregated_profiles = update_control_treatment(shuffled_aggregated_profiles)

# split metadata and morphology feature columns
meta_cols, feat_cols = data_utils.split_meta_and_features(loaded_aggregated_profiles)

# update the Metadata_Pathway column based on the Metadata_treatment column.
# if Metadata_treatment is "DMSO-positive", set Metadata_Pathway to "DMSO-positive".
# if Metadata_treatment is "DMSO-negative", set Metadata_Pathway to "DMSO-negative".
loaded_aggregated_profiles.loc[
    loaded_aggregated_profiles["Metadata_treatment"] == "DMSO-positive",
    "Metadata_Pathway",
] = "DMSO-positive"
loaded_aggregated_profiles.loc[
    loaded_aggregated_profiles["Metadata_treatment"] == "DMSO-negative",
    "Metadata_Pathway",
] = "DMSO-negative"

# store aggregate data profiles as batched
loaded_profiles = {"batch_1": loaded_aggregated_profiles}
shuffled_loaded_profiles = {"batch_1": shuffled_aggregated_profiles}

# display only not shuffled aggregated profiles dmso profiles
print(loaded_profiles["batch_1"].shape)
loaded_profiles["batch_1"].head()


# We save the concatenated profiles of this batch.

# In[5]:


# Creating a for loop to save the profiles
for batch, profile in loaded_profiles.items():
    profile.to_csv(results_dir /
                   f"{batch}_concat_agg_fs.csv", index=False)


# save the shuffled profiles
for batch, profile in shuffled_loaded_profiles.items():
    profile.to_csv(results_dir /
                   f"shuffled_{batch}_concat_agg_fs.csv", index=False)


# In[6]:


loaded_profiles["batch_1"]["Metadata_plate_name"].unique()
