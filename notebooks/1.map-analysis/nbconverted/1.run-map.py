#!/usr/bin/env python

# # 1. Executing Mean Average Precision (mAP)
#
# Mean Average Precision (mAP) is a flexible statistical framework used to measure the **phenotypic activity** of compounds by comparing them to control groups. In this notebook, we utilize high-content screening data, that used the CellPainting assay, to identify potential drug candidates that demonstrate evidence of reversing the effects of cardiac fibrosis. The dataset comprises **image-based profiles at the replicate level (well-level)**.
#
# #### **Controls Used in the Screen**
# To interpret mAP scores, we leverage the following control groups:
# - **Negative control**: Failing CF cells treated with DMSO.
# - **Positive control**: Healthy CF cells treated with DMSO.
#
# #### **Interpreting mAP Scores**
# - **High mAP Scores**:
#   Indicate that wells treated with a specific compound are highly phenotypically distinct compared to the control. This suggests the compound induces a strong and specific phenotypic change.
#
# - **Low mAP Scores**:
#   Indicate that wells treated with a specific compound are phenotypically similar to the control. This suggests the compound has little to no phenotypic effect or a nonspecific one.
#
# #### **Biological Interpretation**
# mAP scores help determine which compounds exhibit phenotypic changes that resemble those of healthy cells, making them potential candidates for reversing the effects of cardiac fibrosis. By comparing the phenotypic activity of compounds to both positive and negative controls, we can prioritize compounds for further validation.
#
# **what is outputted**
# - AP scores generated using both the positive and negative controls
# - mAP scores generated using both the positive and negative controls

# In[1]:


import json
import pathlib
import sys
import warnings

import pandas as pd
from pycytominer.cyto_utils import load_profiles
from tqdm import TqdmWarning

sys.path.append("../../")
from utils import analysis_utils, data_utils, io_utils

# removing warnigns
warnings.filterwarnings("ignore", category=TqdmWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)


# Helper functions

# In[2]:


def update_control_pathways(cell_type:str, treatment:str, pathway:str|None) -> str:
    """ Updates the Metadata Pathway column based on the cell type, treatment, and pathway.

    This function maps the pathway information for each sample based on specific rules:
    - If the cell type is "healthy" and the treatment is "DMSO", the pathway is labeled as "DMSO-positive".
    - If the cell type is "failing" and the treatment is "DMSO", the pathway is labeled as "DMSO-negative".
    - If the treatment is not "DMSO" and the pathway is None or NaN, the pathway is labeled as "No Pathway".
    - Otherwise, the original pathway value is returned.

    Parameters
    ----------
    cell_type : str
        The type of the cell (e.g., "healthy", "failing").
    treatment : str
        The treatment applied to the cell (e.g., "DMSO", "UCD-0159256").
    pathway : str | None
        The pathway associated with the cell, or None/NaN if not available.

    Returns
    -------
    str
        The updated pathway label based on the provided rules.
    """
    # Check if the cell type is "healthy" and treatment is "DMSO"
    if cell_type == "healthy" and treatment == "DMSO":
        return "DMSO-positive"
    # Check if the cell type is "failing" and treatment is "DMSO"
    elif cell_type == "failing" and treatment == "DMSO":
        return "DMSO-negative"
    # Check if the treatment is not "DMSO" and pathway is None/NaN
    elif treatment != "DMSO" and (pathway is None or pd.isna(pathway)):
        return "No Pathway"
    # Return the original pathway if no conditions are met
    return pathway


# This code sets up the necessary file paths and directories required for the notebook, ensuring that input files exist.
# It also creates a results folder if it doesn't already exist to store outputs generated during the analysis.

# In[3]:


# Setting the base data directory and ensure it exists (raises an error if it doesn't)
main_results_dir = pathlib.Path("./results/").resolve(strict=True)
data_dir = pathlib.Path("../data/").resolve(strict=True)
agg_data_dir = (data_dir / "aggregated_profiles").resolve(strict=True)
fs_profiles_paths = list((data_dir / "aggregated_profiles").resolve(strict=True).glob("*.parquet"))

# Setting the metadata directory for updated plate maps and ensure it exists
metadata_dir = pathlib.Path("../data/metadata/updated_platemaps").resolve(strict=True)

# Path to the updated barcode plate map file, ensure it exists
platemap_path = (metadata_dir / "updated_barcode_platemap.csv").resolve(strict=True)

# Path to the configuration file (does not enforce existence check here)
config_path = pathlib.Path("../config.yaml").resolve(strict=True)

# Setting the results directory, resolve the full path, and create it if it doesn't already exist
map_results_dir = pathlib.Path("./results/map_scores").resolve()
map_results_dir.mkdir(exist_ok=True, parents=True)


# Loading in the files

# In[4]:


# loading config and general configs
configs = io_utils.load_config(config_path)
general_configs = configs["general_configs"]
plate_name_lookup = general_configs["plate_name_lookup"]

# loading bar code
barcode = pd.read_csv(platemap_path)


# Since these files have undergone feature selection, it is essential to identify the overlapping feature names to ensure accurate and consistent analysis.

# In[5]:


# finding shared features while deleting duplicate column names
shared_cols = data_utils.find_shared_features(profile_paths=fs_profiles_paths, delete_dups=True)

# saving shared features to a json file
# if the file already exists, it will not be overwritten
if not (main_results_dir / "shared_features.json").exists():
    print("shared_features.json does not exist. Saving shared features to a json file.")
    shared_cols_dict = {}
    shared_cols_dict["shared_features"] = shared_cols
    with open(main_results_dir / "shared_features.json", "w") as shared_file:
        json.dump(shared_cols_dict, shared_file)

# if the file already exists, then we check if the shared features are the same
else:
    with open(main_results_dir / "shared_features.json") as shared_file:
        shared_cols_dict = json.load(shared_file)
    assert shared_cols == shared_cols_dict["shared_features"], "Shared features are not the same"

# total amount of shared columns among all profiles in batch 1
print("Total amount of shared columns among all profiles:")
print(len(shared_cols))


# In this section, the code processes and organizes data by grouping related files and enriching them with additional metadata. Each group is assigned a unique identifier, and the corresponding data files are systematically loaded and prepared. New metadata columns are generated by combining existing information to ensure consistency and clarity. Additional metadata is integrated into the data to provide valuable experimental context, while unique identifiers are added to distinguish the aggregated profiles from different batches.

# In[6]:


# suffix for aggregated profiles
aggregated_file_suffix = "aggregated_post_fs.parquet"

# dictionary to store loaded plate data grouped by batch
loaded_plate_batches = {}
loaded_shuffled_plate_batches = {}

# iterate over unique platemap files and their associated plates
for batch_index, (platemap_filename, associated_plates_df) in enumerate(
    barcode.groupby("platemap_file")
):
    # generate a unique batch ID
    batch_id = f"batch_{batch_index + 1}"

    # load the platemap CSV file
    platemap_path = (metadata_dir / f"{platemap_filename}.csv").resolve(strict=True)
    platemap_data = pd.read_csv(platemap_path)

    # extract all plate names associated with the current platemap
    plate_barcodes = associated_plates_df["plate_barcode"].tolist()

    # list to store all loaded and processed aggregated plates for the current batch
    loaded_aggregated_plates = []
    loaded_shuffled_aggregated_plates = []

    for plate_barcode in plate_barcodes:
        # resolve the file path for the aggregated plate data
        plate_file_path = (
            agg_data_dir / f"{plate_barcode}_{aggregated_file_suffix}"
        ).resolve(strict=True)

        # load the aggregated profile data for the current plate
        aggregated_data = load_profiles(plate_file_path)

        # update loaded data frame with only shared features
        aggregated_data = aggregated_data[shared_cols]

        # add a new column indicating the source plate for each row
        aggregated_data.insert(0,"Metadata_plate_barcode" , plate_barcode)

        # Add plate name
        aggregated_data.insert(1, "Metadata_plate_name", aggregated_data["Metadata_plate_barcode"].map(plate_name_lookup["batch_1"]))

        # Update Metadata_Pathway column
        aggregated_data["Metadata_Pathway"] = aggregated_data.apply(
            lambda row: update_control_pathways(row["Metadata_cell_type"], row["Metadata_treatment"], row["Metadata_Pathway"]), axis=1
        )

        # append the processed aggregated data for this plate to the batch list
        loaded_aggregated_plates.append(aggregated_data)

        # adding shuffled aggregated profiles
        shuffled_aggregated_data = data_utils.shuffle_features(aggregated_data)

        # append the processed and shuffled aggregated data for this plate to the batch list
        loaded_shuffled_aggregated_plates.append(shuffled_aggregated_data)

    # combine all processed plates for the current batch into a single DataFrame
    combined_aggregated_data = pd.concat(loaded_aggregated_plates).reset_index(drop=True)
    meta_concat, feats_concat = data_utils.split_meta_and_features(combined_aggregated_data)

    # combine all shuffled and processed plates for the current batch into a single DataFrame
    # shuffled_combined_aggregated_data = pd.concat(loaded_shuffled_aggregated_plates).reset_index().rename(columns={"index": "Metadata_old_index"})
    shuffled_combined_aggregated_data = pd.concat(loaded_shuffled_aggregated_plates).reset_index(drop=True)
    meta_concat, feats_concat = data_utils.split_meta_and_features(shuffled_combined_aggregated_data)

    # store the combined DataFrame in the loaded_plate_batches dictionary
    loaded_plate_batches[batch_id] = combined_aggregated_data
    loaded_shuffled_plate_batches[batch_id] = shuffled_combined_aggregated_data


# ## Running mAP only on controls across all plates
#

# In this section, we calculate the mAP (mean Average Precision) scores between controls to assess their quality. Specifically, we aim to evaluate how the negative control compares when using a positive control as a reference, and vice versa. This analysis helps determine whether the controls in the experiment are reliable indicators of quality and consistency. Reliable controls are critical for ensuring the validity of the experiment's results.

# In[7]:


# calculating mAP scores only on with original DMSO profiles
analysis_utils.calculate_dmso_map_batch_profiles(
    batched_profiles=loaded_plate_batches,
    configs=configs,
    outdir_path=map_results_dir,
    shuffled=False,
)

# calculating mAP scores only on with shuffled DMSO profiles
analysis_utils.calculate_dmso_map_batch_profiles(
    batched_profiles=loaded_shuffled_plate_batches,
    configs=configs,
    outdir_path=map_results_dir,
    shuffled=True,
)


# ## Calculating mAP scores on only treatments
#

# In this section, we analyze a high-content screening dataset generated from cell painting experiments, where failing cardiac fibroblasts are treated with multiple compounds. Our goal is to calculate the mean average precision (mAP) by comparing the experimental treatments to two controls: a negative control consisting of DMSO-treated failing cardiac fibroblasts and a positive control consisting of DMSO-treated healthy cardiac fibroblasts.

# In[8]:


# Here we execute mAP pipeline with with the original
analysis_utils.calculate_trt_map_batch_profiles(
    batched_profiles=loaded_plate_batches,
    configs=configs,
    outdir_path=map_results_dir,
    shuffled=False
)

# Here we execute mAP pipeline with with the shuffled dataset
analysis_utils.calculate_trt_map_batch_profiles(
    batched_profiles=loaded_shuffled_plate_batches,
    configs=configs,
    outdir_path=map_results_dir,
    shuffled=True
)
