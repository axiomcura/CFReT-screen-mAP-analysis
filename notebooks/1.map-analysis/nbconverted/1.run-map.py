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
# **what is outputed**
# - AP scores generated using both the postive and negvative controls
# - mAP scores generated using both the postive and negative controls

# In[1]:


import pathlib
import sys
import warnings
from pprint import pprint

import pandas as pd
from copairs import map
from pycytominer.cyto_utils import load_profiles
from tqdm import TqdmWarning

sys.path.append("../../")
from utils import data_utils, io_utils

# removing warnigns
warnings.filterwarnings("ignore", category=TqdmWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)


# Helper functions for this notebook

# In[2]:


def label_control_types(metadata_cell_type, heart_failure_type):
    if metadata_cell_type == "healthy" and heart_failure_type is None:
        return "positive"
    elif (
        metadata_cell_type == "failing"
        and heart_failure_type == "dilated_cardiomyopathy"
    ):
        return "negative"
    else:
        raise ValueError(f"Unknown combination added: {metadata_cell_type} {heart_failure_type}")


# This code sets up the necessary file paths and directories required for the notebook, ensuring that input files exist.
# It also creates a results folder if it doesn't already exist to store outputs generated during the analysis.

# In[3]:


# Setting the base data directory and ensure it exists (raises an error if it doesn't)
data_dir = pathlib.Path("../data/").resolve(strict=True)

# Setting the metadata directory for updated plate maps and ensure it exists
metadata_dir = pathlib.Path("../data/metadata/updated_platemaps").resolve(strict=True)

# Path to the updated barcode plate map file, ensure it exists
platemap_path = (metadata_dir / "updated_barcode_platemap.csv").resolve(strict=True)

# Path to the configuration file (does not enforce existence check here)
config_path = pathlib.Path("../config.yaml").resolve(strict=True)

# Setting the results directory, resolve the full path, and create it if it doesn't already exist
results_dir = pathlib.Path("./results").resolve()
map_results_dir = (results_dir /"map_scores").resolve()
map_results_dir.mkdir(exist_ok=True, parents=True)


# Loading in the files and setting config parameters

# In[4]:


# loading config and general configs
configs = io_utils.load_config(config_path)
general_configs = configs["general_configs"]

# loading bar code
barcode = pd.read_csv(platemap_path)

# setting notebook specific parameters
control_list = [("negative", "DMSO", "failing"), ("positive", "DMSO", "healthy")]


# Since these files have undergone feature selection, it is essential to identify the overlapping feature names to ensure accurate and consistent analysis.

# In[5]:


shared_cols = None
for aggregated_profile in list(data_dir.glob("*.parquet")):
    # read aggreagated profiled and column names
    agg_df = pd.read_parquet(aggregated_profile)
    columns = list(agg_df.columns)

    # Update the shared_columns set
    if shared_cols is None:
        # Initialize shared columns with the first profile's columns, preserving order
        shared_cols = columns
    else:
        # Retain only the columns present in both the current profile and shared columns
        shared_cols = [col for col in shared_cols if col in columns]


# ## Formatting Data

# In this section, the code processes and organizes data by grouping related files and enriching them with additional metadata. Each group is assigned a unique identifier, and the corresponding data files are systematically loaded and prepared. New metadata columns are generated by combining existing information to ensure consistency and clarity. Additional metadata is integrated into the data to provide valuable experimental context, while unique identifiers are added to distinguish the aggregated profiles from different batches.

# In[6]:


# Suffix for aggregated profiles
aggregated_file_suffix = "aggregated_post_fs.parquet"

# Dictionary to store loaded plate data grouped by batch
loaded_plate_batches = {}

# Iterate over unique platemap files and their associated plates
for batch_index, (platemap_filename, associated_plates_df) in enumerate(
    barcode.groupby("platemap_file")
):
    # Generate a unique batch ID
    batch_id = f"batch_{batch_index + 1}"

    # Load the platemap CSV file
    platemap_path = (metadata_dir / f"{platemap_filename}.csv").resolve(strict=True)
    platemap_data = pd.read_csv(platemap_path)

    # Extract all plate names associated with the current platemap
    plate_barcodes = associated_plates_df["plate_barcode"].tolist()

    # List to store all loaded and processed aggregated plates for the current batch
    loaded_aggregated_plates = []

    for plate_barcode in plate_barcodes:
        # Resolve the file path for the aggregated plate data
        plate_file_path = (
            data_dir / f"{plate_barcode}_{aggregated_file_suffix}"
        ).resolve(strict=True)

        # Load the aggregated profile data for the current plate
        aggregated_data = load_profiles(plate_file_path)

        # Update loaded data frame with only shared features
        aggregated_data = aggregated_data[shared_cols]

        # Add a new column indicating the source plate for each row
        aggregated_data.insert(0, "Metadata_plate_barcode", plate_barcode)

        # Append the processed aggregated data for this plate to the batch list
        loaded_aggregated_plates.append(aggregated_data)

    # Combine all processed plates for the current batch into a single DataFrame
    combined_aggregated_data = pd.concat(loaded_aggregated_plates)
    meta_concat, feats_concat = data_utils.split_meta_and_features(
        combined_aggregated_data
    )

    # Store the combined DataFrame in the loaded_plate_batches dictionary
    loaded_plate_batches[batch_id] = combined_aggregated_data


# ## Running mAP only on controls across all plates

# In this section, we calculate the mAP (mean Average Precision) scores between controls to assess their quality. Specifically, we aim to evaluate how the negative control compares when using a positive control as a reference, and vice versa. This analysis helps determine whether the controls in the experiment are reliable indicators of quality and consistency. Reliable controls are critical for ensuring the validity of the experiment's results.

# In[7]:


## Loading configurations
cntrl_copairs_ap_configs = configs["cntrl_copairs_ap_configs"]
cntrl_copairs_map_configs = configs["cntrl_copairs_map_configs"]


# In[8]:


profile = loaded_plate_batches["batch_1"]
dmso_profile = profile.loc[profile["Metadata_treatment"] == "DMSO"]
plate_ids = dmso_profile["Metadata_plate_barcode"].unique().tolist()

# add control type information
# adding control_type information into the data frame
dmso_profile["Metadata_treatment_type"] = "control"
dmso_profile["Metadata_control_type"] = dmso_profile.apply(
    lambda row: label_control_types(
        row["Metadata_cell_type"], row["Metadata_heart_failure_type"]
    ),
    axis=1,
)
dmso_profile = dmso_profile.reset_index().rename(columns={"index": "original_index"})


# In[9]:


# List of control types to evaluate
control_list = ["negative", "positive"]

# Iterate over batches of loaded plate profiles
for batch_id, profile in loaded_plate_batches.items():
    # Filter profiles for DMSO-treated wells
    dmso_profile = profile.loc[profile["Metadata_treatment"] == "DMSO"]

    # Get unique plate IDs for DMSO-treated wells
    plate_ids = dmso_profile["Metadata_plate_barcode"].unique().tolist()

    # Add control type information to the dataframe
    dmso_profile["Metadata_treatment_type"] = "control"  # Tag all rows as control
    dmso_profile["Metadata_control_type"] = dmso_profile.apply(
        lambda row: label_control_types(
            row["Metadata_cell_type"], row["Metadata_heart_failure_type"]
        ),
        axis=1,
    )

    # Reset index and store the original index for reference
    dmso_profile = dmso_profile.reset_index().rename(
        columns={"index": "original_index"}
    )

    # Iterate over control types to use them as references
    for ref_type in control_list:
        print(f"Using '{ref_type}' as reference to calculate mAP")

        ap_scores = []  # Initialize list to store AP scores

        # Iterate over all targeted plate IDs
        for targeted_plate_id in plate_ids:
            # Create a deep copy of the DMSO profile for manipulation
            dmso_profile_w_target_plate = dmso_profile.copy(deep=True)

            # Tag rows corresponding to the targeted plate
            dmso_profile_w_target_plate["Metadata_targeted"] = (
                dmso_profile_w_target_plate["Metadata_plate_barcode"].apply(
                    lambda plate_id: plate_id == targeted_plate_id
                )
            )

            # Initialize reference index for mAP calculation
            # Default to -1 for all wells except targeted reference wells
            dmso_profile_w_target_plate["Metadata_reference_index"] = (
                dmso_profile_w_target_plate.index
            )
            dmso_profile_w_target_plate["Metadata_reference_index"] = (
                dmso_profile_w_target_plate.apply(
                    lambda row: row["Metadata_reference_index"]
                    if row["Metadata_targeted"]
                    and row["Metadata_control_type"] == ref_type
                    else -1,
                    axis=1,
                )
            )

            # Split metadata and feature columns for analysis
            dmso_meta, dmso_feats = data_utils.split_meta_and_features(
                dmso_profile_w_target_plate
            )

            # Compute average precision (AP) scores for the current setup
            dmso_ap_scores = map.average_precision(
                meta=dmso_profile_w_target_plate[dmso_meta],
                feats=dmso_profile_w_target_plate[dmso_feats].values,
                pos_sameby=cntrl_copairs_ap_configs["pos_sameby"],
                pos_diffby=[],
                neg_sameby=[],
                neg_diffby=cntrl_copairs_ap_configs["neg_diffby"],
                batch_size=cntrl_copairs_ap_configs["batch_size"],
                distance=cntrl_copairs_ap_configs["distance"],
            )

            # Append the computed AP scores for this targeted plate
            ap_scores.append(dmso_ap_scores)

        # Concatenate all AP scores into a single dataframe
        dmso_ap_scores = pd.concat(ap_scores)

        # Calculate mean Average Precision (mAP) scores
        dmso_map_scores = map.mean_average_precision(
            dmso_ap_scores,
            sameby=cntrl_copairs_map_configs["same_by"],
            null_size=cntrl_copairs_map_configs["null_size"],
            threshold=cntrl_copairs_map_configs["threshold"],
            seed=general_configs["seed"],
        )

        # Store the computed AP and mAP scores as CSV files
        dmso_ap_scores.to_csv(
            map_results_dir / f"{batch_id}_{ref_type}_ref_dmso_AP_scores.csv"
        )
        dmso_map_scores.to_csv(
            map_results_dir / f"{batch_id}_{ref_type}_ref_dmso_mAP_scores.csv"
        )


# ## Calculating mAP scores on only treatments

# In this section, we analyze a high-content screening dataset generated from cell painting experiments, where failing cardiac fibroblasts are treated with multiple compounds. Our goal is to calculate the mean average precision (mAP) by comparing the experimental treatments to two controls: a negative control consisting of DMSO-treated failing cardiac fibroblasts and a positive control consisting of DMSO-treated healthy cardiac fibroblasts.
#
# We start by preparing the dataset, copying the profiles, and assigning a reference index to ensure proper grouping of non-DMSO treatment replicates. Metadata and feature columns are separated to facilitate the calculation of average precision (AP) scores. To calculate these scores, we define positive pairs as treatments with the same metadata values (e.g., same treatment type) across all plates. Negative pairs, on the other hand, are determined by comparing all DMSO-treated wells across all plates with all other treatments.
#
# Once the AP scores are computed, we aggregate them across all plates for each treatment to derive the mean average precision (mAP) score. This process captures the consistency of treatment performance relative to the controls and allows for a comprehensive evaluation of the dataset. Finally, we save both the AP and mAP scores for each control condition, providing a well-structured dataset for further interpretation and downstream analysis.

# In[10]:


# Load configurations for average precision (AP) and mean average precision (mAP)
trt_copairs_ap_configs = configs["trt_copairs_ap_configs"]
trt_copairs_map_configs = configs["trt_copairs_map_configs"]

# displaying there parameters that were used to execute mAP pipeline
print("AP paramters:")
pprint(trt_copairs_ap_configs)

print("\nmAP paramters:")
pprint(trt_copairs_map_configs)


# In[11]:


# Define control conditions for the analysis
# Each tuple specifies the control type, treatment, and associated cell state
control_list = [("negative", "DMSO", "failing"), ("positive", "DMSO", "healthy")]

# Iterate over each batch of loaded plate profiles
for batch_id, profile in loaded_plate_batches.items():
    # Analyze the profile for each control condition
    for control_type, control_treatment, cell_state in control_list:
        # Create a copy of the profile to preserve the original data
        profile = profile.copy()

        # Assign a default reference index based on the row index
        profile["Metadata_reference_index"] = profile.index

        # Mark all non-control replicates (e.g., treatments not matching the current control)
        profile.loc[
            (profile["Metadata_treatment"] != control_treatment)
            & (profile["Metadata_cell_type"] != cell_state),
            "Metadata_reference_index",
        ] = -1

        # Move the "Metadata_reference_index" column to the beginning for clarity
        profile.insert(
            0, "Metadata_reference_index", profile.pop("Metadata_reference_index")
        )

        # Separate metadata columns from feature columns for downstream calculations
        meta_columns, feature_columns = data_utils.split_meta_and_features(profile)

        # Calculate average precision (AP) for the profile
        # Positive pairs are based on treatments with the same metadata
        # Negative pairs compare all DMSO-treated wells to all treatments
        trt_replicate_aps = map.average_precision(
            meta=profile[meta_columns],
            feats=profile[feature_columns].values,
            pos_sameby=trt_copairs_ap_configs["pos_sameby"],
            pos_diffby=trt_copairs_ap_configs["pos_diffby"],
            neg_sameby=[],
            neg_diffby=trt_copairs_ap_configs["neg_diffby"],
        )

        # Calculating mAP scores for only treatments (no controls)
        # Exclude wells treated with the control treatment (DMSO)
        replicate_aps = trt_replicate_aps.loc[
            trt_replicate_aps["Metadata_treatment"] != control_treatment
        ]

        # Save the calculated AP scores to a file for further analysis
        trt_replicate_aps.to_csv(
            map_results_dir
            / f"{control_type}_control_{cell_state}_{control_treatment}_AP_scores.csv",
            index=False,
        )

        # Calculate mean average precision (mAP) from the AP scores
        trt_replicate_maps = map.mean_average_precision(
            trt_replicate_aps,
            sameby=trt_copairs_map_configs["same_by"],  # Grouping criteria for mAP
            null_size=trt_copairs_map_configs["null_size"],  # Null distribution size
            threshold=trt_copairs_map_configs["threshold"],  # Significance threshold
            seed=general_configs["seed"],  # Seed for reproducibility
        )

        # Save the mAP scores to a file for reporting
        trt_replicate_maps.to_csv(
            map_results_dir
            / f"{control_type}_control_{cell_state}_{control_treatment}_mAP_scores.csv",
            index=False,
        )
