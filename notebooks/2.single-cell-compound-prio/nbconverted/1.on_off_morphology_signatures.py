#!/usr/bin/env python

# # Identifying On and Off Morphological Signatures
#
# In this section, we will identify morphological features that are distinctly different between the reference group and the target group. This analysis aims to generate morphological signatures associated with the different cellular states.
#
# For this analysis:
# - The **reference state** is the negative control: Failing CF Cells treated with DMSO.
# - The **target state** is the positive control: Healthy CF Cells treated with DMSO.
#
# ### On-Morphological Features
# **On-morphological features** refer to morphological characteristics that show significant differences between the reference and target states. These features represent the on-target morphology associated with the target state.
#
# On-morphological features are crucial when developing metrics or models to differentiate the target state from the reference. These features signify cellular morphological changes that are specific to the target state and should be prioritized during metric development.
#
# ### Off-Morphological Features
# **Off-morphological features** refer to morphological characteristics that do not show significant differences between the reference and target states. These features are not associated with the target state and may reflect general cellular morphology unaffected by the treatment or target condition.
#
# These features can be leveraged to:
# - Monitor off-target effects.
# - Identify morphological changes unrelated to the target state.
#
# ### Goal of This Analysis
# The goal is to clearly separate on-target morphological features (those associated with healthy CF cells) from off-target features (those not significantly affected by the target condition). This distinction helps in designing metrics for detecting on-target effects while monitoring and minimizing off-target impacts.

# In[1]:


import json
import pathlib
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import ks_2samp
from statsmodels.stats.multitest import multipletests

sys.path.append("../../")
from utils import data_utils

# ## Helper functions

# In[2]:


def weighted_ks_test(
    reference: pd.DataFrame,
    target: pd.DataFrame,
    p_thresh: float | None = 0.05,
    correction_method: str = "fdr_bh",
) -> tuple[list[str], list[str]]:
    """ Performs a weighted Kolmogorov-Smirnov (KS) test between the target and reference
    datasets for each morphological feature. Adjusts for imbalanced sample sizes by applying
    weights to the cumulative distribution functions (CDFs). Includes multiple testing correction.

    Parameters
    ----------
    reference : pd.DataFrame
        A DataFrame containing the morphological features of the reference dataset (e.g.,
        reference cells). Each column represents a different feature, and each row represents
        a single observation (e.g., a cell).

    target : pd.DataFrame
        A DataFrame containing the morphological features of the target dataset (e.g.,
        desired cell state). Each column represents a different feature, and each row
        represents a single observation (e.g., a cell).

    p_thresh : Optional[float], default=0.05
        The significance threshold for the corrected p-value.

    correction_method : str, default="fdr_bh"
        The method for multiple testing correction. Options include:
        - "fdr_bh" (False Discovery Rate, Benjamini-Hochberg)
        - "bonferroni" (Bonferroni correction)
        Refer to `statsmodels.stats.multitest.multipletests` for other options.

    Returns
    -------
    Tuple[List[str], List[str]]
        - A list of features that are not significantly different between the target
          and reference datasets (off-morphology signatures).
        - A list of features that are significantly different (on-morphology signatures).

    Notes
    -----
    - This implementation uses weights proportional to the inverse of the dataset sizes
      to adjust for imbalances between the reference and target datasets.
    - Multiple testing correction is applied to the computed p-values.
    """

    # Store the KS statistics and raw p-values for each feature
    ks_stats = []
    raw_p_values = []
    feature_names = target.columns.tolist()

    # Iterate through each morphological feature in the dataset
    for morphology_feature in feature_names:
        # Step 1: Calculate weights for both target and reference datasets
        # Weights ensure each dataset contributes equally, regardless of sample size
        target_weights = np.ones(len(target)) / len(target)
        reference_weights = np.ones(len(reference)) / len(reference)

        # Step 2: Sort the values of the feature and their corresponding weights.
        # A required step in constructing cumulative distribution functions (CDFs),
        # which will be used on the next step
        sorted_target_indices = np.argsort(target[morphology_feature].to_numpy())
        sorted_reference_indices = np.argsort(reference[morphology_feature].to_numpy())

        sorted_target_data = target[morphology_feature].to_numpy()[
            sorted_target_indices
        ]
        sorted_reference_data = reference[morphology_feature].to_numpy()[
            sorted_reference_indices
        ]

        sorted_target_weights = target_weights[sorted_target_indices]
        sorted_reference_weights = reference_weights[sorted_reference_indices]

        # Step 3: Compute the weighted cumulative distribution functions (CDFs)
        # Use cumulative sums of sorted weights to calculate the CDFs
        weighted_target_cdf = np.cumsum(sorted_target_weights) / np.sum(
            sorted_target_weights
        )
        weighted_reference_cdf = np.cumsum(sorted_reference_weights) / np.sum(
            sorted_reference_weights
        )

        # Step 4: Find all unique feature values across both datasets
        # Unique values are necessary for interpolating CDFs
        all_values = np.unique(
            np.concatenate([sorted_target_data, sorted_reference_data])
        )

        # Step 5: Interpolate the CDFs at the unique values
        # Ensures the CDFs can be compared directly at the same points
        target_cdf_at_values = np.interp(
            all_values, sorted_target_data, weighted_target_cdf, left=0, right=1
        )
        reference_cdf_at_values = np.interp(
            all_values, sorted_reference_data, weighted_reference_cdf, left=0, right=1
        )

        # Step 6: Compute the KS statistic
        # The KS statistic is the maximum absolute difference between the two CDFs
        ks_stat = np.max(np.abs(target_cdf_at_values - reference_cdf_at_values))
        ks_stats.append(ks_stat)

        # Step 7: Compute the raw p-value using an unweighted KS test
        # The p-value is used to assess statistical significance
        _, p_val = ks_2samp(target[morphology_feature], reference[morphology_feature])
        raw_p_values.append(p_val)

    # Step 8: Apply multiple testing correction to raw p-values
    # This controls for false positives when testing multiple features
    corrected_results = multipletests(
        raw_p_values, alpha=p_thresh, method=correction_method
    )

    # we are only extracting the flags if the feature is significant or not
    # corrected p-values are not used in this implementation
    significant_flags = corrected_results[0]  # Boolean flags for significance

    # Step 9: Categorize features based on corrected p-values
    # Separate features into on-morphology (significant) and off-morphology (non-significant)
    on_morphology_signatures = [
        feature_names[i]
        for i, significant in enumerate(significant_flags)
        if significant
    ]
    off_morphology_signatures = [
        feature_names[i]
        for i, significant in enumerate(significant_flags)
        if not significant
    ]

    return off_morphology_signatures, on_morphology_signatures

def find_shared_features(profile_paths: list[str | pathlib.Path]) -> list[str]:
    """ Find the shared features (columns) between the profiles in the provided list of
    file paths, while retaining the order of features as they appear in the first
    profile.

    This function leverages the schema information from the Parquet files to extract the
    the columns names without loading in the entire dataset. The first profile is used
    as the reference for the order of features. Then, the function iterates through the
    remaining profiles to find the common features. If no common features are found, an
    empty list is returned.

    Parameters
    ----------
    profile_paths : list[str | pathlib.Path]
        A list of file paths pointing to the Parquet profile files.

    Returns
    -------
    list[str]
        A list of features (column names) that are common across all profiles, retaining the order
        from the first profile.

    Raises
    ------
    ValueError
        If any of the profile paths do not point to Parquet files.
    """
    # type checker to check if the file provided are parquet files
    for path in profile_paths:
        if not path.suffix == ".parquet" or path.suffix == ".pq":
            raise ValueError("All profile paths must point to Parquet files.")

    # initialize the shared features to None
    shared_features = None

    # iterate through the profile paths
    for profile_path in profile_paths:
        # Load the metadata of the Parquet file
        parquet_metadata = pq.ParquetFile(profile_path)

        # Extract column names from the schema
        column_names = parquet_metadata.schema.names

        if shared_features is None:
            # Initialize shared features on the first iteration
            shared_features = column_names
        else:
            # Retain only the features that are shared, keeping their order
            shared_features = [name for name in shared_features if name in column_names]

    return shared_features if shared_features else []


# ## Applying weighted KS to the CFReT dataset

# In this section, we load the single-cell, feature-selected profiles of the screen data. All controls are pooled into a single group, establishing the "on" and "off" morphological signatures based on these controls. The reference group is the negative control, comprising failing cardiac fibroblast (CF) cells treated with DMSO. In contrast, the positive control, referred to as the target, consists of healthy CF cells treated with DMSO. This approach is guided by the hypothesis of identifying potential compounds that can reverse failing CF cells to a healthy state.
#
# The morphological signatures will be saved as JSON files. These files will include the morphological features present in both the "on" and "off" signatures, preserving the order of the features. This structure ensures consistency for downstream analyses.

# In[3]:


# setting the directory path of the data
data_dir = pathlib.Path("../data").resolve(strict=True)

# setting path to the single-cell feature selected dataset
selected_features_files = list(data_dir.glob("*sc_feature_selected.parquet"))

# setting path where to save our results
results_dir = pathlib.Path("results").resolve()
results_dir.mkdir(exist_ok=True)

# making sub directories for the results
morph_signatures_dir = (results_dir / "morph_signatures").resolve()
morph_signatures_dir.mkdir(exist_ok=True)


# In[4]:


shared_features = find_shared_features(selected_features_files)

# loading all single-cell profiles and updating it with the shared features
loaded_profiles_df = []
for single_cell_path in selected_features_files:

    # loading in single cell feature selected data
    single_cell_df = pd.read_parquet(single_cell_path)

    # split the features in order to retain Metadata
    # we ignore the morphology features because the profiles will be updated with the shared features
    meta, _ = data_utils.split_meta_and_features(single_cell_df)

    # append the updated profiles to the loaded_profiles_df
    loaded_profiles_df.append(single_cell_df[shared_features])

# Concatenate all the profiles
all_profiles_df = pd.concat(loaded_profiles_df, axis=0)

print(all_profiles_df.shape)
all_profiles_df.head()


# In[5]:


# next is to create both the reference and target datasets
# the reference will be the negative control, failing CF cells with DMSO treatment
# the target will be the positive control, healthy CF cells with DMSO treatment
reference_df = all_profiles_df.loc[(all_profiles_df["Metadata_cell_type"] == "failing") & (all_profiles_df["Metadata_treatment"] == "DMSO")]
target_df = all_profiles_df.loc[(all_profiles_df["Metadata_cell_type"] == "healthy") & (all_profiles_df["Metadata_treatment"] == "DMSO")]

# selecting one of the profiles to just split the metadata and morphology features
# they should be the same for all profiles
_, feats = data_utils.split_meta_and_features(reference_df)

print("Shape of the reference dataset:", reference_df.shape)
print("Shape of the target dataset:", target_df.shape)


# In[6]:


# The WKS test works for the given data because it accounts for imbalanced sample sizes
# between the reference and target datasets by applying weights to the cumulative
# distribution functions (CDFs). This ensures that the comparison of morphological
# features is fair and unbiased, regardless of differences in the number of cells
# in the two datasets.
off_morph_signatures, on_morph_signatures = weighted_ks_test(reference_df[feats], target_df[feats], p_thresh=0.05)


# In[7]:


# now lets save the signatures into a json file that also include the amount of features
signatures = {
    "off_morph_signatures": {
        "count": len(off_morph_signatures),
        "features": off_morph_signatures
    },
    "on_morph_signatures": {
        "count": len(on_morph_signatures),
        "features": on_morph_signatures
    }
}

# save into json file in the results directory
with open(morph_signatures_dir / "morph_signatures.json", "w") as f:
    json.dump(signatures, f, indent=4)
