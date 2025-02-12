#!/usr/bin/env python

# # Applying metrics

# In[1]:


import json
import pathlib
import sys
from itertools import product

import numpy as np
import pandas as pd
from pycytominer.cyto_utils import load_profiles
from scipy.stats import entropy
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import MinMaxScaler

# importing analysis utils
sys.path.append("../../utils")
from utils import data_utils

# ## helper functions

# In[2]:


def scale_data_to_non_negative_dist(
    target_df: pd.DataFrame,
    treated_df: pd.DataFrame,
    metadata: list[str],
    features: list[str],
    nan_handle: str | None = "mean",
    method: str | None = "shift",
):
    """Preprocesses two DataFrames by either applying mean shifting or MinMax scaling.

    This function preprocesses two DataFrames by either applying shifting or MinMax
    scaling to ensure non-negative values. The function first separates metadata and features
    from the target and treated DataFrames. It then handles NaN values by filling with the
    column-wise mean (other options include: "tiny" and "zero"). If the method is
    "shift", the function ensures non-negative values by shifting both DataFrames to a non-negative value.
    If the method is "minmax", the function scales both DataFrames to the [0, 1] range.
    The function returns two DataFrames with metadata and processed features.

    Parameters
    ----------
    target_df : pd.DataFrame
        The target DataFrame.
    treated_df : pd.DataFrame
        The treated DataFrame.
    metadata : list
        List of metadata columns to retain without processing.
    nan_handle : str, optional
        The method to handle NaN values. Default is "mean".
    features : list
        List of feature columns to preprocess.
    method : str, optional
        The preprocessing method, either "shift" or "minmax". Default is "shift".

    Returns
    -------
    tuple
        A tuple of two DataFrames (processed_target_df, processed_treated_df) with metadata and scaled features.
    """
    # Separate metadata and features
    target_metadata = target_df[metadata].copy()
    treated_metadata = treated_df[metadata].copy()

    target_features = target_df[features].copy()
    treated_features = treated_df[features].copy()

    # Handle NaN values by filling with column-wise mean
    if nan_handle == "mean":
        target_features = target_features.fillna(target_features.mean())
        treated_features = treated_features.fillna(treated_features.mean())
    elif nan_handle == "tiny":
        target_features = target_features.fillna(np.finfo(np.float32).tiny)
        treated_features = treated_features.fillna(np.finfo(np.float32).tiny)
    elif nan_handle == "zero":
        target_features = target_features.fillna(0)
        treated_features = treated_features.fillna(0)
    elif nan_handle == "impute":
        # create an imputer object
        imputer = SimpleImputer(strategy="mean")
        target_features = imputer.fit(target_features)
        treated_features = imputer.fit(treated_features)
    else:
        raise ValueError("Invalid nan_handle. Choose either 'mean', 'tiny', 'zero', or 'impute'.")

    # Apply preprocessing method to ensure non-negative values
    # shifting values to ensure non-negative values
    if method == "shift":
        min_target = target_features.min().min()
        min_treated = treated_features.min().min()
        shift_value = max(0, -min(min_target, min_treated))  # Find the shift value
        target_features += shift_value
        treated_features += shift_value
    # MinMax scaling: Scale both DataFrames to [0, 1] range
    elif method == "minmax":
        scaler = MinMaxScaler()
        target_features = pd.DataFrame(
            scaler.fit_transform(target_features),
            columns=features,
            index=target_features.index,
        )
        treated_features = pd.DataFrame(
            scaler.fit_transform(treated_features),
            columns=features,
            index=treated_features.index,
        )

    else:
        raise ValueError("Invalid method. Choose either 'shift' or 'minmax'.")

    # Concatenate metadata and processed features
    processed_target_df = pd.concat([target_metadata, target_features], axis=1)
    processed_treated_df = pd.concat([treated_metadata, treated_features], axis=1)

    return processed_target_df, processed_treated_df


# In[3]:


# setting path for data directory
data_dir = pathlib.Path("../data").resolve(strict=True)
results_dir = pathlib.Path("./results").resolve(strict=True)

# setting path for on and off morphology features
morph_sigs_path = (results_dir / "morph_signatures/morph_signatures.json").resolve(
    strict=True
)

# setting path for metadata containing cluster information
metadata_cluster_path = (results_dir / "cluster/metadata_w_clusters.csv").resolve(
    strict=True
)

# setting single-cell profile paths raise an error if no profiles are found
profile_paths = list(data_dir.glob("*sc_feature_selected.parquet"))
if len(profile_paths) == 0:
    raise FileNotFoundError("Profiles were not found at the given directory")

# setting results director for metric scores
metric_results_dir = (results_dir / "metric").resolve()
metric_results_dir.mkdir(exist_ok=True)


# In[4]:


# loading on and off morphological signatures
# off_morph_signatures: indicates morphological features that are not significantly
# associated with specific cellular state
# on_morph_signatures: indicates morphological features that are  significantly
# associated with specific cellular state
with open(morph_sigs_path) as content:
    on_off_sigs = json.load(content)
on_sigs = on_off_sigs["off_morph_signatures"]["features"]
off_sigs = on_off_sigs["on_morph_signatures"]["features"]

# Load the metadata with cluster information
meta_w_cluster_info_df = pd.read_csv(metadata_cluster_path)

# Establishing the feature space that is shared across all plates
shared_features = data_utils.find_shared_features(profile_paths)

# loading all single-cell profiles and updating it with the shared features
loaded_profiles_df = [
    load_profiles(single_cell_path)[shared_features]
    for single_cell_path in profile_paths
]
# Concatenate all the single_cell profiles and reset index and save original shape
all_profiles_df = pd.concat(loaded_profiles_df, axis=0).reset_index(drop=True)

# split the metadata and feature columns
all_meta, all_feats = data_utils.split_meta_and_features(all_profiles_df)

# updating original single-cell profile dataframe with clustering information
all_profiles_df = meta_w_cluster_info_df.merge(
    all_profiles_df[all_meta], on=all_meta, how="left"
).merge(all_profiles_df[all_feats], left_index=True, right_index=True, how="inner")

# Separate metadata and features
metadata_feats, morph_feats = data_utils.split_meta_and_features(all_profiles_df)

# display
print(all_profiles_df.shape)
all_profiles_df.head()


# Below, we calculate the Kullback-Leibler (KL) divergence to quantify the phenotypic effects of treated cells compared to the diseased state. The diseased state, in this case, is represented by the positive control, consisting of healthy CF cells treated with DMSO.
#
# The KL divergence is computed for two sets of morphological features:
# 1. **Off-target morphological features:** To evaluate unintended effects of the treatment.
# 2. **On-target morphological features:** To assess how closely the treated cells revert to the desired phenotypic state.
#
# These KL divergence scores provide insights into both the efficacy and specificity of the compound treatment.
#
# Finally, the computed scores are saved into a CSV file under the `results/metric` folder for further analysis and reporting.

# In[5]:


# parameters
metadata_treatments = "Metadata_cluster_family"
profile = None
target_name = "DMSO-healthy"
score_method = "mean"

# split the metadata and morphology feature
meta_cols, feat_cols = data_utils.split_meta_and_features(all_profiles_df)

# check if the selected metadata column contains the metadata_treatment that represents the control
if metadata_treatments not in meta_cols:
    raise ValueError(f"{metadata_treatments} is a metadata column that does not exist")

# separate the data to target and treated
target_df = all_profiles_df.loc[all_profiles_df[metadata_treatments] == target_name]
treated_df = all_profiles_df.loc[all_profiles_df[metadata_treatments] != target_name]

# Removing the -1 cluster label
# These clusters are known a "noisy" clusters and are not used in the analysis
target_df = target_df.loc[target_df["Metadata_cluster_label"] != -1]
treated_df = treated_df.loc[treated_df["Metadata_cluster_label"] != -1]

# After the noise removal, we can check the number of clusters in the target and treated data'
# and form combinations between the two datasets to compare
score_results = []
for trt_name, trt_df in treated_df.groupby(metadata_treatments):
    # Generate combinations of cluster labels between target and treated data
    clusters_to_compare = list(
        product(
            target_df["Metadata_cluster_label"].unique().tolist(),
            trt_df["Metadata_cluster_label"].unique().tolist(),
        )
    )

    # Calculate KL divergence between target and treated data for each cluster combination
    for target_cluster, treated_cluster in clusters_to_compare:
        # Filter the target and treated data based on the cluster labels pairs
        target_cluster_df = target_df.loc[
            target_df["Metadata_cluster_label"] == target_cluster
        ]
        treated_cluster_df = trt_df.loc[trt_df["Metadata_cluster_label"] == treated_cluster]

        # Next we need to convert morphological features into a probability distribution
        # Assumptions of KL divergence is that the data is a non-negative, probability
        # distribution and the data is continuous.
        target_cluster_df, treated_cluster_df = scale_data_to_non_negative_dist(
            target_df=target_cluster_df,
            treated_df=treated_cluster_df,
            metadata=meta_cols,
            features=feat_cols,
            nan_handle="mean",
            method="shift",
        )

        # Here we are separating the morphological feature spaces in order to generated
        # two scores for both the on and off morphological signatures
        off_target_cluster_df = target_cluster_df[off_sigs].reset_index(drop=True)
        off_treated_cluster_df = treated_cluster_df[off_sigs].reset_index(drop=True)

        on_target_cluster_df = target_cluster_df[on_sigs].reset_index(drop=True)
        on_treated_cluster_df = treated_cluster_df[on_sigs].reset_index(drop=True)

        # next we aggregate the probability distribution by taking the mean of the probability
        # This is a required step since KL divergence requires both distributions to be
        # the same shape
        off_target_cluster_df = off_target_cluster_df.mean(axis=0)
        off_treated_cluster_df = off_treated_cluster_df.mean(axis=0)
        on_target_cluster_df = on_target_cluster_df.mean(axis=0)
        on_treated_cluster_df = on_treated_cluster_df.mean(axis=0)

        # Next we calculate KL divergence for both on and off morphological signatures

        # Calculate KL divergence for off-morphology signatures
        off_kl_divergence = entropy(
            off_treated_cluster_df.values, off_target_cluster_df.values
        )

        # Calculate KL divergence for on-morphology signatures
        on_kl_divergence = entropy(
            on_treated_cluster_df.values, on_target_cluster_df.values
        )

        # now calculating the combined score by taking the mean of the two scores
        # This help by providing a single score that can be used to compare the two off
        # and on morphological signatures. A good way to summarize the overall difference
        # between the two distributions.
        if score_method == "mean":
            combined_score = np.mean([off_kl_divergence, on_kl_divergence])
        if score_method == "sum":
            combined_score = np.sum([off_kl_divergence, on_kl_divergence])

        # storing results in a dictionary
        results = {
            "treatment_name": trt_name,
            "target_cluster": target_cluster,
            "treated_cluster": treated_cluster,
            "off_kl_divergence": off_kl_divergence,
            "on_kl_divergence": on_kl_divergence,
            "combined_score": combined_score,
        }

        # appending the results to the score_results list
        score_results.append(results)

# converting the results to a dataframe
score_results = pd.DataFrame(score_results)


# In[6]:


score_results.to_csv(metric_results_dir / "kl_divergence_scores.csv", index=False)
