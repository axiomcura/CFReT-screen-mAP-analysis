#!/usr/bin/env python

# In[1]:


import pathlib
import sys

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from pycytominer.cyto_utils import load_profiles
from sklearn.cluster import DBSCAN
from sklearn.metrics.pairwise import distance_metrics

sys.path.append("../../utils")
from utils import data_utils

# ## helper functions

# In[2]:


class MergeError(Exception):
    """Raised when there's a merge error captured with the data"""

    def __init__(self, message):
        # leveraging the Exception class attributes to add message
        super().__init__(message)


def DBSCAN_single_cell_heterogeneity(
    profile: str | pathlib.Path | pd.DataFrame,
    metadata_groupby: str | list[str],
    min_points: int | None = 10,
    max_dist: float | None = 0.4,
    distance_method: str = "cosine",
    n_jobs: int | None = -1,
) -> pd.DataFrame:
    """Perform DBSCAN clustering on single-cell profiles to identify heterogeneity within specific groups.

    Parameters
    ----------
    profile : str, pathlib.Path, or pd.DataFrame
        Input data, which can be a file path to the profiles or a pandas DataFrame.
    metadata_groupby : str or list of str
        Column(s) in the metadata to group cells for clustering.
    min_points : int, optional
        Minimum number of samples in a neighborhood to form a cluster, by default 10.
    max_dist : float, optional
        Maximum distance between two samples to be considered in the same neighborhood, by default 0.4.
    distance_method : str, optional
        Distance metric to use for DBSCAN clustering, by default "cosine".
    n_jobs : int, optional
        Number of parallel jobs to run, by default -1 (use all available cores).

    Returns
    -------
    pd.DataFrame
        DataFrame with original metadata and additional columns for DBSCAN cluster labels.

    Raises
    ------
    TypeError
        If the input parameters have invalid types.
    ValueError
        If specified distance metric or metadata columns are invalid.
    MergeError
        If the number of rows changes after merging the metadata with cluster labels.
    """
    # Type checker
    if isinstance(profile, str) or isinstance(profile, pathlib.Path):
        profile = load_profiles(profile)
    if isinstance(metadata_groupby, str):
        metadata_groupby = [metadata_groupby]
    if not isinstance(metadata_groupby, list) and all(
        [isinstance(elm, str) for elm in metadata_groupby]
    ):
        raise TypeError("'groupby' must a list and each element should be a string")
    if not isinstance(profile, pd.DataFrame):
        raise TypeError("'profile' must be a pandas DataFrame.")
    if not isinstance(min_points, int):
        raise TypeError("'min_points' must be an integer.")
    if not isinstance(max_dist, float):
        raise TypeError("'max_dist' must be an integer.")
    if not isinstance(n_jobs, int):
        raise TypeError("'n_jobs' must be a integer.")
    if not isinstance(distance_method, str):
        raise TypeError("'distance_metric' must be a string.")

    # Loading all distance metrics from sklearn
    loaded_sklearn_distance_metrics = distance_metrics()

    # Selecting distance metric and assign it
    # If the metric does not exist, raise an error
    if distance_method not in loaded_sklearn_distance_metrics:
        raise ValueError(
            f"Invalid distance metric: {distance_method}"
            f"supported distance metrics: {list(loaded_sklearn_distance_metrics.keys())}"
        )

    # Split metadata and features columns
    meta_cols, feat_cols = data_utils.split_meta_and_features(profile)

    # Check if the metadata columns selected for groupby exists
    metadata_groupby_check = list(set(metadata_groupby) - set(meta_cols))
    if not len(metadata_groupby_check) == 0:
        raise ValueError(
            "These metadata features do not exist in the dataset:", meta_cols
        )

    # here we are iteration groups of cell that is dictated by the `metadata_groupby`
    # Each group of cells will go through DBSCAN which will cluster single-cell profiles within the
    # population. Each cluster indicates potential heterogenous states of this specific group at hte single cell level
    metadata_w_cluster_labels_df = []
    for group_name, group_profile in profile.groupby(by=metadata_groupby):
        # Create DBSCAN object and apply it to the profile
        # This will attempt to identify sub population of profiles within the selected group (e.g treatment)
        db_scan = DBSCAN(
            eps=max_dist,
            min_samples=min_points,
            metric=distance_method,
            algorithm="auto",
            n_jobs=n_jobs,
        ).fit(group_profile[feat_cols])

        # Adding clusters in to metadata df
        # Cluster represent a subpopulation of profile indicating difference in phenotypic profiles.
        # Metadata_cluster_family = indicates where the cluster was found
        # Metadata_cluster_label = cluster label (sub population) identify within the family
        group_meta_df = group_profile[meta_cols]
        group_meta_df["Metadata_cluster_family"] = group_name[0]
        group_meta_df["Metadata_cluster_label"] = db_scan.labels_

        # store updated metadata
        metadata_w_cluster_labels_df.append(group_meta_df)

    # concat all the metadata into one
    metadata_w_cluster_labels_df = pd.concat(metadata_w_cluster_labels_df)

    # next merge it to the original metadata in order to preserve order of the meta
    meta_df = profile[meta_cols]
    metadata_w_clusters = meta_df.merge(
        metadata_w_cluster_labels_df, on=meta_cols, how="inner"
    )

    # check if number of rows have changed
    if meta_df.shape[0] != metadata_w_clusters.shape[0]:
        raise MergeError(
            "The metadata and the metadata with cluster labels do not have the same number of rows after merging"
        )

    return metadata_w_clusters


def find_shared_features(profile_paths: list[str | pathlib.Path]) -> list[str]:
    """Find the shared features (columns) between the profiles in the provided list of
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


# In[3]:


# setting paths
data_dir = pathlib.Path("../data").resolve(strict=True)

# setting all single-cell profiles paths
profile_paths = list(data_dir.glob("*sc_feature_selected.parquet"))

# setting path where to save our results
results_dir = pathlib.Path("results").resolve()
results_dir.mkdir(exist_ok=True)

# making sub directories for the results
single_cell_hetero_dir = (results_dir / "cluster").resolve()
single_cell_hetero_dir.mkdir(exist_ok=True)


# In[4]:


shared_features = find_shared_features(profile_paths)

# loading all single-cell profiles and updating it with the shared features
loaded_profiles_df = []
for single_cell_path in profile_paths:
    # loading in single cell feature selected data
    single_cell_df = load_profiles(single_cell_path)

    # append the updated profiles to the loaded_profiles_df
    loaded_profiles_df.append(single_cell_df[shared_features])

# Concatenate all the profiles
all_profiles_df = pd.concat(loaded_profiles_df, axis=0)

print(all_profiles_df.shape)
all_profiles_df.head()


# In[5]:


# Vectorized update using np.where
# renaming the Metadata_treatment for cells that have been treated with DMSO
# if treatment is DMSO and cell type is healthy -> DMSO-healthy
# if treatment is DMSO and cell type is failing -> DMSO-failing
all_profiles_df["Metadata_treatment"] = np.where(
    (all_profiles_df["Metadata_treatment"] == "DMSO")
    & (all_profiles_df["Metadata_cell_type"] == "healthy"),
    "DMSO-healthy",
    np.where(
        (all_profiles_df["Metadata_treatment"] == "DMSO")
        & (all_profiles_df["Metadata_cell_type"] == "failing"),
        "DMSO-failing",
        all_profiles_df["Metadata_treatment"],
    ),
)

print(all_profiles_df.shape)
all_profiles_df.head()


# In[6]:


metadata_w_cluster_labels = DBSCAN_single_cell_heterogeneity(
    profile=all_profiles_df, metadata_groupby=["Metadata_treatment"]
)

metadata_w_cluster_labels.to_csv(
    single_cell_hetero_dir / "metadata_w_clusters.csv", index=False
)


# In[7]:


# Group by Metadata_cluster_family and count unique labels, excluding noise (-1)
unique_clusters = (
   metadata_w_cluster_labels.groupby("Metadata_cluster_family")[
       "Metadata_cluster_label"
   ]
   .apply(lambda labels: len(set(labels) - {-1}) if -1 in labels else len(set(labels)))
   .reset_index(name="Unique_Clusters")
)

# Handle cases where all clusters are noise (-1)
unique_clusters["Unique_Clusters"] = unique_clusters["Unique_Clusters"].apply(
   lambda x: x if x > 0 else 0
)

# Display the result
unique_clusters.to_csv(
   single_cell_hetero_dir / "n_clusters_per_treatment.csv", index=False
)


# In[8]:


# Group by Metadata_cluster_family and Metadata_cluster_label to count single-cell compositions
cluster_counts = (
    metadata_w_cluster_labels.groupby(
        ["Metadata_cluster_family", "Metadata_cluster_label"]
    )
    .size()
    .reset_index(name="Single_Cell_Count")
)

# Handle families with only noise (-1) and ensure they are represented
# We are trying single cells labeled as noise (-1) as zero due to not belonging into a cluster
# The objective here is to count all single-cell profiles to counted in a cluster (not -1)
cluster_counts["Cluster_label"] = cluster_counts["Metadata_cluster_label"].apply(
    lambda x: 0 if x == -1 else 1
)

# Display the result
cluster_counts.to_csv(
    single_cell_hetero_dir / "single_cells_counts_per_cluster.csv", index=False
)
