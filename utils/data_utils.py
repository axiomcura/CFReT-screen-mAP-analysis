"""
In this module, we have functions that has to do with providing utilities for
the loaded profiles.
"""

import pathlib

import pandas as pd
import pyarrow.parquet as pq
from pycytominer.cyto_utils import infer_cp_features


def split_meta_and_features(
    profile: pd.DataFrame,
    compartments: list[str] = ["Nuclei", "Cells", "Cytoplasm"],
    metadata_tag: bool | None = False,
) -> tuple[list[str], list[str]]:
    """Splits metadata and feature column names

    Parameters
    ----------
    profile : pd.DataFrame
        Dataframe containing image-based profile
    compartments : list, optional
        compartments used to generated image-based profiles, by default
        ["Nuclei", "Cells", "Cytoplasm"]
    metadata_tag : Optional[bool], optional
        indicating if the profiles have metadata columns tagged with 'Metadata_'
        , by default False

    Returns
    -------
    tuple[List[str], List[str]]
        Tuple containing metadata and feature column names
    """

    # identify features names
    features_cols = infer_cp_features(profile, compartments=compartments)

    # iteratively search metadata features and retain order if the Metadata tag is not added
    if metadata_tag is False:
        meta_cols = [
            colname
            for colname in profile.columns.tolist()
            if colname not in features_cols
        ]
    else:
        meta_cols = infer_cp_features(profile, metadata=metadata_tag)

    return (meta_cols, features_cols)

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


def shuffle_features(profile: pd.DataFrame, seed: int = 0) -> pd.DataFrame:
    """Shuffle the values in the feature columns of a DataFrame while preserving metadata columns.

    This function separates the metadata and feature columns from the input DataFrame, shuffles
    the values within each feature column independently using a specified random seed, and then
    concatenates the shuffled feature columns back with the metadata columns.

    Parameters:
    ----------
    profile : pd.DataFrame
        The input DataFrame containing both metadata and feature columns.
        Metadata columns are preserved, and feature columns are shuffled.
    seed : int, optional (default=0)
        The random seed for reproducibility. Ensures the same shuffle is applied
        for each column when the function is run with the same seed.

    Returns:
    -------
    pd.DataFrame
        A new DataFrame where feature columns are shuffled, and metadata columns remain unchanged.

    Raises:
    TypeError:
        Raised if a 'profiles' is not a pandas dataframe and when 'seed' is not an integer
    """
    # type checking
    if not isinstance(profile, pd.DataFrame):
        raise TypeError("'profile' must be a pandas dataframe")
    if not isinstance(seed, int):
        raise TypeError("'seed' must be an integer")

    # Split metadata and feature columns
    meta_cols, feat_cols = split_meta_and_features(profile)

    # Select only feature columns for shuffling
    feats_df = profile[feat_cols].copy()

    # Shuffle each feature column independently
    for colname in feats_df.columns:
        feats_df[colname] = (
            feats_df[colname].sample(frac=1, random_state=seed).reset_index(drop=True)
        )

    # Concatenate metadata and shuffled feature columns
    return pd.concat([profile[meta_cols], feats_df], axis=1)

def label_control_types(metadata_cell_type, heart_failure_type):
    """ Label control types based on cell metadata and heart failure type.

    This helper function for Pandas DataFrame operations assigns control labels
    ('positive' or 'negative') based on the combination of `metadata_cell_type`
    and `heart_failure_type`. It raises an error for unknown or invalid combinations.

    Parameters
    ----------
    metadata_cell_type : str
        The cell type metadata. Expected values are:
        - "healthy": Indicates a healthy cell type.
        - "failing": Indicates a failing cell type.

    heart_failure_type : str or None
        The type of heart failure associated with the cell type. Expected values are:
        - None: Used when `metadata_cell_type` is "healthy".
        - "dilated_cardiomyopathy": Used when `metadata_cell_type` is "failing".

    Returns
    -------
    str
        A control label based on the input combination:
        - "positive": For healthy cells with no heart failure type.
        - "negative": For failing cells with "dilated_cardiomyopathy".

    Raises
    ------
    ValueError
        If the combination of `metadata_cell_type` and `heart_failure_type` does not
        match the expected criteria.

    """
    if metadata_cell_type == "healthy" and heart_failure_type is None:
        return "positive"
    elif (
        metadata_cell_type == "failing"
        and heart_failure_type == "dilated_cardiomyopathy"
    ):
        return "negative"
    else:
        raise ValueError(f"Unknown combination added: {metadata_cell_type} {heart_failure_type}")
