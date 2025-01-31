"""
In this module, we have functions that has to do with providing utilities for
the loaded profiles.
"""

import pandas as pd
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
