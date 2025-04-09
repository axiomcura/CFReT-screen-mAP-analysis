"""
This module contains all the functions responsible for performing the CFReT analysis in
this project.
"""

import pathlib

import pandas as pd
from copairs import map
from copairs.matching import assign_reference_index

from .data_utils import label_control_types, split_meta_and_features
from .io_utils import load_config


def calculate_dmso_map_batch_profiles(
    batched_profiles: dict,
    configs: str | pathlib.Path | dict,
    outdir_path: str | pathlib.Path,
    shuffled: bool | None = False,
) -> None:
    """ Calculate average precision (AP) and mean average precision (mAP) scores
    for DMSO-treated profiles across batches of replicate plates.

    This function processes batched phenotypic profiles to compute AP and mAP
    scores for wells treated with DMSO (used as control treatments). It saves
    the results as CSV files in the specified output directory. The function
    supports optional shuffling of data, indicated by a 'shuffled' label in
    the output filenames.

    Parameters
    ----------
    batched_profiles : dict
        A dictionary where keys are batch identifiers and values are `pd.DataFrame`
        objects containing phenotypic profiles with metadata and features.

    configs : str, pathlib.Path, or dict
        Configuration settings for the analysis. If a path is provided, the
        configuration is loaded using the `load_config` function.

    outdir_path : str or pathlib.Path
        The directory where the AP and mAP scores will be saved.

    shuffled : bool, optional
        Indicates whether the data should be shuffled. If True, the output file
        names will include a 'shuffled' label. Default is False.

    Returns
    -------
    None
        The function does not return any value. Results are saved as CSV files
        in the specified output directory.
    """
    # type checking
    if not isinstance(batched_profiles, dict):
        raise TypeError("'batched_profiles' must be a dictionary")
    # load configs if a path is provided and then load it
    if isinstance(configs, (str, pathlib.Path)):
        configs = load_config(configs)
    if not isinstance(configs, dict):
        raise TypeError("'configs' must be a dictionary")
    if not isinstance(outdir_path, (str, pathlib.Path)):
        raise TypeError("'outdir_path' must be either pathlib.Path or str")
    if not isinstance(shuffled, bool):
        raise TypeError("'shuffled' must be a boolean")

    # Load configs
    general_configs = configs["general_configs"]
    cntrl_copairs_ap_configs = configs["dmso_copairs_ap_configs"]
    cntrl_copairs_map_configs = configs["dmso_copairs_map_configs"]

    # setting controls
    # where negative controls indicates failing cells and positive controls indicates
    # healthy cells
    control_list = ["negative", "positive"]
    # setting shuffled labels
    shuffled_label = "original"
    if shuffled:
        shuffled_label = "shuffled"

    # Iterate over batches of loaded plate profiles
    for batch_id, profile in batched_profiles.items():
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
        dmso_profile = dmso_profile.copy().reset_index().rename(
            columns={"index": "original_index"}
        )

        # Iterate over control types to use them as references
        for ref_type in control_list:

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

                # set reference index for the targeted plate
                ref_col = "Metadata_reference_index"
                dmso_profile_w_target_plate = assign_reference_index(
                    df = dmso_profile_w_target_plate,
                    condition= f"Metadata_targeted and Metadata_Pathway == 'DMSO-{ref_type}'",
                    reference_col = ref_col,
                    default_value = -1,
                )

                dmso_profile_w_target_plate["Metadata_reference_control_type"] = ref_type

                # Split metadata and feature columns for analysis
                dmso_meta, dmso_feats = split_meta_and_features(
                    dmso_profile_w_target_plate
                )

                # Compute average precision (AP) scores for the current setup
                dmso_ap_scores = map.average_precision(
                    meta=dmso_profile_w_target_plate[dmso_meta],
                    feats=dmso_profile_w_target_plate[dmso_feats].values,
                    pos_sameby=cntrl_copairs_ap_configs["pos_sameby"] + [ref_col],
                    pos_diffby=[],
                    neg_sameby=cntrl_copairs_ap_configs["neg_sameby"],
                    neg_diffby=cntrl_copairs_ap_configs["neg_diffby"] + [ref_col],
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
            dmso_ap_save_path = (
                outdir_path
                / f"{batch_id}_{shuffled_label}_{ref_type}_ref_dmso_AP_scores.csv"
            )
            dmso_ap_scores.to_csv(dmso_ap_save_path, index=False)
            dmso_map_save_path = (
                outdir_path
                / f"{batch_id}_{shuffled_label}_{ref_type}_ref_dmso_mAP_scores.csv"
            )
            dmso_map_scores.to_csv(dmso_map_save_path, index=False)


def calculate_trt_map_batch_profiles(
    batched_profiles: dict,
    configs: str | pathlib.Path | dict,
    outdir_path: str | pathlib.Path,
    shuffled: bool | None = False,
) -> None:
    """Calculate treatment-specific mean average precision (mAP) and average precision (AP) scores
    for batched phenotypic profiles.

    This function processes batches of phenotypic profiles, calculates average precision (AP)
    and mean average precision (mAP) scores, and saves the results to the specified output directory.
    It supports multiple control conditions to evaluate the phenotypic differences between treated
    and control samples.

    Parameters
    ----------
    batched_profiles : dict
        A dictionary where keys are batch identifiers and values are `pd.DataFrame` objects
        containing phenotypic profiles with metadata and features.

    configs : str, pathlib.Path, or dict
        Configuration settings for the analysis. If a path is provided, the configuration is
        loaded using the `load_config` function.

    outdir_path : str, pathlib.Path
        The directory where the AP and mAP scores will be saved. If None, no files are saved.

    shuffled : bool, optional
        Indicates whether the data should be shuffled. If set to True, a 'shuffled' label
        will be appended to the output file names. Default is False and provides 'original'
        label.

    Returns
    -------
    None
        AP and mAP scores are saved in provided directory path.

    Raises
    ------
    TypeError
        If the input types of `batched_profiles`, `configs`, or `outdir_path` are invalid.
    """

    # type checking
    if not isinstance(batched_profiles, dict):
        raise TypeError("'batched_profiles' must be a dictionary")
    # load configs if a path is provided and then load it
    if isinstance(configs, (str, pathlib.Path)):
        configs = load_config(configs)
    if not isinstance(configs, dict):
        raise TypeError("'configs' must be a dictionary")
    if not isinstance(outdir_path, (str, pathlib.Path)):
        raise TypeError("'outdir_path' must be either pathlib.Path or str")
    if not isinstance(shuffled, bool):
        raise TypeError("'shuffled' must be a boolean")

    # setting in project configs
    general_configs = configs["general_configs"]
    copairs_ap_configs = configs["trt_copairs_ap_configs"]
    copairs_map_configs = configs["trt_copairs_map_configs"]

    # setting shuffled labels
    shuffled_label = "original"
    if shuffled:
        shuffled_label = "shuffled"

    # Define control conditions for the analysis
    # Each tuple specifies the control type, treatment, and associated cell state
    control_list = [("negative", "DMSO", "failing"), ("positive", "DMSO", "healthy")]

    # Iterate over each batch of loaded plate profiles
    for batch_id, profile in batched_profiles.items():
        # Analyze the profile for each control condition
        for control_type, control_treatment, cell_state in control_list:

            # Create a copy of the profile to preserve the original data
            profile = profile.copy()

            # Setting reference index for the control treatment
            ref_col = "Metadata_reference_index"
            profile = assign_reference_index(df = profile,
                                    condition = f"Metadata_Pathway == 'DMSO-{control_type}'",
                                    reference_col = ref_col,
                                    default_value = -1)

            # Separate metadata columns from feature columns for downstream calculations
            meta_columns, feature_columns = split_meta_and_features(profile)

            # Calculate average precision (AP) for the profile
            # Positive pairs are based on treatments with the same metadata
            # Negative pairs compare all DMSO-treated wells to all treatments
            replicate_aps = map.average_precision(
                meta=profile[meta_columns],
                feats=profile[feature_columns].values,
                pos_sameby=copairs_ap_configs["pos_sameby"] + [ref_col],
                pos_diffby=copairs_ap_configs["pos_diffby"],
                neg_sameby=[],
                neg_diffby=copairs_ap_configs["neg_diffby"] + [ref_col],
            )

            # Exclude wells treated with the control treatment (DMSO)
            replicate_aps = replicate_aps.loc[
                replicate_aps["Metadata_treatment"] != control_treatment
            ]

            # Save the calculated AP scores to a file for further analysis
            save_ap_path = (
                outdir_path
                / f"{batch_id}_{shuffled_label}_{control_type}_control_{cell_state}_{control_treatment}_trt_AP_scores.csv"
            )
            replicate_aps.to_csv(
                save_ap_path,
                index=False,
            )

            # Calculate mean average precision (mAP) from the AP scores
            replicate_maps = map.mean_average_precision(
                replicate_aps,
                sameby=copairs_map_configs["same_by"] + [ref_col],  # Grouping criteria for mAP
                null_size=copairs_map_configs["null_size"],  # Null distribution size
                threshold=copairs_map_configs["threshold"],  # Significance threshold
                seed=general_configs["seed"],  # Seed for reproducibility
            )

            # Save the mAP scores to a file for reporting
            save_map_path = (
                outdir_path
                / f"{batch_id}_{shuffled_label}_{control_type}_control_{cell_state}_{control_treatment}_trt_mAP_scores.csv"
            )
            replicate_maps.to_csv(
                save_map_path,
                index=False,
            )
