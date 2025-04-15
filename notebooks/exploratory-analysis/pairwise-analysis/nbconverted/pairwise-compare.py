#!/usr/bin/env python

# In[1]:


import pathlib
import sys

import pandas as pd
from comparators.PearsonsCorrelation import PearsonsCorrelation
from comparison_tools.PairwiseCompareManager import PairwiseCompareManager
from pycytominer.cyto_utils import load_profiles

# loading project utils
sys.path.append("../../../")
from utils.data_utils import split_meta_and_features

# In[2]:


# set data path
data_path = pathlib.Path(
    "../UMAP-aggregated-fs-profiles/results/concat_data/batch_1_concat_agg_fs.csv"
).resolve(strict=True)

# setting output path
output_path = pathlib.Path("./results").resolve()
output_path.mkdir(exist_ok=True)


# In[3]:


# Create the "Metadata_plate_well" column using iloc
agg_profile = load_profiles(data_path)

# split the features:
metadata, features = split_meta_and_features(agg_profile)

# now only select DMSO profiles that are DMSO-positive and DMSO-negative
dmso_profiles = agg_profile.loc[
    (agg_profile["Metadata_treatment"] == "DMSO-positive")
    | (agg_profile["Metadata_treatment"] == "DMSO-negative")
]
dmso_profiles["Metadata_plate_well"] = dmso_profiles[["Metadata_plate_name", "Metadata_Well"]].apply(lambda row: f"{row[0]}_{row[1]}", axis=1)

# create a dataframe only containing pathway information and the the treatments
pathway_df = agg_profile[["Metadata_treatment", "Metadata_Pathway"]]


# ## Applying pairwise-compare to only the controls

# **Calculating Pairwise Correlation Scores for Controls**
#
# In this section, we calculate the pairwise correlation scores for both healthy controls (DMSO-positive) and failing controls (DMSO-negative). The goal is to assess whether batch effects are present within the same control groups. Specifically:
#
# - **Healthy Controls (DMSO-positive):** We compare healthy controls across all plates to evaluate their similarity.
# - **Failing Controls (DMSO-negative):** We compare failing controls across all plates to evaluate their similarity.
#
# This analysis helps identify potential inconsistencies or batch effects within the same control groups.

# In[4]:


# Comparing all positive controls (healthy cells) cross all plates to see if they are similar
dmso_pos_cntrl_comparer = PairwiseCompareManager(
    _df=dmso_profiles.loc[dmso_profiles["Metadata_treatment"] == "DMSO-positive"],
    _feat_cols=features,
    _different_columns=["Metadata_plate_well"],
    _same_columns=["Metadata_treatment"],
    _comparator=PearsonsCorrelation(),
)

# collecting all pairwise scores
pos_cntrl_pairwise_scores = dmso_pos_cntrl_comparer()


# In[5]:


# Comparing all negative controls (unhealthy cells) cross all plates to see if they are similar
dmso_neg_cntrl_comparer = PairwiseCompareManager(
    _df=dmso_profiles.loc[dmso_profiles["Metadata_treatment"] == "DMSO-negative"],
    _feat_cols=features,
    _different_columns=["Metadata_plate_well"],
    _same_columns=["Metadata_treatment"],
    _comparator=PearsonsCorrelation(),
)

# collecting all pairwise scores
neg_cntrl_pairwise_scores = dmso_neg_cntrl_comparer()


# In[6]:


# concatenate the scores
final_dmso_pairwise_scores = pd.concat(
    [
        pos_cntrl_pairwise_scores[
            [
                "pearsons_correlation",
                "Metadata_treatment__antehoc_group0",
                "Metadata_plate_well__posthoc_group0",
                "Metadata_plate_well__posthoc_group1",
            ]
        ],
        neg_cntrl_pairwise_scores[
            [
                "pearsons_correlation",
                "Metadata_treatment__antehoc_group0",
                "Metadata_plate_well__posthoc_group0",
                "Metadata_plate_well__posthoc_group1",
            ]
        ],
    ]
)

# update the columns names
final_dmso_pairwise_scores.columns = ["pearsons_correlation", "Metadata_treatment", "plate_well_0", "plate_well_1"]

# save to csv file
final_dmso_pairwise_scores.to_csv(
    output_path / "final_pairwise_scores.csv", index=False
)


# ## Calculating pair wise across treatments
#
# In this section of the notebook, we conduct pairwise comparisons across all treatments and specific controls. Two data frames are created:
#
# - **healthy_ref**: This dataset contains pairwise calculations comparing all treated failing cells to the healthy reference.
# - **failing_ref**: This dataset contains pairwise calculations comparing all treated failing cells to the failing reference.

# In[7]:


# calculating pairwise correlation between healthy control and treated failing wells
healthy_ref_trt_pairwise_comparer = PairwiseCompareManager(
    _df=agg_profile.loc[agg_profile["Metadata_treatment"] != "DMSO-negative"],
    _feat_cols=features,
    _different_columns=[
        "Metadata_control_type",
        "Metadata_treatment",
    ],
    _comparator=PearsonsCorrelation(),
)

# calculating pairwise correlation between failing control and treated failing wells
failing_ref_trt_pairwise_comparer = PairwiseCompareManager(
    _df=agg_profile.loc[agg_profile["Metadata_treatment"] != "DMSO-positive"],
    _feat_cols=features,
    _different_columns=["Metadata_control_type", "Metadata_treatment"],
    _comparator=PearsonsCorrelation(),
)

# collecting all pairwise scores
healthy_ref_trt_pairwise_scores = healthy_ref_trt_pairwise_comparer()
failing_ref_trt_pairwise_scores = failing_ref_trt_pairwise_comparer()


# In[8]:


# Select only the relevant columns and add a reference column for healthy controls
health_ref_pairwise_scores = healthy_ref_trt_pairwise_scores[
    ["pearsons_correlation", "Metadata_treatment__antehoc_group1"]
].copy()
health_ref_pairwise_scores["reference"] = "Healthy"

# Select only the relevant columns and add a reference column for failing controls
failing_ref_pairwise_scores = failing_ref_trt_pairwise_scores[
    ["pearsons_correlation", "Metadata_treatment__antehoc_group1"]
].copy()
failing_ref_pairwise_scores["reference"] = "Failing"

# Combine the healthy and failing control dataframes into a single dataframe
final_trt_pairwise_scores = (
    pd.concat([health_ref_pairwise_scores, failing_ref_pairwise_scores])
    .rename(columns={"Metadata_treatment__antehoc_group1": "Metadata_treatment"})
    .reset_index(drop=True)
)

# Merge the combined dataframe with pathway information from pathway_df
# This adds the Metadata_Pathway column to the final dataframe
final_trt_pairwise_scores = final_trt_pairwise_scores.merge(
    pathway_df, how="left", on="Metadata_treatment"
).drop_duplicates()

# Validate the correctness of the pathway information by comparing the merged data
# with the original pathway dictionary (pathway_dict)
final_trt_pairwise_pathways = dict(
    zip(
        final_trt_pairwise_scores["Metadata_treatment"],
        final_trt_pairwise_scores["Metadata_Pathway"],
    )
)
pathway_dict = dict(
    zip(
        pathway_df["Metadata_treatment"],
        pathway_df["Metadata_Pathway"],
    )
)
# Validate the correctness of the pathway information
for treatment, merged_pathway in final_trt_pairwise_pathways.items():
    original_pathway = pathway_dict.get(treatment)

    if original_pathway is None:
        raise KeyError("Key {treatment} not found in pathway_dict")
    elif merged_pathway != original_pathway:
        if pd.isna(merged_pathway) and pd.isna(original_pathway):
            continue
        raise ValueError(
            f"Pathway mismatch for key {treatment}: {merged_pathway} != {original_pathway}"
        )

# If there are NaN values in the pathway column, fill them with "No Pathway"
final_trt_pairwise_scores["Metadata_Pathway"] = final_trt_pairwise_scores[
    "Metadata_Pathway"
].apply(lambda x: x if pd.notna(x) else "No Pathway")

# Save the final dataframe with pairwise scores and pathway information to a CSV file
final_trt_pairwise_scores.to_csv(
    output_path / "final_trt_pairwise_scores.csv", index=False
)
