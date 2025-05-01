#!/usr/bin/env python

# # Generating DMSO consensus profiles

# In[1]:


import pathlib
import sys

import pandas as pd
from pycytominer import consensus

sys.path.append("../../../")
from utils import data_utils

# In[2]:


# setting profile path
concat_profile_path = pathlib.Path("../UMAP-aggregated-fs-profiles/results/concat_data/batch_1_concat_agg_fs.csv").resolve(strict=True)

# setting output path
# output_path = pathlib.Path("results/").resolve(strict=True)


# In[3]:


# load in aggregate profiles
agg_df = pd.read_csv(concat_profile_path)

# update aggregate profiles to only DMSO treated wells
dmso_agg_df = agg_df.loc[
    (agg_df["Metadata_control_type"] == "positive")
| (agg_df["Metadata_control_type"] == "negative")]

# split the metadata and morphology features
dmso_agg_meta, dmso_agg_feats = data_utils.split_meta_and_features(dmso_agg_df)

# display
print("Shape: ", dmso_agg_df.shape)
dmso_agg_df.head()


# In[4]:


consensus_df = consensus(profiles = dmso_agg_df,
                         replicate_columns=["Metadata_plate_barcode", "Metadata_plate_name", "Metadata_treatment"],
                         operation="median",
                         features=dmso_agg_feats,
)
