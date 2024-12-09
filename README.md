# CFReT-screen-mAP-analysis

## About

Cardiac fibrosis, a condition associated with significant morbidity and mortality, contributes to the global burden of cardiovascular diseases, which remain the leading cause of death worldwide, accounting for approximately 31% of all deaths.
This pathological process is characterized by the excessive deposition of extracellular matrix (ECM) components, especially collagen, resulting in the formation of scar tissue.
While ECM production is a vital part of the body’s response to injury, its overproduction disrupts tissue homeostasis and leads to structural changes that impair organ function.
In the heart, for example, excessive ECM accumulation causes stiffening of the myocardial tissue, compromising its ability to contract and relax efficiently.
This can lead to reduced cardiac output, impaired electrical conductivity, and an increased risk of heart failure or fatal arrhythmias.
Given its global impact and contribution to high mortality rates, there is an urgent need to develop advanced diagnostic tools and targeted therapeutic strategies to better understand, manage, and treat cardiac fibrosis.

Recently, we applied image-based profiling with the [CellPainting](https://www.nature.com/articles/nprot.2016.105)(CP) assay to investigate whether single-cell image-based profiles could be leveraged to identify failing heart cells (cardiac fibrosis), healthy heart cells, and heart cells treated with `drug_x`.
Our study demonstrated that by using image-based profiles with machine learning models, we could accurately distinguish between single cells derived from either heart failure or non-heart failure patients.
Additionally, we identified morphological features that were critical for the model’s ability to differentiate between heart states.
These findings suggest that specific cellular structures are strong indicators of heart cell health or failure.
The processing methods and results of this analysis can be accessed [here](https://github.com/WayScience/cellpainting_predicts_cardiac_fibrosis).

Building upon these results, we advanced our study to the next level by leveraging high-content screening (HCS) with CP assay. 
Specifically, we applied a library of 550 small-molecule compounds to identify those with the potential to reverse the effects of cardiac fibrosis.
This approach not only aims to pinpoint candidate compounds for therapeutic intervention but also seeks to understand the biological mechanisms underlying the reversal of fibrosis.

## Analytical approach

In this notebook, we utilize a statistical framework called [**Mean Average Precision (mAP)**](https://pmc.ncbi.nlm.nih.gov/articles/PMC11014546/) to evaluate and analyze image-based profiles generated through high-content screening (HCS) using the CP assay.
mAP is a robust and versatile metric designed to quantify the ability of a query profile to retrieve related profiles based on phenotypic similarity.

This framework enables us to assess two key aspects:

- **Phenotypic activity:** How distinct perturbations (e.g., compounds or treatments) are from controls.
- **Phenotypic consistency:** How similar profiles are within biologically related groups (e.g., compounds with similar mechanisms of action).

By leveraging mAP, we aim to identify potential compounds that can reverse the effects of cardiac fibroblasts, which are central to the pathology of cardiac fibrosis.
Beyond identifying these candidate compounds, mAP also helps us uncover the biological mechanisms driving this reversal, providing insights into the cellular and molecular processes that contribute to phenotypic recovery.
