# Datasets and computation of scores

# WID_SMK_cpgs (source("1-analysis-pipeline/2-output/WID_SMK_cpgs.Rdata")) are used in the WID_SMK function to compute mean scores.
# Use is res <- WID_SMK(beta)

# The following datasets were used:
# 1. EUTOPS FORECEE dataset (deposited in EGA), discovery
# 2. NSHD Buccal and Blood matched sample methylation, validation (data not publicly available due to restricted use and confidentiality)
# 3. Vaping dataset by Richmond et al. 2021 - available upon request (original authors)
# 4. GEO moist snuff tobacco dataset
# 5. TCGA-LUAC And -LUSC (accessed via TCGA Biolinks package)
# 6. GEO CIS Progression and Regression dataset (GEOquery)

# Example code for download and preprocessing of datasets 4-6 is shown in scripts 4a-c.