# imporant: cd into the project folder before running the code :)
# + make sure you have R installed.

# setup: install packages if not available
Rscript 0-source/install-libraries.R

# # EWAS: run using RAW data (update directories - if you do not have these, skip these steps)
# Rscript 1-analsis-pipeline/1-delta-beta/1a-buccal-delta-beta.R
# Rscript 1-analsis-pipeline/1-delta-beta/1b-cervical-delta-beta.R
# Rscript 1-analsis-pipeline/1-delta-beta/1c-blood-delta-beta.R

# Identify clusters
Rscript 1-analysis-pipeline/2-define-subsets.R

# Gene set enrichment for clusters
Rscript 1-analysis-pipeline/3-gene-sets-and-enrichment.R

# # Get additional datasets + perform QC: NOTE - update directories; these steps take long and output files are computed so do not need to be run
# Rscript 1-analysis-pipeline/4-datasets/4a-snuff-tobacco.R
# Rscript 1-analysis-pipeline/4-datasets/4b-tcga.R
# Rscript 1-analysis-pipeline/4-datasets/4c-cis-progression.R
# Rscript 1-analysis-pipeline/4-datasets/4d-cerv-cancer.R

# Get coefficients for IC correction
Rscript 1-analysis-pipeline/5-ic-correction.R

# # Delta-beta for e-cigarette users: NOTE - update directories - if you do not have these, skip these steps
# Rscript 1-analysis-pipeline/6-ecigarette-db.R

# # Run Expression analysis: NOTE - update directories to expression data - if you do not have these, skip these steps
# Rscript 1-analysis-pipeline/8-expression.R

# Prepare figures/figure panels (code in code chunks in Rmd; eval = F by default)
Rscript -e "rmarkdown::render('2-markdown/WID_SMK.Rmd', param=list(args=F))"

# Prepare supplementary figures/panels/tables (code in code chunks in Rmd; eval = F by default)
Rscript -e "rmarkdown::render('2-markdown/Supplementary_WID_SMK.Rmd', param=list(args=F))"