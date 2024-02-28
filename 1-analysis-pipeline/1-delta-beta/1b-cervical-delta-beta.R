# Step 1. Prepare data for computing delta Beta in epithelial and immune fraction in buccal, cervical, and blood samples.
# Control samples for which smoking history is not unknown are utilised.
# Recommend this script is run in command line/terminal (absolute paths are required)
# To compute delta-beta in Rproject, relative paths can be used
# Cervical -----------------------------------
# EGAD IDs: EGAD00010002079, EGAD00010002232, tbc

# Load in datasets (after eutopsQC) and find intersect of beta (if saved separately)
cat('Starting EWAS: buccal\n')

library(here)
here::i_am("1-analysis-pipeline/1-delta-beta/1b-cervical-delta-beta.R")

# add links to beta matrices:
# beta_3c <- '<path-to-beta1>'
beta_3c <- '~/Dropbox/data/3c/beta_merged.Rdata'
# beta_3cval <- '<path-to-beta2>'
beta_3cval <- '~/Dropbox/data/3c-ext-validation/beta_merged.Rdata'
# beta_brca_ds1 <- '<path-to-beta3>'
beta_brca_ds1 <- '~/Dropbox/data/brca-ds1/beta_merged.Rdata'

cat("Load data ...")
load("1-analysis-pipeline/0-data/data_int.Rdata")
data <- data_int |> 
  dplyr::filter(dataset == "discovery (cervical)")

load(beta_3c)
beta_tmp1 <- beta_merged[,colnames(beta_merged) %in% rownames(data)]

load(beta_3cval)
beta_tmp2 <- beta_merged[,colnames(beta_merged) %in% rownames(data)]

load(beta_brca_ds1)
beta_tmp3 <- beta_merged[,colnames(beta_merged) %in% rownames(data)]

# intersect cpgs
intersect <- intersect(rownames(beta_tmp1), intersect(rownames(beta_tmp2), rownames(beta_tmp3)))
beta <- cbind(beta_tmp1[intersect,],
              beta_tmp2[intersect,],
              beta_tmp3[intersect,])

# match sample names
intersect <- intersect(rownames(data), colnames(beta))
beta <- beta[,intersect]
data <- data[intersect,]
if(!identical(colnames(beta), rownames(data))){
  stop("Beta and pheno names not overlapping")
}

rm(beta_tmp1, beta_tmp2, beta_tmp3, intersect, data_int);invisible(gc())
cat("done\n")
cat("Start delta-beta estimation ...\n")


# Set up delta-beta
source("0-source/estimateDeltaBeta.R")
dir <- "1-analysis-pipeline/1-delta-beta/1b-cervical/"

estimateDeltaBeta(beta = beta,
                  pheno = data,
                  output = dir,
                  typevar = "smoking_history",
                  adjustment = c("age", "ic"),
                  base = "Non-smoker")

cat("\n\nCervical EWAS done.\n")
