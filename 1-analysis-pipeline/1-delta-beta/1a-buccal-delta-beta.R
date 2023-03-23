# Step 1. Prepare data for computing delta Beta in epithelial and immune fraction in buccal, cervical, and blood samples.
# Control samples for which smoking history is not unknown are utilised.
# Recommend this script is run in command line/terminal (absolute paths are required)
# To compute delta-beta in Rproject, relative paths can be used
# Buccal -----------------------------------
# EGAD IDs: EGAD00010002080, EGAD00010002231

# Load in datasets (after eutopsQC) and find intersect of beta (if saved separately)
setwd("~/Dropbox/index-dev/SMK-index/16-manuscript/WID_SMK_code/")
cat("Load data ...")
load("1-analysis-pipeline/0-data/data_int.Rdata")
data <- data_int |> 
  dplyr::filter(dataset == "discovery (buccal)")

load("~/Dropbox/data/3c-buccal/beta_merged.Rdata")
beta_tmp1 <- beta_merged[,colnames(beta_merged) %in% rownames(data)]

load("~/Dropbox/data/4c-buccal/beta_merged.Rdata")
beta_tmp2 <- beta_merged[,colnames(beta_merged) %in% rownames(data)]

load("~/Dropbox/data/brca-ds1/beta_merged.Rdata")
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
dir <- "~/Dropbox/index-dev/SMK-index/16-manuscript/WID_SMK_code/1-analysis-pipeline/1-delta-beta/1a-buccal/"
estimateDeltaBeta(beta = beta,
                  pheno = data,
                  output = dir,
                  typevar = "smoking_history",
                  adjustment = c("age", "ic"),
                  base = "Non-smoker")
