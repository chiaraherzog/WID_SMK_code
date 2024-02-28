# Step 1. Prepare data for computing delta Beta in epithelial and immune fraction in buccal, cervical, and blood samples.
# Control samples for which smoking history is not unknown are utilised.
# Recommend this script is run in command line/terminal (absolute paths are required)
# To compute delta-beta in Rproject, relative paths can be used
# Blood -----------------------------------
# EGAD IDs: tbc

# Load in datasets (after eutopsQC) and find intersect of beta (if saved separately)
cat('Starting EWAS: blood\n')

library(here)
here::i_am("1-analysis-pipeline/1-delta-beta/1c-blood-delta-beta.R")

# add links to beta matrices:
beta_3c <- '<path-to-beta1>'
beta_3c <- '~/Dropbox/data/3c-blood/beta_merged.Rdata'
# beta_brca_ds1 <- '<path-to-beta3>'
beta_brca_ds1 <- '<path-to-beta3>'

cat("Load data ...")
load("1-analysis-pipeline/0-data/data_int.Rdata")
data <- data_int |> 
  dplyr::filter(dataset == "discovery (blood)")

load(beta_3c)
beta_tmp1 <- beta_merged[,colnames(beta_merged) %in% rownames(data)]

load(beta_brca_ds1)
beta_tmp2 <- beta_merged[,colnames(beta_merged) %in% rownames(data)]

# intersect cpgs
intersect <- intersect(rownames(beta_tmp1), rownames(beta_tmp2))
beta <- cbind(beta_tmp1[intersect,],
              beta_tmp2[intersect,])

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
dir <- "1-analysis-pipeline/1-delta-beta/1c-blood/"

estimateDeltaBeta(beta = beta,
                  pheno = data,
                  output = dir,
                  typevar = "smoking_history",
                  adjustment = c("age", "ic"),
                  base = "Non-smoker",
                  names = c("db_lymphoid", "db_myeloid"))

cat("\n\nBlood EWAS done.\n")
