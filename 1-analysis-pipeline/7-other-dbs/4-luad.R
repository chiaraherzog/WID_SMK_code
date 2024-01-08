# db for cancer (luad)
library(dplyr)

# Load in datasets (after eutopsQC) and find intersect of beta (if saved separately)
setwd("~/Dropbox/index-dev/SMK-index/16-manuscript/WID_SMK_code/")
cat("Load data ...")

load("1-analysis-pipeline/4-datasets/1-output/tcga.Rdata")
data <- tcga_data |> 
  dplyr::filter(smoking_history != "Unknown" & project == "TCGA-LUAD")

load("~/Documents/Work/data.nosync/tcga_sola/TCGA-LUAD/beta.Rdata")
beta <- beta[,match(rownames(data), colnames(beta))]
identical(colnames(beta), rownames(data))

load("~/Dropbox/index-dev/SMK-index/16-manuscript/WID_SMK_code/1-analysis-pipeline/2-output/WID_SMK_cpgs.Rdata")

beta <- beta[rownames(beta) %in% WID_SMK_cpgs$cg,]
cat("done\n")
cat("Start delta-beta estimation ...\n")

# Set up delta-beta
source("0-source/estimateDeltaBeta.R")
dir <- "~/Dropbox/index-dev/SMK-index/16-manuscript/WID_SMK_code/1-analysis-pipeline/7-other-dbs/4-output/"
estimateDeltaBeta(beta = beta,
                  pheno = data,
                  output = dir,
                  typevar = "type",
                  adjustment = c("age", "ic", "smoking_history"),
                  base = "Control",
                  plot = F)
