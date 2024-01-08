# db for smokeless
library(dplyr)

# Load in datasets (after eutopsQC) and find intersect of beta (if saved separately)
setwd("~/Dropbox/index-dev/SMK-index/16-manuscript/WID_SMK_code/")
cat("Load data ...")

load("1-analysis-pipeline/4-datasets/1-output/snuff_tobacco.Rdata")

data <- pheno |> 
  dplyr::mutate(sampletype = "saliva",
                smoking_history = case_when(smoking_current == "Control" ~ "Never smoker",
                                            smoking_current == "Cigarette Smoker" ~ "Smoker")) |> 
  dplyr::filter(!is.na(smoking_history))

load("~/Documents/Work/data.nosync/GEO/GSE94876/beta_named.Rdata")
beta <- beta[,match(rownames(data), colnames(beta))]
identical(colnames(beta), rownames(data))

load("~/Dropbox/index-dev/SMK-index/16-manuscript/WID_SMK_code/1-analysis-pipeline/2-output/WID_SMK_cpgs.Rdata")

beta <- beta[rownames(beta) %in% WID_SMK_cpgs$cg,]
cat("done\n")
cat("Start delta-beta estimation ...\n")

# Set up delta-beta
source("0-source/estimateDeltaBeta.R")
dir <- "~/Dropbox/index-dev/SMK-index/16-manuscript/WID_SMK_code/1-analysis-pipeline/7-other-dbs/3-output/"
estimateDeltaBeta(beta = beta,
                  pheno = data,
                  output = dir,
                  typevar = "smoking_history",
                  adjustment = c("age", "ic"),
                  base = "Never smoker",
                  plot = F)
