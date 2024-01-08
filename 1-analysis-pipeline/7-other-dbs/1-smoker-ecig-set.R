# db for smoking (ecig dataset)
library(dplyr)

# Load in datasets (after eutopsQC) and find intersect of beta (if saved separately)
setwd("~/Dropbox/index-dev/SMK-index/16-manuscript/WID_SMK_code/")
cat("Load data ...")

load("~/Dropbox/index-dev/SMK-index/15-manuscript-analysis/3-output/vaping-pheno.Rdata")

data <- pheno |> 
  dplyr::mutate(sampletype = "saliva",
                dataset = "vaping set",
                smoking_history = case_when(type == "Non Smokers" ~ "Never smoker",
                                            type == "Smokers" ~ "Smoker")) |> 
  dplyr::filter(!is.na(smoking_history))

load("~/Documents/Work/data.nosync/vaping/norm.beta.RData")
beta <- norm.beta[,match(data$basename, colnames(norm.beta))]
identical(colnames(beta), as.character(data$basename))
rm(norm.beta);gc()

lowmi <- read.csv("0-source/lowMI_names.csv", header = F)
notonepicv2 <- read.csv("0-source/notonEPICv2_names.csv", header = F)
load("~/Dropbox/index-dev/SMK-index/16-manuscript/WID_SMK_code/1-analysis-pipeline/2-output/WID_SMK_cpgs.Rdata")

beta <- beta[rownames(beta) %in% WID_SMK_cpgs$cg,]
cat("done\n")
cat("Start delta-beta estimation ...\n")


# Set up delta-beta
source("0-source/estimateDeltaBeta.R")
dir <- "~/Dropbox/index-dev/SMK-index/16-manuscript/WID_SMK_code/1-analysis-pipeline/7-other-dbs/1-output/"
estimateDeltaBeta(beta = beta,
                  pheno = data,
                  output = dir,
                  typevar = "smoking_history",
                  adjustment = c("age", "ic"),
                  base = "Never smoker",
                  plot = F)
