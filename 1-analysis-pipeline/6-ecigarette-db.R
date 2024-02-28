# db for e-cigarette use

cat("Run delta-beta for e-cigarette users ----------------------- \n")

library(dplyr)

# Load in datasets (after eutopsQC) and find intersect of beta (if saved separately)
cat("Load data ...")

# pheno <- '<set path to private pheno>'
pheno <- "~/Dropbox/index-dev/SMK-index/15-manuscript-analysis/3-output/vaping-pheno.Rdata"
load(pheno)

data <- pheno |> 
  dplyr::mutate(sampletype = "saliva",
                dataset = "vaping set",
                smoking_history = case_when(type == "Non Smokers" ~ "Never smoker",
                                            type == "Vapers" ~ "e-cigarette user")) |> 
  dplyr::filter(!is.na(smoking_history))

# beta <- '<set path to beta>'
beta <- "~/Documents/Work/data.nosync/vaping/norm.beta.RData"
load(beta)

beta <- norm.beta[,match(data$basename, colnames(norm.beta))]
identical(colnames(beta), as.character(data$basename))
rm(norm.beta);gc()

lowmi <- read.csv("0-source/lowMI_names.csv", header = F)
notonepicv2 <- read.csv("0-source/notonEPICv2_names.csv", header = F)

beta <- beta[!rownames(beta) %in% c(lowmi$V1, notonepicv2$V1),]
cat("done\n")
cat("Start delta-beta estimation ...\n")

# Set up delta-beta
source("0-source/estimateDeltaBeta.R")
dir <- "1-analysis-pipeline/6-output/"
estimateDeltaBeta(beta = beta,
                  pheno = data,
                  output = dir,
                  typevar = "smoking_history",
                  adjustment = c("age", "ic"),
                  base = "Never smoker")
