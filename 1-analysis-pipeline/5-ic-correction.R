# IC correction: raw data for models
library(dplyr)

load("1-analysis-pipeline/0-data/data_plot.Rdata")

# Correction for buccal
# Subset data
data <- data |> 
  dplyr::filter(dataset == "discovery set" & sampletype == "buccal" & !is.na(smoking_history)) |> 
  droplevels()

# Scores
scores <- c("WID_SMK_epithelial_hypoM",
            "WID_SMK_immune_hypoM",
            "WID_SMK_distal_epithelial_hyperM",
            "WID_SMK_proximal_epithelial_hyperM",
            "WID_SMK450_epithelial_hypoM",
            "WID_SMK450_immune_hypoM",
            "WID_SMK450_distal_epithelial_hyperM",
            "WID_SMK450_proximal_epithelial_hyperM",
            "cg05575921", "smoking_mrs")

models_nsmk <- list()
models_exsmk <- list()
models_smk <- list()

for (i in scores){
  m <- lm(data[data$smoking_history=="Never smoker",][[i]] ~ ic, data[data$smoking_history=="Never smoker",])
  models_nsmk[[i]] <- m
  
  m <- lm(data[data$smoking_history=="Ex-smoker",][[i]] ~ ic, data[data$smoking_history=="Ex-smoker",])
  models_exsmk[[i]] <- m
  
  m <- lm(data[data$smoking_history=="Smoker",][[i]] ~ ic, data[data$smoking_history=="Smoker",])
  models_smk[[i]] <- m
}

save(models_nsmk, models_exsmk, models_smk, file = "1-analysis-pipeline/5-output/coef.Rdata")

