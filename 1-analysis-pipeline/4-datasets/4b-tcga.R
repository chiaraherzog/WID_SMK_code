# TCGA Dataset -----------------------

suppressPackageStartupMessages(library(TCGAbiolinks))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(EpiDISH))

# samples
source("0-source/WID_SMK.R")
source("0-source/smoking_mrs.R")

projects <- c("TCGA-LUAD", "TCGA-LUSC")
dir <- "<dir>" # path to directory to save files in

# Access data:

for (id in projects){
  cat("Starting ", id, " (", which(projects==id), "/", length(projects), ").\n", sep = "")
  
  suppressMessages(query <- GDCquery(project = id,
                                     data.category = "DNA Methylation",
                                     platform = "Illumina Human Methylation 450",
                                     sample.type = c("Primary Tumor",
                                                     "Solid Tissue Normal"),
                                     data.type = "Methylation Beta Value"))
  dat <- query$results[[1]]
  query_pheno <- GDCquery_clinic(project = id,
                                 type = "clinical")
  ind <- match(dat$cases.submitter_id, query_pheno$submitter_id)
  tcga_pheno <- query_pheno[ind,]
  # save(tcga_pheno, dat, file = paste0(dir, "/pheno_", substr(id, 6, nchar(id)), ".Rdata")) # can save just in case
  
  cat("Starting download (", as.character(Sys.time()), ").\n", sep = "")
  setwd(paste0(dir, id))
  GDCdownload(query, method = "client", files.per.chunk = 10) # RIP internet
  suppressMessages(data <- GDCprepare(query))
  beta <- assays(data)[[1]]
  beta <- na.omit(beta)
  cat("Done downloading, saving file.\n")
  save(beta, file = paste0(dir, "/beta_", substr(id, 6, nchar(id)), ".Rdata")) # Save just in case need to recompute
  
  # delete files
  unlink("GDCdata/", recursive = T)
  suppressMessages(file.remove(c("gdc-client_v1.6.1_OSX_x64.zip",
                                 "gdc_manifest.txt",
                                 "gdc-client",
                                 "gdc_client_configuration.dtt")))
  
  # Format and compute
  
  # Link submitter ids between beta/phenos
  dat <- dat[match(substr(colnames(beta), 1, 16), dat$sample.submitter_id),]
  tcga_pheno <- tcga_pheno[match(substr(dat$sample.submitter_id, 1, 12), tcga_pheno$submitter_id),]
  rownames(tcga_pheno) <- colnames(beta)
  tcga_pheno$barcode1 <- colnames(beta)
  tcga_pheno$barcode2 <- dat$sample.submitter_id
  tcga_pheno$sample_type <- dat$sample_type
  
  lookup <- c(age = "age_at_index",
              stage_ajcc = "ajcc_pathologic_stage",
              tissue = "tissue_or_organ_of_origin",
              barcode3 = "submitter_id")
  
  cols <- c("project", "age", "type", "sample_type", "stage_ajcc",
            "smoking_frequency", "tobacco_smoking_status", "years_smoked", "pack_years_smoked", "tissue",
            "barcode1", "barcode2", "barcode3", "vital_status", "time_to_followup") 
  
  pheno <- tcga_pheno |> 
    mutate(type = case_when(sample_type %in% c("Solid Tissue Normal") ~ "Control",
                            sample_type %in% c("Primary Tumor") ~ "Cancer"),
           time_to_followup = case_when(is.na(days_to_death) ~ days_to_last_follow_up,
                                          !is.na(days_to_death) ~ days_to_death)) |> 
    rename(any_of(lookup)) %>%
    dplyr::select(any_of(cols))
  
  # Add tobacco history
  query_tobacco <- GDCquery(
    project = id,
    data.category = "Clinical",
    data.type = "Clinical Supplement", 
    data.format = "BCR Biotab"
  )
  
  GDCdownload(query_tobacco)
  clinical <- GDCprepare(query_tobacco)
  patient <- as.data.frame(clinical[[grep("clinical_patient", names(clinical))]])

  # Annotation (as here https://github.com/PoisonAlien/maftools/issues/410)
  # 1 = never smoker
  # 2 = current smoker
  # 3 = former smoker, more than 15 years ago
  # 4 = former smoker, less than 15 years age
  # 5 = former smoker
  
  # append smoking info
  ind <- match(pheno$barcode3, patient$bcr_patient_barcode)
  patient <- patient[ind,] |> 
    dplyr::mutate(smoking_history = case_when(tobacco_smoking_history_indicator == "1" ~ "Never smoker",
                                              tobacco_smoking_history_indicator == "2" ~ "Smoker",
                                              tobacco_smoking_history_indicator %in% c("3", "4", "5") ~ "Ex-smoker",
                                              tobacco_smoking_history_indicator %in% c("[Not Available]", "[Unknown]") ~"Unknown"))
  pheno$smoking_history <- patient$smoking_history
  
  # compute parameters
  out.l <- epidish(beta.m = beta,
                   ref.m = centEpiFibIC.m,
                   method = "RPC")$estF
  pheno$ic <- out.l[,3]
  
  # hepidish
  frac.m <- hepidish(beta.m = beta,
                     ref1.m = centEpiFibIC.m,
                     ref2.m = centBloodSub.m, 
                     h.CT.idx = 3,
                     method = 'RPC',
                     maxit = 500)
  pheno <- cbind(pheno,
                 hepidish_Epi=frac.m[,1],
                 hepidish_Fib=frac.m[,2],
                 hepidish_B=frac.m[,3],
                 hepidish_NK=frac.m[,4],
                 hepidish_CD4T=frac.m[,5],
                 hepidish_CD8T=frac.m[,6],
                 hepidish_Mono=frac.m[,7],
                 hepidish_Neutro=frac.m[,8],
                 hepidish_Eosino=frac.m[,9])
  
  # WID_SMK
  res <- WID_SMK(beta)
  pheno <- cbind(pheno, res)
  
  # cg
  pheno$cg05575921 <- as.numeric(beta["cg05575921",])
  
  # smoking mrs
  pheno$smoking_mrs <- smoking_mrs(beta)
  
  setwd(here())
  save(pheno, file = paste0("1-analysis-pipeline/4-datasets/1-output/pheno_", substr(id, 6, nchar(id)), ".Rdata"))
}
  
# Merge the two sets
load("1-analysis-pipeline/4-datasets/1-output/pheno_LUAD.Rdata")
luad <- pheno
load("1-analysis-pipeline/4-datasets/1-output/pheno_LUSC.Rdata")
lusc <- pheno

tcga_data <- rbind(luad, lusc)
save(tcga_data, file = "1-analysis-pipeline/4-datasets/1-output/tcga.Rdata")

# Remove intermediary files
unlink("1-analysis-pipeline/4-datasets/1-output/pheno_LUSC.Rdata")
unlink("1-analysis-pipeline/4-datasets/1-output/pheno_LUAD.Rdata")

# Access expression data:
for (id in project){
  query <- GDCquery(
    project = id,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification", 
    workflow.type = "STAR - Counts"
  )
  GDCdownload(query = query)
  data <- GDCprepare(query = query)
  save(data, file = file.path('~/Documents/Work/data.nosync/tcga_sola/', id, '/expression/data_expression_full.Rdata'))
  
}
