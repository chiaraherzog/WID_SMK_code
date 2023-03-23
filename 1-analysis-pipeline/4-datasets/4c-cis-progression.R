# CIS progression dataset -----------------------
library(GEOquery)
library(ChAMP)

options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")

setwd("<dir>") # raw local directory of your choice

# Download pheno
gse <- getGEO("GSE108123")
dir.create("GSE108123")
pheno <- pData(gse$GSE108123_series_matrix.txt.gz)
setwd("GSE108123/") # create a folder

# Download raw data
list.files()
getGEOSuppFiles("GSE108123", filter_regex = "GSE108123_matrix_signal_intensities.txt.gz",
                fetch_files = T)
gunzip("GSE108123_matrix_signal_intensities.txt")
file.copy(from = "GSE108123_matrix_signal_intensities.txt",
          to = "./")
unlink("GSE108123/", recursive = T)

# Preprocess raw data from signal intensities
readLines("GSE108123_matrix_signal_intensities.txt", n = 1)
Mset <- minfi::readGEORawFile(filename = "GSE108123_matrix_signal_intensities.txt",
                       sep = "\t",
                       Uname = "Unmethylated Signal",
                       Mname = "Methylated Signal",
                       array = "IlluminaHumanMethylation450k",
                       annotation = "ilmn12.hg19",
                       showProgress = TRUE)

# QC/pipeline. Define global thresholds
INTENSITY_THRESHOLD <- 9.5     # minimum median intensity required
DETECTION_P_THRESHOLD <- 0.01  # maximum detection p-value
FAILED_PROBE_THESHOLD <- 0.1   # maximum proportion of failed probes per sample
qc <- getQC(Mset)
plotQC(qc) # some low intensity samples but still above 9.5 threshold
# Filter any samples with median (un)methylated intensity less than threshold
low_intensity_samples <- rownames(qc)[qc$mMed<INTENSITY_THRESHOLD | qc$uMed<INTENSITY_THRESHOLD] # no samples with intensity below threshold

# Read in detP
colnames <- strsplit(readLines("GSE108123_matrix_signal_intensities.txt", n = 1), "\t")[[1]]
select <- sort(grep("Detection Pval",colnames))
detP <- data.table::fread("GSE108123_matrix_signal_intensities.txt",
              sep = "\t",
              select = select)

# Filter samples with too many failed probes
failed_samples <- colnames(detP)[colSums(detP>DETECTION_P_THRESHOLD) > (nrow(detP) * FAILED_PROBE_THESHOLD)] # No samples failing QC
rm(detP, colnames, select, qc);gc()

# Convert to RatioSet and then beta
RSet <- ratioConvert(Mset, what = "both", keepCN = TRUE)
rm(Mset);gc()
beta <- getBeta(RSet)
beta <- na.omit(beta)
rm(RSet);gc()
densityPlot(beta) # Plot densities

# ChAMP Normalisation
norm <- champ.norm(beta, arraytype = "450k", method = "BMIQ", cores = 5)
beta <- norm
densityPlot(beta)

# Additional QC steps for visualisation can be conducted (eg PCA, UMAP)

# Combine with pheno data
library(here)
setwd(here())

# FORMAT
# Rename columns
ind <- match(colnames(beta), pheno$title)
colnames(beta) <- rownames(pheno)
identical(rownames(pheno), colnames(beta))

# Information on progression not great in GEO annotation - use previously published supplementary files from the research group

# Match names to published file in Supplementary https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7611527/
dat <- readxl::read_xlsx("0-source/EMS128586-supplement-Table_S1.xlsx", sheet = 1) |> 
  janitor::clean_names() |> 
  dplyr::filter(methylation == "TRUE") |> 
  dplyr::mutate(id = gsub("X", "", sample_number_meth))

ind <- match(pheno$title, dat$id)
dat <- dat[ind,]

# get biopsy date from supplementary table 1 in https://doi.org/10.1038/s41591-018-0323-0 
dat2 <- readxl::read_xlsx("0-source/41591_2018_323_MOESM1_ESM.xlsx", sheet = 1, skip =1) |> 
  janitor::clean_names() |> 
  dplyr::mutate(patient_id = as.character(readr::parse_number(stringr::str_split(patient_number, "_", simplify = T)[,1])),
                biopsy_id = stringr::str_split(patient_number, "_", simplify = T)[,2]) |>
  dplyr::rename(type = outcome) |> 
  dplyr::filter(methylation == "YES")

dat <- dat |> 
  dplyr::left_join(dat2, by = c("patient_number", "biopsy_site")) |> 
  dplyr::rename(GEO_id = id) |> 
  as.data.frame()

rownames(dat) <- rownames(pheno)
pheno <- dat
# identical(rownames(dat), colnames(beta))

# compute parameters
library(EpiDISH)
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
source("0-source/WID_SMK.R")
res <- WID_SMK(beta)
pheno <- cbind(pheno, res)

# cg
pheno$cg05575921 <- as.numeric(beta["cg05575921",])

# smoking mrs
source("0-source/smoking_mrs.R")
pheno$smoking_mrs <- smoking_mrs(beta)

pheno <- pheno |>
  dplyr::mutate(type = ifelse(type == "Controls", "Control", type),
                type = factor(type, levels = c("Control", "Regression", "Progression")),
                age = as.numeric(age_at_bronchoscopy),
                sampletype = "lung",
                sex = gender,
                smoking_history = case_when(smoking_status == "Former" ~ "Ex-smoker",
                                            smoking_status == "Current" ~ "Smoker",
                                            smoking_status == "Unknown" ~ "Unknown"))

save(pheno, file = "1-analysis-pipeline/4-datasets/1-output/cis_progression.Rdata")

