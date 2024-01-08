# Cervical cancer set -----------------------
library(GEOquery)
library(ChAMP)

options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")

setwd("<dir>") # raw local directory of your choice

# Download pheno
gse <- getGEO("GSE211668")
dir.create("GSE211668")
pheno <- pData(gse$GSE211668_series_matrix.txt.gz)
setwd("GSE211668/") # create a folder

# Download raw data
list.files()
getGEOSuppFiles("GSE211668", filter_regex = "GSE211668_UCL_Cervical_MatrixSignal.txt.gz", makeDirectory = F)
gunzip("GSE211668_UCL_Cervical_MatrixSignal.txt.gz")
file.copy(from = "GSE211668_UCL_Cervical_MatrixSignal.txt.gz",
          to = "./")
unlink("GSE211668/", recursive = T)

# Preprocess raw data from signal intensities
readLines("GSE211668_UCL_Cervical_MatrixSignal.txt", n = 1)
Mset <- minfi::readGEORawFile(filename = "GSE211668_UCL_Cervical_MatrixSignal.txt",
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
colnames <- strsplit(readLines("GSE211668_UCL_Cervical_MatrixSignal.txt", n = 1), "\t")[[1]]
select <- sort(grep("Detection Pval",colnames))
detP <- data.table::fread("GSE211668_UCL_Cervical_MatrixSignal.txt",
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

# rename cols
pheno$title <- gsub("^\\D+","",pheno$title)
pheno$title <- gsub("\\.", "", pheno$title)
ind <- match(pheno$title, colnames(beta)) # one sample missing -> removed sample.
pheno <- pheno[match(colnames(beta), pheno$title),]
identical(colnames(beta), pheno$title)

# Rename colnames of beta
colnames(beta) <- rownames(pheno)

# rename
pheno <- pheno |> 
  dplyr::mutate(type = case_when(`disease state:ch1` == "Normal" ~ "Control",
                                 `disease state:ch1` == "Tumour" ~ "Cancer"),
                id = stringr::str_split(title, "\ ", simplify = T)[,1]) |>
  dplyr::select(type, id, title, description)

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
library(here)
setwd(here())
source("0-source/WID_SMK.R")
res <- WID_SMK(beta)
pheno <- cbind(pheno, res)

# # cg
# pheno$cg05575921 <- as.numeric(beta["cg05575921",])
# 
# # smoking mrs
# source("0-source/smoking_mrs.R")
# pheno$smoking_mrs <- smoking_mrs(beta)

save(pheno, file = "1-analysis-pipeline/4-datasets/1-output/cervical-cancer.Rdata")


load("1-analysis-pipeline/4-datasets/1-output/cervical-cancer.Rdata")
identical(colnames(beta), rownames(pheno))

pheno <- pheno |> 
  dplyr::select(-contains("WID_SMK"))
source("0-source/WID_SMK.R")
res <- WID_SMK(beta)
pheno <- cbind(pheno, res)
