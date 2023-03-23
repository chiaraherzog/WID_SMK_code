# Snuff tobacco dataset -----------------------
library(GEOquery)
library(ChAMP)

options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")

setwd("<dir>") # raw local directory of your choice

# Download pheno
gse <- getGEO("GSE94876")
dir.create("GSE94876")
pheno <- pData(gse$GSE94876_series_matrix.txt.gz)
setwd("GSE94876/") # create a folder

# Download raw data
list.files()
getGEOSuppFiles("GSE94876", filter_regex = "GSE94876_matrix-signal-intensities.txt.gz",
                fetch_files = T)
gunzip("GSE94876/GSE94876_matrix-signal-intensities.txt.gz")
file.copy(from = "GSE94876/GSE94876_matrix-signal-intensities.txt",
          to = "./")
unlink("GSE94876/", recursive = T)

# Preprocess raw data from signal intensities
readLines("GSE94876_matrix-signal-intensities.txt", n = 1)
Mset <- readGEORawFile(filename = "GSE94876_matrix-signal-intensities.txt",
                       sep = "\t",
                       Uname = "Signal_A",
                       Mname = "Signal_B",
                       array = "IlluminaHumanMethylation450k",
                       annotation = "ilmn12.hg19",
                       showProgress = TRUE,
                       row.names = 2)

# QC/pipeline. Define global thresholds
INTENSITY_THRESHOLD <- 9.5     # minimum median intensity required
DETECTION_P_THRESHOLD <- 0.01  # maximum detection p-value
FAILED_PROBE_THESHOLD <- 0.1   # maximum proportion of failed probes per sample
qc <- getQC(Mset)
plotQC(qc) # QC looks ok
# Filter any samples with median (un)methylated intensity less than threshold
low_intensity_samples <- rownames(qc)[qc$mMed<INTENSITY_THRESHOLD | qc$uMed<INTENSITY_THRESHOLD] # no samples with intensity below threshold

# Read in detP
colnames <- strsplit(readLines("GSE94876_matrix-signal-intensities.txt", n = 1), "\t")[[1]]
select <- sort(grep("Detection.Pval",colnames))
detP <- data.table::fread("GSE94876_matrix-signal-intensities.txt",
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
densityPlot(beta) # Plot densities - doesn't look too bad

# ChAMP Normalisation
norm <- champ.norm(beta, arraytype = "450k", method = "BMIQ", cores = 5)
beta <- norm
densityPlot(beta) # looks better now

# Additional QC steps for visualisation can be conducted (eg PCA, UMAP)
library(here)
setwd(here())

# FORMAT
# Match beta to pheno - IDs are buried in colnames
colnms <- substr(colnames(beta), 2, 9)
ind <- match(pheno$title, colnms)
beta <- beta[,ind]
colnames(beta) <- rownames(pheno)
pheno$age <- as.numeric(pheno$`age:ch1`)
pheno$smoking_current <- ifelse(as.character(pheno$`class:ch1`) == "Non-Tobacco User", "Control", as.character(pheno$`class:ch1`))
pheno$sex <- "male"
pheno$sampletype <- "buccal"

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
  dplyr::select(sex, age, smoking_current, sampletype, ic, grep("hepidish", colnames(pheno)), grep("WID_SMK", colnames(pheno)),
                cg05575921, smoking_mrs)

save(pheno, file = "1-analysis-pipeline/4-datasets/1-output/snuff_tobacco.Rdata")

