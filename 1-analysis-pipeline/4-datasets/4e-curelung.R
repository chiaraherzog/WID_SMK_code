# CURELUNG dataset -----------------------
library(GEOquery)
library(ChAMP)
library(minfi)

options(timeout = max(300, getOption("timeout")))
options(download.file.method.GEOquery = "wget")

setwd("<dir>") # raw local directory of your choice

# Download pheno and beta
gse <- getGEO("GSE39279")
pheno <- pData(gse[[1]])
beta <- exprs(gse[[1]])

# Check beta and normalise
densityPlot(beta)
beta <- na.omit(beta)
norm <- champ.norm(beta, arraytype = "450k", method = "BMIQ", cores = 4)
densityPlot(norm)
beta <- norm
identical(colnames(beta), rownames(pheno))
save(beta, file = "beta.Rdata")

# load("~/Dropbox/index-dev/SMK-index/13-tissue-comparison_all/8-output/pheno-curelung.Rdata")

# Format pheno 
library(dplyr)
pheno_new <- pheno |> 
  dplyr::mutate(cancer = ifelse(grepl("adeno", source_name_ch1), "LUAD", "LUSC"),
                age = as.numeric(`age:ch1`),
                stage_T = case_when(grepl("unknown", `tnm:ch1`) ~ NA,
                                    grepl("T1M1", `tnm:ch1`) ~ 1,
                                    TRUE ~ as.numeric(substr(`tnm:ch1`, 2,2))),
                stage_N = case_when(grepl("unknown|T1M1", `tnm:ch1`) ~ NA,
                                    TRUE ~ as.numeric(substr(`tnm:ch1`, 4,4))),
                stage_M = case_when(grepl("unknown|T1M1|Mx", `tnm:ch1`) ~ NA,
                                    TRUE ~ as.numeric(substr(`tnm:ch1`, 6,6))),
                stage = `Stage:ch1`,
                smoking_current = ifelse(pheno$`smoker:ch1` == "yes", "Current smoker",
                                         ifelse(pheno$`smoker:ch1` == "no", "Non-smoker", NA)),
                recurrence = `recurrence:ch1`,
                recurrence_time = as.numeric(`time_rec (years):ch1`),
                sex = gsub("gender: ", "", characteristics_ch1.8),
                adjuvant = ifelse(`adjuvant chemotherapy:ch1` == "yes", "yes",
                                  ifelse(`adjuvant chemotherapy:ch1` == "no", "no", NA))
                ) |> 
  dplyr::select(age, sex, smoking_current, cancer, stage, stage_T, stage_N, stage_M, recurrence, recurrence_time, adjuvant)

identical(rownames(pheno_new), colnames(beta))

pheno <- pheno_new

# compute age estimate
library(WIDclocks)
res <- WIDclocks::WID_clocks(beta)
pheno <- cbind(pheno, as.data.frame(res))

# compute ic estimates
library(EpiDISH)
out.l <- epidish(beta.m = beta,
                 ref.m = centEpiFibFatIC.m,
                 method = "RPC")$estF
pheno <- cbind(pheno, out.l)

# hepidish
frac.m <- hepidish(beta.m = beta,
                   ref1.m = centEpiFibFatIC.m,
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

# compute WID smoking scores
library(WID.smk)
res <- WID_SMK(beta)
pheno <- cbind(pheno, res)

# EpiSmoker
library(EpiSmokEr)
res <- epismoker(beta, samplesheet = pheno)
pheno$smoking_history_est <- res$PredictedSmokingStatus # most predicted to be smokers? likely due to miscalibration of package for adipose.
# table(pheno$smoking_history_est, pheno$smoking_current)

# Principal components
sds <- rowSds(beta)
names(sds) <- rownames(beta)
sds <- sds[order(sds, decreasing = T)][1:10000]
tmp <- beta[names(sds),]
pc <- prcomp(t(tmp))
umap <- uwot::umap(t(tmp))

pheno$pc1 <- pc$x[,1]
pheno$pc2 <- pc$x[,2]

pheno$umap1 <- umap[,1]
pheno$umap2 <- umap[,2]

save(pheno, file = "~/Dropbox/index-dev/SMK-index/16-manuscript/WID_SMK_code/1-analysis-pipeline/4-datasets/pheno_curelung.Rdata")


pheno |> 
  ggplot(aes(x = IC,
             y = WID_SMK450_immune_hypoM,
             colour = smoking_current)) +
  geom_point() +
  geom_smooth(method = 'lm')
