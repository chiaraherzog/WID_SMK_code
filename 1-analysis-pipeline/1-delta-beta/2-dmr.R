library(here)
library(fs)
library(DMRcate)
library(dplyr)
library(tidyr)

# load scripts
source(here("0-source", "coerce_numeric.R"))
source(here("0-source", "saveDMRs.R"))
source(here("0-source", "plotDMRs.R"))

# set params/paths
here::i_am("1-analysis-pipeline/1-delta-beta/2-dmr.R")
db_path <- fs::path_expand("~/Dropbox")

# Data
cat("Load data ...")
load(here("1-analysis-pipeline/0-data/data_int.Rdata"))
data <- data_int |> 
  dplyr::filter(dataset == "discovery (buccal)")

load(file.path(db_path, "/data/3c-buccal/beta_merged.Rdata"))
beta_tmp1 <- beta_merged[,colnames(beta_merged) %in% rownames(data)]

load(file.path(db_path, "/data/4c-buccal/beta_merged.Rdata"))
beta_tmp2 <- beta_merged[,colnames(beta_merged) %in% rownames(data)]

load(file.path(db_path, "/data/brca-ds1/beta_merged.Rdata"))
beta_tmp3 <- beta_merged[,colnames(beta_merged) %in% rownames(data)]

# intersect cpgs
intersect <- intersect(rownames(beta_tmp1), intersect(rownames(beta_tmp2), rownames(beta_tmp3)))
beta <- cbind(beta_tmp1[intersect,],
              beta_tmp2[intersect,],
              beta_tmp3[intersect,])

# match sample names
intersect <- intersect(rownames(data), colnames(beta))
beta <- beta[,intersect]
data <- data[intersect,]
if(!identical(colnames(beta), rownames(data))){
  stop("Beta and pheno names not overlapping")
}

rm(beta_tmp1, beta_tmp2, beta_tmp3, intersect, data_int);invisible(gc())
cat("done\n")

# Remove lower MI CpG loci
lowmi <- read.csv("0-source/lowMI_names.csv", header = F)
notonepicv2 <- read.csv("0-source/notonEPICv2_names.csv", header = F)
beta <- beta[!rownames(beta) %in% c(lowmi$V1, notonepicv2$V1),]

# set params for DMRcate
type <- data$smoking_history
age <- data$age
ic <- data$ic
beta <- coerce_numeric(beta)

# Set up model matrix
design <- model.matrix(~type + age + ic)

# annotate
anno <- cpg.annotate("array", beta, arraytype = "EPIC",
                     analysis.type="differential",
                     coef = 2,
                     design=design)
dmrc <- dmrcate(anno)
devtools::install_version("dbplyr", version = "2.3.4")
packageVersion('dbplyr')
save(dmrc, file = here("1-analysis-pipeline", "2-output", "dmrc_buccal.Rdata"))

ranges <- extractRanges(dmrc, genome = "hg19")
ranges <- ranges[order(ranges$Stouffer),]
save(ranges, file = here("1-analysis-pipeline",  "1-delta-beta", "2-output",  "ranges_buccal.Rdata"))

pheno <- data
pheno$type <- ifelse(data$smoking_history == 'Non-smoker', 'Control', 'Smoker')

df <- saveDMRs(ranges= ranges,
               beta = beta,
               pheno = pheno,
               file = here("1-analysis-pipeline",  "1-delta-beta", "2-output", 'dmr_list'))

ranges_filtered <- ranges[match(paste0(df$seqnames, df$start, df$end),
                                paste0(seqnames(ranges), start(ranges), end(ranges))),]

pdf(file = here("1-analysis-pipeline", "1-delta-beta", "2-output", "dmr_buccal_discovery.pdf"),
    width = 9, height = 5)
plotDMRs(ranges = ranges_filtered,
         beta = beta,
         pheno = pheno,
         type = "buccal",
         n = 100)
dev.off()


# Cervical -----------------------
cat("Load data ...")
load("1-analysis-pipeline/0-data/data_int.Rdata")
data <- data_int |> 
  dplyr::filter(dataset == "discovery (cervical)")

load("~/Dropbox/data/3c/beta_merged.Rdata")
beta_tmp1 <- beta_merged[,colnames(beta_merged) %in% rownames(data)]

load("~/Dropbox/data/3c-ext-validation/beta_merged.Rdata")
beta_tmp2 <- beta_merged[,colnames(beta_merged) %in% rownames(data)]

load("~/Dropbox/data/brca-ds1/beta_merged.Rdata")
beta_tmp3 <- beta_merged[,colnames(beta_merged) %in% rownames(data)]

# intersect cpgs
intersect <- intersect(rownames(beta_tmp1), intersect(rownames(beta_tmp2), rownames(beta_tmp3)))
beta <- cbind(beta_tmp1[intersect,],
              beta_tmp2[intersect,],
              beta_tmp3[intersect,])

# match sample names
intersect <- intersect(rownames(data), colnames(beta))
beta <- beta[,intersect]
data <- data[intersect,]
if(!identical(colnames(beta), rownames(data))){
  stop("Beta and pheno names not overlapping")
}

rm(beta_tmp1, beta_tmp2, beta_tmp3, intersect, data_int);invisible(gc())
cat("done\n")

# Remove lower MI CpG loci
lowmi <- read.csv("0-source/lowMI_names.csv", header = F)
notonepicv2 <- read.csv("0-source/notonEPICv2_names.csv", header = F)
beta <- beta[!rownames(beta) %in% c(lowmi$V1, notonepicv2$V1),]

# set params for DMRcate
type <- data$smoking_history
age <- data$age
ic <- data$ic
beta <- coerce_numeric(beta)

# Set up model matrix
design <- model.matrix(~type + age + ic)

# annotate
anno <- cpg.annotate("array", beta, arraytype = "EPIC",
                     analysis.type="differential",
                     coef = 2,
                     design=design)
dmrc <- dmrcate(anno)
save(dmrc, file = here("1-analysis-pipeline","1-delta-beta", "2-output", "dmrc_cervical.Rdata"))
ranges <- extractRanges(dmrc, genome = "hg19")
ranges <- ranges[order(ranges$Stouffer),]
save(ranges, file = here("1-analysis-pipeline",  "1-delta-beta", "2-output",  "ranges_cervical.Rdata"))

pheno <- data
pheno$type <- ifelse(data$smoking_history == 'Non-smoker', 'Control', 'Smoker')

df <- saveDMRs(ranges= ranges,
               beta = beta,
               pheno = pheno,
               file = here("1-analysis-pipeline",  "1-delta-beta", "2-output", 'dmr_list_cervical'))

ranges_filtered <- ranges[match(paste0(df$seqnames, df$start, df$end),
                                paste0(seqnames(ranges), start(ranges), end(ranges))),]

pdf(file = here("1-analysis-pipeline", "1-delta-beta", "2-output", "dmr_cervical_discovery.pdf"),
    width = 9, height = 5)
plotDMRs(ranges = ranges_filtered,
         beta = beta,
         pheno = pheno,
         type = "cervical",
         n = 100)
dev.off()

# Blood -----------------------
cat("Load data ...")
load("1-analysis-pipeline/0-data/data_int.Rdata")
data <- data_int |> 
  dplyr::filter(dataset == "discovery (blood)")

load("~/Dropbox/data/3c-blood/beta_merged.Rdata")
beta_tmp1 <- beta_merged[,colnames(beta_merged) %in% rownames(data)]

load("~/Dropbox/data/brca-ds1/beta_merged.Rdata")
beta_tmp2 <- beta_merged[,colnames(beta_merged) %in% rownames(data)]

# intersect cpgs
intersect <- intersect(rownames(beta_tmp1), rownames(beta_tmp2))
beta <- cbind(beta_tmp1[intersect,],
              beta_tmp2[intersect,])

# match sample names
intersect <- intersect(rownames(data), colnames(beta))
beta <- beta[,intersect]
data <- data[intersect,]
if(!identical(colnames(beta), rownames(data))){
  stop("Beta and pheno names not overlapping")
}

rm(beta_tmp1, beta_tmp2, beta_tmp3, intersect, data_int);invisible(gc())
cat("done\n")

# set params for DMRcate
type <- data$smoking_history
age <- data$age
ic <- data$ic
beta <- coerce_numeric(beta)

# Set up model matrix
design <- model.matrix(~type + age + ic)

# annotate
anno <- cpg.annotate("array", beta, arraytype = "EPIC",
                     analysis.type="differential",
                     coef = 2,
                     design=design)
dmrc <- dmrcate(anno)
save(dmrc, file = here("1-analysis-pipeline","1-delta-beta", "2-output", "dmrc_blood.Rdata"))
ranges <- extractRanges(dmrc, genome = "hg19")
ranges <- ranges[order(ranges$Stouffer),]
save(ranges, file = here("1-analysis-pipeline",  "1-delta-beta", "2-output",  "ranges_blood.Rdata"))

pheno <- data
pheno$type <- ifelse(data$smoking_history == 'Non-smoker', 'Control', 'Smoker')

df <- saveDMRs(ranges= ranges,
               beta = beta,
               pheno = pheno,
               file = here("1-analysis-pipeline",  "1-delta-beta", "2-output", 'dmr_list_blood'))

ranges_filtered <- ranges[match(paste0(df$seqnames, df$start, df$end),
                                paste0(seqnames(ranges), start(ranges), end(ranges))),]

pdf(file = here("1-analysis-pipeline", "1-delta-beta", "2-output", "dmr_blood_discovery.pdf"),
    width = 9, height = 5)
plotDMRs(ranges = ranges_filtered,
         beta = beta,
         pheno = pheno,
         type = "buccal",
         n = length(ranges))
dev.off()



# Visualise -----------------
library(circlize)
buccal <- readRDS(here("1-analysis-pipeline",  "1-delta-beta", "2-output", 'dmr_list_buccal.Rds')) |> 
  dplyr::rename(value = diff) |> 
  dplyr::select(seqnames, start, end, value)

cervical <- readRDS(here("1-analysis-pipeline",  "1-delta-beta", "2-output", 'dmr_list_cervical.Rds')) |> 
  dplyr::rename(value = diff) |> 
  dplyr::select(seqnames, start, end, value)

blood <- readRDS(here("1-analysis-pipeline",  "1-delta-beta", "2-output", 'dmr_list_blood.Rds')) |> 
  dplyr::rename(value = diff) |> 
  dplyr::select(seqnames, start, end, value)

circos.initializeWithIdeogram(plotType = c("ideogram", "labels"),chromosome.index = paste0("chr", 1:22))
circos.genomicDensity(buccal, col = c("#FF000080"), track.height = 0.2)
circos.genomicDensity(cervical, col = c("#0000FF80"), track.height = 0.2)
circos.genomicDensity(blood, col = c("#f395f2"), track.height = 0.2)
circos.clear()


list <- list(buccal, cervical, blood)
circos.initializeWithIdeogram(plotType = c("ideogram", "labels"),chromosome.index = paste0("chr", 1:22))
circos.genomicDensity(list, col = c("#FF000080", "#0000FF80", "#f395f2"), track.height = 0.3)
circos.genomicDensity(cervical, col = c("#0000FF80"), track.height = 0.2)
circos.genomicDensity(blood, col = c("#f395f2"), track.height = 0.2)
circos.clear()
