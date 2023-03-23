# 2. Define subsets of CpGs
library(here)
library(dplyr)

# Load in data and keep CpGs with a padj < 0.05 in at least one set
db_buccal <- readRDS("1-analysis-pipeline/1-delta-beta/1a-buccal/delta-beta.Rds") |> 
  as.data.frame() |> 
  mutate(padj = p.adjust(type),
         sig = ifelse(padj < 0.05, "yes", "no")) |> 
  tibble::rownames_to_column("cg")

db_cervical <- readRDS("1-analysis-pipeline/1-delta-beta/1b-cervical/delta-beta.Rds") |> 
  as.data.frame() |> 
  mutate(padj = p.adjust(type),
         sig = ifelse(padj < 0.05, "yes", "no")) |> 
  tibble::rownames_to_column("cg")

db_blood <- readRDS("1-analysis-pipeline/1-delta-beta/1c-blood/delta-beta.Rds") |> 
  as.data.frame() |> 
  mutate(padj = p.adjust(type),
         sig = ifelse(padj < 0.05, "yes", "no")) |> 
  tibble::rownames_to_column("cg")

intersect <- intersect(db_buccal$cg, intersect(db_cervical$cg, db_blood$cg))

# Significant dataframe
sig <- unique(c(db_blood[db_blood$sig=="yes",]$cg,
                db_buccal[db_buccal$sig=="yes",]$cg,
                db_cervical[db_cervical$sig=="yes",]$cg))
sig <- sig[sig %in% intersect]

rownames(db_buccal) <- db_buccal$cg
rownames(db_blood) <- db_blood$cg
rownames(db_cervical) <- db_cervical$cg

# Create dataframe for heatmap (Figure 2)
db <- data.frame(buccal_epi = db_buccal[sig,]$db_epithelial,
                 buccal_imm = db_buccal[sig,]$db_immune,
                 buccal_p = db_buccal[sig,]$type,
                 buccal_padj =db_buccal[sig,]$padj,
                 cervical_epi = db_cervical[sig,]$db_epithelial,
                 cervical_imm = db_cervical[sig,]$db_immune,
                 cervical_p = db_cervical[sig,]$type,
                 cervical_padj =db_cervical[sig,]$padj,
                 blood_mye = db_blood[sig,]$db_myeloid,
                 blood_lym = db_blood[sig,]$db_lymphoid,
                 blood_p = db_blood[sig,]$type,
                 blood_padj =db_blood[sig,]$padj,
                 cg = sig)

# Format the matrix for plotting and identification of subsets
mat <- db[match(sig, db$cg),] |> 
  dplyr::select(grep("epi|imm|mye|lym", colnames(db))) |> 
  as.matrix()
rownames(mat) <- sig

# Save for plotting
save(mat, file = "1-analysis-pipeline/2-output/plot_matrix.Rdata")

# Create subsets for computation
set.seed(2957)
u <- uwot::umap(mat)
u <- as.data.frame(u)
rownames(u) <- sig

subsets <- u |> 
  as.data.frame() |> 
  tibble::rownames_to_column("cg") |> 
  dplyr::rename(umap1 = V1,
                umap2 = V2) |> 
  mutate(set = case_when(umap2 > 1 & umap1 > -0.5 & umap1 < 5 ~ "set3",
                         umap1 > 5 ~ "set4",
                         umap2 < -1.5 ~ "set2",
                         umap1 < -0.5 & umap2 > -1.5 ~ "set1"))

# Save subsets:
# 1. Rename (informed based on data after plotting)
subsets <- subsets |>
  dplyr::mutate(set = case_when(set == "set1" ~ "epithelial_hypoM",
                                set == "set2" ~ "immune_hypoM",
                                set == "set3" ~ "distal_epithelial_hyperM",
                                set == "set4" ~ "proximal_epithelial_hyperM"))


# 2. Remove unreliable probes (low MI) and those which are not found on EPIC version 2.
lowmi <- read.csv("0-source/lowMI_names.csv", header = F)
notonepicv2 <- read.csv("0-source/notonEPICv2_names.csv", header = F)

subsets <- subsets |> 
  dplyr::filter(!cg %in% c(lowmi$V1, notonepicv2$V1)) # Removing some CpGs

# 3. Annotation
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19, lociNames = subsets$cg)
# identical(subsets$cg, rownames(anno))
subsets$chr <- anno$chr
subsets$pos <- anno$pos
subsets$strand <- anno$strand
subsets$on_450k_array <- ifelse(anno$Methyl450_Loci=="TRUE", "yes", "no")

# annotate genes
gene <- stringr::str_split(anno$UCSC_RefGene_Name, ";")
gene2 <- character(nrow(subsets))
for (g in 1:length(gene)){
  gene2[g] <- paste0(unique(gene[[g]]), collapse = ";")
}
subsets$gene <- gene2
rm(gene, gene2);gc()

# Annotate regions
region <- stringr::str_split(anno$UCSC_RefGene_Group, ";")
subsets$region <- character(nrow(subsets))
for (g in 1:length(region)){
  subsets$region[g] <- paste0(unique(region[[g]]), collapse = ";")
}

# Annotate relation to island
subsets$relation_to_island <- anno$Relation_to_Island
WID_SMK_cpgs <- subsets

# Save
save(WID_SMK_cpgs, file = "1-analysis-pipeline/2-output/WID_SMK_cpgs.Rdata")