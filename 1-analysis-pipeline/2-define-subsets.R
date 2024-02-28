# 2. Define subsets of CpGs

cat("Defining subsets of CpGs (clustering) ----------------------- \n")

library(here)
library(dplyr)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# Unreliable or non-v2 probes
lowmi <- read.csv("0-source/lowMI_names.csv", header = F)
notonepicv2 <- read.csv("0-source/notonEPICv2_names.csv", header = F)

# Load in data and keep CpGs with a padj < 0.05 in at least one set
db_buccal <- readRDS("1-analysis-pipeline/1-delta-beta/1a-buccal/delta-beta.Rds") |> 
  as.data.frame() |> 
  tibble::rownames_to_column("cg") |>
  dplyr::filter(!cg %in% c(lowmi$V1, notonepicv2$V1)) |> 
  dplyr::mutate(padj = p.adjust(type),
                sig = ifelse(padj < 0.05, "yes", "no")) 

table(db_buccal$sig)
sig <- db_buccal[db_buccal$sig != 'yes',]
sig[which.min(sig$type),]
format(sig[which.min(sig$type),]$type, scientific = T)


db_cervical <- readRDS("1-analysis-pipeline/1-delta-beta/1b-cervical/delta-beta.Rds") |> 
  as.data.frame() |> 
  tibble::rownames_to_column("cg") |>
  dplyr::filter(!cg %in% c(lowmi$V1, notonepicv2$V1)) |> 
  dplyr::mutate(padj = p.adjust(type, method = 'bonferroni'),
         sig = ifelse(padj < 0.05, "yes", "no"))

table(db_cervical$sig)
sig <- db_cervical[db_cervical$sig != 'yes',]
sig[which.min(sig$type),]
format(sig[which.min(sig$type),]$type, scientific = T)

db_blood <- readRDS("1-analysis-pipeline/1-delta-beta/1c-blood/delta-beta.Rds") |> 
  as.data.frame() |> 
  tibble::rownames_to_column("cg") |>
  dplyr::filter(!cg %in% c(lowmi$V1, notonepicv2$V1)) |> 
  dplyr::mutate(padj = p.adjust(type, method = 'bonferroni'),
         sig = ifelse(padj < 0.05, "yes", "no"))
table(db_blood$sig)
sig <- db_blood[db_blood$sig != 'yes',]
sig[which.min(sig$type),]
format(sig[which.min(sig$type),]$type, scientific = T)

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

# PCA
library(factoextra)
x <- FactoMineR::PCA(mat)
# factoextra::fviz_contrib(x, choice = 'var',axes = 1)
# factoextra::fviz_contrib(x, choice = 'var',axes = 2)
# factoextra::fviz_contrib(x, choice = 'var',axes = 3)
# factoextra::fviz_contrib(x, choice = 'var',axes = 4)
# factoextra::fviz_contrib(x, choice = 'var',axes = 5)
# factoextra::fviz_contrib(x, choice = 'var',axes = 6)
# factoextra::fviz_pca_biplot(x)

# HCLUST - additional clustering method
cols2 <- MetBrewer::met.brewer("Hiroshige", n = 9)

# ComplexHeatmap::Heatmap(t(mat),
#                         clustering_distance_columns = 'manhattan',
#                         clustering_method_columns = 'ward.D',
#                         column_km = 4,show_column_names = F,
#                         col = circlize::colorRamp2(breaks = seq(from = -0.2, to = 0.2,
#                                                                         length.out = 9),
#                                                            colors = rev(cols2)))


# Create subsets for computation
set.seed(2957)
u <- uwot::umap(mat)
u <- as.data.frame(u)
rownames(u) <- sig

library(ggplot2)

u |> 
  ggplot(aes(x = V1, y = V2)) +
  geom_point()

subsets <- u |> 
  dplyr::rename(umap1 = V1,
                umap2 = V2) |> 
  mutate(set = case_when(umap1 > 0 & umap1 < 5 ~ "set3",
                         umap1 > 5 ~ "set4",
                         umap2 < -2 ~ "set2",
                         umap1 < -0 & umap2 > -2 ~ "set1"))

# Save subsets:
# 1. Rename (informed based on data after plotting)
subsets <- subsets |>
  dplyr::mutate(set = case_when(set == "set1" ~ "epithelial_hypoM",
                                set == "set2" ~ "immune_hypoM",
                                set == "set3" ~ "distal_epithelial_hyperM",
                                set == "set4" ~ "proximal_epithelial_hyperM")) |> 
  tibble::rownames_to_column("cg")

# 2. Annotation
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19, lociNames = subsets$cg)
identical(subsets$cg, rownames(anno))
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
