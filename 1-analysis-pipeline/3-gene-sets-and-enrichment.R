# 3 Pathway enrichment and gene sets -------------------

cat("Enriching subsets of CpGs ----------------------- \n")

library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(dplyr)
load("1-analysis-pipeline/2-output/WID_SMK_cpgs.Rdata")

table(WID_SMK_cpgs$on_450k_array)
WID_SMK_cpgs <- WID_SMK_cpgs |> 
  dplyr::filter(chr != "chrX")

# 1. Extract gene sets
epi_hypo <- unique(as.character(stringr::str_split(WID_SMK_cpgs[WID_SMK_cpgs$set == "epithelial_hypoM",]$gene, ";", simplify = T)))
imm_hypo <- unique(as.character(stringr::str_split(WID_SMK_cpgs[WID_SMK_cpgs$set == "immune_hypoM",]$gene, ";", simplify = T)))
dist_epi_hyper <- unique(as.character(stringr::str_split(WID_SMK_cpgs[WID_SMK_cpgs$set == "distal_epithelial_hyperM",]$gene, ";", simplify = T)))
prox_epi_hyper <- unique(as.character(stringr::str_split(WID_SMK_cpgs[WID_SMK_cpgs$set == "proximal_epithelial_hyperM",]$gene, ";", simplify = T)))
epi_hypo <- epi_hypo[epi_hypo!=""]
imm_hypo <- imm_hypo[imm_hypo!=""]
dist_epi_hyper <- dist_epi_hyper[dist_epi_hyper != ""]
prox_epi_hyper <- prox_epi_hyper[prox_epi_hyper != ""]
sets <- list("epithelial\nhypoM" = epi_hypo, "immune\nhypoM" = imm_hypo,
             "distal\nepithelial\nhyperM" = dist_epi_hyper, "proximal\nepithelial\nhyperM" = prox_epi_hyper)
save(sets, file = "1-analysis-pipeline/3-output/venn-sets.Rdata")

# 2. Pathway enrichment
# Using "exclusive" genes (i.e. not in another set)

# Gene universe: genes that are associated in genes at intersect of delta betas
db_buccal <- readRDS("1-analysis-pipeline/1-delta-beta/1a-buccal/delta-beta.Rds")
db_blood <- readRDS("1-analysis-pipeline/1-delta-beta/1c-blood/delta-beta.Rds")
db_cervical <- readRDS("1-analysis-pipeline/1-delta-beta/1b-cervical/delta-beta.Rds")
cpgs <- intersect(rownames(db_buccal), intersect(rownames(db_blood), rownames(db_cervical)))
rm(db_buccal, db_blood, db_cervical);gc()

genes <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19, lociNames = cpgs) |> 
  as.data.frame() |> 
  dplyr::filter(!chr %in% c("chrX", "chrY"))
genes <- unique(as.character(stringr::str_split(genes$UCSC_RefGene_Name, ";", simplify = T)))
genes <- genes[genes !=""]

# ENTREZIDS
universegenes <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# 2.1. epi_hypo ---------------------------
epi_hypo_bp <- clusterProfiler::enrichGO(gene = epi_hypo[!epi_hypo %in% c(imm_hypo, dist_epi_hyper, prox_epi_hyper)],
                                           universe = genes,
                                           OrgDb = org.Hs.eg.db,
                                           keyType = 'SYMBOL',
                                           ont = "BP",
                                         pAdjustMethod = "BH")
as.data.frame(epi_hypo_bp) # Enrichment for detoxification

epi_hypo_mf <- clusterProfiler::enrichGO(gene = epi_hypo[!epi_hypo %in% c(imm_hypo, dist_epi_hyper, prox_epi_hyper)],
                                         universe = genes,
                                         OrgDb = org.Hs.eg.db,
                                         keyType = 'SYMBOL',
                                         ont = "MF",
                                         pAdjustMethod = "BH")
as.data.frame(epi_hypo_mf) # glucuronosyltransferase activity

epi_hypo_cc <- clusterProfiler::enrichGO(gene = epi_hypo[!epi_hypo %in% c(imm_hypo, dist_epi_hyper, prox_epi_hyper)],
                                         universe = genes,
                                         OrgDb = org.Hs.eg.db,
                                         keyType = 'SYMBOL',
                                         ont = "CC",
                                         pAdjustMethod = "BH")
as.data.frame(epi_hypo_cc) # Enrichment for actin cytoskeleton, stress fibers

tmp <- bitr(epi_hypo, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
epi_hypo_react_all <- enrichPathway(tmp$ENTREZID, organism = "human",
                                    universe = universegenes$ENTREZID,
                                    pvalueCutoff = 0.2,
                                    pAdjustMethod = "BH", 
                                    readable = T,
                                    minGSSize = 3)
x <- as.data.frame(epi_hypo_react_all) # ADME
save(epi_hypo_bp, epi_hypo_cc, epi_hypo_mf, epi_hypo_react_all, file = "1-analysis-pipeline/3-output/epi_hypo_GO.Rdata")

# 2.1. imm_hypo---------------------------
imm_hypo_bp <- clusterProfiler::enrichGO(gene = imm_hypo[!imm_hypo %in% c(epi_hypo, dist_epi_hyper, prox_epi_hyper)],
                                         universe = genes,
                                         OrgDb = org.Hs.eg.db,
                                         keyType = 'SYMBOL',
                                         ont = "BP",
                                         pAdjustMethod = "BH")
as.data.frame(imm_hypo_bp) # Enrichment for morphogenesis and development, differentiation


imm_hypo_mf <- clusterProfiler::enrichGO(gene = imm_hypo[!imm_hypo %in% c(epi_hypo, dist_epi_hyper, prox_epi_hyper)],
                                         universe = genes,
                                         OrgDb = org.Hs.eg.db,
                                         keyType = 'SYMBOL',
                                         ont = "MF",
                                         pAdjustMethod = "BH")
as.data.frame(imm_hypo_mf) # No enrichment

imm_hypo_cc <- clusterProfiler::enrichGO(gene = imm_hypo[!imm_hypo %in% c(epi_hypo, dist_epi_hyper, prox_epi_hyper)],
                                         universe = genes,
                                         OrgDb = org.Hs.eg.db,
                                         keyType = 'SYMBOL',
                                         ont = "CC",
                                         pAdjustMethod = "BH")
as.data.frame(imm_hypo_cc) # No enrichment

tmp <- bitr(imm_hypo, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
imm_hypo_react_all <- enrichPathway(tmp$ENTREZID, organism = "human",
                                    universe = universegenes$ENTREZID,
                                    pvalueCutoff = 0.2,
                                    pAdjustMethod = "BH", 
                                    readable = T,
                                    minGSSize = 3)
x <- as.data.frame(imm_hypo_react_all)
enrichplot::cnetplot(imm_hypo_react_all) # Thrombin signaling
save(imm_hypo_bp, imm_hypo_react_all, file = "1-analysis-pipeline/3-output/imm_hypo_GO.Rdata")

# 2.3. dist_epi_hyper ---------------------------
dist_epi_hyper_bp <- clusterProfiler::enrichGO(gene = dist_epi_hyper[!dist_epi_hyper %in% c(epi_hypo, imm_hypo, prox_epi_hyper)],
                                         universe = genes,
                                         OrgDb = org.Hs.eg.db,
                                         keyType = 'SYMBOL',
                                         ont = "BP",
                                         pAdjustMethod = "BH")
as.data.frame(dist_epi_hyper_bp) # No enrichment

dist_epi_hyper_mf <- clusterProfiler::enrichGO(gene = dist_epi_hyper[!dist_epi_hyper %in% c(epi_hypo, imm_hypo, prox_epi_hyper)],
                                         universe = genes,
                                         OrgDb = org.Hs.eg.db,
                                         keyType = 'SYMBOL',
                                         ont = "MF",
                                         pAdjustMethod = "BH")
as.data.frame(dist_epi_hyper_mf) # No enrichment

dist_epi_hyper_cc <- clusterProfiler::enrichGO(gene = dist_epi_hyper[!dist_epi_hyper %in% c(epi_hypo, imm_hypo, prox_epi_hyper)],
                                         universe = genes,
                                         OrgDb = org.Hs.eg.db,
                                         keyType = 'SYMBOL',
                                         ont = "CC",
                                         pAdjustMethod = "BH")
as.data.frame(dist_epi_hyper_cc) # No enrichment

tmp <- bitr(dist_epi_hyper, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
dist_epi_hyper_react_all <- enrichPathway(tmp$ENTREZID, organism = "human",
                                          universe = universegenes$ENTREZID,
                                          pvalueCutoff = 0.2,
                                          pAdjustMethod = "BH", 
                                          readable = T,
                                          minGSSize = 3)
x <- as.data.frame(dist_epi_hyper_react_all)
enrichplot::cnetplot(dist_epi_hyper_react_all) # Nuclear receptors, PPARa, lipid metabolism

save(dist_epi_hyper_bp, dist_epi_hyper_mf, dist_epi_hyper_react_all, file = "1-analysis-pipeline/3-output/dist_epi_hyper_GO.Rdata")

# 2.4. prox_epi_hyper ---------------------------
prox_epi_hyper_bp <- clusterProfiler::enrichGO(gene = prox_epi_hyper[!prox_epi_hyper %in% c(epi_hypo, imm_hypo, dist_epi_hyper)],
                                               universe = genes,
                                               OrgDb = org.Hs.eg.db,
                                               keyType = 'SYMBOL',
                                               ont = "BP",
                                               pAdjustMethod = "BH")
as.data.frame(prox_epi_hyper_bp) # Wnt signaling

prox_epi_hyper_mf <- clusterProfiler::enrichGO(gene = prox_epi_hyper[!prox_epi_hyper %in% c(epi_hypo, imm_hypo, dist_epi_hyper)],
                                               universe = genes,
                                               OrgDb = org.Hs.eg.db,
                                               keyType = 'SYMBOL',
                                               ont = "MF",
                                               pAdjustMethod = "BH")
as.data.frame(prox_epi_hyper_mf) # polyU RNA binding, pyrimidine

prox_epi_hyper_cc <- clusterProfiler::enrichGO(gene = prox_epi_hyper[!prox_epi_hyper %in% c(epi_hypo, imm_hypo, dist_epi_hyper)],
                                               universe = genes,
                                               OrgDb = org.Hs.eg.db,
                                               keyType = 'SYMBOL',
                                               ont = "CC",
                                               pAdjustMethod = "BH")
as.data.frame(prox_epi_hyper_cc) # No enrichment

tmp <- bitr(prox_epi_hyper, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
prox_epi_hyper_react_all <- enrichPathway(tmp$ENTREZID, organism = "human",
                                      universe = universegenes$ENTREZID,
                                      pvalueCutoff = 0.2,
                                      pAdjustMethod = "BH", 
                                      readable = T,
                                      minGSSize = 3)
x <- as.data.frame(prox_epi_hyper_react_all)
enrichplot::cnetplot(prox_epi_hyper_react_all)

save(prox_epi_hyper_mf, prox_epi_hyper_react_all, file = "1-analysis-pipeline/3-output/prox_epi_hyper_GO.Rdata")
rm(list=ls())

# Write tables
files <- c("1-analysis-pipeline/3-output/epi_hypo_GO.Rdata",
           "1-analysis-pipeline/3-output/imm_hypo_GO.Rdata",
           "1-analysis-pipeline/3-output/dist_epi_hyper_GO.Rdata",
           "1-analysis-pipeline/3-output/prox_epi_hyper_GO.Rdata")

for (i in files){
  load(i)
  
  filename = case_when(grepl("epi_hypo", i) ~ "epithelial_hypoM",
                       grepl("imm_hypo", i) ~ "immune_hypoM",
                       grepl("dist_epi_hyper", i) ~ "distal_epithelial_hyperM",
                       grepl("prox_epi_hyper", i) ~ "proximal_epithelial_hyperM")
  vars <- ls()
  vars <- vars[grepl("hyper|hypo", vars)]
    
  for (v in vars){
    name = case_when(grepl("bp", v) ~ "GO_BP",
                     grepl("mf", v) ~ "GO_MF",
                     grepl("cc", v) ~ "GO_CC",
                     grepl("react", v) ~ "REACT_PA")
    
    append = ifelse (v == vars[1], F, T)
    xlsx::write.xlsx2(get(v), file = paste0("1-analysis-pipeline/3-output/", filename, ".xlsx"),
                      col.names = T, row.names = F,
                      append = append,
                      sheetName = name)
    
    if(v == vars[length(vars)]){
      rm <- ls()
      rmfiles <- c(rm[!rm %in% c("files", "i")], "rm", "rmfiles")
      rm(list = rmfiles)
    }
  }
}

rm(list=ls())

# Identify PCGT cpgs -----------------------------------------
# Working from Lee et al. 2006 paper based on promoter occupancy of genes (single/triple)

library(readxl)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

# genes on array
data <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19, what = "everything")
genes <- data$UCSC_RefGene_Name
tss <- data$UCSC_RefGene_Group
pb <- txtProgressBar(min = 0, max = length(genes), style = 3)
for (i in 1:length(genes)){
  genes[i] <- strsplit(genes[i], ";")[[1]][1]
  tss[i] <- strsplit(tss[i], ";")[[1]][1]
  setTxtProgressBar(pb, i)
}

data$gene <- genes # 26650 unique genes
data$tss <- tss

# Single occupancy from paper ---------------------

dat <- read_xls("0-source/mmc10.xls", skip = 1, col_types = "text")
dat <- dat[(dat$`Suz12 Occupancy`==1 | dat$`Eed Occupancy` == 1 | dat$`H3K27me3 Occupancy` == 1),]

pcgt_single <- intersect(data$gene, dat$`Gene Name`) # 1343 genes
cpgs_single <- data[data$gene %in% pcgt_single,]
nrow(cpgs_single) #57840 CpGs

# Filtering by TSS200
cpgs_single_tss200 <- cpgs_single[cpgs_single$tss == "TSS200",] # only 4625 left
save(pcgt_single, cpgs_single, cpgs_single_tss200, file = "0-source/single_occupancy.Rdata")

rm(list=ls());gc()