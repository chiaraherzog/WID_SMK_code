# Enrichment

cat("CpG subset enrichment for e-cigarette users ----------------------- \n")

# Load with smk cpgs
load("1-analysis-pipeline/2-output/WID_SMK_cpgs.Rdata")

# Helper function for dbs
load_and_name <- function(path = "",
                          set = "",
                          type2 = ""){
  tmp <- readRDS(path) |> 
    as.data.frame() |> 
    tibble::rownames_to_column("cg") |> 
    dplyr::filter(cg %in% WID_SMK_cpgs$cg) |> 
    dplyr::mutate(set = set,
                  type2 = type2)
  
  return(tmp)
}

# Load data
db_orig <- load_and_name(path = "1-analysis-pipeline/1-delta-beta/1a-buccal/delta-beta.Rds",
                         set = "discovery",
                         type2 = "smoker (discovery set)")

db_ecig <- load_and_name(path = "1-analysis-pipeline/6-output/delta-beta.Rds",
                         set = "ecig",
                         type2 = "e-cigarette user")

dat <- plyr::rbind.fill(db_orig, db_ecig) |> 
  tidyr::pivot_wider(id_cols = "cg",
                     names_from = c("set", "type2"),
                     values_from = c("db_immune", "db_epithelial")) |> 
  tidyr::drop_na() |> 
  dplyr::left_join(dplyr::select(WID_SMK_cpgs, cg, set, chr, gene)) |> 
  dplyr::filter(!chr %in% c("chrX", "chrY"))

# Keep those that are in the same direction
same_epi_hypoM <- dat |> 
  dplyr::filter(set == "epithelial_hypoM") |> 
  dplyr::filter((`db_epithelial_discovery_smoker (discovery set)` < 0 & `db_epithelial_ecig_e-cigarette user` < 0) | (`db_epithelial_discovery_smoker (discovery set)` > 0 & `db_epithelial_ecig_e-cigarette user` > 0))

same_prox <- dat |> 
  dplyr::filter(set == "proximal_epithelial_hyperM") |> 
  dplyr::filter((`db_epithelial_discovery_smoker (discovery set)` < 0 & `db_epithelial_ecig_e-cigarette user` < 0) | (`db_epithelial_discovery_smoker (discovery set)` > 0 & `db_epithelial_ecig_e-cigarette user` > 0))

same_dist <- dat |> 
  dplyr::filter(set == "distal_epithelial_hyperM") |> 
  dplyr::filter((`db_epithelial_discovery_smoker (discovery set)` < 0 & `db_epithelial_ecig_e-cigarette user` < 0) | (`db_epithelial_discovery_smoker (discovery set)` > 0 & `db_epithelial_ecig_e-cigarette user` > 0))

# Enrichments --------------
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

# 1. Extract gene sets
epi_hypo <- unique(same_epi_hypoM$gene)
dist_epi_hyper <- unique(same_dist$gene)
prox_epi_hyper <-  unique(same_prox$gene)
epi_hypo <- epi_hypo[epi_hypo!=""]
dist_epi_hyper <- dist_epi_hyper[dist_epi_hyper != ""]
prox_epi_hyper <- prox_epi_hyper[prox_epi_hyper != ""]

# 2. Pathway enrichment
# Gene universe: genes that are associated in genes at intersect of delta betas
db_buccal <- readRDS("1-analysis-pipeline/1-delta-beta/1a-buccal/delta-beta.Rds")
db_ecig <- readRDS("1-analysis-pipeline/6-output/delta-beta.Rds")
cpgs <- intersect(rownames(db_buccal), rownames(db_ecig))
rm(db_buccal, db_ecig);gc()
genes <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19, lociNames = cpgs) |> 
  as.data.frame() |> 
  dplyr::filter(!chr %in% c("chrX", "chrY"))
genes <- unique(as.character(stringr::str_split(genes$UCSC_RefGene_Name, ";", simplify = T)))
genes <- genes[genes !=""]

# ENTREZIDS
universegenes <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

# epi_hypo ---------------------------
epi_hypo_bp <- clusterProfiler::enrichGO(gene = epi_hypo,
                                         universe = genes,
                                         OrgDb = org.Hs.eg.db,
                                         keyType = 'SYMBOL',
                                         ont = "BP",
                                         pAdjustMethod = "BH")
as.data.frame(epi_hypo_bp) # Enrichment for detoxification

epi_hypo_mf <- clusterProfiler::enrichGO(gene = epi_hypo,
                                         universe = genes,
                                         OrgDb = org.Hs.eg.db,
                                         keyType = 'SYMBOL',
                                         ont = "MF",
                                         pAdjustMethod = "BH")
as.data.frame(epi_hypo_mf) # antioxidant activity

epi_hypo_cc <- clusterProfiler::enrichGO(gene = epi_hypo,
                                         universe = genes,
                                         OrgDb = org.Hs.eg.db,
                                         keyType = 'SYMBOL',
                                         ont = "CC",
                                         pAdjustMethod = "BH")
as.data.frame(epi_hypo_cc) # no enrichment

tmp <- bitr(epi_hypo, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
epi_hypo_react_all <- enrichPathway(tmp$ENTREZID, organism = "human",
                                    universe = universegenes$ENTREZID,
                                    pvalueCutoff = 0.2,
                                    pAdjustMethod = "BH", 
                                    readable = T,
                                    minGSSize = 3)
x <- as.data.frame(epi_hypo_react_all) # cellular response to chemical stress
save(epi_hypo_bp, epi_hypo_cc, epi_hypo_mf, epi_hypo_react_all, file = "1-analysis-pipeline/7-output/epi_hypo_GO.Rdata")

# prox_epi_hyper ---------------------------
prox_epi_hyper_bp <- clusterProfiler::enrichGO(gene = prox_epi_hyper,
                                         universe = genes,
                                         OrgDb = org.Hs.eg.db,
                                         keyType = 'SYMBOL',
                                         ont = "BP",
                                         pAdjustMethod = "BH")
as.data.frame(prox_epi_hyper_bp) # none

prox_epi_hyper_mf <- clusterProfiler::enrichGO(gene = prox_epi_hyper,
                                         universe = genes,
                                         OrgDb = org.Hs.eg.db,
                                         keyType = 'SYMBOL',
                                         ont = "MF",
                                         pAdjustMethod = "BH")
as.data.frame(prox_epi_hyper_mf) # polypyrimidine binding

prox_epi_hyper_cc <- clusterProfiler::enrichGO(gene = prox_epi_hyper,
                                         universe = genes,
                                         OrgDb = org.Hs.eg.db,
                                         keyType = 'SYMBOL',
                                         ont = "CC",
                                         pAdjustMethod = "BH")
as.data.frame(prox_epi_hyper_cc) # none

tmp <- bitr(prox_epi_hyper, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
prox_epi_hyper_react_all <- enrichPathway(tmp$ENTREZID, organism = "human",
                                    universe = universegenes$ENTREZID,
                                    pvalueCutoff = 0.2,
                                    pAdjustMethod = "BH", 
                                    readable = T,
                                    minGSSize = 3)
x <- as.data.frame(prox_epi_hyper_react_all) # cellular response to chemical stress
enrichplot::cnetplot(prox_epi_hyper_react_all)
enrichplot::cnetplot(epi_hypo_react_all)

save(prox_epi_hyper_bp, prox_epi_hyper_cc, prox_epi_hyper_mf, prox_epi_hyper_react_all, file = "1-analysis-pipeline/7-output/prox_epi_hyper_GO.Rdata")

# Write tables
files <- c("1-analysis-pipeline/7-output/epi_hypo_GO.Rdata",
           "1-analysis-pipeline/7-output/prox_epi_hyper_GO.Rdata")

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
    xlsx::write.xlsx2(get(v), file = paste0("1-analysis-pipeline/7-output/", filename, ".xlsx"),
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

