cat("Expression analysis (TCGA-LUAD/LUSC) ----------------------- \n")

library(here)
library(SummarizedExperiment)
library(purrr)
library(broom)

here::i_am("1-analysis-pipeline/8-expression.R")
# path_dat <- '<path to your expression data>'
path_dat <- '~/Documents/Work/data.nosync/tcga_sola/'
load(here("1-analysis-pipeline", '2-output', 'WID_SMK_cpgs.Rdata'))

# Luad
id <- 'TCGA-LUAD'
load(file.path(path_dat, id, '/expression/data_expression_full.Rdata'))
exp_luad <- assay(data[,,'unstranded'])
load(file.path(path_dat, id, '/beta.Rdata'))
beta_luad <- beta[rownames(beta) %in% WID_SMK_cpgs$cg,]

intersect(substr(colnames(beta_luad), 1, 16),
          substr(colnames(exp_luad), 1, 16))

# LUSC
id <- 'TCGA-LUSC'
load(file.path(path_dat, id, '/expression/data_expression_full.Rdata'))
exp_lusc <- assay(data[,,'tpm_unstrand'])
load(file.path(path_dat, id, '/beta.Rdata'))
beta_lusc <- beta[rownames(beta) %in% WID_SMK_cpgs$cg,]

# Merge expression files
exp <- cbind(exp_luad, exp_lusc)

nrow(beta_luad)
nrow(beta_lusc)

intersect <- intersect(rownames(beta_luad), rownames(beta_lusc))
beta <- cbind(beta_luad[intersect,], beta_lusc[intersect,])

# match files
intersect <- intersect(substr(colnames(beta), 1, 16),
                       substr(colnames(exp), 1, 16))

exp <- exp[,match(intersect, substr(colnames(exp), 1, 16))]
beta <- beta[,match(intersect, substr(colnames(beta), 1, 16))]
identical(substr(colnames(beta), 1, 16),
          substr(colnames(exp), 1, 16))

# Filter genes which have 0 variance
exp_filtered <- exp[apply(exp, 1, function(x) var(x, na.rm = T) != 0),]

# Filter genes which have mean < 20 counts
exp_filtered <- exp_filtered[apply(exp_filtered, 1, function(x) mean(x) < 20),]
nrow(exp_filtered)

# Correlation: log format
plot(beta[1,], exp[1,])
plot(beta[1,], log2(exp[1,]+0.000001))
rm(beta_luad, beta_lusc, exp_luad, exp_lusc, data);gc()

# # Overall scores
# load("1-analysis-pipeline/4-datasets/1-output/tcga.Rdata")
# 
# ind <- match(substr(colnames(exp_filtered), 1, 16), substr(rownames(tcga_data), 1, 16))
# tcga_data <- tcga_data[ind,]  
# identical(substr(colnames(exp_filtered), 1, 16), substr(rownames(tcga_data), 1, 16))  
# 
# scores <- colnames(tcga_data)[grepl("WID_SMK450", colnames(tcga_data))]
# 
# pb = txtProgressBar(min = 0, max = length(scores), style = 3, width = 50) 
# 
# out_scores <- lapply(1:length(scores), function(x){
#   setTxtProgressBar(pb,x)
#   apply(exp_filtered, 1, function(y){
#     cor.test(tcga_data[,scores[x]], log2(y+0.000001))
#   })
# })
# 
# corr_result <- dplyr::bind_rows(lapply(out_scores, function(x) { map_dfr(x, tidy, .id = "ensg") }), .id = 'score') |> 
#   dplyr::select(c('score', 'ensg', 'estimate', 'p.value', 'conf.low', 'conf.high'))
# 
# 
# 
# corr_result2 <- data.table::rbindlist(lapply(out_scores, function(x) { map_dfr(x, tidy, .id = "ensg") }),idcol = 'score')
# 
# corr_result <- corr_result |> 
#   dplyr::mutate(score = dplyr::case_when(score == '1' ~ "WID_SMK450_distal_epithelial_hyperM",
#                                score == '2' ~ "WID_SMK450_epithelial_hypoM",
#                                score == '3' ~ "WID_SMK450_immune_hypoM",
#                                score == '4' ~ "WID_SMK450_proximal_epithelial_hyperM"))
# library(ggplot2)
# 
# corr_result |> 
#   ggplot(aes(x = p.value)) +
#   geom_histogram() +
#   facet_wrap(~score)
# 
# sig <- corr_result |> 
#   dplyr::mutate(padj = p.adjust(p.value)) |> 
#   dplyr::filter(padj < 0.05)
# 
# table(sig$score)
# table(sig$score, sig$ensg)
# 
# sig |> 
#   ggplot(aes(x = estimate)) +
#   geom_histogram() +
#   facet_wrap(~score)
# 
# # Enrich
# library(clusterProfiler)
# library(org.Hs.eg.db)
# distal_hypoM <- bitr(gsub("\\..*", "", sig[sig$score=='WID_SMK450_distal_epithelial_hyperM',]$ensg), fromType="ENSEMBL", toType="SYMBOL", OrgDb="org.Hs.eg.db")
# 
# distal_hyperM <- clusterProfiler::enrichGO(gene = gsub("\\..*", "", sig[sig$score=='WID_SMK450_distal_epithelial_hyperM',]$ensg),
#                                          OrgDb = org.Hs.eg.db,
#                                          keyType = 'ENSEMBL',
#                                          ont = "BP",
#                                          pAdjustMethod = "BH")
# epi_hypoM <- clusterProfiler::enrichGO(gene = unique(gsub("\\..*", "", sig[sig$score=='WID_SMK450_epithelial_hypoM',]$ensg)),
#                                            OrgDb = org.Hs.eg.db,
#                                            keyType = 'ENSEMBL',
#                                            ont = "BP",
#                                            pAdjustMethod = "BH")
# immune_hypoM <- clusterProfiler::enrichGO(gene = unique(gsub("\\..*", "", sig[sig$score=='WID_SMK450_immune_hypoM',]$ensg)),
#                                        OrgDb = org.Hs.eg.db,
#                                        keyType = 'ENSEMBL',
#                                        ont = "BP",
#                                        pAdjustMethod = "BH")
# 
# prox_hypoM <- clusterProfiler::enrichGO(gene = unique(gsub("\\..*", "", sig[sig$score=='WID_SMK450_proximal_epithelial_hyperM',]$ensg)),
#                                           OrgDb = org.Hs.eg.db,
#                                           keyType = 'ENSEMBL',
#                                           ont = "BP",
#                                           pAdjustMethod = "BH")
# 
# x<-as.data.frame(prox_hypoM)
# rm(out_scores, corr_result);gc()



# # Individual
# 
# df <- data.frame(matrix(nrow = nrow(beta)*nrow(exp_filtered), ncol = 4))
# colnames(df) <-  c("cg", 'ensg', 'estimate', 'p')
# df$cg <- rep(rownames(beta), each = nrow(exp_filtered))
# df$ensg <- rep(rownames(exp_filtered), times = nrow(beta))
# table(df$cg)
# table(df$ensg)
# 
# rownames <- rownames(beta)
# expnames <- rownames(exp_filtered)
# 
# pb = txtProgressBar(min = 0, max = nrow(beta), style = 3, width = 50) 
# 
# for (i in rownames){
#   setTxtProgressBar(pb,which(rownames==i))
#   for (j in expnames){
#     tmp <- cor.test(beta[i,], log2(exp_filtered[j,]+0.000001))
#     df[df$cg == i & df$ensg == j,]$p <- tmp$p.value
#     df[df$cg == i & df$ensg == j,]$estimate <- tmp$estimate
#   }
# }
# 
# 
# out <- lapply(1:nrow(beta), function(x){
#   setTxtProgressBar(pb,x)
#   apply(exp_filtered, 1, function(y){
#    cor.test(beta[x,], log2(y+0.000001))
#   })
# })
# 
# library(broom)
# library(purrr)
# pb = txtProgressBar(min = 0, max = nrow(beta), style = 3, width = 50) 
# 
# corr_result <- dplyr::bind_rows(lapply(out, function(x) { map_dfr(x, tidy, .id = "ensg") }),.id = 'cg') |> 
#   dplyr::select(c('cg', 'ensg', 'estimate', 'p.value', 'conf.low', 'conf.high'))
# 
# head(corr_result)
# 
# hist(corr_result$p.value)
# sum(corr_result$p.value<0.05)
# sum(p.adjust(corr_result$p.value)<0.05)
# 
# x <- corr_result[p.adjust(corr_result$p.value)<0.05,]
# 
# 
# 
# 
# 
# # Individual
# zs
# 
# 
# pb = txtProgressBar(min = 0, max = nrow(beta), style = 3, width = 50) 
# 
# out <- lapply(1:nrow(beta), function(x){
#   setTxtProgressBar(pb,x)
#   apply(exp_filtered, 1, function(y){
#    cor.test(beta[x,], log2(y+0.000001))
#   })
# })
# 
# library(broom)
# library(purrr)
# pb = txtProgressBar(min = 0, max = nrow(beta), style = 3, width = 50) 
# 
# corr_result <- dplyr::bind_rows(lapply(out, function(x) { map_dfr(x, tidy, .id = "ensg") }),.id = 'cg') |> 
#   dplyr::select(c('cg', 'ensg', 'estimate', 'p.value', 'conf.low', 'conf.high'))
# 
# head(corr_result)
# 
# hist(corr_result$p.value)
# sum(corr_result$p.value<0.05)
# sum(p.adjust(corr_result$p.value)<0.05)
# 
# x <- corr_result[p.adjust(corr_result$p.value)<0.05,]

# Keep CpGs in tss
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19, lociNames = rownames(beta))
# anno <- anno[grepl("TSS200|TSS1500", anno$UCSC_RefGene_Group),] # 52 loci only

genes <- stringr::str_split(anno$UCSC_RefGene_Name, ";")
regions <- stringr::str_split(anno$UCSC_RefGene_Group, ";")

names(genes) <- anno$Name
names(regions) <- anno$Name

for (i in anno$Name){
  
  # find gene
  if(length(genes[[i]]) == 1 && genes[[i]] != ""){
    
    for (x in 1:length(genes[[i]])){
      ensg <- tryCatch(clusterProfiler::bitr(genes[[i]][x], fromType="SYMBOL", toType="ENSEMBL", OrgDb="org.Hs.eg.db"),
                   error = function(e) {
                     df <- data.frame(cg = i,
                                      ensg = NA,
                                      symbol = genes[[i]][x],
                                      region = regions[[i]][x],
                                      cor = NA, 
                                      p = NA)
                   })
  
      # if more than 1 ensg mapped, need to do each one separately
      for(j in ensg$ENSEMBL){
      
        # if only one ensg, do only one
        if(sum(grepl(j, rownames(exp))>1)){
          e <- rowSums(exp[grepl(j, rownames(exp)),])
        } else {
          e <- exp[grepl(j, rownames(exp)),]
        }
      
        # Correlation
        cor <- ifelse(length(e) != 0, cor(beta[i,], log2(e+1)), NA)
        p <- ifelse(length(e) != 0, cor.test(beta[i,], log2(e+1))$p.value, NA)
      
        # save
        if(j == ensg$ENSEMBL[1]){
          df <- data.frame(cg = i,
                           ensg = j,
                           symbol = genes[[i]][x],
                           region = regions[[i]][x],
                           cor = cor,
                           p = p)
        } else {
          tmp <- data.frame(cg = i,
                            ensg = j,
                            symbol = genes[[i]][x],
                            region = regions[[i]][x],
                            cor = cor,
                            p = p)
          df <- rbind(df, tmp)
        }
      } 
    }
    } else {
    df <- data.frame(cg = i,
                     ensg = NA,
                     symbol = NA,
                     cor = NA, 
                     p = NA)
    }

  if(i == anno$Name[1]){
      df_full <- df
    } else {
      df_full <- dplyr::bind_rows(df_full, df)
    }
} 


hist(df_full$p)
df <- df_full |> 
  dplyr::filter(!is.na(cor)) |> 
  dplyr::distinct() |> 
  dplyr::mutate(padj = p.adjust(p))

library(ggplot2)

df |> 
  ggplot(aes(x = cor))+
  geom_histogram() +
  facet_wrap(~region)

df |> 
  dplyr::filter(padj < 0.05) |> 
  ggplot(aes(x = cor))+
  geom_histogram() +
  facet_wrap(~region)

df |> 
  dplyr::filter(padj < 0.05) |> 
  ggplot(aes(x = cor))+
  geom_histogram() +
  facet_wrap(~region)

df |> 
  dplyr::filter(padj < 0.05) |> 
  ggplot(aes(x = cor))+
  geom_density(aes(fill = region),
               alpha = 0.3,
               bw = 0.1)

df |> 
  ggplot(aes(x = p)) +
  geom_histogram(bins = 7) +
  facet_wrap(~region)

save(df, file = "1-analysis-pipeline/8-output/expression.Rdata")

df |> 
  writexl::write_xlsx(path = here("2-markdown", "2-tables", "Suppl_Table_7_expression.xlsx"))
