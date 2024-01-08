saveDMRs <- function(ranges, beta, pheno,
                     file){
  
  library(stringr)
  library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  
  r <- as.data.frame(ranges)
  mean_beta <- data.frame(matrix(nrow = nrow(r),
                                 ncol = nrow(pheno)))
  rownames(mean_beta) <- paste0(r$seqnames, ";", r$start, "-", r$end)
  colnames(mean_beta) <- rownames(pheno)
  
  data <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19, lociNames = rownames(beta))
  data$gene <- str_split(data$UCSC_RefGene_Name, ";", simplify = TRUE)[,1]
  
  chr <- character(nrow(r))
  start <- numeric(nrow(r))
  end <- numeric(nrow(r))
  
  pb <- txtProgressBar(min = 1, max = nrow(r), style = 3)
  
  for (i in 1:nrow(r)){
    setTxtProgressBar(pb, i)
    chr <- r$seqnames[i]
    start <- r$start[i]
    end <- r$end[i]
    c <- data[data$chr==chr & data$pos >= start & data$pos <= end,]
    tmp <- beta[match(rownames(c), rownames(beta)),]
    mean_beta[i,] <- apply(tmp, 2, mean)
  } 
  
  # type
  # ctrl <- apply(mean_beta[,pheno$type=="Control"], 1, mean)
  # mean_beta <- mean_beta[ctrl<0.2,]
  
  ctrl <- apply(mean_beta[,pheno$type=="Control"], 1, mean)
  case <- apply(mean_beta[,pheno$type != "Control"], 1, mean)
  
  df <- r[match(rownames(mean_beta),
                paste0(r$seqnames, ";", r$start, "-", r$end)),]
  df <- as.data.frame(df)
  rownames(df) <- NULL
  
  df$ucsc_gene <- ""
  for (i in 1:nrow(df)){
    c <- data[data$chr==df$seqnames[i] & data$pos >= df$start[i] & data$pos <= df$end[i],]
    df$ucsc_gene[i] <-   paste0(unique(c$gene), collapse = ",")
  }
  
  df <- df |> 
    dplyr::mutate(ctrl = round(ctrl,2),
                  case = round(case, 2),
                  diff = round(case-ctrl,2)) %>%
    dplyr::select(seqnames, start, end, width, no.cpgs, overlapping.genes, ctrl, case, diff) 
  
  saveRDS(df, file = paste0(file, ".Rds"))
  write.table(df, file = paste0(file, ".csv"), sep = ",",
              row.names = FALSE)
  return(df)
}
