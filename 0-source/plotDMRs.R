plotDMRs <- function(ranges, beta, pheno, type = "cervical", n = 100){
  
  if(length(find.package('MetBrewer', quiet = T)) == 0){
    cat('Installing MetBrewer')
    devtools::install_github('BlakeRMills/MetBrewer'); library(MetBrewer, character.only = T)
  }
  
  library(ggplot2)
  library(patchwork)
  
  
  # Load in data + cols
  data <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19,
                        lociNames = rownames(beta))
  cols <- MetBrewer::met.brewer("Hiroshige", type = "continuous")
  
  plots <- length(1:n)
  
  if(type == "blood"){
    pheno$ic <- as.numeric(pheno$hepidish_Eosino + pheno$hepidish_Mono + pheno$hepidish_Neutro)
  }
  
  for (i in 1:plots){
    # Ranges
    min <- start(ranges)[i]
    max <- end(ranges)[i]
    chr <- as.character(seqnames(ranges)[i])
    c <- data[data$chr==chr & data$pos >= min & data$pos <= max,]
    
    # Data
    tmp <- beta[match(rownames(c), rownames(beta)),]
    tmp_pheno <- cbind(pheno, t(tmp))
    tmp_pheno <- tmp_pheno |>
      tidyr::pivot_longer(rownames(tmp), names_to = "cpg", values_to = "beta") |> 
      dplyr::mutate(type2 = case_when(type == "Control" ~ "Control",
                               type != "Control" ~ "Case")) |> 
      dplyr::mutate(type2 = factor(type2, levels = c("Control", "Case")))
    
    ind <- match(tmp_pheno$cpg, rownames(c))
    tmp_pheno$pos <- c$pos[ind]
    t <- tmp_pheno |> 
      dplyr::group_by(type2, cpg) |> 
      dplyr::reframe(mean = mean(beta),
                     pos = pos) |> 
      dplyr::distinct() 
    
    # Plot linked ranges
    (a <- t %>%
        ggplot(aes(x = pos,
                   y = mean,
                   colour = type2)) +
        geom_point() +
        geom_line() +
        ggtitle(paste0(chr, ":", min, "-", max, ",\n",
                       ranges$overlapping.genes[i], ", p=", signif(ranges[i]$Stouffer,2))) +
        theme(legend.position="top") +
        ylim(0, 1) +
        scale_colour_manual(values = cols[c(7, 3)],
                            name = "") +
        xlab(paste0(unique(c$chr))) +
        ylab(paste0("mean beta at ", nrow(c), " CpGs")))
    
    
    if(type %in% c("buccal", "cervical")){
      tmpsum <- apply(tmp, 2, mean)
      pheno_tmp2 <- cbind(pheno, 
                          mean = tmpsum)
      
      b <- pheno_tmp2 |> 
        ggplot(aes(x = ic,
                   y = mean,
                   colour = type)) +
        geom_point(size = 0.8) +
        geom_smooth(method = "lm") +
        scale_colour_manual(values = c(cols[c(7,3)])) +
        theme(legend.position = "top") +
        ylab(paste0("mean beta at ", nrow(c), " CpGs"))
      
      plot <- a | b
      print(plot)
    } else {
      tmpsum <- apply(tmp, 2, mean)
      pheno_tmp2 <- cbind(pheno, 
                          mean = tmpsum)
      
      c <- pheno_tmp2 |> 
        ggplot(aes(x = ic,
                   y = mean,
                   colour = type)) +
        geom_point(size = 0.8) +
        geom_smooth(method = "lm") +
        scale_colour_manual(values = c(cols[c(7,3)])) +
        theme(legend.position = "top") +
        ylab(paste0("mean beta at ", nrow(c), " CpGs")) +
        xlab("neutrophil fraction")
      
      plot <- a | b 
      print(plot)
    }
  }
}
