# Delta beta 

estimateDeltaBeta <- function(beta, pheno, output, adjustment = c("age"),
                              typevar = "type", base = "Control",
                              names = c("db_epithelial", "db_immune")){
  
  # Adjustment needs a minimum of type and ic
  # if not type naming, need to mention what the "base" level should be
  
  cat('Beginning delta-beta estimation script...\n\n')
  cat('Current time:',as.character(Sys.time()),'\n')
  cat('cwd:',getwd(),'\n\n')
  
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(tibble))
  suppressPackageStartupMessages(library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(stringr))
  
  # Check directories and files exist
  if(!exists("beta")){
    stop('beta not present')
  }
  
  if(!exists("pheno")){
    stop('pheno not present')
  }
  
  if(!dir.exists(output)){
    dir.create(output)
    cat('output directory has been created\n')
  }
  
  ind <- match(colnames(beta), rownames(pheno))
  pheno <- pheno[ind,]
  
  cat('Statistical adjustment will be made for:\n')
  
  adjustment <- c("type", "ic", adjustment)
  adjustment <- unique(adjustment[adjustment!=""])
  for(i in adjustment){cat(i,'\n')}
  cat('\n')
  
  # Set factor levels correctly
  pheno <- pheno |> 
    dplyr::mutate(type = case_when(get(typevar) == base ~ "Control",
                                   get(typevar) != base ~ "Case"),
                  type = factor(type, levels = c("Control", "Case")),
                  ic = as.numeric(ic)) |> 
    dplyr::select(all_of(adjustment))
  
  if("age" %in% adjustment){
    pheno |> 
      mutate(age = as.numeric(age))
  }
  
  # fit linear models ----
  cat('Estimating delta beta...\n')
  db <- matrix(NA, nrow=nrow(beta),ncol=2+ncol(pheno))
  rownames(db) <- rownames(beta)
  colnames(db) <- c('db_epithelial','db_immune',colnames(pheno))
  
  ic.control <- pheno$ic[which(pheno$type=='Control')]
  ic.case <- pheno$ic[which(pheno$type!='Control')]
  
  pB <- txtProgressBar(min=1,max=nrow(beta), width =50L, style = 3)
  for (i in 1:nrow(beta)){
    setTxtProgressBar(pB, i)
    
    ldat <- data.frame(beta=as.numeric(beta[i,]),
                       pheno)
    
    lfit <- lm(beta ~ ., data=ldat)
    db[i,3:ncol(db)] <- summary(lfit)$coefficients[2:(ncol(pheno)+1),4]
    
    beta.control <- as.numeric(beta[i,which(pheno$type=='Control')])
    beta.case <- as.numeric(beta[i,which(pheno$type!='Control')])
    
    dat.control <- data.frame(beta=beta.control, ic=ic.control)
    dat.case <- data.frame(beta=beta.case, ic=ic.case)
    
    fit.control <- lm(beta ~ ic, data=dat.control)
    fit.case <- lm(beta ~ ic, data=dat.case)
    
    delta_beta_epithelial <- round(fit.case$coefficients[1]-fit.control$coefficients[1],digits=3)
    delta_beta_immune <- round(-fit.control$coefficients[1] - fit.control$coefficients[2]
                               +fit.case$coefficients[1] + fit.case$coefficients[2],digits=3)
    db[i,1] <- delta_beta_epithelial
    db[i,2] <- delta_beta_immune
    
  }
  close(pB)
  cat('done\n\n')

  # plot top 100 epithelial delta-betas -----  
  cols <- c("#83A58C", "#DF8A5A")
  names_sub <- stringr::str_split(names, "_", simplify = T)[,2]
  
  cat('Plotting top 100 ', names_sub[1], ' CpGs...')
  
  tmp <- db |> 
    as.data.frame() |> 
    arrange(desc(abs(db_epithelial))) |> 
    dplyr::slice(1:100) |> 
    tibble::rownames_to_column("cg") |> 
    dplyr::rename(p_type = type,
                  p_ic = ic)
  
  if("age" %in% adjustment){
    tmp <- tmp |> 
      dplyr::rename(p_age = age)
  }
  
  beta_tmp <- data.frame(t(beta[tmp$cg,match(rownames(pheno), colnames(beta))]))
  # identical(colnames(beta_tmp), rownames(tmp))
  # identical(rownames(beta_tmp), rownames(pheno))
  
  pheno_tmp <- cbind(pheno, beta_tmp) |> 
    tidyr::pivot_longer(all_of(tmp$cg),
                        names_to = "cg",
                        values_to = "beta") |> 
    dplyr::full_join(tmp) 
  
  # Annotate gene
  anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19, lociNames = tmp$cg) |> 
    data.frame() |> 
    mutate(gene = stringr::str_split(UCSC_RefGene_Name, ";", simplify = T)[,1]) |> 
    tibble::rownames_to_column("cg") |> 
    select(cg, gene)
  
  pheno_tmp <- pheno_tmp |> 
    full_join(anno) |> 
    mutate(cg_label = case_when(gene != "" ~ paste0(cg, " (", gene,
                                                        ")\ndb = ", signif(db_epithelial, 3),
                                                        ",\np = ", signif(p_type, 3)),
                                    TRUE ~ paste0(cg,
                                                  "\ndb = ", signif(db_epithelial, 3),
                                                  ",\np = ", signif(p_type, 3))))
  
  pdf(file=paste(output,"/top-", names_sub[1], "-db.pdf",sep=""), width=4,height=4.5)
  
  for (i in 1:100){
    plot <- pheno_tmp |>
      dplyr::filter(cg == tmp$cg[i]) |> 
      ggplot(aes(x = ic,
                 y = beta,
                 colour = type)) +
      geom_point(size = 0.75) +
      geom_smooth(method = "lm",
                  se = F) +
      facet_wrap(~ cg_label) +
      scale_colour_manual(values = cols,
                          name = "") +
      theme_bw() +
      theme(strip.text = element_text(hjust = 0))
    
    if(i == 1){
      plot <- plot +
        theme(legend.position = "top")
    } else {
      plot <- plot +
        theme(legend.position = "none")
    }
    
    print(plot)
  }
    
  dev.off()
  
  cat(' done\n')
  
  # plot top 100 immune/lymphoid delta-betas -----  
  
  cat('Plotting top 100 ', names_sub[2], ' CpGs...')
  
  tmp <- db |> 
    as.data.frame() |> 
    arrange(desc(abs(db_immune))) |> 
    dplyr::slice(1:100) |> 
    tibble::rownames_to_column("cg") |> 
    dplyr::rename(p_type = type,
                  p_ic = ic)
  beta_tmp <- data.frame(t(beta[tmp$cg,match(rownames(pheno), colnames(beta))]))
  # identical(colnames(beta_tmp), tmp$cg)
  # identical(rownames(beta_tmp), rownames(pheno))
  
  pheno_tmp <- cbind(pheno, beta_tmp) |> 
    tidyr::pivot_longer(all_of(tmp$cg),
                        names_to = "cg",
                        values_to = "beta") |> 
    dplyr::full_join(tmp, by = "cg") 
  
  # Annotate gene
  anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19, lociNames = tmp$cg) |> 
    data.frame() |> 
    mutate(gene = stringr::str_split(UCSC_RefGene_Name, ";", simplify = T)[,1]) |> 
    tibble::rownames_to_column("cg") |> 
    select(cg, gene)
  
  pheno_tmp <- pheno_tmp |> 
    full_join(anno) |> 
    mutate(cg_label = case_when(gene != "" ~ paste0(cg, " (", gene,
                                                        ")\ndb = ", signif(db_immune, 3),
                                                        ",\np = ", signif(p_type, 3)),
                                    TRUE ~ paste0(cg,
                                                  "\ndb = ", signif(db_immune, 3),
                                                  ",\np = ", signif(p_type, 3))))
  
  pdf(file=paste(output,"/top-", names_sub[2], "-db.pdf",sep=""), width=4,height=4.5)
  
  for (i in 1:100){
    plot <- pheno_tmp |>
      dplyr::filter(cg == tmp$cg[i]) |> 
      ggplot(aes(x = ic,
                 y = beta,
                 colour = type)) +
      geom_point(size = 0.75) +
      geom_smooth(method = "lm",
                  se = F) +
      facet_wrap(~ cg_label) +
      scale_colour_manual(values = cols,
                          name = "") +
      theme_bw() +
      theme(strip.text = element_text(hjust = 0))
    
    if(i == 1){
      plot <- plot +
        theme(legend.position = "top")
    } else {
      plot <- plot +
        theme(legend.position = "none")
    }
    
    print(plot)
  }
  
  dev.off()
  
  names
  
  colnames(db) <- c(names,colnames(pheno))
  saveRDS(db, file=paste(output,"/delta-beta.Rds",sep=""))
  
  cat(' done\n\n')
  cat('Session info:\n\n')
  sessionInfo()
}
