# Comparison of AUCS for different indices and types
smallAUCplot <- function(data, indices, types, typevar = "type", base = "Control",
                    cols = c("red", "blue"),
                    reorder = F,
                    reorder.index = F,
                    reorder.index.levels = NULL,
                    reverse_base = F,
                    labels = NULL,
                    direction = "no_dir",
                    ylim = c(0, 1),
                    aspect.ratio = 1,
                    alpha = 0.1,
                    print.annotation = F){
  
  library(pROC)
  
  df <- data.frame(matrix(nrow = length(indices)*length(types),
                          ncol = 7))
  colnames(df) <- c("index", "AUC", "lo", "hi", "type", "annotation", "sig")
  df$index <- indices
  
  df$type <- rep(types, each = length(indices))
  df$type <- factor(df$type, levels = types)
  
  for (t in 1:length(types)){
    for (i in 1:length(indices)){
      
      tmp <- data |> 
        dplyr::filter(UQ(sym(typevar)) %in% c(base, types[t])) |> 
        droplevels()
      
      
      if(reverse_base == T){
        tmp[[typevar]] <- factor(tmp[[typevar]], levels = c(types[t], base))
      }
      
      if(direction == "no_dir"){
        roc <- pROC::roc(tmp[[typevar]], tmp[[indices[i]]])
      } else if (direction == "strict") {
        if(indices[i] %in% c("WID_SMK450_epithelial_hypoM",
                             "WID_SMK450_immune_hypoM")){
          roc <- pROC::roc(tmp[[typevar]], tmp[[indices[i]]], direction = ">")
        } else {
          roc <- pROC::roc(tmp[[typevar]], tmp[[indices[i]]], direction = "<")
        }
      } else {
        roc <- pROC::roc(tmp[[typevar]], tmp[[indices[i]]], direction = "<")
        
      }
      ci <- ci(roc)
      
      df[df$index==indices[i] & df$type==types[t],]$AUC <- as.numeric(ci)[2]
      df[df$index==indices[i] & df$type==types[t],]$lo <- as.numeric(ci)[1]
      df[df$index==indices[i] & df$type==types[t],]$hi <- as.numeric(ci)[3]
      
      if(print.annotation!=F){
        df[df$index==indices[i] & df$type==types[t],]$annotation <- as.character(paste0("AUC=",
                                                                                        signif(as.numeric(ci)[2], 2),
                                                                                        "\n(", signif(as.numeric(ci)[1],2),
                                                                                        "-",
                                                                                        signif(as.numeric(ci)[3], 2),
                                                                                        ")"))
        
        df[df$index==indices[i] & df$type==types[t],]$sig <- case_when((as.numeric(ci)[2] > 0.5 & as.numeric(ci)[1] > 0.5) | (as.numeric(ci)[2] < 0.5 & as.numeric(ci)[3] < 0.5) ~ "yes",
                                                                       TRUE ~ "no")
        df[df$index==indices[i] & df$type==types[t],]$annotation <- ifelse(df[df$index==indices[i] & df$type==types[t],]$sig  == "yes", paste0(df[df$index==indices[i] & df$type==types[t],]$annotation, "*"), paste0(df[df$index==indices[i] & df$type==types[t],]$annotation))
      } else {
        df[df$index==indices[i] & df$type==types[t],]$annotation <- ""
      }
    } 
  }
  
  if(reverse_base == T){
    df$type <- rep(paste0(base, "\nversus ", types), each = length(indices))
    df$type <- factor(df$type, levels = paste0(base, "\nversus ", types))
  }
  
  if(reorder == T){
    
    
    p <- df |> 
      ggplot(aes(x = index,
                 y = AUC,
                 group = type)) +
      geom_point(size = 3,
                 aes(colour = type)) +
      geom_errorbar(aes(ymin = lo,
                        ymax = hi,
                        colour = type),
                    width = 0) +
      geom_ribbon(aes(ymin = lo,
                      ymax = hi,
                      fill = type),
                  alpha = alpha) +
      geom_text(aes(y = 0.38,
                    label = annotation,
                    colour = type,
                    fontface = ifelse(sig == "yes", 2, 1)),
                size = 3,
                nudge_y = case_when(df$type == types[1] & length(types) == 2 ~ -0.05,
                                    df$type == types[2] & length(types) == 2 ~ +0.05,
                                    df$type == types[1] & length(types) == 4 ~ 0.2,
                                    df$type == types[2] & length(types) == 4 ~ 0.08,
                                    df$type == types[3] & length(types) == 4 ~ -0.08,
                                    df$type == types[4] & length(types) == 4 ~ -0.2,
                                    reverse_base == T & grepl(types[1], df$type) ~ 0.1,
                                    reverse_base == T & grepl(types[2], df$type) ~ -0.1),
                show.legend = F) +
      geom_hline(yintercept = 0.5,
                 linetype = "dotted") +
      theme_bw() +
      ylim(ylim) +
      scale_colour_manual(name = "",
                          values = cols,
                          aesthetics = c("colour", "fill")) +
      theme(legend.position = "top",
            aspect.ratio = aspect.ratio,
            axis.text.x = element_text(angle = 60,
                                       hjust = 1)) +
      xlab("") +
      scale_x_discrete(breaks = indices, labels = labels)
    
  } else if(reorder.index == T){
    
    p <- df |> 
      dplyr::mutate(index = factor(index, levels = reorder.index.levels)) |> 
      ggplot(aes(x = index,
                 y = AUC,
                 group = type)) +
      geom_point(size = 3,
                 aes(colour = type)) +
      geom_errorbar(aes(ymin = lo,
                        ymax = hi,
                        colour = type),
                    width = 0) +
      geom_ribbon(aes(ymin = lo,
                      ymax = hi,
                      fill = type),
                  alpha = alpha) +
      geom_text(aes(y = 0.38,
                    label = annotation,
                    colour = type,
                    fontface = ifelse(sig == "yes", 2, 1)),
                size = 3,
                nudge_y = case_when(df$type == types[1] & length(types) == 2 ~ -0.05,
                                    df$type == types[2] & length(types) == 2 ~ +0.05,
                                    df$type == types[1] & length(types) == 4 ~ 0.2,
                                    df$type == types[2] & length(types) == 4 ~ 0.08,
                                    df$type == types[3] & length(types) == 4 ~ -0.08,
                                    df$type == types[4] & length(types) == 4 ~ -0.2,
                                    reverse_base == T & grepl(types[1], df$type) ~ 0.1,
                                    reverse_base == T & grepl(types[2], df$type) ~ -0.1),
                show.legend = F) +
      geom_hline(yintercept = 0.5,
                 linetype = "dotted") +
      theme_bw() +
      ylim(ylim) +
      scale_colour_manual(name = "",
                          values = cols,
                          aesthetics = c("colour", "fill")) +
      theme(legend.position = "top",
            aspect.ratio = aspect.ratio,
            axis.text.x = element_text(angle = 60,
                                       hjust = 1)) +
      xlab("") +
      scale_x_discrete(breaks = indices, labels = labels)
    
  } else {
    
    if(is.null(labels)){
      labels <- indices
    }
    p <-  df |> 
      ggplot(aes(x = index,
                 y = AUC,
                 group = type)) +
      geom_point(size = 3,
                 aes(colour = type)) +
      geom_errorbar(aes(ymin = lo,
                        ymax = hi,
                        colour = type),
                    width = 0) +
      geom_ribbon(aes(ymin = lo,
                      ymax = hi,
                      fill = type),
                  alpha = alpha) +
      geom_hline(yintercept = 0.5,
                 linetype = "dotted") +
      geom_text(aes(y = 0.38,
                    label = annotation,
                    colour = type,
                    fontface = ifelse(sig == "yes", 2, 1)),
                size = 3,
                nudge_y = case_when(df$type == types[1] & length(types) == 2 ~ -0.05,
                                    df$type == types[2] & length(types) == 2 ~ 0.05,
                                    df$type == types[1] & length(types) == 4 ~ 0.2,
                                    df$type == types[2] & length(types) == 4 ~ 0.08,
                                    df$type == types[3] & length(types) == 4 ~ -0.08,
                                    df$type == types[4] & length(types) == 4 ~ -0.2,
                                    reverse_base == T & grepl(types[1], df$type) ~ 0.1,
                                    reverse_base == T & grepl(types[2], df$type) ~ -0.1),
                show.legend = F) +
      theme_bw() +
      ylim(ylim) +
      scale_colour_manual(name = "",
                          values = cols,
                          aesthetics = c("colour", "fill")) +
      theme(legend.position = "top",
            aspect.ratio = aspect.ratio,
            axis.text.x = element_text(angle = 60,
                                       hjust = 1)) +
      xlab("")  +
      scale_x_discrete(limits = rev(indices),
                       labels = rev(labels))
  }
  
  return(list(plot = p,
              df = df))
}

