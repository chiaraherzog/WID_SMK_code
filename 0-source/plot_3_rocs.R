plot_3_rocs <- function(type,
                        index1, index2, index3,
                        col1 = "black", col2 = "grey40", col3 = "#E9BA6C",
                        title1 = "index1", title2 = "index2", title3 = "index3",
                        style = "default"){
  
  #--------------------------------#
  # make annotation labels
  auc <- round(as.numeric(roc(type,index1, quiet=T, ci=T)$auc),digits=2)
  cil <- round(as.numeric(roc(type,index1, quiet=T, ci=T)$ci[1]),digits=2)
  ciu <- round(as.numeric(roc(type,index1, quiet=T, ci=T)$ci[3]),digits=2)
  anno1 <- paste('<b>AUC (<i>', title1, '</i>) = ',auc,'</b><br>(95% CI: ', cil,'-',ciu,')',sep='')
  
  auc <- round(as.numeric(roc(type,index2, quiet=T, ci=T)$auc),digits=2)
  cil <- round(as.numeric(roc(type,index2, quiet=T, ci=T)$ci[1]),digits=2)
  ciu <- round(as.numeric(roc(type,index2, quiet=T, ci=T)$ci[3]),digits=2)
  anno2 <- paste('<b>AUC (<i>', title2, '</i>) = ',auc,'</b><br>(95% CI: ', cil,'-',ciu,')',sep='')
  
  auc <- round(as.numeric(roc(type,index3, quiet=T, ci=T)$auc),digits=2)
  cil <- round(as.numeric(roc(type,index3, quiet=T, ci=T)$ci[1]),digits=2)
  ciu <- round(as.numeric(roc(type,index3, quiet=T, ci=T)$ci[3]),digits=2)
  anno3 <- paste('<b>AUC (', title3, ') = ',auc,'</b><br>(95% CI: ', cil,'-',ciu,')',sep='')
  
  #--------------------------------#
  
  roc1 <- roc(type, index1)
  roc2 <- roc(type, index2)
  roc3 <- roc(type, index3)
  title1 <- as.character(title1)
  title2 <- as.character(title2)
  
  p<- ggplot() +
    geom_path(aes(x=1-roc1$specificities,
                  y=(roc1$sensitivities)),
              colour = col1,
              linewidth = 0.7) +
    geom_path(aes(x=1-roc2$specificities,
                  y=(roc2$sensitivities)),
              colour = col2,
              linewidth = 0.7) +
    geom_path(aes(x=1-roc3$specificities,
                  y=(roc3$sensitivities)),
              colour = col3,
              linewidth = 0.7) +
    annotate("segment", x = 0, y = 0,
             xend = 1, yend = 1,
             colour = "gray60") +
    xlab("1 - Specificity") +
    ylab("Sensitivity") +
    theme_minimal() +
    theme(plot.title = element_text(size=10),
          aspect.ratio = 1) +
    geom_textbox(aes(x = 0.95, y = 0.4, label = anno1),
                 hjust = 1, fill = col1, size = 3.5,
                 width = unit(2.3, unit = "in"),
                 colour = "white") +
    geom_textbox(aes(x = 0.95, y = 0.25, label = anno2),
                 hjust = 1, fill = col2, size = 3.5,
                 width = unit(2.3, unit = "in"),
                 colour = "white") +
    geom_textbox(aes(x = 0.95, y = 0.1, label = anno3),
                 hjust = 1, fill = col3, size = 3.5,
                 width = unit(2.3, unit = "in"),
                 colour = "white")
  print(p)
  
}