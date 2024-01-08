surv_function <- function(data, adj = T){
  
  library(ggtext)
  library(survival)
  cols <- MetBrewer::met.brewer("Hokusai1", n = 9)
  
  # generate dataframes
  
  # Overall survival
  os <- data.frame(matrix(nrow = length(unique(data$project))*length(unique(data$index)),
                          ncol = 18))
  colnames(os) <- c("project", "index",
                    "hr_unadj", "cilo_unadj", "cihi_unadj", "p_unadj",
                    "hr_adj1", "cilo_adj1", "cihi_adj1", "p_adj1", 
                    "hr_adj2", "cilo_adj2", "cihi_adj2", "p_adj2",
                    "hr_adj3", "cilo_adj3", "cihi_adj3", "p_adj3")
  os$project <- rep(unique(data$project), each = length(unique(data$index)))
  os$index <- unique(data$index)
  
  # Progression
  pf <- os
  
  # Metastasis or recurrence
  metastasis <- os
  
  # Distant metastasis
  distmet <- os
  
  # Local recurrence
  localrec <- os
  
  # Correlation with stage
  stage <- pf[,1:6]
  colnames(stage) <- c("project", "index", "r_unadj", "p_unadj", "r_adj1", "p_adj1")
  
  for(p in unique(data$project)){
    
    for (i in unique(data$index)){
      
      # OS 
      # unadj
      fit <- coxph(Surv(os_time, os) ~ value_scale, data = data[data$index == i & data$project == p,])
      os[os$project==p & os$index==i,]$hr_unadj <- summary(fit)$conf.int[1]
      os[os$project==p & os$index==i,]$cilo_unadj <- summary(fit)$conf.int[3]
      os[os$project==p & os$index==i,]$cihi_unadj <- summary(fit)$conf.int[4]
      os[os$project==p & os$index==i,]$p_unadj <- summary(fit)$coefficients[1,5]
      
      # adj (minimal)
      fit <- coxph(Surv(os_time, os) ~ value_scale + sex + age + ic, data = data[data$index == i & data$project == p,])
      os[os$project==p & os$index==i,]$hr_adj1 <- summary(fit)$conf.int[1,1]
      os[os$project==p & os$index==i,]$cilo_adj1 <- summary(fit)$conf.int[1,3]
      os[os$project==p & os$index==i,]$cihi_adj1 <- summary(fit)$conf.int[1,4]
      os[os$project==p & os$index==i,]$p_adj1 <- summary(fit)$coefficients[1,5]
      
      # adj (+ stage)
      fit <- coxph(Surv(os_time, os) ~ value_scale + sex + age + ic + stage, data = data[data$index == i & data$project == p,])
      os[os$project==p & os$index==i,]$hr_adj2 <- summary(fit)$conf.int[1,1]
      os[os$project==p & os$index==i,]$cilo_adj2 <- summary(fit)$conf.int[1,3]
      os[os$project==p & os$index==i,]$cihi_adj2 <- summary(fit)$conf.int[1,4]
      os[os$project==p & os$index==i,]$p_adj2 <- summary(fit)$coefficients[1,5]
      
      # adj (+ smoking)
      fit <- coxph(Surv(os_time, os) ~ value_scale + sex + age + ic + smoking_history, data = data[data$index == i & data$project == p,])
      os[os$project==p & os$index==i,]$hr_adj3 <- summary(fit)$conf.int[1,1]
      os[os$project==p & os$index==i,]$cilo_adj3 <- summary(fit)$conf.int[1,3]
      os[os$project==p & os$index==i,]$cihi_adj3 <- summary(fit)$conf.int[1,4]
      os[os$project==p & os$index==i,]$p_adj3 <- summary(fit)$coefficients[1,5]
      
      # PF 
      # unadj
      fit <- coxph(Surv(pfi_time, pfi) ~ value_scale, data = data[data$index == i & data$project == p,])
      pf[pf$project==p & pf$index==i,]$hr_unadj <- summary(fit)$conf.int[1]
      pf[pf$project==p & pf$index==i,]$cilo_unadj <- summary(fit)$conf.int[3]
      pf[pf$project==p & pf$index==i,]$cihi_unadj <- summary(fit)$conf.int[4]
      pf[pf$project==p & pf$index==i,]$p_unadj <- summary(fit)$coefficients[1,5]
      
      # adj (minimal)
      fit <- coxph(Surv(pfi_time, pfi) ~ value_scale + sex + age + ic, data = data[data$index == i & data$project == p,])
      pf[pf$project==p & pf$index==i,]$hr_adj1 <- summary(fit)$conf.int[1,1]
      pf[pf$project==p & pf$index==i,]$cilo_adj1 <- summary(fit)$conf.int[1,3]
      pf[pf$project==p & pf$index==i,]$cihi_adj1 <- summary(fit)$conf.int[1,4]
      pf[pf$project==p & pf$index==i,]$p_adj1 <- summary(fit)$coefficients[1,5]
      
      # adj (+ stage)
      fit <- coxph(Surv(pfi_time, pfi) ~ value_scale + sex + age + ic + stage, data = data[data$index == i & data$project == p,])
      pf[pf$project==p & pf$index==i,]$hr_adj2 <- summary(fit)$conf.int[1,1]
      pf[pf$project==p & pf$index==i,]$cilo_adj2 <- summary(fit)$conf.int[1,3]
      pf[pf$project==p & pf$index==i,]$cihi_adj2 <- summary(fit)$conf.int[1,4]
      pf[pf$project==p & pf$index==i,]$p_adj2 <- summary(fit)$coefficients[1,5]
      
      # adj (+ smoking_history)
      fit <- coxph(Surv(pfi_time, pfi) ~ value_scale + sex + age + ic + smoking_history, data = data[data$index == i & data$project == p,])
      pf[pf$project==p & pf$index==i,]$hr_adj3 <- summary(fit)$conf.int[1,1]
      pf[pf$project==p & pf$index==i,]$cilo_adj3 <- summary(fit)$conf.int[1,3]
      pf[pf$project==p & pf$index==i,]$cihi_adj3 <- summary(fit)$conf.int[1,4]
      pf[pf$project==p & pf$index==i,]$p_adj3 <- summary(fit)$coefficients[1,5]
      
      # metastasis 
      # unadj
      fit <- coxph(Surv(recur_or_met_time, recur_or_met) ~ value_scale, data = data[data$index == i & data$project == p,])
      metastasis[metastasis$project==p & metastasis$index==i,]$hr_unadj <- summary(fit)$conf.int[1]
      metastasis[metastasis$project==p & metastasis$index==i,]$cilo_unadj <- summary(fit)$conf.int[3]
      metastasis[metastasis$project==p & metastasis$index==i,]$cihi_unadj <- summary(fit)$conf.int[4]
      metastasis[metastasis$project==p & metastasis$index==i,]$p_unadj <- summary(fit)$coefficients[1,5]
      
      # adj (minimal)
      fit <- coxph(Surv(recur_or_met_time, recur_or_met) ~ value_scale + sex + age + ic, data = data[data$index == i & data$project == p,])
      metastasis[metastasis$project==p & metastasis$index==i,]$hr_adj1 <- summary(fit)$conf.int[1,1]
      metastasis[metastasis$project==p & metastasis$index==i,]$cilo_adj1 <- summary(fit)$conf.int[1,3]
      metastasis[metastasis$project==p & metastasis$index==i,]$cihi_adj1 <- summary(fit)$conf.int[1,4]
      metastasis[metastasis$project==p & metastasis$index==i,]$p_adj1 <- summary(fit)$coefficients[1,5]
      
      # adj (+ stage)
      fit <- coxph(Surv(recur_or_met_time, recur_or_met) ~ value_scale + sex + age + ic + stage, data = data[data$index == i & data$project == p,])
      metastasis[metastasis$project==p & metastasis$index==i,]$hr_adj2 <- summary(fit)$conf.int[1,1]
      metastasis[metastasis$project==p & metastasis$index==i,]$cilo_adj2 <- summary(fit)$conf.int[1,3]
      metastasis[metastasis$project==p & metastasis$index==i,]$cihi_adj2 <- summary(fit)$conf.int[1,4]
      metastasis[metastasis$project==p & metastasis$index==i,]$p_adj2 <- summary(fit)$coefficients[1,5]
      
      # adj (+ smoking)
      fit <- coxph(Surv(recur_or_met_time, recur_or_met) ~ value_scale + sex + age + ic + smoking_history, data = data[data$index == i & data$project == p,])
      metastasis[metastasis$project==p & metastasis$index==i,]$hr_adj3 <- summary(fit)$conf.int[1,1]
      metastasis[metastasis$project==p & metastasis$index==i,]$cilo_adj3 <- summary(fit)$conf.int[1,3]
      metastasis[metastasis$project==p & metastasis$index==i,]$cihi_adj3 <- summary(fit)$conf.int[1,4]
      metastasis[metastasis$project==p & metastasis$index==i,]$p_adj3 <- summary(fit)$coefficients[1,5]
      
      # distant metastasis 
      # unadj
      fit <- coxph(Surv(distant_met_time, distant_met) ~ value_scale, data = data[data$index == i & data$project == p,])
      distmet[distmet$project==p & distmet$index==i,]$hr_unadj <- summary(fit)$conf.int[1]
      distmet[distmet$project==p & distmet$index==i,]$cilo_unadj <- summary(fit)$conf.int[3]
      distmet[distmet$project==p & distmet$index==i,]$cihi_unadj <- summary(fit)$conf.int[4]
      distmet[distmet$project==p & distmet$index==i,]$p_unadj <- summary(fit)$coefficients[1,5]
      
      # adj (minimal)
      fit <- coxph(Surv(distant_met_time, distant_met) ~ value_scale + sex + age + ic, data = data[data$index == i & data$project == p,])
      distmet[distmet$project==p & distmet$index==i,]$hr_adj1 <- summary(fit)$conf.int[1,1]
      distmet[distmet$project==p & distmet$index==i,]$cilo_adj1 <- summary(fit)$conf.int[1,3]
      distmet[distmet$project==p & distmet$index==i,]$cihi_adj1 <- summary(fit)$conf.int[1,4]
      distmet[distmet$project==p & distmet$index==i,]$p_adj1 <- summary(fit)$coefficients[1,5]
      
      # adj (+ stage)
      fit <- coxph(Surv(distant_met_time, distant_met) ~ value_scale + sex + age + ic + stage, data = data[data$index == i & data$project == p,])
      distmet[distmet$project==p & distmet$index==i,]$hr_adj2 <- summary(fit)$conf.int[1,1]
      distmet[distmet$project==p & distmet$index==i,]$cilo_adj2 <- summary(fit)$conf.int[1,3]
      distmet[distmet$project==p & distmet$index==i,]$cihi_adj2 <- summary(fit)$conf.int[1,4]
      distmet[distmet$project==p & distmet$index==i,]$p_adj2 <- summary(fit)$coefficients[1,5]
      
      # adj (+ smoking)
      fit <- coxph(Surv(distant_met_time, distant_met) ~ value_scale + sex + age + ic + smoking_history, data = data[data$index == i & data$project == p,])
      distmet[distmet$project==p & distmet$index==i,]$hr_adj3 <- summary(fit)$conf.int[1,1]
      distmet[distmet$project==p & distmet$index==i,]$cilo_adj3 <- summary(fit)$conf.int[1,3]
      distmet[distmet$project==p & distmet$index==i,]$cihi_adj3 <- summary(fit)$conf.int[1,4]
      distmet[distmet$project==p & distmet$index==i,]$p_adj3 <- summary(fit)$coefficients[1,5]
      
      # local recurrence 
      # unadj
      fit <- coxph(Surv(regional_rec_time, regional_rec) ~ value_scale, data = data[data$index == i & data$project == p,])
      localrec[localrec$project==p & localrec$index==i,]$hr_unadj <- summary(fit)$conf.int[1]
      localrec[localrec$project==p & localrec$index==i,]$cilo_unadj <- summary(fit)$conf.int[3]
      localrec[localrec$project==p & localrec$index==i,]$cihi_unadj <- summary(fit)$conf.int[4]
      localrec[localrec$project==p & localrec$index==i,]$p_unadj <- summary(fit)$coefficients[1,5]
      
      # adj (minimal)
      fit <- coxph(Surv(regional_rec_time, regional_rec) ~ value_scale + sex + age + ic, data = data[data$index == i & data$project == p,])
      localrec[localrec$project==p & localrec$index==i,]$hr_adj1 <- summary(fit)$conf.int[1,1]
      localrec[localrec$project==p & localrec$index==i,]$cilo_adj1 <- summary(fit)$conf.int[1,3]
      localrec[localrec$project==p & localrec$index==i,]$cihi_adj1 <- summary(fit)$conf.int[1,4]
      localrec[localrec$project==p & localrec$index==i,]$p_adj1 <- summary(fit)$coefficients[1,5]
      
      # adj (+ stage)
      fit <- coxph(Surv(regional_rec_time, regional_rec) ~ value_scale + sex + age + ic + stage, data = data[data$index == i & data$project == p,])
      localrec[localrec$project==p & localrec$index==i,]$hr_adj2 <- summary(fit)$conf.int[1,1]
      localrec[localrec$project==p & localrec$index==i,]$cilo_adj2 <- summary(fit)$conf.int[1,3]
      localrec[localrec$project==p & localrec$index==i,]$cihi_adj2 <- summary(fit)$conf.int[1,4]
      localrec[localrec$project==p & localrec$index==i,]$p_adj2 <- summary(fit)$coefficients[1,5]
      
      # adj (+ smoking)
      fit <- coxph(Surv(regional_rec_time, regional_rec) ~ value_scale + sex + age + ic + smoking_history, data = data[data$index == i & data$project == p,])
      localrec[localrec$project==p & localrec$index==i,]$hr_adj3 <- summary(fit)$conf.int[1,1]
      localrec[localrec$project==p & localrec$index==i,]$cilo_adj3 <- summary(fit)$conf.int[1,3]
      localrec[localrec$project==p & localrec$index==i,]$cihi_adj3 <- summary(fit)$conf.int[1,4]
      localrec[localrec$project==p & localrec$index==i,]$p_adj3 <- summary(fit)$coefficients[1,5]
      
      # Association with stage
      # unadj
      corrtest <- cor.test(data[data$index == i & data$project == p,]$stage,
                           data[data$index == i & data$project == p,]$value)
      stage[stage$project==p & stage$index==i,]$r_unadj <- corrtest$estimate
      stage[stage$project==p & stage$index==i,]$p_unadj <- corrtest$p.value
      
      # Adj 1
      data$value_adj <- NA
      fit <- lm(value ~ sex + age + ic + smoking_history, data = data[data$index == i & data$project == p,])
      data[data$index == i & data$project == p,]$value_adj <- data[data$index == i & data$project == p,]$value - predict(fit, newdata = data[data$index == i & data$project == p,])
      corrtest <- cor.test(data[data$index == i & data$project == p,]$stage,
                           data[data$index == i & data$project == p,]$value_adj)
      stage[stage$project==p & stage$index==i,]$r_adj1 <- corrtest$estimate
      stage[stage$project==p & stage$index==i,]$p_adj1 <- corrtest$p.value
      
      }
  }
  
  
  out1 <- list(os = os,
               pf = pf,
               met_rec = metastasis,
               distant_met = distmet,
               local_rec = localrec,
               stage = stage)
  
  
  # Plots
  plot_os <- out1$os |> 
    tidyr::pivot_longer(hr_unadj:p_adj3,
                        names_to = "metric",
                        values_to = "value") |> 
    tidyr::separate(metric, sep = "_",
                    into = c("metric", "adjustment")) |> 
    tidyr::pivot_wider(id_cols = c(project, index, adjustment),
                       names_from = metric,
                       values_from = value) |> 
    dplyr::mutate(adjustment = case_when(adjustment == "unadj" ~ "<b>Unadjusted</b><br><br>",
                                         adjustment == "adj1" ~ "<b>Model 1</b><br>(sex, age, and<br>immune cell proportion)",
                                         adjustment == "adj2" ~ "<b>Model 2</b><br>Model 1 + stage<br>",
                                         adjustment == "adj3" ~ "<b>Model 3</b><br>Model 1 + smoking history"),
                  adjustment = factor(adjustment, levels = c("<b>Unadjusted</b><br><br>",
                                                             "<b>Model 1</b><br>(sex, age, and<br>immune cell proportion)",
                                                             "<b>Model 2</b><br>Model 1 + stage<br>",
                                                             "<b>Model 3</b><br>Model 1 + smoking history")),
                  sig = ifelse(p < 0.05, "significant association", NA),
                  index = case_when(grepl("WID_SMK450", index) ~ gsub("_", " ", index),
                                    index == "smoking_mrs" ~ "smoking methylation risk score",
                                    index == "cg05575921" ~ "cg05575921"),
                  index = ifelse(grepl("WID SMK450", index), gsub("WID SMK450 ", "", index), index),
                  index = factor(index, levels = rev(c("cg05575921", "smoking methylation risk score",
                                                   "epithelial hypoM", "immune hypoM", "proximal epithelial hyperM", "distal epithelial hyperM")))) |> 
    ggplot(aes(x = index,
               y = hr)) +
    geom_hline(yintercept = 1, linetype = "dashed",
               colour = "grey50") +
    geom_errorbar(aes(ymin = cilo,
                      ymax = cihi),
                  width = 0.15,
                  colour = "grey60") +
    geom_point(shape = 18,
               aes(colour = sig,
                   size = sig)) +
    facet_grid(project~adjustment) +
    scale_colour_manual(values = cols[3],na.value = "grey50",
                        name = '',
                        breaks = c("significant association")) +
    scale_size_manual(values = 4, na.value = 3.5,
                      name = "",
                      breaks = c("significant association")) +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "top",
          axis.title.x = element_markdown(),
          strip.text = element_markdown(hjust = 0)) +
    labs(x = "", y = "<b>log (HR)</b><br>per standard deviation", title = "Overall survival")
  
  plot_pf <- out1$pf |> 
    tidyr::pivot_longer(hr_unadj:p_adj3,
                        names_to = "metric",
                        values_to = "value") |> 
    tidyr::separate(metric, sep = "_",
                    into = c("metric", "adjustment")) |> 
    tidyr::pivot_wider(id_cols = c(project, index, adjustment),
                       names_from = metric,
                       values_from = value) |> 
    dplyr::mutate(adjustment = case_when(adjustment == "unadj" ~ "<b>Unadjusted</b><br><br>",
                                         adjustment == "adj1" ~ "<b>Model 1</b><br>(sex, age, and<br>immune cell proportion)",
                                         adjustment == "adj2" ~ "<b>Model 2</b><br>Model 1 + stage<br>",
                                         adjustment == "adj3" ~ "<b>Model 3</b><br>Model 1 + smoking history"),
                  adjustment = factor(adjustment, levels = c("<b>Unadjusted</b><br><br>", "<b>Model 1</b><br>(sex, age, and<br>immune cell proportion)", "<b>Model 2</b><br>Model 1 + stage<br>", "<b>Model 3</b><br>Model 1 + smoking history")),
                  sig = ifelse(p < 0.05, "significant association", NA),
                  index = case_when(grepl("WID_SMK450", index) ~ gsub("_", " ", index),
                                    index == "smoking_mrs" ~ "smoking methylation risk score",
                                    index == "cg05575921" ~ "cg05575921"),
                  index = ifelse(grepl("WID SMK450", index), gsub("WID SMK450 ", "", index), index),
                  index = factor(index, levels = rev(c("cg05575921", "smoking methylation risk score",
                                                       "epithelial hypoM", "immune hypoM", "proximal epithelial hyperM", "distal epithelial hyperM")))) |> 
    ggplot(aes(x = index,
               y = hr)) +
    geom_hline(yintercept = 1, linetype = "dashed",
               colour = "grey50") +
    geom_errorbar(aes(ymin = cilo,
                      ymax = cihi),
                  width = 0.15,
                  colour = "grey60") +
    geom_point(shape = 18,
               aes(colour = sig,
                   size = sig)) +
    facet_grid(project~adjustment) +
    scale_colour_manual(values = cols[3],na.value = "grey50",
                        name = '',
                        breaks = c("significant association")) +
    scale_size_manual(values = 4, na.value = 3.5,
                      name = "",
                      breaks = c("significant association")) +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "top",
          axis.title.x = element_markdown(),
          strip.text = element_markdown(hjust = 0)) +
    labs(x = "", y = "<b>log (HR)</b><br>per standard deviation", title = "Progression-free interval")
  
  plot_met <- out1$met_rec |> 
    tidyr::pivot_longer(hr_unadj:p_adj3,
                        names_to = "metric",
                        values_to = "value") |> 
    tidyr::separate(metric, sep = "_",
                    into = c("metric", "adjustment")) |> 
    tidyr::pivot_wider(id_cols = c(project, index, adjustment),
                       names_from = metric,
                       values_from = value) |> 
    dplyr::mutate(adjustment = case_when(adjustment == "unadj" ~ "<b>Unadjusted</b><br><br>",
                                         adjustment == "adj1" ~ "<b>Model 1</b><br>(sex, age, and<br>immune cell proportion)",
                                         adjustment == "adj2" ~ "<b>Model 2</b><br>Model 1 + stage<br>",
                                         adjustment == "adj3" ~ "<b>Model 3</b><br>Model 1 + smoking history"),
                  adjustment = factor(adjustment, levels = c("<b>Unadjusted</b><br><br>", "<b>Model 1</b><br>(sex, age, and<br>immune cell proportion)", "<b>Model 2</b><br>Model 1 + stage<br>", "<b>Model 3</b><br>Model 1 + smoking history")),
                  sig = ifelse(p < 0.05, "significant association", NA),
                  index = case_when(grepl("WID_SMK450", index) ~ gsub("_", " ", index),
                                    index == "smoking_mrs" ~ "smoking methylation risk score",
                                    index == "cg05575921" ~ "cg05575921"),
                  index = ifelse(grepl("WID SMK450", index), gsub("WID SMK450 ", "", index), index),
                  index = factor(index, levels = rev(c("cg05575921", "smoking methylation risk score",
                                                       "epithelial hypoM", "immune hypoM", "proximal epithelial hyperM", "distal epithelial hyperM")))) |> 
    ggplot(aes(x = index,
               y = hr)) +
    geom_hline(yintercept = 1, linetype = "dashed",
               colour = "grey50") +
    geom_errorbar(aes(ymin = cilo,
                      ymax = cihi),
                  width = 0.15,
                  colour = "grey60") +
    geom_point(shape = 18,
               aes(colour = sig,
                   size = sig)) +
    facet_grid(project~adjustment) +
    scale_colour_manual(values = cols[3],na.value = "grey50",
                        name = '',
                        breaks = c("significant association")) +
    scale_size_manual(values = 4, na.value = 3.5,
                      name = "",
                      breaks = c("significant association")) +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "top",
          axis.title.x = element_markdown(),
          strip.text = element_markdown(hjust = 0)) +
    labs(x = "", y = "<b>log (HR)</b><br>per standard deviation", title = "Metastasis or recurrence")
  
  plot_distmet <- out1$distant_met |> 
    tidyr::pivot_longer(hr_unadj:p_adj3,
                        names_to = "metric",
                        values_to = "value") |> 
    tidyr::separate(metric, sep = "_",
                    into = c("metric", "adjustment")) |> 
    tidyr::pivot_wider(id_cols = c(project, index, adjustment),
                       names_from = metric,
                       values_from = value) |> 
    dplyr::mutate(adjustment = case_when(adjustment == "unadj" ~ "<b>Unadjusted</b><br><br>",
                                         adjustment == "adj1" ~ "<b>Model 1</b><br>(sex, age, and<br>immune cell proportion)",
                                         adjustment == "adj2" ~ "<b>Model 2</b><br>Model 1 + stage<br>",
                                         adjustment == "adj3" ~ "<b>Model 3</b><br>Model 1 + smoking history"),
                  adjustment = factor(adjustment, levels = c("<b>Unadjusted</b><br><br>", "<b>Model 1</b><br>(sex, age, and<br>immune cell proportion)", "<b>Model 2</b><br>Model 1 + stage<br>", "<b>Model 3</b><br>Model 1 + smoking history")),
                  sig = ifelse(p < 0.05, "significant association", NA),
                  index = case_when(grepl("WID_SMK450", index) ~ gsub("_", " ", index),
                                    index == "smoking_mrs" ~ "smoking methylation risk score",
                                    index == "cg05575921" ~ "cg05575921"),
                  index = ifelse(grepl("WID SMK450", index), gsub("WID SMK450 ", "", index), index),
                  index = factor(index, levels = rev(c("cg05575921", "smoking methylation risk score",
                                                       "epithelial hypoM", "immune hypoM", "proximal epithelial hyperM", "distal epithelial hyperM")))) |> 
    ggplot(aes(x = index,
               y = hr)) +
    geom_hline(yintercept = 1, linetype = "dashed",
               colour = "grey50") +
    geom_errorbar(aes(ymin = cilo,
                      ymax = cihi),
                  width = 0.15,
                  colour = "grey60") +
    geom_point(shape = 18,
               aes(colour = sig,
                   size = sig)) +
    facet_grid(project~adjustment) +
    scale_colour_manual(values = cols[3],na.value = "grey50",
                        name = '',
                        breaks = c("significant association")) +
    scale_size_manual(values = 4, na.value = 3.5,
                      name = "",
                      breaks = c("significant association")) +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "top",
          axis.title.x = element_markdown(),
          strip.text = element_markdown(hjust = 0)) +
    labs(x = "", y = "<b>log (HR)</b><br>per standard deviation", title = "Distant metastasis")
  
  plot_localrec <- out1$local_rec |> 
    tidyr::pivot_longer(hr_unadj:p_adj3,
                        names_to = "metric",
                        values_to = "value") |> 
    tidyr::separate(metric, sep = "_",
                    into = c("metric", "adjustment")) |> 
    tidyr::pivot_wider(id_cols = c(project, index, adjustment),
                       names_from = metric,
                       values_from = value) |> 
    dplyr::mutate(adjustment = case_when(adjustment == "unadj" ~ "<b>Unadjusted</b><br><br>",
                                         adjustment == "adj1" ~ "<b>Model 1</b><br>(sex, age, and<br>immune cell proportion)",
                                         adjustment == "adj2" ~ "<b>Model 2</b><br>Model 1 + stage<br>",
                                         adjustment == "adj3" ~ "<b>Model 3</b><br>Model 1 + smoking history"),
                  adjustment = factor(adjustment, levels = c("<b>Unadjusted</b><br><br>", "<b>Model 1</b><br>(sex, age, and<br>immune cell proportion)", "<b>Model 2</b><br>Model 1 + stage<br>", "<b>Model 3</b><br>Model 1 + smoking history")),
                  sig = ifelse(p < 0.05, "significant association", NA),
                  index = case_when(grepl("WID_SMK450", index) ~ gsub("_", " ", index),
                                    index == "smoking_mrs" ~ "smoking methylation risk score",
                                    index == "cg05575921" ~ "cg05575921"),
                  index = ifelse(grepl("WID SMK450", index), gsub("WID SMK450 ", "", index), index),
                  index = factor(index, levels = rev(c("cg05575921", "smoking methylation risk score",
                                                       "epithelial hypoM", "immune hypoM", "proximal epithelial hyperM", "distal epithelial hyperM")))) |> 
    ggplot(aes(x = index,
               y = hr)) +
    geom_hline(yintercept = 1, linetype = "dashed",
               colour = "grey50") +
    geom_errorbar(aes(ymin = cilo,
                      ymax = cihi),
                  width = 0.15,
                  colour = "grey60") +
    geom_point(shape = 18,
               aes(colour = sig,
                   size = sig)) +
    facet_grid(project~adjustment) +
    scale_colour_manual(values = cols[3],na.value = "grey50",
                        name = '',
                        breaks = c("significant association")) +
    scale_size_manual(values = 4, na.value = 3.5,
                      name = "",
                      breaks = c("significant association")) +
    coord_flip() +
    theme_bw() +
    theme(legend.position = "top",
          axis.title.x = element_markdown(),
          strip.text = element_markdown(hjust = 0)) +
    labs(x = "", y = "<b>log (HR)</b><br>per standard deviation", title = "Regional recurrence")
  
  plot_stage <- out1$stage |> 
    tidyr::pivot_longer(r_unadj:p_adj1,
                        names_to = "metric",
                        values_to = "value") |> 
    tidyr::separate(metric, sep = "_",
                    into = c("metric", "adjustment")) |> 
    tidyr::pivot_wider(id_cols = c(project, index, adjustment),
                       names_from = metric,
                       values_from = value) |> 
    dplyr::mutate(adjustment = case_when(adjustment == "unadj" ~ "<b>Unadjusted</b><br><br><br>",
                                         adjustment == "adj1" ~ "<b>Adjusted</b><br>(sex, age, and<br>immune cell proportion,<br>smoking history)"),
                  adjustment = factor(adjustment, levels = c("<b>Unadjusted</b><br><br><br>", "<b>Adjusted</b><br>(sex, age, and<br>immune cell proportion,<br>smoking history)")),
                  sig = ifelse(p < 0.05, "significant association", NA),
                  index = case_when(grepl("WID_SMK450", index) ~ gsub("_", " ", index),
                                    index == "smoking_mrs" ~ "smoking methylation risk score",
                                    index == "cg05575921" ~ "cg05575921"),
                  index = ifelse(grepl("WID SMK450", index), gsub("WID SMK450 ", "", index), index),
                  index = factor(index, levels = rev(c("cg05575921", "smoking methylation risk score",
                                                       "epithelial hypoM", "immune hypoM", "proximal epithelial hyperM", "distal epithelial hyperM")))) |> 
    ggplot(aes(x = index,
               y = r)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey50") +
    geom_point(shape = 18,
               aes(colour = sig,
                   size = sig)) +
    facet_grid(project~adjustment) +
    scale_colour_manual(values = cols[3],na.value = "grey50",
                        name = '',
                        breaks = c("significant association")) +
    scale_size_manual(values = 4, na.value = 3.5,
                      name = "",
                      breaks = c("significant association")) +
    coord_flip() +
    theme_bw() +
    ylim(c(-1, 1)) +
    theme(legend.position = "top",
          axis.title.x = element_markdown(),
          strip.text = element_markdown(hjust = 0)) +
    labs(x = "", y = "<b>R</b>", title = "Correlation with stage")
  
  
  out2 <- list(os = plot_os,
               pf = plot_pf,
               met = plot_met,
               distant = plot_distmet,
               local = plot_localrec,
               stage = plot_stage)
  
  return(list(data = out1, plots = out2))
  
  
}