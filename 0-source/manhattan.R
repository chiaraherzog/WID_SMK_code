manhattan <- function(db,
                      label_cutoff = 13,
                      sig_level = 5e-8,
                      sample = "buccal"){
  
  require(IlluminaHumanMethylationEPICanno.ilm10b2.hg19) 
  require(ggtext)
  require(ggrepel)
  cols <- MetBrewer::met.brewer("Hokusai1", n = 9)
  
  load(here("1-analysis-pipeline/2-output/WID_SMK_cpgs.Rdata"))
  
  anno <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19, lociNames = db$cg)
  db$chr <- anno$chr
  db$pos <- anno$pos
  db <- db |> # Rename X to 23 for interim
    dplyr::mutate(chr = gsub("chr", "", chr),
                  chr = ifelse(chr=="X", 23, chr),
                  chr = as.numeric(chr))
  
  # Add cumulative numbers
  cumulative <- db |> 
    group_by(chr) |> 
    summarise(max_bp = max(pos)) |> 
    mutate(bp_add = lag(cumsum(as.numeric(max_bp)), default = 0)) |> 
    select(chr, bp_add) |> 
    ungroup()
  
  chrs <- sort(unique(db$chr))
  
  tmp <- db |> 
    left_join(cumulative, by = "chr") |> 
    mutate(bp_cum = as.numeric(pos) + bp_add,
           label = ifelse(-log10(type)>label_cutoff, cg, ""),
           colour = case_when(cg %in% WID_SMK_cpgs[WID_SMK_cpgs$set=="proximal_epithelial_hyperM",]$cg ~ "proximal epithelial hyperM",
                              cg %in% WID_SMK_cpgs[WID_SMK_cpgs$set=="distal_epithelial_hyperM",]$cg ~ "distal epithelial hyperM",
                              cg %in% WID_SMK_cpgs[WID_SMK_cpgs$set=="immune_hypoM",]$cg ~ "immune hypoM",
                              cg %in% WID_SMK_cpgs[WID_SMK_cpgs$set=="epithelial_hypoM",]$cg ~ "epithelial hypoM",
                              chr %in% chrs[seq(1, length(chrs), 2)] & !cg %in% WID_SMK_cpgs$cg ~ "chrs1",
                              chr %in% chrs[seq(0, length(chrs), 2)] & !cg %in% WID_SMK_cpgs$cg ~ "chrs2"),
           colour = factor(colour, levels = c("chrs1",
                                              "chrs2",
                                              "epithelial hypoM",
                                              "immune hypoM",
                                              "distal epithelial hyperM",
                                              "proximal epithelial hyperM")))
  
  axis_set <- tmp |>  
    group_by(chr) |> 
    summarize(center = mean(bp_cum)) |> 
    ungroup() |> 
    dplyr::mutate(chr = ifelse(chr==23, "X", chr))
  
  ylim <- tmp |> 
    filter(type == min(type)) |> 
    mutate(ylim = abs(floor(log10(type))) + 2) |> 
    pull(ylim)
  
  manhplot <- tmp |> 
    ggplot(
      aes(x = bp_cum, y = -log10(type), 
          size = -log10(type),
          colour = colour
      ),) +
    geom_hline(yintercept = -log10(sig_level), color = "grey40", linetype = "dashed") +
    geom_point() +
    scale_x_continuous(label = axis_set$chr, breaks = axis_set$center) +
    scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
    scale_color_manual(values = c("black", "grey60",
                                         cols[c(9, 5, 3, 1)]),
                       name = "",
                       labels = c("",
                                  "",
                                  "epithelial hypoM",
                                  "immune hypoM",
                                  "distal epithelial hyperM",
                                  "proximal epithelial hyperM")) +
    scale_size_continuous(range = c(0.25,1.5)) +
    labs(x = NULL, 
         y = "-log<sub>10</sub>(p)",
         title = paste0("<span style = 'font-size:14pt; font-family:Helvetica;'><b>Smoking EWAS </b> in ", sample, " samples</span>")) + 
    theme_minimal() +
    theme( 
      legend.position = "top",
      panel.grid.major.x = element_blank(),
      plot.title = element_markdown(),
      panel.grid.minor.x = element_blank(),
      axis.title.y = element_markdown(),
      axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
    ) +
    ggrepel::geom_text_repel(aes(label = label),
                             min.segment.length = 0,
                             nudge_x = 100,
                             nudge_y = 2,
                             force = 2,
                             size = 2.6) +
    guides(size = guide_legend(show_guide = FALSE),
           colour = guide_legend(nrow = 1))
  
  return(manhplot)
}
