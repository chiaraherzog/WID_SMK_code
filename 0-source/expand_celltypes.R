expand_celltypes <- function(dat){
  
  tmp <- dat |> 
    tidyr::pivot_longer(c(hepidish_Epi,
                   hepidish_Fib,
                   hepidish_B,
                   hepidish_NK,
                   hepidish_CD4T,
                   hepidish_CD8T,
                   hepidish_Mono,
                   hepidish_Neutro,
                   hepidish_Eosino),
                 names_to = "celltype",
                 values_to = "cell_prop") |> 
    dplyr::mutate(celltype = case_when(celltype == "hepidish_Epi" ~ "Epithelial cells",
                                celltype == "hepidish_Fib" ~ "Fibroblasts",
                                celltype == "hepidish_B" ~ "B cells",
                                celltype == "hepidish_NK" ~ "NK cells",
                                celltype == "hepidish_CD4T" ~ "CD4 T cells",
                                celltype == "hepidish_CD8T" ~ "CD8 T cells",
                                celltype == "hepidish_Mono" ~ "Monocytes",
                                celltype == "hepidish_Neutro" ~ "Neutrophils",
                                celltype == "hepidish_Eosino" ~ "Eosinophils")) %>%
    dplyr::mutate(celltype = factor(celltype, levels = c("Epithelial cells", "Fibroblasts", "Monocytes", "Neutrophils", "Eosinophils", "B cells", "NK cells", "CD4 T cells", "CD8 T cells")))
  
  return(tmp)
}
