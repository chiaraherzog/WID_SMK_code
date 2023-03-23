smoking_mrs <- function(beta){
  weights <- readxl::read_xlsx("0-source/brenner_smoking_mrs.xlsx") |> 
    janitor::clean_names() |> 
    dplyr::select(cp_gs, weight)
  
  intersect <- intersect(rownames(beta), weights$cp_gs)
  B <- beta[intersect,]
  weights <- weights[match(intersect, weights$cp_gs),]$weight
  mu <- rowMeans(B)
  sigma <- matrixStats::rowSds(B)
  score <- (1/length(weights))* colSums(weights*((B-mu)/sigma))
  return(score)
}