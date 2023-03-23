# Enrichment function (PCGT)

enrich <- function(A, B, universe){
  df <- data.frame(matrix(nrow = length(A),
                          ncol = 5))
  colnames(df) <- c("set", "or", "ci_lo", "ci_hi", "p")
  
  for(i in 1:length(A)){
    df$set[i] <- names(A)[i]
    
    # 			      in set		not in set
    # in pcgt       a       b
    # not pcgt      c       d
    
    a <- length(intersect(A[[i]], B))
    b <- length(B[!B %in% A[[i]]])
    c <- length(A[[i]][!A[[i]] %in% B])
    d <- length(universe[!universe %in% c(A[[i]], B)])
    
    x <- matrix(c(a, c, b, d), ncol = 2)
    t <- fisher.test(x)
    
    df$or[i] <- t$estimate
    df$ci_lo[i] <- t$conf.int[1]
    df$ci_hi[i] <- t$conf.int[2]
    df$p[i] <- t$p.value
  }
  return(df)
}