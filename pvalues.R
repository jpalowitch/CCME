pvalFun <- function (B, nSet = 1:n) {
  
  nodesToB <- intersect(which(colSums(adjMat[B, ]) > 0), nSet)
  m <- length(B)
  
  if (!binary) {
    
    sB_out <- strengths[B]
    sB_in  <- strengths[nodesToB]
    dB_out <- degrees[B]
    dB_in  <- degrees[nodesToB]
    means <- sB_in * sum(sB_out) / sT    
    
    vars <- varFun(dB_in, sB_in, dB_out, sB_out, theta, dT, sT)
    stats <- rowSums(adjMat[nodesToB, B, drop = FALSE])
    pvalsToB <- pnorm(stats, means, sqrt(vars), lower.tail = FALSE)
    
    
  } else {
    
    stats <- rowSums(adjMat[nSet, B])
    pvalsToB <- pbinom(stats - 1, 
                       degrees, 
                       sum(degrees[B]) / dT, 
                       lower.tail = FALSE)
    
  }
  
  pvals <- rep(1, length(nSet))
  pvals[match(nodesToB, nSet)] <- pvalsToB
  return(pvals)
  
}
