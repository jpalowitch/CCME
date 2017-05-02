varpar <- function () {
  
  
  if(generalOutput)
    cat("--Calculating variance..\n")
  
  # Making sparse version of adjMat for varCalc
  if (class(adjMat) != "dgCMatrix") {
    sparseSummary <- summary(Matrix(adjMat, sparse = TRUE))
  } else {
    sparseSummary <- summary(adjMat)
  }
  
  # Variance calc
  suv <- strengths[sparseSummary[,1]] * strengths[sparseSummary[,2]] / sT
  duv <- degrees[sparseSummary[,1]] * degrees[sparseSummary[,2]] / dT
  duv[duv > 1] <- 1
  SSo <- sum((sparseSummary[,3] - suv/duv)^2)
  SSe <- sum((suv/duv)^2)
  theta <- SSo/SSe
  
  if(generalOutput)
    cat("--Variance parameter is",theta,"\n")
  
  return(theta)
  
}