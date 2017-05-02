library(Matrix)
library(plyr)
library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
library(igraph)

sourceCpp("new_funs.cpp")

CCME <- function (edge_list, 
                  alpha = 0.05, 
                  binary = FALSE,
                  n.samples = NULL,
                  OL_thres = 0.90,
                  updateOutput = FALSE, 
                  loopOutput = TRUE, 
                  generalOutput = TRUE,
                  throwInitial = FALSE,
                  fastInitial = FALSE,
                  zoverlap = FALSE) {
  if (FALSE) {
    alpha = 0.05
    binary = FALSE
    n.samples = NULL
    OL_thres = 0.90
    updateOutput = FALSE
    loopOutput = TRUE 
    generalOutput = TRUE
    throwInitial = FALSE
    fastInitial = FALSE
    zoverlap = FALSE
  }

  
  # Operation Functions --------------------------------------------------------
  
  source("auxiliary.R", local = TRUE)
  source("overlap.R", local = TRUE)
  source("varpar.R", local = TRUE)
  source("initialize.R", local = TRUE)
  source("pvalues.R", local = TRUE)
  source("extract.R", local = TRUE)

  # Setup Calculations ---------------------------------------------------------
  if (generalOutput)
    cat("Performing basic calculations...\n")
  
  n <- max(c(edge_list[ , 1], edge_list[ , 2]))

  if (binary) {
    adjMat <- sparseMatrix(i = edge_list[ , 1],
                           j = edge_list[ , 2],
                           dims = c(n, n),
                           symmetric = TRUE)
  } else {
    adjMat <- sparseMatrix(i = edge_list[ , 1],
                           j = edge_list[ , 2],
                           dims = c(n, n),
                           x = edge_list[ , 3],
                           symmetric = TRUE)
  }
  
  degrees <- colSums(adjMat > 0)
  
  if (!binary) {
    
    degrees <- colSums(adjMat > 0)
    strengths <- colSums(adjMat)
    
    sT <- sum(strengths)
    dT <- sum(degrees)
    
    theta <- varpar()
    
  }
  
  # Sampling initial sets ------------------------------------------------------
  
  initRes <- initialize(n.samples)
  inSets <- initRes$inSets
  inNodes <- initRes$inNodes
  rm(initRes)
  gc()
  
  # Returning if none significant ----------------------------------------------
  
  if (length(inSets) == 0) {
    returnList = list("communities" = NULL,
                      "background" = 1:n,
                      "initial.sets" = NULL, 
                      "final.sets" = NULL,
                      "strengths" = strengths,
                      "degrees" = degrees,
                      "theta" = theta,
                      "set_info" = NULL)
    return(returnList)
  }

  # Extractions ----------------------------------------------------------------
  
  if(generalOutput)
    cat("Beginning extractions...\n")
  
  unionC <- unionInitial <- integer(0)
  extractRes <- lapply(seq_along(inNodes), extract)
  comms <- lapply(extractRes, function (L) L$comm)
  update_info <- lapply(extractRes, function (L) L$update_info)
  
  #  Clean-up and return -------------------------------------------------------
  
  # Removing blanks and trivial sets
  nonNullIndxs <- which(unlist(lapply(comms, length)) > 1)
  comms0 <- comms[nonNullIndxs]
  filtered_comms <- list(NULL)
  if (length(comms0) > 0) {
    OLfilt <- filter_overlap(comms0, tau = OL_thres, Zoverlap = TRUE)
    filtered_comms <- OLfilt$final_comms
  } else {
    OLfilt <- NULL
  }
  
  background <- setdiff(1:n, unique(unlist(filtered_comms)))
  
  returnList <- list("communities" = filtered_comms,
                     "background" = background,
                     "initial.sets" = inSets, 
                     "final.sets" = comms,
                     "communities_before_OLfilt" = comms0,
                     "OLfilt" = OLfilt,
                     "update_info" = update_info,
                     "nonNullIndxs" = nonNullIndxs)
  
  return(returnList)
  
}

