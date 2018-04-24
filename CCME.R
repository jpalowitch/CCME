library(Matrix)
library(plyr)
library(Rcpp)
library(RcppArmadillo)
library(microbenchmark)
library(igraph)

#sourceCpp("new_funs.cpp")

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
                  zoverlap = FALSE,
                  bipartite = FALSE) {
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
  
  if (bipartite && !binary) {
    message("Bipartite version not available for weighted networks.")
    stop()
  }

  
  # Operation Functions --------------------------------------------------------
  
  symdiff <- function (s1, s2) {
    return(union(setdiff(s1, s2), setdiff(s2, s1)))
  }
  
  jaccard <- function (s1, s2) {
    return(length(symdiff(s1, s2)) / length(union(s1, s2)))
  }
  
  bh_reject <- function (pvals, alpha) {
    pvals_adj <- length(pvals) * pvals / rank(pvals)
    if (sum(pvals_adj <= alpha) > 0) {
      thres <- max(pvals[pvals_adj <= alpha])
      return(which(pvals <= thres))
    } else {
      return(integer(0))
    }
  }
  
  get_set_z2 <- function (B) {
    stat <- sum(edge_list$weight[edge_list$node1 %in% B & edge_list$node2 %in% B])
    B_s <- strengths[B]
    B_d <- degrees[B]
    return(get_set_z(stat, B_s, B_d, theta, dT, sT))
  }
  
  
  filter_overlap <- function (comms, tau, Zoverlap = zoverlap) {
    
    K <- length(comms)
    if (Zoverlap) {
      scores <- unlist(lapply(comms, get_set_z2))
    } else {
      scores <- unlist(lapply(comms, length))
    }
    
    
    jaccard_mat0 <- matrix(0, K, K)
    for (i in 1:K) {
      for (j in 1:K) {
        jaccard_mat0[i, j] <- length(intersect(comms[[i]], comms[[j]])) / 
          length(comms[[i]])
      }
    }
    
    jaccard_mat <- jaccard_mat0
    diag(jaccard_mat) <- 0
    max_jacc <- max(jaccard_mat)
    deleted_comms <- integer(0)
    
    while (max_jacc > tau) {
      
      inds <- which(jaccard_mat == max_jacc, arr.ind = TRUE)[1, ]
      
      # keep comm with larger score
      delete_comm <- inds[which.min(c(scores[inds[1]], scores[inds[2]]))]
      jaccard_mat[delete_comm, ] <- 0
      jaccard_mat[, delete_comm] <- 0
      deleted_comms <- c(deleted_comms, delete_comm)
      max_jacc <- max(jaccard_mat)
      
    }
    
    kept_comms <- setdiff(1:K, deleted_comms)
    
    return(list("final_comms" = comms[kept_comms],
                "kept_comms" = kept_comms))
    
  }
  
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
  
  halfUpdate <- function (S1, side) {
    pvals <- pvalFun(S, side)
    S2 <- side[bh_reject(pvals)]
    if (length(S1) == 0 || length(S2) == 0) {
      return(integer(0))
    } else {
      return(c(S1, S2))
    }
  }  
  
  initializeBinary <- function (n.samples=NULL) {
  
    if (is.null(n.samples)) {
      n.samples <- n
    }
    
    if (bipartite) {
      inSets <- rep(list(NULL), length(side1) + length(side2))
      for (u in side1) {
        B <- which(adjMat[u, ] > 0)
        inSets[[side1 == u]] <- halfUpdate(B, side1)
      }
      for (v in side2) {
        A <- which(adjMat[v, ] > 0)
        inSets[[side2 == v]] <- halfUpdate(A, side2)
      }
      inNodes <- c(side1, side2)
    } else {
      for (u in 1:n) {
        inSets[[u]] <- which(adjMat[u, ] > 0)
      }
      inNodes <- 1:n
    }    
    return(list("inSets" = inSets, "inNodes" = inNodes))
  }

  initializeWeighted <- function (n.samples = NULL) {
    
    if (is.null(n.samples)) n.samples <- n
    
    if (!fastInitial) {
      
      # Making expanded edge_list
      lower_edge_list <- edge_list[ , c(2, 1, 3)]
      names(lower_edge_list) <- c("node1", "node2", "weight")
      edge_list_ex <- rbind(edge_list, lower_edge_list)
      edge_list_ex <- edge_list_ex[order(edge_list_ex$node1), ]
      rm(lower_edge_list)
      
      big_deg_nodes <- seq_along(degrees)[degrees > 1]
      functional_n.samples <- min(n.samples, length(big_deg_nodes))
      nodeSet <- sort(sample(big_deg_nodes, 
                             functional_n.samples, 
                             replace = FALSE))
      sample_from <- which(edge_list_ex$node1 %in% nodeSet)
      
      
      # Creating normalized prob weight
      str_ex <- strengths[edge_list_ex$node1] * strengths[edge_list_ex$node2] /
        sT
      deg_ex <- degrees[edge_list_ex$node1] * degrees[edge_list_ex$node2] / 
        dT
      deg_ex[deg_ex > 1] <- 1
      mean_ex <- str_ex / deg_ex
      var_ex <- mean_ex^2 * theta
      norm_wt <- (edge_list_ex$weight - mean_ex) / sqrt(var_ex)
      norm_wt <- pmax(norm_wt, 0)
      
      # Making sampling data frame
      dataToSample <- data.frame("fromNodes" = edge_list_ex[ , 1],
                                 "toNodes" = edge_list_ex[ , 2], 
                                 "probs"   = norm_wt)
      rownames(dataToSample) <- NULL
      trackerVec <- match(dataToSample$fromNodes, 
                          sort(unique(dataToSample$fromNodes)))
      dataToSample <- data.frame(dataToSample,
                                 "tracker" = trackerVec)
      rm(trackerVec)
      b_indxs <- which(diff(c(0, dataToSample$tracker)) > 0)
      nSets <- length(b_indxs)
      from_degs <- degrees[dataToSample$fromNodes[b_indxs]]
      probs <- dataToSample$probs
      nEdges <- length(probs)
      chosen <- rep(0, nEdges)
      e_indxs <- c(b_indxs[2:nSets] - 1, nEdges)
      
      # Sampling
      inc_pos <- 0
      increment <- nSets / 10
      if (generalOutput) {
        cat("Sampling candidate initial sets:\n")
        cat("--0%")
      }
      
      for (i in 1:nSets) {
        if (floor(i / increment) > inc_pos) {
          inc_pos <- inc_pos + 1
          if (generalOutput)
            cat(paste0("--", inc_pos * 10, "%"))
        }
        indxs_i <- b_indxs[i]:e_indxs[i]
        if (sum(probs[indxs_i] > 0)) {
          if (length(indxs_i) == 1) {
            chosen_i <- indxs_i
          } else {
            chosen_i <- sample(indxs_i,
                               from_degs[i],
                               TRUE, prob = probs[indxs_i])
            chosen_i <- unique(chosen_i)
          }
          chosen[chosen_i] <- 1
        }
      }
      
      if (generalOutput)
        cat("\n")
      
      # Making reduced data frame
      samples <- dataToSample[chosen == 1, ]
      inNodes <- unique(samples$fromNodes)
      
      # Removing singleton sets
      singleton_initializers <- inNodes[table(samples$fromNodes) == 1]
      samples <- samples[!samples$fromNodes %in% singleton_initializers, ]
      inNodes <- setdiff(inNodes, singleton_initializers)
      
      if (generalOutput)
        cat("Calculating set test statistics:\n")
      # Get stats
      nSets2 <- length(inNodes)
      b_indxs2 <- which(diff(c(0, samples$tracker)) > 0)
      e_indxs2 <- c(b_indxs2[2:nSets2] - 1, nrow(samples))
      set_stats <- get_set_stats(b_indxs2, e_indxs2, samples$toNodes,
                                 edge_list$node1, edge_list$node2,
                                 edge_list$weight,
                                 generalOutput = generalOutput)
      if (generalOutput)
        cat("\n")
      
      # Extracting sets into list
      inSets <- split(samples$toNodes, 
                      samples$fromNodes)
      
      # Make master list of the set + the statistic
      inSets_plus_stats <- rep(list(rep(list(NULL), 3)), length(inSets))
      for (i in seq_along(set_stats)) {
        inSets_plus_stats[[i]][[1]] <- inSets[[i]]
        inSets_plus_stats[[i]][[2]] <- set_stats[i]
        inSets_plus_stats[[i]][[3]] <- i
      }
      
      # Calculating p-values
      report_vec <- seq(1, nSets2, length.out = 11)
      if (generalOutput)
        cat("Calculating set p-values:\n")
      # p-value function for lapply (utilizes Rcpp)
      get_set_pvalue_R <- function (inSet_stat) {
        report_score <- inSet_stat[[3]] - report_vec
        last_report <- suppressWarnings(
          min(report_score[report_score >= 0])
        )
        if (generalOutput && last_report < 1) {
          cat(paste0("--", (sum(report_score >= 0) - 1) * 10, "%"))
        }
        return (get_set_pvalue(inSet_stat[[2]],
                               strengths[inSet_stat[[1]]],
                               degrees[inSet_stat[[1]]],
                               theta, dT, sT))
      }
      
      pvals <- unlist(lapply(inSets_plus_stats, get_set_pvalue_R))
      if (generalOutput)
        cat("\n")
      
      # Filtering
      rejections <- bh_reject(pvals, alpha)  
      inNodes <- inNodes[rejections]
      inSets <- inSets[rejections]
      
    } else {
      
      G <- graph.edgelist(as.matrix(edge_list[ , 1:2]), directed = FALSE)
      E(G)$weight <- edge_list[ , 3]
      
      ig_results1 <- cluster_walktrap(G)
      ig_results2 <- cluster_fast_greedy(G)
      
      inSets1 <- lapply(unique(ig_results1$membership), 
                        function (i) which(ig_results1$membership == i))
      inSets2 <- lapply(unique(ig_results2$membership), 
                        function (i) which(ig_results2$membership == i))
      inSets <- c(inSets1, inSets2)
      
      inSets <- inSets[unlist(lapply(inSets, length)) > 1]
      inNodes <- unlist(lapply(inSets, function (c) c[1]))
      
    }
    
    return(list("inSets" = inSets, "inNodes" = inNodes))
  }
  
  pvalFun <- function (B, nSet = 1:n) {
    
    nodesToB <- intersect(which(colSums(adjMat[B, ]) > 0), nSet)
    m <- length(B)
    dB_out <- degrees[B]
    dB_in  <- degrees[nodesToB]
    
    if (!binary) {
      
      sB_out <- strengths[B]
      sB_in  <- strengths[nodesToB]
      means <- sB_in * sum(sB_out) / sT    
      
      vars <- varFun(dB_in, sB_in, dB_out, sB_out, theta, dT, sT)
      stats <- rowSums(adjMat[nodesToB, B, drop = FALSE])
      pvalsToB <- pnorm(stats, means, sqrt(vars), lower.tail = FALSE)
      
      
    } else {
      
      stats <- rowSums(adjMat[nodesToB, B])
      pvalsToB <- pbinom(stats - 1, dB_in, sum(dB_out) / dT, lower.tail = FALSE)
      
    }
    
    pvals <- rep(1, length(nSet))
    pvals[match(nodesToB, nSet)] <- pvalsToB
    return(pvals)
    
  }
  
  bipartitePvalFun <-  function (B) {
    pvals <- rep(1, length(B))
    B1 <- B[B %in% side1]
    B2 <- B[B %in% side2]
    if (length(B1) == 0 || length(B2) == 0) {
        return(pvals)
    }
    B1pvals <- pvalFun(B1, side2)
    B2pvals <- pvalFun(B2, side1)    
    pvals[B %in% side1] <- B1pvals
    pvals[B %in% side2] <- B2pvals
    return(pvals)
  }
  
  extract <- function (i) {
    
    if (!fastInitial) {
      doTheNode <- !inNodes[i] %in% unionC
      
      if (throwInitial) {
        doTheNode <- doTheNode && !inNodes[i] %in% unionInitial
        if (inNodes[i] %in% unionInitial) {
          cat("###\nFound a node in unionInitial\n")
          cat("!inNodes[i] %in% unionInitial", !inNodes[i] %in% unionInitial, "\n")
          cat("!inNodes[i] %in% unionC", !inNodes[i] %in% unionC, "\n")
          cat("doTheNode now", doTheNode, "\n###\n")
        }
      }
      
    } else {
      doTheNode <- TRUE
    }
    
    
    if (doTheNode) {
      
      if (!fastInitial) {
        B0 <- c(inNodes[i], inSets[[i]])
      } else {
        B0 <- inSets[[i]]
      }
      unionInitial <<- union(unionInitial, B0)
      
      if (loopOutput) {
        cat("#-------------------------------------\n")
        cat(paste0("Initial set ", i,
                   " of ", length(inSets), " is size ", length(B0), ".\n"))
        cat("length(unionInitial) =", length(unionInitial), "\n")
        cat("length(unionC) =", length(unionC), "\n")
      }
      
      # Updates ----------------------------------------------------------------
      
      B_old <- B0
      B_new <- 1:n
      chain <- list(B_old)
      consec_jaccards <- NULL
      found_cycle <- found_break <- NULL
      mean_jaccards <- NULL ## add all these to update_info
      itCount <- 0
      cycledSets <- NULL
      
      if (loopOutput)
        cat("Beginning Updates\n")
      
      while (length(B_new) > 1) {
        
        itCount <- itCount + 1
        if (!bipartite) {
          pvals <- pvalFun(B_old, nSet = 1:n)
        } else {
          pvals <- bipartitePvalFun(B_old)
        }
        B_new <- bh_reject(pvals, alpha)
        
        consec_jaccard <- jaccard(B_new, B_old)
        jaccards <- unlist(lapply(chain, function (B) jaccard(B, B_new)))
        found_cycle <- c(found_cycle, FALSE)
        found_break <- c(found_break, FALSE)
        consec_jaccards <- c(consec_jaccards, consec_jaccard)
        mean_jaccards <- c(mean_jaccards, mean(jaccards))
        
        if (updateOutput) {
          cat(paste0("Update ", itCount, 
                     " is size ", length(B_new), ", ",
                     "jaccard to last is ", round(consec_jaccard, 3), ", ",
                     "mean jaccard along chain is ", round(mean(jaccards), 3), "\n", sep=""))
        }
        
        # Checking for cycles (4.4.1 in paper)
        if (jaccard(B_new, B_old) > 0) { # Otherwise loop will end naturally
          
          jaccards <- unlist(lapply(chain, function (B) jaccard(B, B_new)))
          
          if (sum(jaccards == 0) > 0) { # Cycle has been found
            
            found_cycle[itCount] <- TRUE
            
            if (updateOutput)
              cat("---- Cycle found")
            
            Start <- max(which(jaccards == 0))
            cycle_chain <- chain[Start:length(chain)]
            
            # Checking for cycle break (4.4.1a)
            seq_pair_jaccards <- rep(0, length(cycle_chain) - 1)
            for (j in seq_along(seq_pair_jaccards)) {
              seq_pair_jaccards[j] <- jaccard(chain[[Start + j - 1]], chain[[Start + j]])
            }
            if (sum(seq_pair_jaccards == 1) > 0) {# then break needed
              if (updateOutput)
                cat(" ---- Break found\n")
              found_break[itCount] <- TRUE
              B_new <- NULL
              break
            }
            
            # Create conglomerate set (and check, 4.4.1b)
            B_J <- unique(unlist(cycle_chain))
            B_new <- B_J
            B_J_check <- unlist(lapply(chain, function (B) jaccard(B_J, B)))
            if (sum(B_J_check == 0) > 0) {
              if (updateOutput)
                cat(" ---- Old cycle\n")
              break
            } else {
              if (updateOutput)
                cat(" ---- New cycle\n")
            }
            
          } # From checking jaccards to cycle_chain
          
        } else { # From checking B_new to B_old; if B_new = B_old: 
          break
        }
        
        B_old <- B_new
        chain <- c(chain, list(B_new))
        
      } # From Updates
      
      unionC <<- union(B_new, unionC)
      if (loopOutput)
        cat(paste0(length(unionC), " vertices in communities.\n"))
      
      update_info <- list("mean_jaccards" = mean_jaccards,
                          "consec_jaccards" = consec_jaccards,
                          "found_cycle" = found_cycle,
                          "found_break" = found_break,
                          "didNode" = doTheNode)
      
    } else {
      
      B_new <- integer(0)
      update_info <- list("didNode" = doTheNode)
      
    }
    
    
    return(list("comm" = B_new, "update_info" = update_info))
    
  }
  
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
    
    strengths <- colSums(adjMat)
    
    sT <- sum(strengths)
    dT <- sum(degrees)
    
    theta <- varpar()
    
  }
  
  if (bipartite) {
    side1 <- unique(sort(edge_list[ , 1]))
    side2 <- unique(sort(edge_list[ , 2]))
    if (length(intersect(side1, side2)) > 0) {
      message("Error: bipartite sides are not disjoint.")
      stop()
    }
  }
  
  # Sampling initial sets ------------------------------------------------------
  
  if (binary) {
    initRes <- initializeBinary(n.samples)
  } else {
    initRes <- initializeWeighted(n.samples)
  }
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

