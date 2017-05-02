initialize <- function (n.samples = NULL) {
  
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
    cat("Sampling candidate initial sets:\n")
    cat("--0%")
    
    for (i in 1:nSets) {
      if (floor(i / increment) > inc_pos) {
        inc_pos <- inc_pos + 1
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
    
    cat("\n")
    
    # Making reduced data frame
    samples <- dataToSample[chosen == 1, ]
    inNodes <- unique(samples$fromNodes)
    
    # Removing singleton sets
    singleton_initializers <- inNodes[table(samples$fromNodes) == 1]
    samples <- samples[!samples$fromNodes %in% singleton_initializers, ]
    inNodes <- setdiff(inNodes, singleton_initializers)
    
    cat("Calculating set test statistics:\n")
    # Get stats
    nSets2 <- length(inNodes)
    b_indxs2 <- which(diff(c(0, samples$tracker)) > 0)
    e_indxs2 <- c(b_indxs2[2:nSets2] - 1, nrow(samples))
    set_stats <- get_set_stats(b_indxs2, e_indxs2, samples$toNodes,
                               edge_list$node1, edge_list$node2,
                               edge_list$weight)
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
    cat("Calculating set p-values:\n")
    # p-value function for lapply (utilizes Rcpp)
    get_set_pvalue_R <- function (inSet_stat) {
      report_score <- inSet_stat[[3]] - report_vec
      last_report <- suppressWarnings(
        min(report_score[report_score >= 0])
      )
      if (last_report < 1) {
        cat(paste0("--", (sum(report_score >= 0) - 1) * 10, "%"))
      }
      return (get_set_pvalue(inSet_stat[[2]],
                             strengths[inSet_stat[[1]]],
                             degrees[inSet_stat[[1]]],
                             theta, dT, sT))
    }
    
    pvals <- unlist(lapply(inSets_plus_stats, get_set_pvalue_R))
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
}