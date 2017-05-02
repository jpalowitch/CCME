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