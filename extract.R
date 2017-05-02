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
    
    cat("Beginning Updates\n")
    
    while (length(B_new) > 1) {
      
      itCount <- itCount + 1
      pvals <- pvalFun(B_old, nSet = 1:n)
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