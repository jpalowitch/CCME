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