MakeBinaryBipartiteNetwork <- function () {
    n1 <- 50
    n2 <- 50
    P <- matrix(0.2, 5, 5)
    diag(P) <- 0.8
    comm.sizes1 <- comm.sizes2 <- rep(10, 5)
    nodes <- sample(n1 + n2)
    side1 <- head(nodes, n1)
    side2 <- tail(nodes, n2)
    comms1 <- split(side1, ceiling(seq_along(side1) / 10))
    comms2 <- split(side2, ceiling(seq_along(side2) / 10))        
    edge.lists <- list()
    comms <- list()
    for (i in 1:5) {        
        for (j in 1:5) {
            commi <- comms1[[i]]
            commj <- comms2[[j]]
            nedges <- rbinom(1, comm.sizes1[i] * comm.sizes2[j], P[i, j])
            side1.draws <- sample(rep(commi, comm.sizes2[j]), nedges,
                                  replace=FALSE)
            side2.draws <- sample(rep(commj, comm.sizes1[i]), nedges,
                                  replace=FALSE)
            edge.lists <- c(edge.lists, list(cbind(side1.draws, side2.draws)))
        }
    }
    edgelist <- do.call(rbind, edge.lists)
    return(list(edgelist=edgelist, comms=comms))
}   
