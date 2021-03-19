threshold_mat <- function(m, t) {
  mt <- m
  mt[which(abs(m) < t)] <- 0.0
  return(mt)
}

DesignMatrixAlt <- function(X, rango = .95) {
  MAT <- matrix(0, dim(X)[1], dim(X)[2])
  MAT2 <- c()
  expa <- c()
  for (p in 1:ncol(X)) {
    x <- X[, p]
    MAT[, p] <- x#spikeSlabGAM::lin(x)
    z <- spikeSlabGAM::sm(x,
      K = 20, spline.degree = 3, diff.ord = 2,
      rankZ = rango, centerBase = T, centerx = x,
      decomposition = "ortho"
    )
    MAT2 <- cbind(MAT2, z)
    expa <- c(expa, rep((p - 1), ncol(z)))
  }
  dj <- as.vector(table(expa))
  return(list(Xlin = MAT, Xnlin = MAT2, dj = dj))
}
