#' Simulates data from a Dirichlet-Multinomial regression model.
#'
#' \code{gendata} returns the data simulated according to input parameters
#'
#' @param J number of taxa
#' @param n number of samples
#' @param P number of continuous covariates
#' @param Q number of discrete covariates
#' @param beta_min minimum absolute value of the regression parameters
#' @param beta_max maximum absolute value of the regression parameters
#' @param n_rel_taxa number of relevant taxa
#' @param n_rel_x number of relevant continuous covariates
#' @param n_rel_z number of relevant discrete covariates
#' @param n_rel_inter vector of length four
#' \enumerate{
#'   \item number of continuous-discrete interactions
#'   \item number of discrete-discrete interactions
#'   \item number of linear continuous-continuous interactions
#'   \item number of nonlinear continuous-continuous interactions
#' }
#' @param corrx the correlation between continuous covariates
#' @param corrz the correlation between discrete covariates
#' @param theta0 the overdispersion parameter
#' @param WH the weak heredity flag
#' @param nli whether to include nonlinear interactions between continuous covariates
#' @param randint whether should be included a random intercept
#' @param repmeas number of repeated measurements for each sample

#' @return Y: (count matrix) \eqn{n \times J}
#' @return X: (design matrix) \eqn{n \times P} (continuous variables)
#' @return Z: (design matrix) \eqn{n \times Q} (discrete variables)
# #' @return Xl: linear  = Xilp, Xnl = Xinlp, dj = dj,
#' @return intercept: simulated intercept \eqn{J-}dimensional vector
#' @return thetax: simulated coefficient matrix \eqn{P \times J}
#' @return thetaz: simulated coefficient matrix \eqn{Q \times J}
#' @return inter_xz: simulated coefficient matrix (continuous-discrete
#' interactions) \eqn{P \times Q \times J}
#' @return inter_xx: simulated coefficient matrix (continuous-continuous linear interactions)
#' \eqn{P \times P \times J}
#' @return inter_xx_nl: simulated coefficient matrix (continuous-continuous nonlinear interactions)
#' \eqn{n \times P \times J}
#' @return nlind: indicator for nonlinear interactions
#' @return inter_zz: simulated coefficient matrix (discrete-discrete
#' interactions) \eqn{Q \times Q \times J}
#' @return RI: simulated random intercept \eqn{n-}dimensional vector
#' @return GI: group indicator for repeated measurements
#' @return ZETA: log linear predictor
#' @return part1, part2: subject-specific coefficients

#' @examples
#' gendata(10, 100, 10, 5, 0.75, 1.5, 2, 2, 2, c(2, 2, 3, 3), 0.4, 0.4, 0.01, FALSE, TRUE)
#' @export
#'

gendata <- function(J = 10, n = 100, P = 10, Q = 5, beta_min = 0.75, beta_max = 1.5,
                    n_rel_taxa = 2, n_rel_x = 2, n_rel_z = 2,
                    n_rel_inter = c(2, 2, 3, 3), corrx = 0.4,
                    corrz = 0.4, theta0 = 0.01, WH = FALSE, nli = TRUE,
                    randint = FALSE, repmeas = 3) {

  # construct continuous covariates
  if (corrx != 0) {
    Sigmax <- matrix(1, P, P)
    Sigmax <- corrx^abs(row(Sigmax) - col(Sigmax))
  } else {
    Sigmax <- diag(1, P, P)
  }
  X <- scale(MASS::mvrnorm(n = n, mu = rep(0, P), Sigma = Sigmax))

  # construct discrete covariates
  if (corrz != 0) {
    Sigmaz <- matrix(1, Q, Q)
    Sigmaz <- corrz^abs(row(Sigmaz) - col(Sigmaz))
  } else {
    Sigmaz <- diag(1, Q, Q)
  }
  Z <- scale(MASS::mvrnorm(n = n, mu = rep(0.0, Q), Sigma = Sigmaz))
  Z[which(Z >= 0.0)] <- 1
  Z[which(Z < 0.0)] <- 0

  # main effects X
  thetaX <- matrix(0, J, P)
  st <- 0
  low_side <- beta_min
  high_side <- beta_max
  if (n_rel_taxa != 1) {
    # warning if the lengths don't match
    coef <- suppressWarnings(seq(low_side, high_side, len = n_rel_taxa) *
      c(1, -1))
  } else {
    coef <- (low_side + high_side) / 2
  }
  coef_g <- rep(1.0, len = n_rel_x)
  for (ii in 1:n_rel_x) {
    # overlap species
    thetaX[(st:(st + n_rel_taxa - 1)) %% J + 1, 3 * ii - 2] <- coef_g[ii] *
      sample(coef)[((ii - 1):(ii + n_rel_taxa - 2)) %% n_rel_taxa + 1]
    st <- st + 1
  }

  thetaX <- t(thetaX)

  # main effects Z
  thetaZ <- matrix(0, J, Q)
  st <- 0
  low_side <- beta_min
  high_side <- beta_max
  if (n_rel_taxa != 1) {
    # warning if the lengths don't match
    coef <- suppressWarnings(seq(low_side * 1.5, high_side * 1.5, len = n_rel_taxa) * c(1, -1))
  } else {
    coef <- (low_side + high_side) / 2
  }
  coef_g <- rep(1.0, len = n_rel_z)
  for (ii in 1:n_rel_z) {
    # overlap species
    thetaZ[(st:(st + n_rel_taxa - 1)) %% J + 1, 3 * ii - 2] <- coef_g[ii] * sample(coef)[((ii - 1):(ii + n_rel_taxa - 2)) %% n_rel_taxa + 1]
    st <- st + 1
  }

  thetaZ <- t(thetaZ)

  inter_xz <- array(0, dim = c(P, Q, J))
  # interactions xz
  for (j in 1:J) {
    count <- 0
    for (p in 1:P) {
      for (q in 1:Q) {
        if ((thetaX[p, j] != 0) & (thetaZ[q, j] != 0)) {
          if (count <= n_rel_inter[1]) {
            inter_xz[p, q, j] <- thetaX[p, j] * thetaZ[q, j]#sample(c(thetaX[p, j] * thetaZ[q, j], 0), 1)
            count <- count + 1
          }
        }
        if (WH == T) {
          if ((thetaX[p, j] != 0) | (thetaX[q, j] != 0)) {
            if (count <= n_rel_inter[1]) {
              inter_xz[p, q, j] <- thetaX[p, j] * thetaZ[q, j]
              count <- count + 1
            }
          }
        }
      }
    }
  }

  inter_zz <- array(0, dim = c(Q, Q, J))
  # interactions zz
  for (j in 1:J) {
    count <- 0
    for (p in 1:Q) {
      for (q in 1:Q) {
        if (q > p) {
          if ((thetaZ[p, j] != 0) & (thetaZ[q, j] != 0)) {
            if (count <= n_rel_inter[2]) {
              inter_zz[p, q, j] <- thetaZ[p, j] * thetaZ[q, j]
              count <- count + 1
            }
          }
          if (WH == T) {
            if ((thetaZ[p, j] != 0) | (thetaZ[q, j] != 0)) {
              if (count <= n_rel_inter[2]) {
                inter_zz[p, q, j] <- c(thetaZ[p, j] * thetaZ[q, j])
                count <- count + 1
              }
            }
          }
        }
      }
    }
  }

  inter_xx <- array(0, dim = c(P, P, J))
  # interactions xx
  for (j in 1:J) {
    count <- 0
    for (p in 1:P) {
      for (q in 1:P) {
        if (q > p) {
          if ((thetaX[p, j] != 0) & (thetaX[q, j] != 0)) {
            if (count <= n_rel_inter[3]) {
              inter_xx[p, q, j] <- thetaX[p, j] * thetaX[q, j]
              count <- count + 1
            }
          }
          if (WH == T) {
            if ((thetaX[p, j] != 0) | (thetaX[q, j] != 0)) {
              if (count <= n_rel_inter[3]) {
                if (stats::runif(1) < .5) {
                  inter_xx[p, q, j] <- thetaX[p, j] * thetaX[q, j]
                }
                count <- count + 1
              }
            }
          }
        }
      }
    }
  }

  inter_xx_nl <- array(0, dim = c(n, P, J))
  nlind <- array(0, dim = c(P, P, J))
  if (nli == TRUE) {
    for (j in 1:J) {
      count <- 0
      for (p in 1:P) {
        for (q in 1:P) {
          if (q > p) {
            if (inter_xx[p, q, j] != 0) {
              if (count <= n_rel_inter[4]) {
                inter_xx_nl[, p, j] <- 1.5 * (X[, q]^2 - 1)
                count <- count + 1
                nlind[p, q, j] <- 1
              }
            }
          }
        }
      }
    }
  }

  basis <- DesignMatrixAlt(X)
  Xilp <- basis$Xlin
  Xinlp <- basis$Xnlin
  dj <- basis$dj

  random_intercept <- NULL
  groupindex <- NULL

  ZETA <- matrix(0, n, J)
  intercept <- stats::runif(J, -2.3, 2.3)

  subspecxz <- array(0, dim = c(n, P, J))
  subspeczz <- array(0, dim = c(n, Q, J))
  subspecxx <- array(0, dim = c(n, P, J))
  part1 <- array(0, dim = c(n, P, J))
  part2 <- array(0, dim = c(n, Q, J))
  SS1 <- matrix(0, n, J)
  SS2 <- matrix(0, n, J)

  so <- 1.0

  # construction of linear predictor
  for (j in 1:J) {
    subspecxz[, , j] <- threshold_mat(t(inter_xz[, , j] %*% t(Z)) +
      t(inter_xx[, , j] %*% t(X)) + inter_xx_nl[, , j], so)
    subspeczz[, , j] <- threshold_mat(t(inter_zz[, , j] %*% t(Z)), so)
    part1[, , j] <- t(thetaX[, j] * matrix(1, P, n)) + subspecxz[, , j]
    part2[, , j] <- t(thetaZ[, j] * matrix(1, Q, n)) + subspeczz[, , j]
    SS1[, j] <- as.matrix(apply(X * part1[, , j], 1, sum))
    SS2[, j] <- as.matrix(apply(Z * part2[, , j], 1, sum))
    ZETA[, j] <- intercept[j] + SS1[, j] + SS2[, j]
  }
  # cat(numb / den, "\n")

  n_reads_min <- 1000
  n_reads_max <- n_reads_min * 2
  theta0 <- theta0
  
  if (randint == TRUE) {
    ZETA <- ZETA[rep(1:nrow(ZETA), times = rep(repmeas, n)), ]
    X <- X[rep(1:nrow(X), times = rep(repmeas, n)), ]
    Z <- Z[rep(1:nrow(Z), times = rep(repmeas, n)), ]
    ri <- rnorm(n, 0, 1)
    rie <- ri[rep(1:length(ri), times = rep(repmeas, n))]
    for(ii in 1:nrow(ZETA)){
      ZETA[ii,] <- ZETA[ii,] + rie[ii]
    }
    labelX <- rep(1:n, times = rep(repmeas, n))
    rownames(X) <- as.character(labelX)
  }
  
  
  Y <- matrix(0, nrow(ZETA), J)
  phi <- matrix(0, nrow(Y), J)
  ct0 <- sample(n_reads_min:n_reads_max, nrow(Y), replace = T)
  
  # sampling counts
  for (i in 1:nrow(Y)) {
    thisrow <- as.vector(exp(ZETA[i, ]))
    phi[i, ] <- thisrow / sum(thisrow)
    Y[i, ] <- dirmult::simPop(J = 1, n = ct0[i], pi = phi[i, ], theta = theta0)$data[1, ]
  }
  
  

  return(list(
    Y = Y, X = X, Z = Z, Xl = Xilp, Xnl = Xinlp, dj = dj,
    intercept = intercept, thetax = thetaX, thetaz = thetaZ,
    inter_xz = inter_xz, inter_xx = inter_xx,
    inter_xx_nl = inter_xx_nl, nlind = nlind, inter_zz = inter_zz,
    RI = random_intercept, GI = groupindex,
    ZETA = ZETA, part1 = part1, part2 = part2 # SS1 = SS1, SS2 = SS2
  ))
} # closes the function
