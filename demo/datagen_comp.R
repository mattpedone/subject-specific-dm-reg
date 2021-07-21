  #' simulate data from a Dirichlet-Multinomial regression model
  #'
  #' @note Requires the dirmult and MASS packages
  #' @note Adattata dal codice di Raffaele
  #'
  #' @param J: number of species
  #' @param n: the number of samples
  #' @param P: number of continuous covariates
  #' @param Q: number of discrete covariates
  #' @param beta_min: minimum absolute value of the regression parameters
  #' @param beta_max: maximum absolute value of the regression parameters
  #' @param n_relevant_taxa: number of relevant species
  #' @param n_relevant_x: number of relevant continuous covariates
  #' @param n_relevant_z: number of relevant discrete covariates
  #' @param n_relevant_inter: vector of length four. number of relevant 
  #' interactions. first element is the number of continuous-discrete interactions,
  #' second element is the number of discrete-discrete interactions, third
  #' element is the number of linear continuous-continuous interactions, fourth
  #' element is the number of nonlinear continuous-continuous interactions
  #' @param corrx: the correlation between continuous covariates
  #' @param corrz: the correlation between discrete covariates
  #' @param theta0: the dispersion parameter
  #' @param WH: the weak heredity flag
  #' @param nli: wheter to include nonlinear interactions between continuous covariates
  #' @param f: interactions between continuous covariates are quadratic (default), or sin \code{f = "sin"} sinusoidal function is adopted
  #' @param seed: for reproducibility. default is \code{seed = 121}
  
  #'
  #' @return Y: (count matrix) rows: n samples, columns: J species
  #' @return X: (design matrix) n * P (continuous variables)
  #' @return Z: (design matrix) n * Q (discrete variables)
  #' @return Xl: linear matrix for peNMIG
  #' @return Xnl: nonlinear matrix for peNMIG
  #' @return dj: vector of length \code{P = ncol(X)}. Indicates into how many 
  #' columns have been expanded each X column. 
  #' @return intercept: simulated intercept vector 
  #' @return thetax: simulated coefficient matrix P * J
  #' @return thetaz: simulated coefficient matrix Q * J
  #' @return inter_xz: simulated coefficient matrix (continuous-discrete 
  #' interactions) P x Q x J 
  #' @return inter_xx: simulated coefficient matrix (continuous-continuous
  #' interactions - linear part) P x P x J 
  #' @return inter_xx_nl: simulated coefficient matrix (continuous-continuous
  #' interactions - nonlinear part) n x P x J 
  #' @return nlind: number of non null continuous-continuous interactions
  #' @return inter_zz: simulated coefficient matrix (discrete-discrete 
  #' interactions) Q x Q x J 
  #' @return thetaX_exp: analogous coefficient matrix for comparison methods
  #' @return subspecxz, subspeczz subspecxx: subject specific effects
  #' @return part1, part2, SS1, SS2: addendi del predittore lineare
  #' @return ZETA: linear predictor
  #' @return thetaX_exp: analogous design matrix for comparison methods
  #'
  #' @export
  #' 
  
gendata_simple_confr <- function(J = 10, n = 100, P = 5, Q = 5, beta_min = 0.75, 
                                 beta_max = 2.00, n_relevant_taxa = 2, 
                                 n_relevant_x = 2, n_relevant_z = 2, n_relevant_inter = c(2, 2, 3, 3), 
                                 corrx = 0, corrz = 0, theta0 = 0.01, WH = FALSE, 
                                 nli = T, f = "quadratic", seed = 121) {
  
  set.seed(seed)
  
  threshold_mat <- function(m, t){
    mt <- m
    mt[which(abs(m)<t)] <- 0.0
    return(mt)
  }
  
  #n_relevant_inter <- c(2, 2, 3, 3) in strong heredity 
  #n_relevant_inter <- c(2, 2, 5, 5) in weak heredity 
  # function definition
  threshold_mat <- function(m, t) {
    mt <- m
    mt[which(abs(m) < t)] <- 0.0
    return(mt)
  }
  
  #construct continuous covariates
  if (corrx != 0) {
    Sigmax = matrix(1, P, P)
    Sigmax = corrx^abs(row(Sigmax) - col(Sigmax))
  } else {
    Sigmax = diag(1, P, P)
  }
  X <- scale(MASS::mvrnorm(n = n, mu = rep(0, P), Sigma = Sigmax))
  
  #construct discrete covariates
  if (corrz != 0) {
    Sigmaz = matrix(1, Q, Q)
    Sigmaz = corrz^abs(row(Sigmaz) - col(Sigmaz))
  } else {
    Sigmaz = diag(1, Q, Q)
  }
  Z <- scale(MASS::mvrnorm(n = n, mu = rep(0.5, Q), Sigma = Sigmaz))
  Z[which(Z >= 0.0)] <- 1; Z[which(Z < 0.0)] <- 0
  
  # main effects X
  thetaX <- matrix(0, P, J)
  st = 0
  low_side = beta_min
  high_side = beta_max
  if (n_relevant_x != 1) {
    # warning if the lengths don't match
    coef = suppressWarnings(seq(low_side, high_side, len = n_relevant_x) *                              
                              c(1,-1))
  } else{
    coef = (low_side + high_side) / 2
  }
  coef_g = rep(1.0, len = n_relevant_taxa)
  for (ii in 1:n_relevant_taxa) {
    # overlap species
    thetaX[(st:(st + n_relevant_x - 1)) %% P + 1, 3 * ii - 2] = coef_g[ii] *
      sample(coef)[((ii - 1):(ii + n_relevant_x - 2)) %% n_relevant_x + 1]
    st = st + 1
  }
  
  #thetaX <- t(thetaX)
  
  # main effects Z
  thetaZ <- matrix(0, Q, J)
  st = 0
  low_side = beta_min
  high_side = beta_max
  if (n_relevant_z != 1) {
    # warning if the lengths don't match
    coef = suppressWarnings(seq(low_side * 1.5, high_side * 1.5, len = n_relevant_z) * c(1,-1))
  } else{
    coef = (low_side + high_side) / 2
  }
  coef_g = rep(1.0, len = n_relevant_taxa)
  for (ii in 1:n_relevant_taxa) {
    # overlap species
    thetaZ[(st:(st + n_relevant_z - 1)) %% Q + 1, 3 * ii - 2] = coef_g[ii] * sample(coef)[((ii - 1):(ii + n_relevant_z - 2)) %% n_relevant_z + 1]
    st = st + 1
  }
  
  inter_xz <- array(0, dim = c(P, Q, J))
  #interactions xz
  inter_xz[sample(P, 1), sample(Q, 1), sample(J, 1)] <- runif(1,1,2)*sample(c(-1,1), 1)
  inter_xz[sample(P, 1), sample(Q, 1), sample(J, 1)] <- runif(1,1,2)*sample(c(-1,1), 1)
  inter_xz[sample(P, 1), sample(Q, 1), sample(J, 1)] <- runif(1,1,2)*sample(c(-1,1), 1)
  
  inter_zz <- array(0, dim = c(Q, Q, J))
  #interactions zz
  inter_zz[sample(Q, 1), sample(Q, 1), sample(J, 1)] <- runif(1,1,2)*sample(c(-1,1), 1)
  inter_zz[sample(Q, 1), sample(Q, 1), sample(J, 1)] <- runif(1,1,2)*sample(c(-1,1), 1)
  inter_zz[sample(Q, 1), sample(Q, 1), sample(J, 1)] <- runif(1,1,2)*sample(c(-1,1), 1)
  inter_zz[sample(Q, 1), sample(Q, 1), sample(J, 1)] <- runif(1,1,2)*sample(c(-1,1), 1)
  inter_zz[sample(Q, 1), sample(Q, 1), sample(J, 1)] <- runif(1,1,2)*sample(c(-1,1), 1)
  
  for(j in 1:J) {
    count = 0
    for (p in 1:Q) {
      for (q in 1:Q) {
        if (q <= p) {
          inter_zz[p, q, j] <- 0
        }
      }
    }
  }
  
  inter_xx <- array(0, dim = c(P, P, J))
  #interactions xx
  inter_xx[sample(P, 1), sample(P, 1), sample(J, 1)] <- runif(1,1,2)*sample(c(-1,1), 1)
  inter_xx[sample(P, 1), sample(P, 1), sample(J, 1)] <- runif(1,1,2)*sample(c(-1,1), 1)
  inter_xx[sample(P, 1), sample(P, 1), sample(J, 1)] <- runif(1,1,2)*sample(c(-1,1), 1)
  inter_xx[sample(P, 1), sample(P, 1), sample(J, 1)] <- runif(1,1,2)*sample(c(-1,1), 1)
  inter_xx[sample(P, 1), sample(P, 1), sample(J, 1)] <- runif(1,1,2)*sample(c(-1,1), 1)
  
  for(j in 1:J) {
    count = 0
    for (p in 1:P) {
      for (q in 1:P) {
        if (q <= p) {
          inter_xx[p, q, j] <- 0
        }
      }
    }
  }
  
  inter_xx_nl <- array(0, dim = c(n, P, J))
  nlind <- array(0, dim = c(P, P, J))
  
  if(nli == TRUE){
    #if(f == "quadratic"){
    for (j in 1:J) {
      count = 0
      for (p in 1:P) {
        for(q in 1:P){
          if(q > p){
            if (inter_xx[p, q, j] != 0) {
              if (count < n_relevant_inter[4]) {
                inter_xx_nl[, p, j] = 1.5*(X[, q]^2-1)
                count = count + 1
                nlind[p, q, j] <- 1
              }
            }
          }
        }
      }
    }
    #}
  }
  
  thetaX_exp <- matrix(NA, P+Q+(P*Q)+(Q*(Q-1))/2 + (P*(P-1))/2, J)
  thetaX_exp[c(1:P),c(1:J)] <- thetaX
  thetaX_exp[c((P + 1):(P + Q)),c(1:J)] <- thetaZ
  
  for(j in 1:J){
    for(p in 1:P){
      for(q in 1:Q){
        thetaX_exp[P+Q+(p-1)*P+(q-1)+1,j] <- inter_xz[p, q, j]
        #cat(P+Q+p+(q-1), "\n")
      }
    }
  }
  
  # I am sure there is a better way ZZ
  pis <- c()
  posp <- c()
  posq <- c()
  for(p in 1:(Q - 1)){
    for(q in (p+1):Q){
      pis <- c(pis, 1)
      posp[sum(cumsum(pis))-length(pis)+1] <- p
      posq[sum(cumsum(pis))-length(pis)+1] <- q
    }
  }
  
  posp <- posp[is.na(posp)!=1]
  posq <- posq[is.na(posq)!=1]
  pis <- cumsum(pis)
  
  count <- 1
  for(j in 1:J){
    for(i in 1:length(pis)){
      thetaX_exp[P+Q+(P*Q)+pis[i],j] <- inter_zz[posp[i], posq[i], j]
    }
  }
  
  # I am sure there is a better way XX
  pis <- c()
  posp <- c()
  posq <- c()
  for(p in 1:(P - 1)){
    for(q in (p+1):P){
      pis <- c(pis, 1)
      posp[sum(cumsum(pis))-length(pis)+1] <- p
      posq[sum(cumsum(pis))-length(pis)+1] <- q
    }
  }
  
  posp <- posp[is.na(posp)!=1]
  posq <- posq[is.na(posq)!=1]
  pis <- cumsum(pis)
  
  count <- 1
  for(j in 1:J){
    for(i in 1:length(pis)){
      thetaX_exp[P+Q+(P*Q)+(Q*(Q-1))/2+pis[i],j] <- inter_xx[posp[i], posq[i], j]
    }
  }
  
  DesignMatrixAlt <- function(X, rango = .95){
    MAT <- matrix(0, dim(X)[1], dim(X)[2])
    MAT2 <- c()
    expa <- c()
    for(p in 1:ncol(X)){
      x <- X[,p]
      #z0 <- spikeSlabGAM::lin(x)
      MAT[,p] <- spikeSlabGAM::lin(x)
      #MAT <- cbind(MAT, z) #Scheipl's version
      z <- spikeSlabGAM::sm(x, K = 20, spline.degree = 3, diff.ord = 2,
                            rankZ = rango, centerBase = T, centerx = x, 
                            decomposition = "ortho")
      #z <- orthoDesign(psBasis(x)$X, rankZ = .9)
      MAT2 <- cbind(MAT2, z)
      expa <- c(expa, rep((p-1), ncol(z)))
      #Design[[p]] <- MAT
    }
    dj <- as.vector(table(expa))
    return(list(Xlin = MAT, Xnlin = MAT2, dj=dj))
  }
  
  basis <- DesignMatrixAlt(X)
  Xilp <- basis$Xlin
  Xinlp <- basis$Xnlin
  dj <- basis$dj
  
  
  ZETA <- matrix(0, n, J)
  intercept = runif(J, -2.3, 2.3)
  
  subspecxz <- array(0, dim=c(n, P, J))
  subspeczz <- array(0, dim=c(n, Q, J))
  subspecxx <- array(0, dim=c(n, P, J))
  part1 <- array(0, dim=c(n, P, J))
  part2 <- array(0, dim=c(n, Q, J))
  SS1 <- matrix(0, n, J)
  SS2 <- matrix(0, n, J)
  
  numb <- den <- 0
  so <- 1.5
  
  #construction of linear predictor
  for (j in 1:J) {
    numb <- numb + sum(abs(inter_xz[,,j] %*% t(Z))>so) + sum(abs(inter_zz[,,j] %*% t(Z))>so) + 
      sum(abs(inter_xx[,,j] %*% t(X))>so)
    den <- den + sum(abs(inter_xz[,,j] %*% t(Z)) != 0.0) + sum(abs(inter_zz[,,j] %*% t(Z)) != 0.0) +
      sum(abs(inter_xx[,,j] %*% t(X))>so)
    subspecxz[,,j] <- threshold_mat(t(inter_xz[,,j] %*% t(Z)) + 
                                      t(inter_xx[,,j] %*% t(X)) + inter_xx_nl[,,j], so)
    subspeczz[,,j] <- threshold_mat(t(inter_zz[,,j] %*% t(Z)), so)
    #subspecxx[,,j] <- threshold_mat(t(inter_xx[,,j] %*% t(X)), so)
    part1[,,j] <- t(thetaX[, j]*matrix(1, P, n)) + subspecxz[,,j]# + subspecxx[,,j]
    part2[,,j] <- t(thetaZ[, j]*matrix(1, Q, n)) + subspeczz[,,j]
    SS1[,j] <- as.matrix(apply(X*part1[,,j], 1, sum))
    SS2[,j] <- as.matrix(apply(Z*part2[,,j], 1, sum))
    ZETA[, j] <- intercept[j] + SS1[,j] + SS2[,j]
  }
  cat(numb / den, "\n")
  
  n_reads_min <- 1000
  n_reads_max <- n_reads_min * 2
  theta0 = theta0
  phi = matrix(0, n, J)
  ct0 = sample(n_reads_min:n_reads_max, n, replace = T)
  
  Y <- matrix(0, n, J)
  #Ypred <- matrix(0, n, J)
  #sampling counts
  for (i in 1:n) {
    thisrow = as.vector(exp(ZETA[i,]))
    phi[i, ] = thisrow / sum(thisrow)
    Y[i, ] = dirmult::simPop(J = 1, n = ct0[i], pi = phi[i,], theta = theta0)$data[1,]
    #Ypred[i, ] = dirmult::simPop(J = 1, n = ct0[i], pi = phi[i,], theta = theta0)$data[1,]
  }
  
  #Pi <- P; Qu <- Q;
  XX <- cbind(X, Z)
  
  Xex = XX
  for (i in 1:(P)) {
    for(j in (P+1):(P+Q)){
      Xex = cbind(Xex, XX[,i]*XX[,j])
    }
  }
  for (i in (P+1):(P+Q)) {
    if(i<(P+Q)){
      for(j in (i+1):(P+Q)){
        Xex = cbind(Xex, XX[,i]*XX[,j])
      }
    }
  }
  for (i in (1):(P)) {
    if(i<(P)){
      for(j in (i+1):(P)){
        Xex = cbind(Xex, XX[,i]*XX[,j])
      }
    }
  }
  XX=Xex
  
  #TRUTH <- list(thetaX=thetaX, thetaZ=thetaZ)
  return(list(Y = Y, X = X, Z=Z, Xl = Xilp, Xnl = Xinlp, dj = dj, intercept=intercept, 
              thetax=thetaX, thetaz=thetaZ, inter_xz=inter_xz, inter_xx = inter_xx, 
              inter_xx_nl = inter_xx_nl, nlind=nlind, inter_zz=inter_zz, 
              thetaX_exp = thetaX_exp, subspecxz = subspecxz, subspeczz = subspeczz, 
              subspecxx = subspecxx, part1=part1, part2=part2, SS1=SS1, SS2=SS2, 
              ZETA = ZETA, DesMat = XX))
  #Ypred = Ypred, ME_exp=thetaX_exp, 
}

  #gendata_simple_confr <- function(X, Z) {
  #  
  #  P = ncol(X)
  #  Q = ncol(Z)
  #  
  #  XX <- cbind(X, Z)
  #  
  #  Xex = XX
  #  for (i in 1:(P)) {
  #    for(j in (P+1):(P+Q)){
  #      Xex = cbind(Xex, XX[,i]*XX[,j])
  #    }
  #  }
  #  for (i in (P+1):(P+Q)) {
  #    if(i<(P+Q)){
  #      for(j in (i+1):(P+Q)){
  #        Xex = cbind(Xex, XX[,i]*XX[,j])
  #      }
  #    }
  #  }
  #  for (i in (1):(P)) {
  #    if(i<(P)){
  #      for(j in (i+1):(P)){
  #        Xex = cbind(Xex, XX[,i]*XX[,j])
  #      }
  #    }
  #  }
  #  XX=Xex
  #  
  #  #TRUTH <- list(thetaX=thetaX, thetaZ=thetaZ)
  #  return(XX)
  #  #Ypred = Ypred, ME_exp=thetaX_exp, 
  #}
  
  
  
  gendata_simple_confr_old <- function(J = 10, n = 100, P = 5, Q = 5, beta_min = 0.75, 
                                   beta_max = 2.00, n_relevant_taxa = 2, 
                                   n_relevant_x = 2, n_relevant_z = 2, n_relevant_inter = c(2, 2, 3, 3), 
                                   corrx = 0, corrz = 0, theta0 = 0.01, WH = FALSE, 
                                   nli = T, f = "quadratic", seed = 121) {
    
    threshold_mat <- function(m, t){
      mt <- m
      mt[which(abs(m)<t)] <- 0.0
      return(mt)
    }
    
    #n_relevant_inter <- c(2, 2, 3, 3) in strong heredity 
    #n_relevant_inter <- c(2, 2, 5, 5) in weak heredity 
    # function definition
    threshold_mat <- function(m, t) {
      mt <- m
      mt[which(abs(m) < t)] <- 0.0
      return(mt)
    }
    
    #construct continuous covariates
    if (corrx != 0) {
      Sigmax = matrix(1, P, P)
      Sigmax = corrx^abs(row(Sigmax) - col(Sigmax))
    } else {
      Sigmax = diag(1, P, P)
    }
    X <- scale(MASS::mvrnorm(n = n, mu = rep(0, P), Sigma = Sigmax))
    
    #construct discrete covariates
    if (corrz != 0) {
      Sigmaz = matrix(1, Q, Q)
      Sigmaz = corrz^abs(row(Sigmaz) - col(Sigmaz))
    } else {
      Sigmaz = diag(1, Q, Q)
    }
    Z <- scale(MASS::mvrnorm(n = n, mu = rep(0.5, Q), Sigma = Sigmaz))
    Z[which(Z >= 0.0)] <- 1; Z[which(Z < 0.0)] <- 0
    
    # main effects X
    thetaX <- matrix(0, P, J)
    st = 0
    low_side = beta_min
    high_side = beta_max
    if (n_relevant_x != 1) {
      # warning if the lengths don't match
      coef = suppressWarnings(seq(low_side, high_side, len = n_relevant_x) *                              
                                c(1,-1))
    } else{
      coef = (low_side + high_side) / 2
    }
    coef_g = rep(1.0, len = n_relevant_x)
    for (ii in 1:n_relevant_x) {
      # overlap species
      thetaX[(st:(st + n_relevant_x - 1)) %% P + 1, 3 * ii - 2] = coef_g[ii] *
        sample(coef)[((ii - 1):(ii + n_relevant_x - 2)) %% n_relevant_x + 1]
      st = st + 1
    }
    
    # main effects Z
    thetaZ <- matrix(0, Q, J)
    st = 0
    low_side = beta_min
    high_side = beta_max
    if (n_relevant_z != 1) {
      # warning if the lengths don't match
      coef = suppressWarnings(seq(low_side * 1.5, high_side * 1.5, len = n_relevant_z) * c(1,-1))
    } else{
      coef = (low_side + high_side) / 2
    }
    coef_g = rep(1.0, len = n_relevant_z)
    for (ii in 1:n_relevant_z) {
      # overlap species
      thetaZ[(st:(st + n_relevant_z - 1)) %% Q + 1, 3 * ii - 2] = coef_g[ii] * sample(coef)[((ii - 1):(ii + n_relevant_z - 2)) %% n_relevant_z + 1]
      st = st + 1
    }
    
    inter_xz <- array(0, dim = c(P, Q, J))
    #interactions xz
    for (j in 1:J) {
      count = 0
      for (p in 1:P) {
        for (q in 1:Q) {
          if ((thetaX[p, j] != 0) & (thetaZ[q, j] != 0)) {
            if (count < n_relevant_inter[1]) {
              inter_xz[p, q, j] = sample(c(thetaX[p, j] * thetaZ[q, j], 0), 1)
              count = count + 1
            }
          }
          if (WH == T) {
            if ((thetaX[p, j] != 0) | (thetaX[q, j] != 0)) {
              if (count < n_relevant_inter[1]) {
                inter_xz[p, q, j] = thetaX[p, j] * thetaZ[q, j]
                count = count + 1
              }
            }
          }
        }
      }
    }
    
    inter_zz <- array(0, dim = c(Q, Q, J))
    #interactions zz
    for (j in 1:J) {
      count = 0
      for (p in 1:Q) {
        for (q in 1:Q) {
          if (q > p) {
            if ((thetaZ[p, j] != 0) & (thetaZ[q, j] != 0)) {
              if (count < n_relevant_inter[2]) {
                inter_zz[p, q, j] = thetaZ[p, j] * thetaZ[q, j]
                count = count + 1
              }
            }
            if (WH == T) {
              if ((thetaZ[p, j] != 0) | (thetaZ[q, j] != 0)) {
                if (count < n_relevant_inter[2]) {
                  inter_zz[p, q, j] = c(thetaZ[p, j] * thetaZ[q, j])
                  count = count + 1
                }
              }
            }
          }
        }
      }
    }
    
    inter_xx <- array(0, dim = c(P, P, J))
    #interactions xx
    for(j in 1:J) {
      count = 0
      for (p in 1:P) {
        for (q in 1:P) {
          if (q > p) {
            if ((thetaX[p, j] != 0) & (thetaX[q, j] != 0)) {
              if (count < n_relevant_inter[3]) {
                inter_xx[p, q, j] = thetaX[p, j] * thetaX[q, j]
                count = count + 1
              }
            }
            if (WH == T) {
              if ((thetaX[p, j] != 0) | (thetaX[q, j] != 0)) {
                if (count < n_relevant_inter[3]) {
                  if(runif(1)<.5){
                    #inter_xx[p, q, j] = sample(c(thetaX[p, j] * thetaX[q, j], 0), 1)
                    inter_xx[p, q, j] = thetaX[p, j] * thetaX[q, j]
                  }#sample(c(thetaZ[p, j] * thetaZ[q, j], 0), 1)
                  count = count + 1
                }
              }
            }
          }
        }
      }
    }
    
    inter_xx_nl <- array(0, dim = c(n, P, J))
    nlind <- array(0, dim = c(P, P, J))
    if(nli == TRUE){
      if(f == "quadratic"){
        for (j in 1:J) {
          count = 0
          for (p in 1:P) {
            for(q in 1:P){
              if(q > p){
                if ((thetaX[p, j] != 0) & (thetaX[q, j] != 0)) {
                  if (count < n_relevant_inter[4]) {
                    inter_xx_nl[, p, j] = 1.5*(X[, q]^2-1)
                    count = count + 1
                    nlind[p, q, j] <- 1
                  }
                }
              }
            }
          }
          if (WH == T) {
            for (p in 1:P) {
              for(q in 1:P){
                if(q > p){
                  if ((thetaX[p, j] != 0) | (thetaX[q, j] != 0)) {
                    if (count < n_relevant_inter[4]) {
                      if(runif(1)<.5){
                        inter_xx_nl[, p, j] = 1.5*(X[, q]^2-1)
                      }#sample(c(thetaZ[p, j] * thetaZ[q, j], 0), 1)
                      count = count + 1
                      nlind[p, q, j] <- 1
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    
    thetaX_exp <- matrix(NA, P+Q+(P*Q)+(Q*(Q-1))/2 + (P*(P-1))/2, J)
    thetaX_exp[c(1:P),c(1:J)] <- thetaX
    thetaX_exp[c((P + 1):(P + Q)),c(1:J)] <- thetaZ
    
    for(j in 1:J){
      for(p in 1:P){
        for(q in 1:Q){
          thetaX_exp[P+Q+(p-1)*P+(q-1)+1,j] <- inter_xz[p, q, j]
          #cat(P+Q+p+(q-1), "\n")
        }
      }
    }
    
    # I am sure there is a better way ZZ
    pis <- c()
    posp <- c()
    posq <- c()
    for(p in 1:(Q - 1)){
      for(q in (p+1):Q){
        pis <- c(pis, 1)
        posp[sum(cumsum(pis))-length(pis)+1] <- p
        posq[sum(cumsum(pis))-length(pis)+1] <- q
      }
    }
    
    posp <- posp[is.na(posp)!=1]
    posq <- posq[is.na(posq)!=1]
    pis <- cumsum(pis)
    
    count <- 1
    for(j in 1:J){
      for(i in 1:length(pis)){
        thetaX_exp[P+Q+(P*Q)+pis[i],j] <- inter_zz[posp[i], posq[i], j]
      }
    }
    
    # I am sure there is a better way XX
    pis <- c()
    posp <- c()
    posq <- c()
    for(p in 1:(P - 1)){
      for(q in (p+1):P){
        pis <- c(pis, 1)
        posp[sum(cumsum(pis))-length(pis)+1] <- p
        posq[sum(cumsum(pis))-length(pis)+1] <- q
      }
    }
    
    posp <- posp[is.na(posp)!=1]
    posq <- posq[is.na(posq)!=1]
    pis <- cumsum(pis)
    
    count <- 1
    for(j in 1:J){
      for(i in 1:length(pis)){
        thetaX_exp[P+Q+(P*Q)+(Q*(Q-1))/2+pis[i],j] <- inter_xx[posp[i], posq[i], j]
      }
    }
    
    DesignMatrixAlt <- function(X, rango = .95){
      MAT <- matrix(0, dim(X)[1], dim(X)[2])
      MAT2 <- c()
      expa <- c()
      for(p in 1:ncol(X)){
        x <- X[,p]
        #z0 <- spikeSlabGAM::lin(x)
        MAT[,p] <- spikeSlabGAM::lin(x)
        #MAT <- cbind(MAT, z) #Scheipl's version
        z <- spikeSlabGAM::sm(x, K = 20, spline.degree = 3, diff.ord = 2,
                              rankZ = rango, centerBase = T, centerx = x, 
                              decomposition = "ortho")
        #z <- orthoDesign(psBasis(x)$X, rankZ = .9)
        MAT2 <- cbind(MAT2, z)
        expa <- c(expa, rep((p-1), ncol(z)))
        #Design[[p]] <- MAT
      }
      dj <- as.vector(table(expa))
      return(list(Xlin = MAT, Xnlin = MAT2, dj=dj))
    }
    
    basis <- DesignMatrixAlt(X)
    Xilp <- basis$Xlin
    Xinlp <- basis$Xnlin
    dj <- basis$dj
    
    
    ZETA <- matrix(0, n, J)
    intercept = runif(J, -2.3, 2.3)
    
    subspecxz <- array(0, dim=c(n, P, J))
    subspeczz <- array(0, dim=c(n, Q, J))
    subspecxx <- array(0, dim=c(n, P, J))
    part1 <- array(0, dim=c(n, P, J))
    part2 <- array(0, dim=c(n, Q, J))
    SS1 <- matrix(0, n, J)
    SS2 <- matrix(0, n, J)
    
    numb <- den <- 0
    so <- 2.0
    
    #construction of linear predictor
    for (j in 1:J) {
      numb <- numb + sum(abs(inter_xz[,,j] %*% t(Z))>so) + sum(abs(inter_zz[,,j] %*% t(Z))>so) + 
        sum(abs(inter_xx[,,j] %*% t(X))>so)
      den <- den + sum(abs(inter_xz[,,j] %*% t(Z)) != 0.0) + sum(abs(inter_zz[,,j] %*% t(Z)) != 0.0) +
        sum(abs(inter_xx[,,j] %*% t(X))>so)
      subspecxz[,,j] <- threshold_mat(t(inter_xz[,,j] %*% t(Z)) + 
                                        t(inter_xx[,,j] %*% t(X)) + inter_xx_nl[,,j], so)
      subspeczz[,,j] <- threshold_mat(t(inter_zz[,,j] %*% t(Z)), so)
      #subspecxx[,,j] <- threshold_mat(t(inter_xx[,,j] %*% t(X)), so)
      part1[,,j] <- t(thetaX[, j]*matrix(1, P, n)) + subspecxz[,,j]# + subspecxx[,,j]
      part2[,,j] <- t(thetaZ[, j]*matrix(1, Q, n)) + subspeczz[,,j]
      SS1[,j] <- as.matrix(apply(X*part1[,,j], 1, sum))
      SS2[,j] <- as.matrix(apply(Z*part2[,,j], 1, sum))
      ZETA[, j] <- intercept[j] + SS1[,j] + SS2[,j]
    }
    cat(numb / den, "\n")
    
    n_reads_min <- 1000
    n_reads_max <- n_reads_min * 2
    theta0 = theta0
    phi = matrix(0, n, J)
    ct0 = sample(n_reads_min:n_reads_max, n, replace = T)
    
    Y <- matrix(0, n, J)
    #Ypred <- matrix(0, n, J)
    #sampling counts
    for (i in 1:n) {
      thisrow = as.vector(exp(ZETA[i,]))
      phi[i, ] = thisrow / sum(thisrow)
      Y[i, ] = dirmult::simPop(J = 1, n = ct0[i], pi = phi[i,], theta = theta0)$data[1,]
      #Ypred[i, ] = dirmult::simPop(J = 1, n = ct0[i], pi = phi[i,], theta = theta0)$data[1,]
    }
    
    #Pi <- P; Qu <- Q;
    XX <- cbind(X, Z)
    
    Xex = XX
    for (i in 1:(P)) {
      for(j in (P+1):(P+Q)){
        Xex = cbind(Xex, XX[,i]*XX[,j])
      }
    }
    for (i in (P+1):(P+Q)) {
      if(i<(P+Q)){
        for(j in (i+1):(P+Q)){
          Xex = cbind(Xex, XX[,i]*XX[,j])
        }
      }
    }
    for (i in (1):(P)) {
      if(i<(P)){
        for(j in (i+1):(P)){
          Xex = cbind(Xex, XX[,i]*XX[,j])
        }
      }
    }
    XX=Xex
    
    #TRUTH <- list(thetaX=thetaX, thetaZ=thetaZ)
    return(list(Y = Y, X = X, Z=Z, Xl = Xilp, Xnl = Xinlp, dj = dj, intercept=intercept, 
                thetax=thetaX, thetaz=thetaZ, inter_xz=inter_xz, inter_xx = inter_xx, 
                inter_xx_nl = inter_xx_nl, nlind=nlind, inter_zz=inter_zz, 
                thetaX_exp = thetaX_exp, subspecxz = subspecxz, subspeczz = subspeczz, 
                subspecxx = subspecxx, part1=part1, part2=part2, SS1=SS1, SS2=SS2, 
                ZETA = ZETA, DesMat = XX))
    #Ypred = Ypred, ME_exp=thetaX_exp, 
  }