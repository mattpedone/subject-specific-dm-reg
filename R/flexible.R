rm(list=ls())
library(Compositional)
Rcpp::sourceCpp("src/sampler.cpp")
source("R/utils.R")
source("R/wrapfn.R")

RMSEK <- RMSEB <- rep(0, 30)
ADK <- ADB <- rep(0, 30)

for(k in 1:30){
gendata_simple_confr <- function(J = 10, n = 100, P = 5, Q = 5, beta_min = 0.75, 
                                 beta_max = 1.00, n_relevant_taxa = 2, 
                                 n_relevant_x = 3, n_relevant_z = 3, n_relevant_inter = c(2, 2, 3, 3), 
                                 corrx = .4, corrz = .4, theta0 = 0.1, WH = T, 
                                 nli = T, f = "quadratic", seed = sample(1:1000, 1)) {
  
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

rmse <- function(x, y) {
  if (length(x) != length(y)) {
    stop(sprintf("x and y are different lengths : %d, %d",
                 length(x), length(y)))
  }
  
  residuals <- x-y
  res2      <- residuals^2
  res2.mean <- sum(res2) / length(x)
  rms       <- sqrt(res2.mean)
  
  return(rms)
  
}

aitdist <- function(x, y){
  if (length(x) != length(y)) {
    stop(sprintf("x and y are different lengths : %d, %d",
                 length(x), length(y)))
  }
  if(any(x<0.00000001)){ 
    x[which(x<0.00000001)] <- 0.00000001
  }
  if(any(y<0.00000001)){ 
    y[which(y<0.00000001)] <- 0.00000001
  }
  gmx <- FSA::geomean(x, zneg.rm = TRUE)
  gmy <- FSA::geomean(y, zneg.rm = TRUE)
  ad <- sqrt(sum((log(x/gmx)-log(y/gmy))^2))
  return(ad)
}

data <- gendata_simple_confr()
#data$DesMat
idx <- c(76:100)
y <- data$Y[-idx,]
yp <- data$Y[idx,]
y <- y / rowSums(y)
ypprop <- yp/rowSums(yp)
x <- data$X[-idx,]
z <- data$Z[-idx,]
dm <- data$DesMat[-idx,]
xp <- data$DesMat[idx,]
mod <- aknnreg.tune(y, dm, a = c(0.4, 0.6), k = 2:4, nfolds = 5)
mod <- aknn.reg(xp, y, dm, a = mod$kl.alpha, k = mod$kl.k, apostasi = "euclidean")
mod <- (mod[[1]])
rmseaknn <- c()
for(i in 1:length(idx)){
  rmseaknn[i] <- rmse(as.vector(mod[[1]][i,]), as.vector(ypprop[i,]))
}
#rmseaknn
#mean(rmseaknn)



out <- ssdm(
  YY = y, XX = x, ZZ = z, Niter = 10000, burn = 5000, 
  penmig_lin = c(5, 25, 0.00025, 1.00, 1.00),
  penmig_nl = c(5, 25, 0.00025, 1.00, 1.00),
  thin = 10, randint = FALSE, hereditariety = 2, init = TRUE,
  init_fd = 0.1, conlnl = TRUE
)

EFF=(100-50)/10#(Niter-burn)/thin

#devo ricostruire i loggamma con le covariate del validation set
loggamma <- out$loglinpred[1:25,,]
phi <- c()
Ypred <- matrix(0, length(idx), ncol(y))
prop <- array(0, dim=c(length(idx), ncol(y), EFF))
ct0 = sample(1000:2000, length(idx), rep = T)
for(eff in 1:EFF){
  for (i in 1:length(idx)) {
    thisrow = as.vector(exp(loggamma[i, ,eff]))
    phi = thisrow / sum(thisrow)
    Ypred[i,] = dirmult::simPop(J = 1, n = ct0[i], pi = phi, theta = 0.01)$data[1,]
    prop[i, , eff] <- Ypred[i,]/rowSums(Ypred)[i]
  }
}

prop <- apply(prop, c(1, 2), mean)

rmsessdm <- c()
for(i in 1:length(idx)){
  rmsessdm[i] <- rmse(as.vector(prop[i,]), as.vector(ypprop[i,]))
}

RMSEK[k] <- sum(rmseaknn)
RMSEB[k] <- sum(rmsessdm)
work <- 0
for(i in 1:length(idx)){
  work <- work +  aitdist(as.vector(mod[[1]][i,]), as.vector(ypprop[i,]))
}
ADK[k] <- work/length(idx)

work <- 0
for(i in 1:length(idx)){
  work <- work +  aitdist(as.vector(prop[i,]), as.vector(ypprop[i,]))
}
ADB[k] <- work/length(idx)
}

mean(RMSEK)
sd(RMSEK)
mean(RMSEB)
sd(RMSEB)
mean(ADK)
sd(ADK)
mean(ADB)
sd(ADB)
#rmsessdm
#mean(rmsessdm)
#prop contiene la distribuzione a posteriori
#ne prendo la media e ci calcolo rmse

#avrmsek <- apply(RMSEK, 2, mean)
#sdrmsek <- sqrt(apply(RMSEK, 2, var))
#avrmseb <- apply(RMSEB, 2, mean)
#sdrmseb <- sqrt(apply(RMSEB, 2, var))
#
#data <- data.frame(rmse = c(avrmsek, avrmseb), method = c(rep("akNN", 25), rep("SSDM", 25)))
#p <- ggplot2::ggplot(data, aes(x=method, y=rmse, color = method)) + 
#  geom_boxplot() + scale_color_manual(values=c("#E69F00", "#56B4E9")) 
#p
#
#lab <- c("aknn", " ", "ssdm", "")
#sc1 <- round(c(mean(avrmsek), sd(avrmsek), mean(avrmseb), sd(avrmseb)), digits = 4)
#tab <- rbind(lab, sc1)
#xtable::xtable(tab)

