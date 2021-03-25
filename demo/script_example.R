rm(list = ls())

devtools::load_all()

J <- 20
n <- 50
P <- 5
Q <- 5

Niter <- 300000
burn <- Niter / 2
thin <- 10

Eff <- (Niter - burn) / thin

data <- gendata(J, n, P, Q, .25, .75, 2, 2, 2, c(1, 1, 1, 1),
  corrx = 0.4, corrz = 0.4, nli = T, WH = F,
  randint = T
)

system.time(out <- ssdm(
  YY = data$Y, XX = data$X, ZZ = data$Z, Niter = Niter, burn = burn,
  penmig_lin = c(5, 25, 0.00025, 1.00, 1.00),
  penmig_nl = c(5, 25, 0.00025, 1.00, 1.00),
  thin = thin, randint = FALSE, hereditariety = 2, init = TRUE,
  init_fd = 0.1, conlnl = TRUE
))

###### effetti principali ######
par(mfrow = c(4, 2))

intercept_post <- apply(out$muposterior, 2, mean)

thetaXpost <- round(apply(out$thetaxposterior, c(1, 2), mean), 2)
thetaZpost <- round(apply(out$thetazposterior, c(1, 2), mean), 2)
ppithetaX <- apply((out$thetaxposterior != 0), c(1, 2), mean)
ppithetaZ <- apply((out$thetazposterior != 0), c(1, 2), mean)

thetaXpostM <- thetaXpost
thetaXpostM[which(ppithetaX < .5)] <- 0.0 # this is the median model
betaincluded1 <- rep(0, length((as.vector(thetaXpostM))))
betaincluded1[which(as.vector(ppithetaX) > .5)] <- 1

plot(betaincluded1,
  type = "h", ylim = c(0, 1), xlab = expression(paste("index for ", theta[x])),
  ylab = expression(paste("PPI ", theta[x]))
)

thetaZpostM <- thetaZpost
thetaZpostM[which(ppithetaZ < .5)] <- 0.0 # this is the median model
betaincluded2 <- rep(0, length((as.vector(thetaZpostM))))
betaincluded2[which(as.vector(ppithetaZ) > .5)] <- 1

plot(betaincluded2,
  type = "h", ylim = c(0, 1), xlab = expression(paste("index for ", theta[z])),
  ylab = expression(paste("PPI ", theta[z]))
)

bipostxz <- array(0, dim = c(P, Q, J, Eff))

for (ss in 1:Eff) {
  for (jj in 1:J) {
    for (qq in 1:Q) {
      for (pp in 1:P) {
        hh <- qq + ((pp - 1) * Q)
        if ((out$thetaxposterior[pp, jj, ss] == 0.0) | (out$thetazposterior[qq, jj, ss] == 0)) {
          bipostxz[pp, qq, jj, ss] <- NA
        } else {
          bipostxz[pp, qq, jj, ss] <- out$bxzpost[jj, hh, ss]
        }
      }
    }
  }
  # ll <- ll + 1
}

bixzmean <- array(0, dim = c(P, Q, J))
bixzincluded <- array(0, dim = c(P, Q, J))
quant <- matrix(0, nrow = (P * Q * J), ncol = 3)
ll <- 1
for (jj in 1:J) {
  for (qq in 1:Q) {
    for (pp in 1:P) {
      quant[ll, ] <- quantile(bipostxz[pp, qq, jj, ], probs = c(0.05, .5, 0.95), na.rm = T)
      bixzmean[pp, qq, jj] <- mean(bipostxz[pp, qq, jj, ], na.rm = T)
      bixzincluded[pp, qq, jj] <- as.numeric(!dplyr::between(0, quant[ll, 1], quant[ll, 3]))
      ll <- ll + 1
    }
  }
}

bixzincluded[is.na(bixzincluded)] <- 0
col <- rep("orange", length(c(bixzincluded)))
col[which(c(bixzincluded) != 0)] <- "blue"
# col <- c("orange", "blue")
mycol <- c(col[c(bixzincluded)])
plotrix::plotCI(1:(P * Q * J), quant[, 2],
  li = quant[, 1], ui = quant[, 3],
  scol = col, xlab = "Interaction Index",
  ylab = expression(b[xz]), xaxt = "n", ylim = c(-2.0, 2.0)
)
axis(1, at = seq(0, (P * Q * J), by = 100), las = 2)
abline(h = 0, col = "red", lty = 2)
plot(c(bixzincluded),
  type = "h", ylim = c(0, 1), xlab = expression(paste("index for ", b[xz])),
  ylab = expression(paste("Inclusion ", b[xz]))
)

# interazioni xx

betaincluded_l <- round(apply(simplify2array(out$xi0posterior) == 1, c(1, 2, 3), mean), 5)

plot(betaincluded_l,
  type = "h", ylim = c(0, 1), xlab = expression(paste("index for ", alpha^{
    0
  })),
  ylab = expression(paste("PPI ", alpha^{
    0
  }))
)

betaincluded_nl <- round(apply(simplify2array(out$xistarposterior) == 1, c(1, 2, 3), mean), 5)

plot(betaincluded_nl,
  type = "h", ylim = c(0, 1), xlab = expression(paste("index for ", alpha, "*")),
  ylab = expression(paste("PPI ", alpha, "*"))
)

# interazioni zz
bipostzz <- array(0, dim = c(Q, Q, J, Eff))

for (ss in 1:Eff) {
  for (jj in 1:J) {
    for (qq in 1:Q) {
      for (pp in 1:Q) {
        hh <- qq + ((pp - 1) * Q)
        if ((out$thetazposterior[pp, jj, ss] == 0.0) | (out$thetazposterior[qq, jj, ss] == 0)) {
          bipostzz[pp, qq, jj, ss] <- NA
        } else {
          bipostzz[pp, qq, jj, ss] <- out$bzzposterior[jj, hh, ss]
        }
      }
    }
  }
  # ll <- ll + 1
}

bizzmean <- array(0, dim = c(Q, Q, J))
bizzincluded <- array(0, dim = c(Q, Q, J))
quant2 <- matrix(0, nrow = (Q * Q * J), ncol = 3)
ll <- 1
for (jj in 1:J) {
  for (qq in 1:Q) {
    for (pp in 1:Q) {
      quant2[ll, ] <- quantile(bipostzz[pp, qq, jj, ], probs = c(0.05, .5, 0.95), na.rm = T)
      bizzmean[pp, qq, jj] <- mean(bipostzz[pp, qq, jj, ], na.rm = T)
      bizzincluded[pp, qq, jj] <- as.numeric(!dplyr::between(0, quant2[ll, 1], quant2[ll, 3]))
      ll <- ll + 1
    }
  }
}

bizzincluded[is.na(bizzincluded)] <- 0
col <- rep("orange", length(c(bizzincluded)))
col[which(c(bizzincluded) != 0)] <- "blue"
# col <- c("orange", "blue")
mycol <- c(col[c(bizzincluded)])
plotrix::plotCI(1:(Q * Q * J), quant2[, 2],
  li = quant2[, 1], ui = quant2[, 3],
  scol = col, xlab = "Interaction Index",
  ylab = expression(b[zz]), xaxt = "n", ylim = c(-2.0, 2.0)
)
axis(1, at = seq(0, (Q * Q * J), by = 100), las = 2)
abline(h = 0, col = "red", lty = 2)
plot(c(bizzincluded),
  type = "h", ylim = c(0, 1), xlab = expression(paste("index for ", b[zz])),
  ylab = expression(paste("Inclusion ", b[zz]))
)