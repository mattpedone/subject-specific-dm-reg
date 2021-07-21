rm(list=ls())
library(gdata)
#source("R/datagenerator.R")
source("demo/datagen_comp.R")

n <- 100
P <- 5
Q <- 5
K <- 3
add <- floor(runif(K, 0, 1000))
scenarij <- c(10, 50, 100, 150)

sim <- matrix(0, 8, 3)
idx <- 1

lab <- c("tpr", "fpr", "mcc")
calculate_statistics <- function(tbl) {
  tbl <- tbl/1000
  tpr <- tbl[2,2]/sum(tbl[,2])
  fpr <- tbl[2,1]/sum(tbl[,1])
  num <- (tbl[2,2] * tbl[1,1]) - (tbl[2,1] * tbl[1,2])
  den <- (sum(tbl[2,])) * (sum(tbl[,2])) * (sum(tbl[,1])) * (sum(tbl[1,]))
  mcc <- num/sqrt(den)
  res <- c(tpr, fpr, mcc)
  return(res)
}

##### Raffaele et al ##### 
source("demo/raffaele/code/wrapper.R")
source(file.path("demo/raffaele/code", "helper_functions.R"))
system("cd demo/raffaele/code; make")
executable_location = "demo/raffaele/code/dmbvs.x"
save_prefix = "simulation"

res_pp_R <- matrix(0, K, 3)
colnames(res_pp_R) <- lab

output_R <- list()
for(scen in 1:length(scenarij)){
# prepare and check data
for(k in 1:K){
  cat("replica", k, "raff")
  
  SimData <- gendata_simple_confr(scenarij[scen], n, P, Q, .75, 1.5, 3, 2, 2, c(2, 2, 3, 3),
                                  corrx = 0.4, corrz = 0.4, nli = T, WH = T, 
                                  seed = 121 + add[k])
  YY = as.matrix(SimData$Y)
  XX = scale(as.matrix(SimData$DesMat), center = T, scale = T)
  
  #GG = 150001L; thin = 10L; burn = 100001L; # good defaults, in this case
  GG = 101L; thin = 1L; burn = 51L; # fast for testing
  bb_alpha = 0.02; bb_beta = 2 - bb_alpha
  proposal_alpha = 0.5; proposal_beta = 0.5
  slab_variance = 10; intercept_variance = 10
  
  
  # run the algorithm
  results = dmbvs(XX = XX, YY = YY, intercept_variance = intercept_variance,
                  slab_variance = slab_variance, bb_alpha = bb_alpha,
                  bb_beta = bb_beta, GG = GG, thin = thin, burn = burn,
                  init_beta = "warmstart", init_alpha = "warmstart",
                  proposal_alpha = proposal_alpha, proposal_beta = proposal_beta,
                  exec = executable_location, selection_type = "ss",
                  #output_location = ".")
                  output_location = "demo/simcompraff")
  
  output_R[[k]] <- results
  
  params = data.frame(GG, burn, thin, intercept_variance,
                      slab_variance, bb_alpha, bb_beta,
                      proposal_alpha, proposal_beta)
  
  pp_pred <- colMeans((results$beta != 0) + 0)
  pp_pred <- as.numeric(pp_pred > .5)
  pp_pred <- factor(pp_pred, levels = 0:1)
  
  pp_act <- as.vector(SimData$thetaX_exp)
  pp_act[which(pp_act != 0)] <- 1
  
  tbl <- table(pp_pred, pp_act)
  res_pp_R[k,] <- calculate_statistics(tbl)
}

sim[idx,] <- apply(res_pp_R, 2, mean)
idx <- idx+1
sim[idx,] <- sqrt(apply(res_pp_R, 2, var))
idx <- idx+1
}