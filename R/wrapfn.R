#' an R wrapper to Rcpp function sampler
#'
#' MCMC sampler to perform bayesian variable selection, hierarchical
#' Dirichlet-Multinomial regression
#'
#' @param YY (count matrix) rows n samples, columns: J species
#' @param XX (design matrix) n x P (continuous variables)
#' @param ZZ (design matrix) n * Q (discrete variables)
#' @param Niter number of MCMC iteration
#' @param burn number of MCMC iteration discarded due to burn-in period
#' @param thin thinning, retain one iteration each \code{thin}
#' @param upsv flag to update slab variance (default is \code{FALSE})
#' @param hereditariety hereditariety assumption on interactions. If \code{2}
#' strong heredity is assumed, if \code{1} weak heredity is assumed, if 
#' \code{0} no assumption on heredity condition is made. 
#' Default is \code{TRUE}
#' @param conlnl it is a constrain imposed on nonlinear continuous continuous 
#' interactions. Default is \code{TRUE}. A nonlinear continuous continuous interaction
#' may be included if and only if the correspondent linear is included.
#' @param randint random intercept is considered for different observation from
#' the same patient. Default is \code{TRUE} (per ora FALSE)
#' @param hyper_theta_x vector of dimension two. The first element is the a 
#' parameter of Beta function after eq 12 manuscript and default is 
#' \code{a = 0.02}. The second element is the variance of the slab. It is not 
#' updated and default is \code{tau_slab = 10}.
#' @param hyper_theta_z vector of dimension two. The first element is the a 
#' parameter of Beta function after eq 12 manuscript and default is 
#' \code{a = 0.2}. The second element is the variance of the slab. It is not 
#' updated and default is \code{tau_slab = 10}.
#' @param penmig_lin vector of length 5. 
#' These are the hyperparameters for 
#' peNMIG on linear continuous continuous interactions. The first two are
#' \code{a_tau, b_tau} that are the parameters for Inverse Gamma on 
#' \code{alpha}'s variance. The third is \code{v_0}. The remaining are 
#' \code{a_omega, b_omega}.
#' 
#' Default is \code{penmig_lin = c(5.0, 25.0, 0.0025, 1.0, 1.0)}
#' @param penmig_nl vector of length 5. 
#' These are the hyperparameters for 
#' peNMIG on nonlinear continuous continuous interactions. The first two are
#' \code{a_tau, b_tau} that are the parameters for Inverse Gamma on 
#' \code{alpha}'s variance. The third is \code{v_0}. The remaining are 
#' \code{a_omega, b_omega}. 
#' 
#' Default is \code{penmig_lin = c(5.0, 25.0, 0.0025, 1.0, 1.0)}
#' @param tx vector of length 2. Hyperparameters for Uniform prior on 
#' first threshold function's parameter. Default is \code{tx = c(0, 0.25)}
#' @param tz vector of length 2. Hyperparameters for Uniform prior on 
#' second threshold function's parameter. Default is \code{tz = c(0, 0.25)}
#' @param prior_int vector of length 2. Hyperparameters for Normal prior on model's intercept. 
#' Default is \code{prior_int = c(0, 10)}
#' @param init whether main effects coefficients shold be initialized with marginal correlations. Default is \code{TRUE}
#' @param init_fd false discovery rate adopted as threshold if \code{init = TRUE}

#' @section Details
#' math stuff to be inserted here
 
#' @return thetaxposterior 
#' samples from posterior distribution of \eqn{ \theta_x }
#' @return thetazposterior samples from posterior distribution of \eqn{ \theta_z }
#' @return muposterior samples from posterior distribution of \eqn{ \mu }
#' @return thrxposterior samples from posterior distribution of \eqn{t_x}
#' @return thrzposterior samples from posterior distribution of \eqn{t_z}
#' @return loglinpred log linear predictor at each MCMC iteration
#' @return loglik loglikelihood at each MCMC iteration
#' @return alpha0posterior samples from posterior distribution of \eqn{\alpha^0}
#' @return alphastarposterior samples from posterior distribution of \eqn{\alpha^{\star}}
#' @return eta0posterior samples from posterior distribution of \eqn{\eta^{0}}
#' @return etastarposterior samples from posterior distribution of \eqn{\eta^{\star}}
#' @return xi0posterior samples from posterior distribution of \eqn{\dot{\xi}^{0}}
#' @return xistarposterior samples from posterior distribution of \eqn{\dot{\xi}^{\star}}
#' @return bxzposterior samples from posterior distribution of \eqn{b^{xz}}
#' @return bzzposterior samples from posterior distribution of \eqn{b^{zz}}
#' @return betaxz samples from posterior distribution of \eqn{\beta(x,z)}
#' @return betaz samples from posterior distribution of \eqn{\beta(z)}
#' @return rint samples from posterior distribution of \eqn{\iota}
#' @return varrint samples from posterior distribution of \eqn{\sigma_{\iota}}
 
#' @export
 
ssdm <-  function(
  YY, XX, ZZ,  
  Niter, burn, thin, 
  upsv = FALSE, hereditariety = 2, conlnl = TRUE, randint = TRUE,
  hyper_theta_x = c(0.05, 10), hyper_theta_z = c(0.05, 10),
  penmig_lin = c(5.0, 25.0, 0.00025, 1.0, 1.0),
  penmig_nl = c(5.0, 25.0, 0.00025, 1.0, 1.0),
  tx = c(0, 1.0), tz = c(0, 1.0),
  prior_int = c(0, 10),
  init = TRUE,
  init_fd = 0.1){
  
  # data dimensions
  numero_categorie = ncol(YY)
  numero_oss = nrow(XX)
  numero_var_cont = ncol(XX)
  numero_var_disc = ncol(ZZ)
  
  # argument checking
  if(hyper_theta_x[1] >= 1.0 | hyper_theta_z[1] >= 1.0){
    warning("Beta-Bernoulli prior is symmetric or left skewed: 
            right skewed enforces sparsity")
    }
  #if(penmig_lin[4] >= penmig_lin[5] | penmig_nl[4] >= penmig_nl[5]){warning("Beta-Bernoulli prior for peNMIG is symmetric or left skewed: right skewed enforces sparsity")}
  if(all(XX[,1] == 1)){
    stop("XX appears to have an intercept column")
    }
  if(all(ZZ[,1] == 1)){
    stop("ZZ appears to have an intercept column")
    }
  if((numero_oss != nrow(YY)) | (numero_oss != nrow(ZZ)) | (nrow(YY) != nrow(ZZ))){
    stop("dimensions of design matrix and YY do not match")
    }
  
  if(burn > Niter){
    stop("burnin greater than the number of iterations")
    }
  if(any(hyper_theta_x[2] < 0, hyper_theta_z[2] < 0, 
         hyper_theta_x[1] < 0, hyper_theta_z[1] < 0)){
    stop("check hyperparameter values (hyper_theta_x or hyper_theta_z)")
  }
  if(any(penmig_lin[1] < 0, penmig_lin[2] < 0, penmig_lin[3] < 0, 
         penmig_lin[4] < 0, penmig_lin[5] < 0)){
    stop("check hyperparameter values (penmig_lin)")
  }
  if(any(penmig_nl[1] < 0, penmig_nl[2] < 0, penmig_nl[3] < 0, 
         penmig_nl[4] < 0, penmig_nl[5] < 0)){
    stop("check hyperparameter values (penmig_nl)")
  }
  if((any(tx[1] < 0, tx[2] < 0, tz[1] < 0, tz[2] < 0)) | 
     (any(tx[1] > tx[2], tz[1] > tz[2]))){
    stop("check thresholds hyperparameter values")
  }
  if(prior_int[2] <= 0){
    stop("intercept variance value is not valid")
  }
  if((hereditariety != 2) & (hereditariety != 1) & (hereditariety != 0)){
    stop("option for hereditariety is not valid")
  }
  
  X <- XX
  basis <- DesignMatrixAlt(X)
  Xilp <- basis$Xlin
  Xinlp <- basis$Xnlin
  dj <- basis$dj
  
  if(randint == FALSE){
    if(init == TRUE){
      Jc = ncol(YY)
      Qc = ncol(XX)
      
      cormat = matrix(0, Jc, Qc)
      pmat = matrix(0, Jc, Qc)
      
      for(rr in 1:Jc){
        for(cc in 1:Qc){
          pmat[rr, cc] = stats::cor.test(XX[, cc], YY[, rr], method = "spearman",
                                         exact = F)$p.value
          cormat[rr, cc] = stats::cor(XX[, cc], YY[, rr], method = "spearman")
        }
      }
      # defaults to 0.2 false discovery rate
      pm = matrix((stats::p.adjust(c(pmat), method = "fdr") <= init_fd) + 0, Jc, Qc)
      betmat = cormat * pm
      betmat = t(betmat)
      
      Qc = ncol(ZZ)
      
      cormat = matrix(0, Jc, Qc)
      pmat = matrix(0, Jc, Qc)
      
      for(rr in 1:Jc){
        for(cc in 1:Qc){
          pmat[rr, cc] = stats::cor.test(ZZ[, cc], YY[, rr], method = "spearman",
                                         exact = F)$p.value
          cormat[rr, cc] = stats::cor(ZZ[, cc], YY[, rr], method = "spearman")
        }
      }
      # defaults to 0.2 false discovery rate
      pm = matrix((stats::p.adjust(c(pmat), method = "fdr") <= init_fd) + 0, Jc, Qc)
      betmat2 = cormat * pm
      betmat2 = t(betmat2)
    } else{
      Jc = ncol(YY)
      Qc = ncol(XX)
      betmat = rep(0, Jc*Qc)
      
      Qc = ncol(ZZ)
      betmat2 = rep(0, Jc*Qc)
    }
    
    output <- sampler(YY = YY, XX = XX, ZZ = ZZ, XXl = Xilp, XXT = Xinlp, dj = dj, 
                      Niter = Niter, burn = burn, thin = thin, 
                      hyper_theta_x = hyper_theta_x, hyper_theta_z = hyper_theta_z, 
                      penmig_lin = penmig_lin, penmig_nl = penmig_nl, tx = tx, 
                      tz = tz, prior_int = prior_int, upsv = upsv, 
                      hereditariety = hereditariety, conlnl = conlnl,
                      theta_init = as.vector(betmat), theta_init2 = as.vector(betmat2))
  } else {
    k <- 0
    grpind <- c()
    for(i in 1:length(table(rownames(X)))){
      grpind <- c(grpind, rep(k, table(rownames(X))[i]))
      k <- k + 1
    }
    
    if(init == TRUE){
      Jc = ncol(YY)
      Qc = ncol(XX)
      
      cormat = matrix(0, Jc, Qc)
      pmat = matrix(0, Jc, Qc)
      
      for(rr in 1:Jc){
        for(cc in 1:Qc){
          pmat[rr, cc] = stats::cor.test(XX[, cc], YY[, rr], method = "spearman",
                                         exact = F)$p.value
          cormat[rr, cc] = stats::cor(XX[, cc], YY[, rr], method = "spearman")
        }
      }
      # defaults to 0.2 false discovery rate
      pm = matrix((stats::p.adjust(c(pmat), method = "fdr") <= init_fd) + 0, Jc, Qc)
      betmat = cormat * pm
      betmat = t(betmat)
      
      Qc = ncol(ZZ)
      
      cormat = matrix(0, Jc, Qc)
      pmat = matrix(0, Jc, Qc)
      
      for(rr in 1:Jc){
        for(cc in 1:Qc){
          pmat[rr, cc] = stats::cor.test(ZZ[, cc], YY[, rr], method = "spearman",
                                         exact = F)$p.value
          cormat[rr, cc] = stats::cor(ZZ[, cc], YY[, rr], method = "spearman")
        }
      }
      # defaults to 0.2 false discovery rate
      pm = matrix((stats::p.adjust(c(pmat), method = "fdr") <= init_fd) + 0, Jc, Qc)
      betmat2 = cormat * pm
      betmat2 = t(betmat2)
    } else {
      Jc = ncol(YY)
      Qc = ncol(XX)
      betmat = runif(Jc*Qc, -1, 1)
      
      Qc = ncol(ZZ)
      betmat2 = runif(Jc*Qc, -1, 1)
    }
    
    output <- sampler_randint(YY = YY, XX = XX, ZZ = ZZ, XXl = Xilp, XXT = Xinlp, dj = dj, 
                      Niter = Niter, burn = burn, thin = thin, 
                      hyper_theta_x = hyper_theta_x, hyper_theta_z = hyper_theta_z, 
                      penmig_lin = penmig_lin, penmig_nl = penmig_nl, tx = tx, 
                      tz = tz, prior_int = prior_int, upsv = upsv, 
                      hereditariety = hereditariety, conlnl = conlnl, 
                      grplabel = grpind, ngroup = grpind[length(grpind)]+1, 
                      theta_init = as.vector(betmat), theta_init2 = as.vector(betmat2))
  }
  return(output)
}
