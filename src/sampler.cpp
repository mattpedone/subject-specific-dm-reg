#include <RcppArmadillo.h>
#include <math.h>
#include <Rcpp.h>
#include <cstdlib>

//[[Rcpp::depends(RcppArmadillo)]]

double inner(arma::vec x,
             arma::vec y) {
  arma::mat ip = x.t() * y ;
  return(ip(0)) ;
}

double threshold(double val, double soglia){
  if(abs(val) < soglia)
    val = 0.0;
  return val;
}

arma::mat threshold_mat(arma::mat M, double t, double sub = 0.0){
  const int rows = M.n_rows;
  const int cols = M.n_cols;
  arma::mat Mt(rows, cols);
  
  for(int i = 0; i< rows; i++){
    for(int j = 0; j< cols; j++){
      if((std::abs(M(i, j)/1.0)) < t){
        Mt(i, j) = sub;
      } else {
        Mt(i, j) = M(i,j);
      }
    }
  }
  return Mt;
}

arma::vec myrowSums(arma::mat X){
  
  //ritorna le somme sulle righe
  int nRows = X.n_rows;
  arma::vec out(nRows);
  
  for(int i = 0; i < nRows; i++){
    out(i) = sum(X.row(i));
  }
  return out;
}

arma::colvec calculate_gamma_ri(arma::mat XX, arma::mat ZZ, arma::mat Xl,
                                arma::mat Xnl, arma::vec &alpha, arma::vec rand_int,
                                arma::vec theta,
                                arma::vec theta2, arma::mat bXZ, arma::mat bXXl,
                                arma::mat bXXnl,
                                double soglia, arma::mat bZZ, double soglia2,
                                int jj){
  
  //calcola la j-esima colonna del predittore lineare (subject-specific)
  int PP = XX.n_cols;
  int QQ = ZZ.n_cols;
  int n = XX.n_rows;
  int JJ = alpha.size();
  int pp, qq, hh;
  
  arma::vec thetaX(PP);
  arma::vec thetaZ(QQ);
  
  for(pp = 0; pp < PP; pp++){
    hh = jj + pp * JJ;
    thetaX(pp) = theta(hh);
  }
  
  for(qq = 0; qq < QQ; qq++){
    hh = jj + qq * JJ;
    thetaZ(qq) = theta2(hh);
  }
  
  arma::colvec ZETAj(n);
  arma::mat temp1(n, PP);
  arma::mat temp2(n, QQ);
  
  arma::mat v1(1, n);
  v1.fill(1.0);
  
  arma::mat matricesoglia(n, PP);
  matricesoglia = trans(bXZ * ZZ.t()) + trans(bXXl * Xl.t()) + trans(bXXnl * Xnl.t());
  temp1 = trans(thetaX * v1) + threshold_mat(matricesoglia, soglia);
  temp2 = trans(thetaZ * v1) + threshold_mat(trans(bZZ * ZZ.t()), soglia2);
  ZETAj = trans(alpha(jj)*v1) + rand_int + myrowSums(XX % temp1) + myrowSums(ZZ % temp2);
  
  return(ZETAj);
}

arma::rowvec uplam(arma::rowvec lambda, arma::vec beta, int pp, double tau=1,
                   double sigma=1){
  
  //adaptated from https://rdrr.io/cran/horseshoe/src/R/horseshoe.R
  int  qq, hh;
  
  int QQ = beta.size();
  arma::vec etas(QQ), etas2(QQ), upsi(QQ), ub(QQ), temp(QQ), Fub(QQ), up(QQ);
  
  arma::vec work_lambda(QQ);
  for(qq = 0; qq < QQ; qq ++){
    hh = qq + pp * QQ;
    work_lambda(qq) = lambda(hh);
  }
  
  for(qq = 0; qq < QQ; qq ++){
    etas(qq) = 1.0/pow(work_lambda(qq), 2.0);
    upsi(qq) = ::Rf_runif(0, 1.0/(1+etas(qq)));
    ub(qq) = (1-upsi(qq))/upsi(qq);
    temp(qq) = pow(beta(qq), 2.0)/(2.0*sigma*pow(tau, 2.0));/*tempps*/
    Fub(qq) = 1 - exp(-temp(qq)*ub(qq));/*Fub*/
    if(Fub(qq) < .0001){
      Fub(qq) = .0001;
    }
    up(qq) = ::Rf_runif(0.0, Fub(qq));/*up*/
    etas2(qq) = -log(1-up(qq))/temp(qq);
    work_lambda(qq) = 1.0/pow(etas2(qq), .5);
    
    if(work_lambda(qq) < .0001){
      work_lambda(qq) = .0001;
    }
  }
  
  for(qq = 0; qq < QQ; qq ++){
    hh = qq + pp * QQ;
    lambda(hh) = work_lambda(qq);
  }
  
  return lambda;
}

double uptau(arma::vec lambda, arma::vec beta, int qq, int PP, double tau,
             bool truncated = false){
  
  double sigma = 1;
  double tempt;
  arma::vec lambda1(PP);
  int hh;
  for(int pp = 0; pp < PP;pp++){
    hh = pp + qq * PP;
    lambda1(pp) = lambda(hh);
  }
  if(truncated == true){
    tempt=arma::sum(pow(beta/lambda1,2))/(2*sigma);
    double et = 1.0/pow(tau, 2.0);
    double utau= ::Rf_runif(0.0, 1/(1+et));
    double ubt = std::min((1-utau)/utau, pow(PP, 2.0));
    double Fubt1 = ::Rf_pgamma(1.0,(((double)PP)+1)/2.0, (1/tempt), 1, 0);
    double Fubt = ::Rf_pgamma(ubt,(((double)PP)+1)/2.0, (1/tempt), 1, 0);
    if(Fubt1>Fubt){
      Fubt1 = std::max(Fubt1 - (Fubt1-Fubt+.0001), 0.0);
    }
    double ut = ::Rf_runif(Fubt1,Fubt);
    
    et = ::Rf_qgamma(ut,(((double)PP)+1)/2.0,(1/tempt), 1, 0);
    if(et < .0001){
      et = .0001;
    }
    tau = 1/sqrt(et);
    if(tau < .0001){
      tau = .0001;
    }
    if(tau > 1.0001){
      tau = 1.0001;
    }
    
    if(std::isnan(tau)){
      tau = 0.1;
    }
    
  } else {
    tempt=arma::sum(pow(beta/lambda1,2))/(2*sigma);
    double et = 1.0/pow(tau, 2);
    double utau= ::Rf_runif(0.0, 1/(1+et));
    double ubt = (1-utau)/utau;
    double Fubt = ::Rf_pgamma(ubt,(((double)PP)+1)/2.0, (1/tempt), 1, 0);
    Fubt = std::max(Fubt, .0001);
    double ut = ::Rf_runif(0.0,Fubt);
    et = ::Rf_qgamma(ut,(((double)PP)+1)/2.0,(1/tempt), 1, 0);
    if(et < .0001){
      et = .0001;
    }
    tau = 1/sqrt(et);
    if(tau < .0001){
      tau = .0001;
    }
  }
  
  return tau;
  
}//closes uptau() function

arma::colvec calculate_gamma(arma::mat XX, arma::mat ZZ, arma::mat Xl,
                             arma::mat Xnl, arma::vec &alpha, arma::vec theta,
                             arma::vec theta2, arma::mat bXZ, arma::mat bXXl,
                             arma::mat bXXnl,
                             double soglia, arma::mat bZZ, double soglia2,
                             int jj){
  
  //calcola la j-esima colonna del predittore lineare (subject-specific)
  int PP = XX.n_cols;
  int QQ = ZZ.n_cols;
  int n = XX.n_rows;
  int JJ = alpha.size();
  int pp, qq, hh;
  
  arma::vec thetaX(PP);
  arma::vec thetaZ(QQ);
  
  for(pp = 0; pp < PP; pp++){
    hh = jj + pp * JJ;
    thetaX(pp) = theta(hh);
  }
  
  for(qq = 0; qq < QQ; qq++){
    hh = jj + qq * JJ;
    thetaZ(qq) = theta2(hh);
  }
  
  arma::colvec ZETAj(n);
  arma::mat temp1(n, PP);
  arma::mat temp2(n, QQ);
  
  arma::mat v1(1, n);
  v1.fill(1.0);
  
  arma::mat matricesoglia(n, PP);
  matricesoglia = trans(bXZ * ZZ.t()) + trans(bXXl * Xl.t()) + trans(bXXnl * Xnl.t());
  temp1 = trans(thetaX * v1) + threshold_mat(matricesoglia, soglia);
  temp2 = trans(thetaZ * v1) + threshold_mat(trans(bZZ * ZZ.t()), soglia2);
  ZETAj = trans(alpha(jj)*v1) + myrowSums(XX % temp1) + myrowSums(ZZ % temp2);
  
  return(ZETAj);
}

double adap_prop(double curr_var, int PP, int JJ, int itr){
  
  //adaptive proposal for the intercepts
  double dd = (double)JJ * (double)PP * (double)itr;
  
  // ensures bounded convergence
  int safeguard = ::Rf_rbinom(1.0, .05);
  double usually = pow(2.38, 2.0) * curr_var/dd;
  double unusually = pow(0.1, 2.0)/dd;
  
  double prop_var = (1 - safeguard) * ::Rf_rnorm(0.0, pow(usually, 0.5)) +
    safeguard * ::Rf_rnorm(0.0, pow(unusually, 0.5));
  
  return(prop_var);
}

// valuto la prior beta mixed binomial
double lprior_bb(double betajk, int sjk, double sig_bejk, double mu_bejk,
                 double aa_hp, double bb_hp){
  //https://github.com/duncanwadsworth/dmbvs
  
  // calculate additional beta factor
  double aa = sjk + aa_hp;
  double bb = 1 - sjk + bb_hp;
  double lbeta_factor = ::Rf_lbeta(aa, bb) - ::Rf_lbeta(aa_hp, bb_hp);
  
  // piece-wise function
  if(sjk == 0){
    return(log(1.0 - exp(lbeta_factor)));
  }
  else{
    return(lbeta_factor - 0.5 * log(2.0 * M_PI * sig_bejk) -
           1.0/(2.0 * sig_bejk) * pow(betajk - mu_bejk, 2.0));
  }
}

double online_mean(int iteration, double last_mean, double curr_obs){
  //https://github.com/duncanwadsworth/dmbvs
  
  // don't fall off by one
  int n_vals = iteration++;
  
  // assume that iteration > 0 since there's a burnin proposal
  double curr_mean = last_mean + (curr_obs - last_mean)/(double)n_vals;
  
  return(curr_mean);
}

double online_var(int iteration, double last_mean, double last_var,
                  double curr_mean, double curr_obs){
  //https://github.com/duncanwadsworth/dmbvs
  
  // note that iteration == n - 1 and that (n - 1) * last_var == last_ss
  // assume that iteration > 0 since there's a burnin proposal
  double curr_ss = (double)iteration * last_var + (curr_obs - last_mean) *
    (curr_obs - curr_mean);
  
  return(curr_ss/(double)iteration);
}

void update_theta(arma::mat XX, arma::mat ZZ, arma::mat SS, arma::mat &loggamma,
                  arma::vec &alpha, arma::vec &theta_temp,
                  arma::vec &theta_temp2, arma::vec theta, arma::vec theta2,
                  arma::mat &bi_xz, double soglia, arma::mat &bi_zz,
                  double soglia2, arma::vec &inclusion_indicator,
                  arma::vec &inclusion_indicator2, arma::vec &prop_per_theta,
                  arma::vec &prop_per_theta2, arma::mat mu_theta,
                  arma::mat mu_theta2, arma::mat sigma_theta,
                  arma::mat sigma_theta2, double aa_hp, double bb_hp,
                  double aa_hp2, double bb_hp2, int jj, bool primo = true){
  
  int PP = XX.n_cols;
  int QQ = ZZ.n_cols;
  int JJ = SS.n_cols;
  int n = XX.n_rows;
  
  double sigma_prop;
  arma::vec loggamma_p(n);
  
  int hh;
  
  double lnu, ln_acp;
  double log_full_theta, log_full_theta_p;
  
  double u;
  u = arma::randu();
  
  if(primo == true){
    
    int pp;
    arma::vec theta_p = theta_temp;
    
    //scelgo una covariata
    u = u*(PP/1.0);
    pp = Rcpp::as<int>(Rcpp::wrap(std::floor(u)));
    
    hh = jj + pp * JJ;
    
    if(inclusion_indicator(hh)==1){
      
      //calcolo il denominatore del MH ratio "allo stato corrente"
      log_full_theta = 0.0;
      log_full_theta = log_full_theta - arma::sum(lgamma(exp(loggamma.col(jj))))
        + arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)))
        + lprior_bb(theta_temp(hh), inclusion_indicator(hh), sigma_theta(pp, jj), mu_theta(pp, jj), aa_hp, bb_hp);
        
        //campiono da adaptive proposal
        sigma_prop = prop_per_theta(hh);
        theta_p(hh) = theta_temp(hh)+ adap_prop(sigma_prop, PP, JJ, 1);
        
        //calcolo predittore lineare con valore proposti
        loggamma_p = loggamma.col(jj) - theta_temp(hh)*XX.col(pp) + theta_p(hh)*XX.col(pp);
        
        //calcolo il numeratore del MH ratio "allo stato \star"
        log_full_theta_p = 0.0;
        log_full_theta_p = log_full_theta_p - arma::sum(lgamma(exp(loggamma_p)))
          + arma::sum(exp(loggamma_p)%log(SS.col(jj)))
          + lprior_bb(theta_p(hh), inclusion_indicator(hh), sigma_theta(pp, jj), mu_theta(pp, jj), aa_hp, bb_hp);
          
          ln_acp = log_full_theta_p - log_full_theta;
          lnu = log(arma::randu());
          
          //MH ratio
          if(lnu < ln_acp){
            theta_temp(hh) = theta_p(hh);
            loggamma.col(jj)= loggamma_p;
          }
    }//chiude if dell indicatore
    
  } else {//SECONDA PARTE
    
    int qq;
    arma::vec theta_p = theta_temp2;
    
    u = u*(QQ/1.0);
    qq = Rcpp::as<int>(Rcpp::wrap(std::floor(u)));
    
    hh = jj + qq * JJ;
    
    if(inclusion_indicator2(hh)==1){
      
      log_full_theta = 0.0;
      
      log_full_theta = log_full_theta - arma::sum(lgamma(exp(loggamma.col(jj))))
        + arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)))
        + lprior_bb(theta_temp2(hh), inclusion_indicator2(hh), sigma_theta2(qq, jj), mu_theta2(qq, jj), aa_hp2, bb_hp2);
        
        sigma_prop = prop_per_theta2(hh);
        theta_p(hh) = theta_temp2(hh)+ adap_prop(sigma_prop, QQ, JJ, 1);
        
        loggamma_p = loggamma.col(jj) - theta_temp2(hh)*ZZ.col(qq) + theta_p(hh)*ZZ.col(qq);
        
        log_full_theta_p = 0.0;
        log_full_theta_p = log_full_theta_p - arma::sum(lgamma(exp(loggamma_p)))
          + arma::sum(exp(loggamma_p)%log(SS.col(jj)))
          + lprior_bb(theta_p(hh), inclusion_indicator2(hh), sigma_theta2(qq, jj), mu_theta2(qq, jj), aa_hp2, bb_hp2);
          
          ln_acp =log_full_theta_p -log_full_theta;
          
          lnu = log(arma::randu());
          
          if(lnu < ln_acp){
            theta_temp2(hh) = theta_p(hh);
            loggamma.col(jj)= loggamma_p;
          }
    }//chiude if dell indicatore
  }//chiude seconda parte
}//chiude la function

double update_soglia(arma::mat XX, arma::mat ZZ, arma::mat Xl, arma::mat Xnl, arma::mat SS,
                     arma::mat &loggamma, arma::vec &alpha,
                     arma::vec theta, arma::vec theta2, arma::mat &bi_xz, arma::cube bXXl,
                     arma::cube bXXnl, double soglia, arma::mat &bi_zz, double soglia2,
                     arma::vec &acc_soglia_flag, double lb, double ub,
                     int pos = 0){
  if(pos == 0){
    
    int ii, jj;
    int JJ = SS.n_cols;
    int n = XX.n_rows;
    
    int PP = XX.n_cols;
    int QQ = ZZ.n_cols;
    
    double soglia_p;
    
    double lnu, ln_acp;
    
    arma::mat loggamma_p(n, JJ);
    
    double log_full_soglia, log_full_soglia_p;
    
    log_full_soglia = 0.0;
    for(jj = 0; jj < JJ; jj ++){
      for(ii = 0; ii < n; ii++){
        log_full_soglia = log_full_soglia - lgamma(exp(loggamma(ii, jj)));
        log_full_soglia = log_full_soglia + exp(loggamma(ii, jj))*log(SS(ii, jj));
      }
    }
    
    log_full_soglia = log_full_soglia + ::Rf_dunif(soglia, lb, ub, 1);
    
    soglia_p = soglia + ::Rf_rnorm(0.0, .01);
    if(soglia_p >= ub){
      soglia_p = ub - (soglia_p - ub);
    }
    if(soglia_p <= lb){
      soglia_p = lb + (lb - soglia_p);
    }
    
    arma::mat bXZ(PP, QQ);
    arma::mat bZZ(QQ, QQ);
    
    for(jj = 0; jj < JJ; jj++){
      for(int qq2 = 0; qq2 < QQ; qq2++){
        for(int pp2 = 0; pp2 < PP; pp2++){
          int hh2 = qq2+((pp2)*QQ);
          bXZ(pp2, qq2)=bi_xz(jj, hh2);
        }
        
        for(int pp2 = 0; pp2 < QQ; pp2++){
          int hh2 = qq2 + ((pp2) * QQ);
          bZZ(pp2, qq2) = bi_zz(jj, hh2);
        }
      }
      
      loggamma_p.col(jj)=calculate_gamma(XX, ZZ, Xl, Xnl, alpha, theta, theta2, bXZ,
                     bXXl.slice(jj), bXXnl.slice(jj),soglia_p, bZZ, soglia2, jj);
    }
    
    log_full_soglia_p = 0.0;
    for(jj = 0; jj < JJ; jj ++){
      for(ii = 0; ii <n; ii++){
        log_full_soglia_p = log_full_soglia_p - lgamma(exp(loggamma(ii, jj)));
        log_full_soglia_p = log_full_soglia_p + exp(loggamma(ii, jj))*log(SS(ii, jj));
      }
    }
    
    log_full_soglia_p = log_full_soglia_p + ::Rf_dunif(soglia_p, lb, ub, 1);
    
    
    ln_acp = log_full_soglia_p - log_full_soglia;
    lnu = log(arma::randu());
    
    if(lnu < ln_acp){
      soglia = soglia_p;
      acc_soglia_flag(pos) = 1;
      for(jj = 0; jj < JJ; jj++){
        loggamma.col(jj)= loggamma_p.col(jj);
      }
    }
  } else {
    
    int ii, jj;
    int JJ = SS.n_cols;
    int n = XX.n_rows;
    
    int PP = XX.n_cols;
    int QQ = ZZ.n_cols;
    
    double soglia_p;
    
    double lnu, ln_acp;
    
    arma::mat loggamma_p(n, JJ);
    
    double log_full_soglia, log_full_soglia_p;
    
    log_full_soglia = 0.0;
    for(jj = 0; jj < JJ; jj ++){
      for(ii = 0; ii < n; ii++){
        log_full_soglia = log_full_soglia - lgamma(exp(loggamma(ii, jj)));
        log_full_soglia = log_full_soglia + exp(loggamma(ii, jj))*log(SS(ii, jj));
      }
    }
    
    log_full_soglia = log_full_soglia + ::Rf_dunif(soglia2, lb, ub, 1);
    
    soglia_p = soglia2 + ::Rf_rnorm(0.0, .01);
    if(soglia_p >= ub){
      soglia_p = ub - (soglia_p - ub);
    }
    if(soglia_p <= lb){
      soglia_p = lb + (lb - soglia_p);
    }
    
    arma::mat bXZ(PP, QQ);
    arma::mat bZZ(QQ, QQ);
    
    for(jj = 0; jj < JJ; jj++){
      for(int qq2 = 0; qq2 < QQ; qq2++){
        for(int pp2 = 0; pp2 < PP; pp2++){
          int hh2 = qq2+((pp2)*QQ);
          bXZ(pp2, qq2)=bi_xz(jj, hh2);
        }
        
        for(int pp2 = 0; pp2 < QQ; pp2++){
          int hh2 = qq2 + ((pp2) * QQ);
          bZZ(pp2, qq2) = bi_zz(jj, hh2);
        }
      }
      
      loggamma_p.col(jj)=calculate_gamma(XX, ZZ, Xl, Xnl, alpha, theta, theta2, bXZ,
                     bXXl.slice(jj), bXXnl.slice(jj), soglia, bZZ, soglia_p, jj);
    }
    
    log_full_soglia_p = 0.0;
    for(jj = 0; jj < JJ; jj ++){
      for(ii = 0; ii <n; ii++){
        log_full_soglia_p = log_full_soglia_p - lgamma(exp(loggamma(ii, jj)));
        log_full_soglia_p = log_full_soglia_p + exp(loggamma(ii, jj))*log(SS(ii, jj));
      }
    }
    
    log_full_soglia_p = log_full_soglia_p + ::Rf_dunif(soglia_p, lb, ub, 1);
    
    
    ln_acp = log_full_soglia_p - log_full_soglia;
    lnu = log(arma::randu());
    
    if(lnu < ln_acp){
      soglia = soglia_p;
      acc_soglia_flag(pos) = 1;
      for(jj = 0; jj < JJ; jj++){
        loggamma.col(jj)= loggamma_p.col(jj);
      }
    }
  }
  return soglia;
}

void update_bi(arma::mat XX, arma::mat ZZ, arma::mat XXl, arma::mat XXT,
               arma::mat SS, arma::mat &loggamma,
               arma::mat lambda_hs, arma::vec tau_hs, arma::vec &alpha,
               arma::vec theta, arma::vec theta2, arma::mat &bi_xz,
               arma::mat &bi_zz, arma::mat &bi_p_mat, arma::mat bXXl,
               arma::mat bXXnl, double soglia, double soglia2,
               arma::mat &prop_per_bi, arma::mat &acc_bi_flag,
               double mu_hs, int pp, int qq, int jj, int hered){
  
  int PP = XX.n_cols;
  int QQ = ZZ.n_cols;
  int JJ = SS.n_cols;
  int n = ZZ.n_rows;
  
  arma::vec loggamma_p(n);
  
  int ll, hh;
  
  ll = qq + jj * QQ;//indice per tau
  hh = qq + pp * QQ;//indice per bi e lambda
  
  //indici hereditariety
  int hh_theta = jj + pp * JJ;
  int ll_theta = jj + qq * JJ;
  
  arma::mat bXZ(PP, QQ);
  
  for(int qq2 = 0; qq2 < QQ; qq2++){
    for(int pp2 = 0; pp2 < PP; pp2++){
      int hh2 = qq2 + ((pp2) * QQ);
      bXZ(pp2, qq2) = bi_xz(jj, hh2);
    }
  }
  
  arma::mat bZZ(QQ, QQ);
  
  for(int qq2 = 0; qq2 < QQ; qq2++){
    for(int pp2 = 0; pp2 < QQ; pp2++){
      int hh2 = qq2 + ((pp2) * QQ);
      bZZ(pp2, qq2) = bi_zz(jj, hh2);
    }
  }
  
  //if hereditariety
  if(hered == 2){
    if((theta(hh_theta) == 0) | (theta2(ll_theta) == 0.0)){
      acc_bi_flag(jj, hh) = 0.0;
      bi_p_mat(pp, qq) = 0.0;
      
      loggamma.col(jj) = calculate_gamma(XX, ZZ, XXl, XXT, alpha, theta, theta2, bi_p_mat,
                   bXXl, bXXnl, soglia, bZZ, soglia2, jj);
      
      bi_xz(jj, hh) = bi_p_mat(pp, qq);
    } else {
      
      double lnu, ln_acp;
      
      double log_full_bi, log_full_bi_p;
      
      log_full_bi = 0.0;
      log_full_bi = log_full_bi - arma::sum(lgamma(exp(loggamma.col(jj)))) +
        arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)))
        + ::Rf_dnorm4(bi_xz(jj, hh), mu_hs, lambda_hs(jj, hh)*tau_hs(ll), 1);
      
      double propbipm = prop_per_bi(jj, hh);
      propbipm = adap_prop(propbipm, PP, JJ, QQ);
      
      bi_p_mat(pp, qq) += propbipm;
      
      loggamma.col(jj) = calculate_gamma(XX, ZZ, XXl, XXT, alpha, theta, theta2, bi_p_mat,
                   bXXl, bXXnl, soglia, bZZ, soglia2, jj);
      
      log_full_bi_p = 0.0;
      log_full_bi_p = log_full_bi_p - arma::sum(lgamma(exp(loggamma_p))) +
        arma::sum(exp(loggamma_p)%log(SS.col(jj)))
        + ::Rf_dnorm4(bi_p_mat(pp, qq), mu_hs, lambda_hs(jj, hh)*tau_hs(ll), 1);
      
      ln_acp = log_full_bi_p - log_full_bi;
      
      lnu = log(arma::randu());
      
      if(lnu < ln_acp){
        bi_xz(jj, hh) = bi_p_mat(pp, qq);
        acc_bi_flag(jj, hh) = 1;
        loggamma.col(jj)= loggamma_p;
      } else {
        bi_p_mat(pp, qq) -= propbipm;
      }//chiude ratio
    }
  }
  if(hered == 1){
    if((theta(hh_theta) == 0) & (theta2(ll_theta) == 0.0)){
      acc_bi_flag(jj, hh) = 0.0;
      bi_p_mat(pp, qq) = 0.0;
      
      loggamma.col(jj) = calculate_gamma(XX, ZZ, XXl, XXT, alpha, theta, theta2, bi_p_mat,
                   bXXl, bXXnl, soglia, bZZ, soglia2, jj);
      
      bi_xz(jj, hh) = bi_p_mat(pp, qq);
    } else {
      
      double lnu, ln_acp;
      
      double log_full_bi, log_full_bi_p;
      
      log_full_bi = 0.0;
      log_full_bi = log_full_bi - arma::sum(lgamma(exp(loggamma.col(jj)))) +
        arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)))
        + ::Rf_dnorm4(bi_xz(jj, hh), mu_hs, lambda_hs(jj, hh)*tau_hs(ll), 1);
      
      double propbipm = prop_per_bi(jj, hh);
      propbipm = adap_prop(propbipm, PP, JJ, QQ);
      
      bi_p_mat(pp, qq) += propbipm;
      
      loggamma.col(jj) = calculate_gamma(XX, ZZ, XXl, XXT, alpha, theta, theta2, bi_p_mat,
                   bXXl, bXXnl, soglia, bZZ, soglia2, jj);
      
      log_full_bi_p = 0.0;
      log_full_bi_p = log_full_bi_p - arma::sum(lgamma(exp(loggamma_p))) +
        arma::sum(exp(loggamma_p)%log(SS.col(jj)))
        + ::Rf_dnorm4(bi_p_mat(pp, qq), mu_hs, lambda_hs(jj, hh)*tau_hs(ll), 1);
      
      ln_acp = log_full_bi_p - log_full_bi;
      
      lnu = log(arma::randu());
      
      if(lnu < ln_acp){
        bi_xz(jj, hh) = bi_p_mat(pp, qq);
        acc_bi_flag(jj, hh) = 1;
        loggamma.col(jj)= loggamma_p;
      } else {
        bi_p_mat(pp, qq) -= propbipm;
      }//chiude ratio
    }
  }
  if(hered == 0){
    double lnu, ln_acp;
    
    double log_full_bi, log_full_bi_p;
    
    log_full_bi = 0.0;
    log_full_bi = log_full_bi - arma::sum(lgamma(exp(loggamma.col(jj)))) +
      arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)))
      + ::Rf_dnorm4(bi_xz(jj, hh), mu_hs, lambda_hs(jj, hh)*tau_hs(ll), 1);
    
    double propbipm = prop_per_bi(jj, hh);
    propbipm = adap_prop(propbipm, PP, JJ, QQ);
    
    bi_p_mat(pp, qq) += propbipm;
    
    loggamma.col(jj) = calculate_gamma(XX, ZZ, XXl, XXT, alpha, theta, theta2, bi_p_mat,
                 bXXl, bXXnl, soglia, bZZ, soglia2, jj);
    
    log_full_bi_p = 0.0;
    log_full_bi_p = log_full_bi_p - arma::sum(lgamma(exp(loggamma_p))) +
      arma::sum(exp(loggamma_p)%log(SS.col(jj)))
      + ::Rf_dnorm4(bi_p_mat(pp, qq), mu_hs, lambda_hs(jj, hh)*tau_hs(ll), 1);
    
    ln_acp = log_full_bi_p - log_full_bi;
    
    lnu = log(arma::randu());
    
    if(lnu < ln_acp){
      bi_xz(jj, hh) = bi_p_mat(pp, qq);
      acc_bi_flag(jj, hh) = 1;
      loggamma.col(jj)= loggamma_p;
    } else {
      bi_p_mat(pp, qq) -= propbipm;
    }//chiude ratio
  }
}

void update_bi2(arma::mat XX, arma::mat ZZ, arma::mat SS, arma::mat &loggamma,
                arma::mat lambda_hs2, arma::vec tau_hs2, arma::vec &alpha,
                arma::vec theta, arma::vec theta2, arma::mat &bi_xz,
                arma::mat &bi_zz, arma::mat &bi_p_mat2, double soglia,
                double soglia2, arma::mat &prop_per_bi2, arma::mat &acc_bi_flag2,
                double mu_hs, int pp, int qq, int jj, int hered){
  
  int QQ = ZZ.n_cols;
  int JJ = SS.n_cols;
  int n = ZZ.n_rows;
  
  arma::vec loggamma_p(n);
  
  int ll, hh;
  ll = qq + jj * QQ;//indice per tau
  hh = qq + pp * QQ;//indice per bi e lambda QQ è il numero di colonne
  
  //indici hereditariety
  int hh_theta = jj + pp * JJ;
  int ll_theta = jj + qq * JJ;
  
  arma::mat bZZ(QQ, QQ);
  
  for(int qq2 = 0; qq2 < QQ; qq2++){
    for(int pp2 = 0; pp2 < QQ; pp2++){
      int hh2 = qq2 + ((pp2) * QQ);
      bZZ(pp2, qq2) = bi_zz(jj, hh2);
    }
  }
  
  //if hereditariety
  if(hered == 2){
    if((theta2(hh_theta) == 0) | (theta2(ll_theta) == 0.0)){
      acc_bi_flag2(jj, hh) = 0;
      bi_p_mat2(pp, qq) = 0.0;
      
      loggamma.col(jj) = loggamma.col(jj) - trans(threshold_mat(bZZ(pp, qq) *
        ZZ.col(qq).t(), soglia2)) % ZZ.col(qq) + trans(threshold_mat(bi_p_mat2(pp, qq) *
        ZZ.col(qq).t(), soglia2)) % ZZ.col(qq);
      
      bi_zz(jj, hh) = bi_p_mat2(pp, qq);
    } else {
      double lnu, ln_acp;
      
      double log_full_bi, log_full_bi_p;
      
      log_full_bi = 0.0;
      log_full_bi = log_full_bi - arma::sum(lgamma(exp(loggamma.col(jj)))) +
        arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)))
        + ::Rf_dnorm4(bi_zz(jj, hh), mu_hs, lambda_hs2(jj, hh)*tau_hs2(ll), 1);
      
      double propbipm = prop_per_bi2(jj, hh);
      propbipm = adap_prop(propbipm, QQ, JJ, QQ);
      bi_p_mat2(pp, qq) += propbipm;
      
      loggamma_p = loggamma.col(jj) - trans(threshold_mat(bZZ(pp, qq) *
        ZZ.col(qq).t(), soglia2)) % ZZ.col(qq) + trans(threshold_mat(bi_p_mat2(pp, qq) *
        ZZ.col(qq).t(), soglia2)) % ZZ.col(qq);
      
      log_full_bi_p = 0.0;
      log_full_bi_p = log_full_bi_p - arma::sum(lgamma(exp(loggamma_p))) +
        arma::sum(exp(loggamma_p)%log(SS.col(jj)))
        + ::Rf_dnorm4(bi_p_mat2(pp, qq), mu_hs, lambda_hs2(jj, hh)*tau_hs2(ll), 1);
      
      ln_acp = log_full_bi_p - log_full_bi;
      
      lnu = log(arma::randu());
      
      if(lnu < ln_acp){
        bi_zz(jj, hh) = bi_p_mat2(pp, qq);
        acc_bi_flag2(jj, hh) = 1;
        loggamma.col(jj)= loggamma_p;
      } else {
        bi_p_mat2(pp, qq) -= propbipm;
      }//chiude ratio
    }
  }
  if(hered == 1){
    if((theta2(hh_theta) == 0) & (theta2(ll_theta) == 0.0)){
      acc_bi_flag2(jj, hh) = 0;
      bi_p_mat2(pp, qq) = 0.0;
      
      loggamma.col(jj) = loggamma.col(jj) - trans(threshold_mat(bZZ(pp, qq) *
        ZZ.col(qq).t(), soglia2)) % ZZ.col(qq) + trans(threshold_mat(bi_p_mat2(pp, qq) *
        ZZ.col(qq).t(), soglia2)) % ZZ.col(qq);
      
      bi_zz(jj, hh) = bi_p_mat2(pp, qq);
    } else {
      double lnu, ln_acp;
      
      double log_full_bi, log_full_bi_p;
      
      log_full_bi = 0.0;
      log_full_bi = log_full_bi - arma::sum(lgamma(exp(loggamma.col(jj)))) +
        arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)))
        + ::Rf_dnorm4(bi_zz(jj, hh), mu_hs, lambda_hs2(jj, hh)*tau_hs2(ll), 1);
      
      double propbipm = prop_per_bi2(jj, hh);
      propbipm = adap_prop(propbipm, QQ, JJ, QQ);
      bi_p_mat2(pp, qq) += propbipm;
      
      loggamma_p = loggamma.col(jj) - trans(threshold_mat(bZZ(pp, qq) *
        ZZ.col(qq).t(), soglia2)) % ZZ.col(qq) + trans(threshold_mat(bi_p_mat2(pp, qq) *
        ZZ.col(qq).t(), soglia2)) % ZZ.col(qq);
      
      log_full_bi_p = 0.0;
      log_full_bi_p = log_full_bi_p - arma::sum(lgamma(exp(loggamma_p))) +
        arma::sum(exp(loggamma_p)%log(SS.col(jj)))
        + ::Rf_dnorm4(bi_p_mat2(pp, qq), mu_hs, lambda_hs2(jj, hh)*tau_hs2(ll), 1);
      
      ln_acp = log_full_bi_p - log_full_bi;
      
      lnu = log(arma::randu());
      
      if(lnu < ln_acp){
        bi_zz(jj, hh) = bi_p_mat2(pp, qq);
        acc_bi_flag2(jj, hh) = 1;
        loggamma.col(jj)= loggamma_p;
      } else {
        bi_p_mat2(pp, qq) -= propbipm;
      }//chiude ratio
    }
  }
  if(hered == 0){
    double lnu, ln_acp;
    
    double log_full_bi, log_full_bi_p;
    
    log_full_bi = 0.0;
    log_full_bi = log_full_bi - arma::sum(lgamma(exp(loggamma.col(jj)))) +
      arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)))
      + ::Rf_dnorm4(bi_zz(jj, hh), mu_hs, lambda_hs2(jj, hh)*tau_hs2(ll), 1);
    
    double propbipm = prop_per_bi2(jj, hh);
    propbipm = adap_prop(propbipm, QQ, JJ, QQ);
    bi_p_mat2(pp, qq) += propbipm;
    
    loggamma_p = loggamma.col(jj) - trans(threshold_mat(bZZ(pp, qq) *
      ZZ.col(qq).t(), soglia2)) % ZZ.col(qq) + trans(threshold_mat(bi_p_mat2(pp, qq) *
      ZZ.col(qq).t(), soglia2)) % ZZ.col(qq);
    
    log_full_bi_p = 0.0;
    log_full_bi_p = log_full_bi_p - arma::sum(lgamma(exp(loggamma_p))) +
      arma::sum(exp(loggamma_p)%log(SS.col(jj)))
      + ::Rf_dnorm4(bi_p_mat2(pp, qq), mu_hs, lambda_hs2(jj, hh)*tau_hs2(ll), 1);
    
    ln_acp = log_full_bi_p - log_full_bi;
    
    lnu = log(arma::randu());
    
    if(lnu < ln_acp){
      bi_zz(jj, hh) = bi_p_mat2(pp, qq);
      acc_bi_flag2(jj, hh) = 1;
      loggamma.col(jj)= loggamma_p;
    } else {
      bi_p_mat2(pp, qq) -= propbipm;
    }//chiude ratio
  }
}

void update_alpha(arma::mat XX, arma::mat ZZ, arma::mat SS, arma::mat &loggamma,
                  arma::vec &alpha, arma::vec theta, arma::vec theta2,
                  arma::mat &bi_xz, arma::mat &bi_zz, double soglia,
                  double soglia2, arma::vec &prop_per_alpha, arma::vec &acc_alpha_flag,
                  arma::vec mu_al, arma::vec sig_al, int jj){
  
  int PP = XX.n_cols;
  int QQ = ZZ.n_cols;
  int n = SS.n_rows;
  
  double sig_prop = prop_per_alpha(jj);
  arma::vec loggamma_p(n);
  
  arma::vec alpha_p = alpha;
  double lnu, ln_acp;
  
  arma::mat bXZ(PP, QQ);
  arma::mat bZZ(QQ, QQ);
  
  for(int qq2 = 0; qq2 < QQ; qq2++){
    for(int pp2 = 0; pp2 < PP; pp2++){
      int hh2 = qq2+((pp2)*QQ);
      bXZ(pp2, qq2)=bi_xz(jj, hh2);
    }
    
    for(int pp2 = 0; pp2 < QQ; pp2++){
      int hh2 = qq2 + ((pp2) * QQ);
      bZZ(pp2, qq2) = bi_zz(jj, hh2);
    }
  }
  
  // Prepare the current and proposed full conditional values
  double log_full_alpha, log_full_alpha_p;
  
  // Calculate the full conditional for the current value
  log_full_alpha = 0.0;
  log_full_alpha = log_full_alpha - arma::sum(lgamma(exp(loggamma.col(jj))))
    + arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)))
    + ::Rf_dnorm4(alpha(jj), mu_al(jj), pow(sig_al(jj), .5), 1);
    
    // Propose
    alpha_p(jj)= alpha(jj) + ::Rf_rnorm(0.0, sig_prop);
    
    loggamma_p = loggamma.col(jj) - alpha(jj) + alpha_p(jj);
    
    // Calculate the full conditional for the proposed value
    log_full_alpha_p = 0;
    log_full_alpha_p = log_full_alpha_p - arma::sum(lgamma(exp(loggamma_p)))
      + arma::sum(exp(loggamma_p) % log(SS.col(jj)))
      + ::Rf_dnorm4(alpha_p(jj), mu_al(jj), pow(sig_al(jj), .5), 1);
      
      ln_acp = log_full_alpha_p - log_full_alpha;
      lnu = log(arma::randu());
      
      if(lnu < ln_acp){
        
        acc_alpha_flag(jj) = 1;
        alpha(jj) = alpha_p(jj);
        
        loggamma.col(jj) = loggamma_p;
      }
} // Close function

void swap(arma::mat XX, arma::mat ZZ, arma::mat SS, arma::mat &loggamma,
          arma::vec &alpha, arma::vec &theta_temp, arma::vec &theta_temp2,
          arma::vec theta, arma::vec theta2, arma::mat &bi_xz, double soglia,
          arma::mat &bi_zz, double soglia2, arma::vec &acc_theta_flag,
          arma::vec &acc_theta_flag2, arma::vec &inclusion_indicator,
          arma::vec &inclusion_indicator2, arma::mat mu_theta, arma::mat mu_theta2,
          arma::mat sigma_theta, arma::mat sigma_theta2, double aa_hp,
          double bb_hp, double aa_hp2, double bb_hp2, int jj, bool primo = true){
  
  int PP = XX.n_cols;
  int QQ = ZZ.n_cols;
  int JJ = SS.n_cols;
  int n = XX.n_rows;
  
  double sigma_prop;
  double mu_prop;
  
  arma::vec loggamma_p(n);
  
  double lnu, ln_acp;
  
  double log_full_theta, log_full_theta_p;
  
  arma::mat bXZ(PP, QQ);
  arma::mat bZZ(QQ, QQ);
  
  for(int qq2 = 0; qq2 < QQ; qq2++){
    for(int pp2 = 0; pp2 < PP; pp2++){
      int hh2 = qq2+((pp2)*QQ);
      bXZ(pp2, qq2)=bi_xz(jj, hh2);
    }
    
    for(int pp2 = 0; pp2 < QQ; pp2++){
      int hh2 = qq2 + ((pp2) * QQ);
      bZZ(pp2, qq2) = bi_zz(jj, hh2);
    }
  }
  
  if(primo == true){
    
    int pp;
    arma::vec theta_p = theta_temp;
    
    int inclusion_indicator_p = 0;
    
    for(pp = 0; pp < PP; pp++){
      int hh = jj + pp * JJ;
      
      log_full_theta = 0.0;
      log_full_theta = log_full_theta - arma::sum(lgamma(exp(loggamma.col(jj))))
        + arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)))
        + lprior_bb(theta_temp(hh), inclusion_indicator(hh), sigma_theta(pp, jj),
                    mu_theta(pp, jj), aa_hp, bb_hp);
      
      sigma_prop = pow(sigma_theta(pp,jj), 0.5);
      mu_prop = mu_theta(pp,jj);
      
      if(inclusion_indicator(hh)==0){
        theta_p(hh) = mu_prop + ::Rf_rnorm(0.0, sigma_prop);
        inclusion_indicator_p = 1;
      } else {
        theta_p(hh) = 0.0;
        inclusion_indicator_p = 0;
      }
      
      loggamma_p = loggamma.col(jj) - theta_temp(hh)*XX.col(pp) + theta_p(hh)*XX.col(pp);
      
      log_full_theta_p = 0.0;
      log_full_theta_p = log_full_theta_p - arma::sum(lgamma(exp(loggamma_p)))
        + arma::sum(exp(loggamma_p)%log(SS.col(jj)))
        + lprior_bb(theta_p(hh), inclusion_indicator_p, sigma_theta(pp, jj), mu_theta(pp, jj), aa_hp, bb_hp);
        
        ln_acp =log_full_theta_p -log_full_theta;
        lnu = log(arma::randu());
        
        if(lnu < ln_acp){
          acc_theta_flag(hh) = 1;
          theta_temp(hh)= theta_p(hh);
          inclusion_indicator(hh)=inclusion_indicator_p;
          loggamma.col(jj)=loggamma_p;
        }//chiude if del ratio
    }//chiude il ciclo su pp
  } else {//SECONDA PARTE
    
    int qq;
    arma::vec theta_p = theta_temp2;
    
    int inclusion_indicator_p = 0;
    
    for(qq = 0; qq < QQ; qq++){
      int hh = jj + qq * JJ;
      
      //calcolo denominatore MH ratio "allo stato corrente"
      log_full_theta = 0.0;
      log_full_theta = log_full_theta - arma::sum(lgamma(exp(loggamma.col(jj))))
        + arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)))
        + lprior_bb(theta_temp2(hh), inclusion_indicator2(hh), sigma_theta2(qq, jj), mu_theta2(qq, jj), aa_hp2, bb_hp2);
        
        //propongo un nuovo valore per il beta_qj dalla media e varianza a priori
        sigma_prop = pow(sigma_theta2(qq,jj), 0.5);
        mu_prop = mu_theta2(qq,jj);
        
        //swap and sample
        if(inclusion_indicator2(hh)==0){
          theta_p(hh) = mu_prop + ::Rf_rnorm(0.0, sigma_prop);
          inclusion_indicator_p = 1;
        } else {
          theta_p(hh) = 0.0;
          inclusion_indicator_p = 0;
        }
        
        loggamma_p = loggamma.col(jj) - theta_temp2(hh)*ZZ.col(qq) + theta_p(hh)*ZZ.col(qq);
        
        //calcolo numeratore MH ratio "allo stato \star"
        log_full_theta_p = 0.0;
        log_full_theta_p = log_full_theta_p - arma::sum(lgamma(exp(loggamma_p)))
          + arma::sum(exp(loggamma_p)%log(SS.col(jj)))
          + lprior_bb(theta_p(hh), inclusion_indicator_p, sigma_theta2(qq, jj), mu_theta2(qq, jj), aa_hp2, bb_hp2);
          
          ln_acp = log_full_theta_p -log_full_theta;
          lnu = log(arma::randu());
          
          // MH ratio
          if(lnu < ln_acp){
            acc_theta_flag2(hh) = 1;
            theta_temp2(hh)= theta_p(hh);
            inclusion_indicator2(hh)=inclusion_indicator_p;
            loggamma.col(jj)=loggamma_p;
          }//chiude if del ratio
    }//chiude for su qq
  }
}//chiude la funzione

arma::mat alpha_ksi_prod(arma::mat alpha, arma::mat ksi, arma::vec dj,
                         bool linear){
  
  int PP = alpha.n_rows;
  int KK = ksi.n_cols;
  
  if (linear == true) {
    arma::mat beta(PP, PP);
    beta = alpha % ksi;
    return beta;
  } else {
    arma::mat beta(PP, KK, arma::fill::zeros);
    arma::vec cdj(PP-1);
    cdj = cumsum(dj);
    
    for (int pp2 = 0; pp2 < PP; pp2++) {
      for(int qq2 = pp2 + 1; qq2 < PP; qq2++){
        for(int ii = cdj(qq2-1); ii < cdj(qq2); ii++){
          beta(pp2, ii) = alpha(pp2, qq2) * ksi(pp2, ii);
        }
      }
    }
    return beta;
  }
}

arma::mat update_alpha_xx(arma::mat XX, arma::mat ZZ, arma::mat XXl, arma::mat XXT,
                          arma::mat SS, arma::mat &loggamma, arma::vec &alpha, arma::vec theta,
                          arma::vec theta2, arma::mat &bi_xz, double soglia,
                          arma::mat &bi_zz, double soglia2, arma::vec dj, arma::mat alpha_lin,
                          arma::mat alpha_nl, arma::mat tausq_lin, arma::mat tausq_nl,
                          arma::mat gamma_lin, arma::mat gamma_nl, arma::mat ksi_lin, arma::mat ksi_nl,
                          bool linear, int jj, int hered, bool constraint){
  
  if(linear == true){
    int PP = XX.n_cols;
    int QQ = ZZ.n_cols;
    int KK = XXT.n_cols;
    int JJ = SS.n_cols;
    int n = ZZ.n_rows;
    
    arma::vec loggamma_p(n);
    
    int hh_theta, ll_theta;//indici x hereditariety
    
    arma::mat bXZ(PP, QQ);
    
    for(int qq2 = 0; qq2 < QQ; qq2++){
      for(int pp2 = 0; pp2 < PP; pp2++){
        int hh2 = qq2 + ((pp2) * QQ);
        bXZ(pp2, qq2) = bi_xz(jj, hh2);
      }
    }
    
    arma::mat bZZ(QQ, QQ);
    
    for(int qq2 = 0; qq2 < QQ; qq2++){
      for(int pp2 = 0; pp2 < QQ; pp2++){
        int hh2 = qq2 + ((pp2) * QQ);
        bZZ(pp2, qq2) = bi_zz(jj, hh2);
      }
    }
    
    double mhr_num, mhr_den, mhr_diff, mhr_u;
    
    arma::mat alpha_p(PP, PP);
    alpha_p.fill(0.0);
    
    for(int pp2 = 0; pp2 < PP; pp2++){
      for(int qq2 = pp2 + 1; qq2 < PP; qq2++){
        hh_theta = jj + pp2 * JJ;
        ll_theta = jj + qq2 * JJ;
        //cambia |/& a seconda dell'assunzione sull ereditarietà
        if(hered == 2){
          if((theta(hh_theta) == 0.0) | (theta(ll_theta) == 0.0)){
            alpha_lin(pp2, qq2) = 0.0;
            arma::mat beta_lin_j(PP, PP, arma::fill::zeros);
            arma::mat beta_nl_j(PP, KK, arma::fill::zeros);
            
            beta_lin_j = alpha_ksi_prod(alpha_lin, ksi_lin, dj, true);
            beta_nl_j = alpha_ksi_prod(alpha_nl, ksi_nl, dj, false);
            loggamma_p = calculate_gamma(XX, ZZ, XXl, XXT, alpha, theta, theta2, bXZ,
                                         beta_lin_j, beta_nl_j, soglia, bZZ, soglia2, jj);
            
          } else{
            alpha_p(pp2, qq2) = alpha_lin(pp2, qq2) + ::Rf_rnorm(0, .1);
            
            mhr_den = 0.0;
            mhr_den = mhr_den - arma::sum(lgamma(exp(loggamma.col(jj)))) +
              arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)));
            
            double gammatau;
            gammatau = pow(tausq_lin(pp2, qq2), .5) * gamma_lin(pp2, qq2);
            mhr_den += ::Rf_dnorm4(alpha_lin(pp2, qq2), 0.0, gammatau, 1);
            
            arma::mat beta_lin_j_p(PP, PP, arma::fill::zeros);
            arma::mat beta_nl_j(PP, KK, arma::fill::zeros);
            
            beta_lin_j_p = alpha_ksi_prod(alpha_p, ksi_lin, dj, true);
            beta_nl_j = alpha_ksi_prod(alpha_nl, ksi_nl, dj, false);
            loggamma_p = calculate_gamma(XX, ZZ, XXl, XXT, alpha, theta, theta2, bXZ,
                                         beta_lin_j_p, beta_nl_j, soglia, bZZ, soglia2, jj);
            
            mhr_num = 0.0;
            mhr_num = mhr_num - arma::sum(lgamma(exp(loggamma_p))) +
              arma::sum(exp(loggamma_p)%log(SS.col(jj)));
            
            gammatau = pow(tausq_lin(pp2, qq2), .5) * gamma_lin(pp2, qq2);
            mhr_num += ::Rf_dnorm4(alpha_p(pp2, qq2), 0.0, gammatau, 1);
            
            mhr_diff = mhr_num - mhr_den;
            
            mhr_u = log(arma::randu());
            
            if (mhr_u < mhr_diff) {
              alpha_lin(pp2, qq2) = alpha_p(pp2, qq2);
              loggamma.col(jj) = loggamma_p;
            }//chiude ratio
          }//chiude l else dell ereditarietà
        }
        if(hered == 1){
          if((theta(hh_theta) == 0.0) & (theta(ll_theta) == 0.0)){
            alpha_lin(pp2, qq2) = 0.0;
            arma::mat beta_lin_j(PP, PP, arma::fill::zeros);
            arma::mat beta_nl_j(PP, KK, arma::fill::zeros);
            
            beta_lin_j = alpha_ksi_prod(alpha_lin, ksi_lin, dj, true);
            beta_nl_j = alpha_ksi_prod(alpha_nl, ksi_nl, dj, false);
            loggamma_p = calculate_gamma(XX, ZZ, XXl, XXT, alpha, theta, theta2, bXZ,
                                         beta_lin_j, beta_nl_j, soglia, bZZ, soglia2, jj);
            
          } else{
            alpha_p(pp2, qq2) = alpha_lin(pp2, qq2) + ::Rf_rnorm(0, .1);
            
            mhr_den = 0.0;
            mhr_den = mhr_den - arma::sum(lgamma(exp(loggamma.col(jj)))) +
              arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)));
            
            double gammatau;
            gammatau = pow(tausq_lin(pp2, qq2), .5) * gamma_lin(pp2, qq2);
            mhr_den += ::Rf_dnorm4(alpha_lin(pp2, qq2), 0.0, gammatau, 1);
            
            arma::mat beta_lin_j_p(PP, PP, arma::fill::zeros);
            arma::mat beta_nl_j(PP, KK, arma::fill::zeros);
            
            beta_lin_j_p = alpha_ksi_prod(alpha_p, ksi_lin, dj, true);
            beta_nl_j = alpha_ksi_prod(alpha_nl, ksi_nl, dj, false);
            loggamma_p = calculate_gamma(XX, ZZ, XXl, XXT, alpha, theta, theta2, bXZ,
                                         beta_lin_j_p, beta_nl_j, soglia, bZZ, soglia2, jj);
            
            mhr_num = 0.0;
            mhr_num = mhr_num - arma::sum(lgamma(exp(loggamma_p))) +
              arma::sum(exp(loggamma_p)%log(SS.col(jj)));
            
            gammatau = pow(tausq_lin(pp2, qq2), .5) * gamma_lin(pp2, qq2);
            mhr_num += ::Rf_dnorm4(alpha_p(pp2, qq2), 0.0, gammatau, 1);
            
            mhr_diff = mhr_num - mhr_den;
            
            mhr_u = log(arma::randu());
            
            if (mhr_u < mhr_diff) {
              alpha_lin(pp2, qq2) = alpha_p(pp2, qq2);
              loggamma.col(jj) = loggamma_p;
            }//chiude ratio
          }//chiude l else dell ereditarietà
        }
        if(hered == 0){
          alpha_p(pp2, qq2) = alpha_lin(pp2, qq2) + ::Rf_rnorm(0, .1);
          
          mhr_den = 0.0;
          mhr_den = mhr_den - arma::sum(lgamma(exp(loggamma.col(jj)))) +
            arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)));
          
          double gammatau;
          gammatau = pow(tausq_lin(pp2, qq2), .5) * gamma_lin(pp2, qq2);
          mhr_den += ::Rf_dnorm4(alpha_lin(pp2, qq2), 0.0, gammatau, 1);
          
          arma::mat beta_lin_j_p(PP, PP, arma::fill::zeros);
          arma::mat beta_nl_j(PP, KK, arma::fill::zeros);
          
          beta_lin_j_p = alpha_ksi_prod(alpha_p, ksi_lin, dj, true);
          beta_nl_j = alpha_ksi_prod(alpha_nl, ksi_nl, dj, false);
          loggamma_p = calculate_gamma(XX, ZZ, XXl, XXT, alpha, theta, theta2, bXZ,
                                       beta_lin_j_p, beta_nl_j, soglia, bZZ, soglia2, jj);
          
          mhr_num = 0.0;
          mhr_num = mhr_num - arma::sum(lgamma(exp(loggamma_p))) +
            arma::sum(exp(loggamma_p)%log(SS.col(jj)));
          
          gammatau = pow(tausq_lin(pp2, qq2), .5) * gamma_lin(pp2, qq2);
          mhr_num += ::Rf_dnorm4(alpha_p(pp2, qq2), 0.0, gammatau, 1);
          
          mhr_diff = mhr_num - mhr_den;
          
          mhr_u = log(arma::randu());
          
          if (mhr_u < mhr_diff) {
            alpha_lin(pp2, qq2) = alpha_p(pp2, qq2);
            loggamma.col(jj) = loggamma_p;
          }//chiude ratio
        }//chiude l else dell ereditarietà
      }
    }//
    return alpha_lin;
  } else {
    int PP = XX.n_cols;
    int KK = XXT.n_cols;
    int QQ = ZZ.n_cols;
    int JJ = SS.n_cols;
    int n = ZZ.n_rows;
    
    arma::vec loggamma_p(n);
    
    arma::mat bXZ(PP, QQ);
    
    for(int qq2 = 0; qq2 < QQ; qq2++){
      for(int pp2 = 0; pp2 < PP; pp2++){
        int hh2 = qq2 + ((pp2) * QQ);
        bXZ(pp2, qq2) = bi_xz(jj, hh2);
      }
    }
    
    arma::mat bZZ(QQ, QQ);
    
    for(int qq2 = 0; qq2 < QQ; qq2++){
      for(int pp2 = 0; pp2 < QQ; pp2++){
        int hh2 = qq2 + ((pp2) * QQ);
        bZZ(pp2, qq2) = bi_zz(jj, hh2);
      }
    }
    
    double mhr_num, mhr_den, mhr_diff, mhr_u;
    
    arma::mat alpha_p(PP, PP);
    alpha_p.fill(0.0);
    
    //indici x hereditariety
    int hh_theta;
    int ll_theta;
    
    for(int pp2 = 0; pp2 < PP; pp2++){
      for(int qq2 = pp2 + 1; qq2 < PP; qq2++){
        hh_theta = jj + pp2 * JJ;
        ll_theta = jj + qq2 * JJ;
        //unica condizione sul constraint
        if(constraint == true){
          if(gamma_lin(pp2, qq2) != 1.0){
            alpha_nl(pp2, qq2) = 0.0;
            arma::mat beta_lin_j(PP, PP, arma::fill::zeros);
            arma::mat beta_nl_j(PP, KK, arma::fill::zeros);
            
            beta_lin_j = alpha_ksi_prod(alpha_lin, ksi_lin, dj, true);
            beta_nl_j = alpha_ksi_prod(alpha_nl, ksi_nl, dj, false);
            loggamma_p = calculate_gamma(XX, ZZ, XXl, XXT, alpha, theta, theta2, bXZ,
                                         beta_lin_j, beta_nl_j, soglia, bZZ, soglia2, jj);
          } else{
            alpha_p(pp2, qq2) = alpha_nl(pp2, qq2) + ::Rf_rnorm(0, .1);
            
            mhr_den = 0.0;
            mhr_den = mhr_den - arma::sum(lgamma(exp(loggamma.col(jj)))) +
              arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)));
            
            double gammatau;
            gammatau = pow(tausq_nl(pp2, qq2), .5) * gamma_nl(pp2, qq2);
            mhr_den += ::Rf_dnorm4(alpha_nl(pp2, qq2), 0.0, gammatau, 1);
            arma::mat beta_lin_j(PP, PP, arma::fill::zeros);
            arma::mat beta_nl_j_p(PP, KK, arma::fill::zeros);
            
            beta_lin_j = alpha_ksi_prod(alpha_lin, ksi_lin, dj, true);
            beta_nl_j_p = alpha_ksi_prod(alpha_p, ksi_nl, dj, false);
            loggamma_p = calculate_gamma(XX, ZZ, XXl, XXT, alpha, theta, theta2, bXZ,
                                         beta_lin_j, beta_nl_j_p, soglia, bZZ, soglia2, jj);
            
            mhr_num = 0.0;
            mhr_num = mhr_num - arma::sum(lgamma(exp(loggamma_p))) +
              arma::sum(exp(loggamma_p)%log(SS.col(jj)));
            
            gammatau = pow(tausq_nl(pp2, qq2), .5) * gamma_nl(pp2, qq2);
            mhr_num += ::Rf_dnorm4(alpha_p(pp2, qq2), 0.0, gammatau, 1);
            
            
            mhr_diff = mhr_num - mhr_den;
            
            mhr_u = log(arma::randu());
            
            if (mhr_u < mhr_diff) {
              alpha_nl(pp2, qq2) = alpha_p(pp2, qq2);
              loggamma.col(jj) = loggamma_p;
            }//chiude ratio
          }//chiude else hereditariety
        } else {//chiude if(constraint == true)
          /*if(gamma_lin(pp2, qq2) != 1.0) {
            alpha_nl(pp2, qq2) = 0.0;
            arma::mat beta_lin_j(PP, PP, arma::fill::zeros);
            arma::mat beta_nl_j(PP, KK, arma::fill::zeros);
            
            beta_lin_j = alpha_ksi_prod(alpha_lin, ksi_lin, dj, true);
            beta_nl_j = alpha_ksi_prod(alpha_nl, ksi_nl, dj, false);
            loggamma_p = calculate_gamma(XX, ZZ, XXl, XXT, alpha, theta, theta2, bXZ,
                                         beta_lin_j, beta_nl_j, soglia, bZZ, soglia2, jj);
          } else{*/
            alpha_p(pp2, qq2) = alpha_nl(pp2, qq2) + ::Rf_rnorm(0, .1);
            
            
            mhr_den = 0.0;
            mhr_den = mhr_den - arma::sum(lgamma(exp(loggamma.col(jj)))) +
              arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)));
            
            double gammatau;
            gammatau = pow(tausq_nl(pp2, qq2), .5) * gamma_nl(pp2, qq2);
            mhr_den += ::Rf_dnorm4(alpha_nl(pp2, qq2), 0.0, gammatau, 1);
            
            arma::mat beta_lin_j(PP, PP, arma::fill::zeros);
            arma::mat beta_nl_j_p(PP, KK, arma::fill::zeros);
            
            beta_lin_j = alpha_ksi_prod(alpha_lin, ksi_lin, dj, true);
            beta_nl_j_p = alpha_ksi_prod(alpha_p, ksi_nl, dj, false);
            loggamma_p = calculate_gamma(XX, ZZ, XXl, XXT, alpha, theta, theta2, bXZ,
                                         beta_lin_j, beta_nl_j_p, soglia, bZZ, soglia2, jj);
            
            mhr_num = 0.0;
            mhr_num = mhr_num - arma::sum(lgamma(exp(loggamma_p))) +
              arma::sum(exp(loggamma_p)%log(SS.col(jj)));
            
            gammatau = pow(tausq_nl(pp2, qq2), .5) * gamma_nl(pp2, qq2);
            mhr_num += ::Rf_dnorm4(alpha_p(pp2, qq2), 0.0, gammatau, 1);
            
            mhr_diff = mhr_num - mhr_den;
            
            mhr_u = log(arma::randu());
            
            if (mhr_u < mhr_diff) {
              alpha_nl(pp2, qq2) = alpha_p(pp2, qq2);
              loggamma.col(jj) = loggamma_p;
            }//chiude ratio
          //}//chiude else gammalin =! 1.0
        }//chiude else constraint
      }//chiude qq
    }//chiude pp
    return alpha_nl;
  }//chiude else
}//chiude la funzione update_alpha_xx()

arma::mat update_ksi(arma::mat XX, arma::mat ZZ, arma:: mat XXl, arma::mat XXT, arma::mat SS,
                     arma::mat &loggamma, arma::vec &alpha, arma::vec theta, arma::vec theta2,
                     arma::mat &bi_xz, double soglia, arma::mat &bi_zz, double soglia2,
                     arma::vec dj, arma::mat alpha_lin,
                     arma::mat alpha_nl, arma::mat tausq, arma::mat tausq_nl, arma::mat gamma_lin,
                     arma::mat gamma_nl, arma::mat ksi_lin, arma::mat ksi_nl, arma::mat m_lin,
                     arma::mat m_nl, bool linear, int jj){
  
  if(linear == true){
    int PP = XX.n_cols;
    int KK = XXT.n_cols;
    int QQ = ZZ.n_cols;
    int n = ZZ.n_rows;
    
    arma::vec loggamma_p(n);
    
    arma::mat bXZ(PP, QQ);
    
    for(int qq2 = 0; qq2 < QQ; qq2++){
      for(int pp2 = 0; pp2 < PP; pp2++){
        int hh2 = qq2 + ((pp2) * QQ);
        bXZ(pp2, qq2) = bi_xz(jj, hh2);
      }
    }
    
    arma::mat bZZ(QQ, QQ);
    
    for(int qq2 = 0; qq2 < QQ; qq2++){
      for(int pp2 = 0; pp2 < QQ; pp2++){
        int hh2 = qq2 + ((pp2) * QQ);
        bZZ(pp2, qq2) = bi_zz(jj, hh2);
      }
    }
    
    double mhr_num, mhr_den, mhr_diff, mhr_u;
    
    arma::mat ksi_p(PP, PP, arma::fill::zeros);
    
    for(int pp2 = 0; pp2 < PP; pp2++){
      for(int qq2 = pp2 + 1; qq2 < PP; qq2++){
        ksi_p(pp2, qq2) = ksi_lin(pp2, qq2) + ::Rf_rnorm(0, .1);
        
        mhr_den = 0.0;
        mhr_den = mhr_den - arma::sum(lgamma(exp(loggamma.col(jj)))) +
          arma::sum(exp(loggamma.col(jj)) % log(SS.col(jj)));
        
        mhr_den += ::Rf_dnorm4(ksi_lin(pp2, qq2), m_lin(pp2, qq2), 1.0, 1);
        
        arma::mat beta_lin_j_p(PP, PP, arma::fill::zeros);
        arma::mat beta_nl_j(PP, KK, arma::fill::zeros);
        beta_lin_j_p = alpha_ksi_prod(alpha_lin, ksi_p, dj, true);
        beta_nl_j = alpha_ksi_prod(alpha_nl, ksi_nl, dj, false);
        
        loggamma_p = calculate_gamma(XX, ZZ, XXl, XXT, alpha, theta, theta2, bXZ,
                                     beta_lin_j_p, beta_nl_j, soglia, bZZ, soglia2, jj);
        
        mhr_num = 0.0;
        mhr_num = mhr_num - arma::sum(lgamma(exp(loggamma_p))) +
          arma::sum(exp(loggamma_p)%log(SS.col(jj)));
        
        mhr_num += ::Rf_dnorm4(ksi_p(pp2, qq2), m_lin(pp2, qq2), 1.0, 1);
        
        mhr_diff = mhr_num - mhr_den;
        
        mhr_u = log(arma::randu());
        
        if (mhr_u < mhr_diff) {
          ksi_lin(pp2, qq2) = ksi_p(pp2, qq2);
          loggamma.col(jj) = loggamma_p;
        }//chiude ratio
      }
    }
    return ksi_lin;
  } else {//chiude (linear == true)
    int PP = XX.n_cols;
    int QQ = ZZ.n_cols;
    int KK = XXT.n_cols;
    int n = ZZ.n_rows;
    
    arma::vec loggamma_p(n);
    
    arma::mat bXZ(PP, QQ);
    
    for(int qq2 = 0; qq2 < QQ; qq2++){
      for(int pp2 = 0; pp2 < PP; pp2++){
        int hh2 = qq2 + ((pp2) * QQ);
        bXZ(pp2, qq2) = bi_xz(jj, hh2);
      }
    }
    
    arma::mat bZZ(QQ, QQ);
    
    for(int qq2 = 0; qq2 < QQ; qq2++){
      for(int pp2 = 0; pp2 < QQ; pp2++){
        int hh2 = qq2 + ((pp2) * QQ);
        bZZ(pp2, qq2) = bi_zz(jj, hh2);
      }
    }
    
    double mhr_num, mhr_den, mhr_diff, mhr_u;
    
    arma::mat ksi_p(PP, KK, arma::fill::zeros);
    
    arma::vec cdj(PP-1);
    cdj = cumsum(dj);
    double rs = 0.0;
    for(int pp2 = 0; pp2 < PP; pp2++){
      for(int qq2 = pp2 + 1; qq2 < PP; qq2++){
        rs = ::Rf_rnorm(0, .1);
        for(int ii = cdj(qq2-1); ii < cdj(qq2); ii++){
          ksi_p(pp2, ii) = ksi_nl(pp2, ii) + rs;
        }
      }
    }
    
    mhr_den = 0.0;
    mhr_den = mhr_den - arma::sum(lgamma(exp(loggamma.col(jj)))) +
      arma::sum(exp(loggamma.col(jj)) % log(SS.col(jj)));
    
    for(int pp2 = 0; pp2 < PP; pp2++){
      for(int qq2 = pp2 + 1; qq2 < PP; qq2++){
        for(int ii = cdj(qq2-1); ii < cdj(qq2); ii++){
          mhr_den += ::Rf_dnorm4(ksi_nl(pp2, ii), m_nl(pp2, ii), 1.0, 1);
        }
      }
    }
    
    arma::mat beta_lin_j(PP, PP, arma::fill::zeros);
    arma::mat beta_nl_j_p(PP, KK, arma::fill::zeros);
    beta_lin_j = alpha_ksi_prod(alpha_lin, ksi_lin, dj, true);
    beta_nl_j_p = alpha_ksi_prod(alpha_nl, ksi_p, dj, false);
    loggamma_p = calculate_gamma(XX, ZZ, XXl, XXT, alpha, theta, theta2, bXZ,
                                 beta_lin_j, beta_nl_j_p, soglia, bZZ, soglia2, jj);
    
    mhr_num = 0.0;
    mhr_num = mhr_num - arma::sum(lgamma(exp(loggamma_p))) +
      arma::sum(exp(loggamma_p)%log(SS.col(jj)));
    
    for(int pp2 = 0; pp2 < PP; pp2++){
      for(int qq2 = pp2 + 1; qq2 < PP; qq2++){
        for(int ii = cdj(qq2-1); ii < cdj(qq2); ii++){
          mhr_num += ::Rf_dnorm4(ksi_p(pp2, ii), m_nl(pp2, ii), 1.0, 1);
        }
      }
    }
    
    mhr_diff = mhr_num - mhr_den;
    
    mhr_u = log(arma::randu());
    
    if (mhr_u < mhr_diff) {
      ksi_nl = ksi_p;
      loggamma.col(jj) = loggamma_p;
    }
    return ksi_nl;
  }
}//chiude funzione update_ksi()

arma::mat update_m(arma::mat m, arma::mat ksi){
  
  int PP = m.n_rows;
  int KK = m.n_cols;
  
  double u_m;
  
  for(int pp2 = 0; pp2 < PP; pp2++){
    for(int qq2 = 0; qq2 < KK; qq2++){
      u_m = 1/(1 + exp(-2 * ksi(pp2, qq2))) - ::Rf_runif(0, 1);
      if(u_m > 0){
        m(pp2, qq2) = 1;
      } else {
        m(pp2, qq2) = -1;
      }
    }//closes for loop kk
  }//closes for loop pp
  return m;//OK
}//closes update_m()

Rcpp::List rescale(arma::mat alpha, arma::mat ksi, arma::vec dj, bool linear){
  int PP = alpha.n_rows;
  if (linear == true) {
    double val;
    for(int pp2 = 0; pp2 < PP; pp2++){
      for(int qq2 = pp2 + 1; qq2 < PP; qq2++){
        val = abs(ksi(pp2, qq2));
        ksi(pp2, qq2) = ksi(pp2, qq2)/val;
        
        alpha(pp2, qq2) = alpha(pp2, qq2) * val;
      }
    }
    Rcpp::List rescale_out(2);
    rescale_out[0] = alpha;
    rescale_out[1] = ksi;
    return rescale_out;
  } else {
    arma::mat myval(PP, PP);
    myval.fill(0.0);
    arma::mat fs(PP, PP);
    fs.fill(0.0);
    arma::vec cdj(PP-1);
    cdj = cumsum(dj);
    
    for(int pp2 = 0; pp2 < PP; pp2++){
      for(int qq2 = pp2 + 1; qq2 < PP; qq2++){
        for(int ii = cdj(qq2-1); ii < cdj(qq2); ii++){
          myval(pp2, qq2) += abs(ksi(pp2, ii));
        }
        fs(pp2, qq2) = dj(pp2)/myval(pp2, qq2);
      }
    }
    
    for(int pp2 = 0; pp2 < PP; pp2++){
      for(int qq2 = pp2 + 1; qq2 < PP; qq2++){
        for(int ii = cdj(qq2-1); ii < cdj(qq2); ii++){
          ksi(pp2, ii) = fs(pp2, qq2) * ksi(pp2, ii);
        }
        alpha(pp2, qq2) = pow(fs(pp2, qq2), -1.0) * alpha(pp2, qq2);
      }
    }
    
    Rcpp::List rescale_out(2);
    rescale_out[0] = alpha;
    rescale_out[1] = ksi;
    return rescale_out;
  }
}// closes rescale() function

arma::mat update_tau(arma::mat tausq, arma::mat alpha, double atau, double btau,
                     arma::mat gamma, int jj){
  
  int PP = tausq.n_cols;
  double a = atau + .5;
  double b;
  
  for(int pp2 = 0; pp2 < PP; pp2++){
    for(int qq2 = pp2 + 1; qq2 < PP; qq2++){
      b = 1/(btau + std::pow(alpha(pp2, qq2), 2.0)/(2.0 * gamma(pp2, qq2)));
      tausq(pp2, qq2) = 1.0/::Rf_rgamma(a, b);
    }
  }
  return tausq;
}//closes update_tau()



arma::vec update_omega(arma::mat gamma, arma::vec omega, double v0, double aomega,
                       double bomega){
  int PP = gamma.n_cols;
  for(int pp2 = 0; pp2 < PP; pp2++){
    int included = 0;
    int excluded = 0;
    for(int qq2 = pp2 + 1; qq2 < PP; qq2++){
      included += (gamma(pp2, qq2) == 1.0);
      excluded += (gamma(pp2, qq2) == v0);
    }
    omega(pp2) = ::Rf_rbeta(aomega + included, bomega + excluded);
  }
  return omega;
}//closes update_omega()


arma::mat update_gamma(arma::mat alpha, arma::mat tausq, arma::mat gamma,
                       arma::vec omega, double v0, int jj, int JJ, arma::vec theta,
                       int hered){
  
  int PP = alpha.n_cols;
  double A;
  arma::mat p1(PP, PP);
  p1.fill(0);
  for(int pp2 = 0; pp2 < PP; pp2++){
    double c1 = log(omega(pp2)/(1-omega(pp2))) + log(v0) * .5;
    double c2 = (1 - v0)/(2 * v0);
    for(int qq2 = pp2 + 1; qq2 < PP; qq2++){
      double c3 = std::pow(alpha(pp2, qq2), 2.0)/tausq(pp2, qq2);
      A = exp( c1  + c2 * c3);
      if(!ISNAN(A)){
        if(A > 10000){
          p1(pp2, qq2) = 1;
          gamma(pp2, qq2) = 1;
        } else {
          if(A < 0.00001) {
            p1(pp2, qq2) = 0;
            gamma(pp2, qq2) = v0;
          } else {
            p1(pp2, qq2) =  A / (1 + A);
            gamma(pp2, qq2) = (::Rf_runif(0, 1) < p1(pp2, qq2)) ? 1.0 : v0;
          }
        }
      } // end if(A==NA)
    }//closes for loop qq2
  }//closes for loop pp2
  for(int pp2 = 0; pp2 < PP; pp2++){
    for(int qq2 = pp2 + 1; qq2 < PP; qq2++){
      int hh_theta = jj + pp2 * JJ;
      int ll_theta = jj + qq2 * JJ;
      //cambia |/& a seconda dell'assunzione sull ereditarietà
      if(hered == 2){
        if((theta(hh_theta) == 0.0) | (theta(ll_theta) == 0.0)){
          gamma(pp2, qq2) = v0;
        }
      }
      if(hered == 1){
        if((theta(hh_theta) == 0.0) & (theta(ll_theta) == 0.0)){
          gamma(pp2, qq2) = v0;
        }
      }
    }
  }
  return gamma;
}//closes UpdateGamma()

double logliksimpleC(arma::mat Y, arma::mat zeta){
  
  arma::vec zeta_rowsums(zeta.n_rows);
  arma::vec m(Y.n_rows);
  zeta_rowsums = myrowSums(zeta);
  m = myrowSums(Y);
  
  double ll = sum((lgamma(m + 1) + myrowSums(lgamma(Y + zeta)) + lgamma(zeta_rowsums)) -
                  (myrowSums(lgamma(Y + 1)) + myrowSums(lgamma(zeta)) + lgamma(zeta_rowsums + m)));
  
  return ll;
  
}

double update_soglia_ri(arma::mat XX, arma::mat ZZ, arma::mat Xl, arma::mat Xnl, arma::mat SS,
                        arma::mat &loggamma, arma::vec &alpha, arma::vec ran_int_exp,
                        arma::vec theta, arma::vec theta2, arma::mat &bi_xz, arma::cube bXXl,
                        arma::cube bXXnl, double soglia, arma::mat &bi_zz, double soglia2,
                        arma::vec &acc_soglia_flag, double lb, double ub,
                        int pos = 0){
  if(pos == 0){
    
    int ii, jj;
    int JJ = SS.n_cols;
    int n = XX.n_rows;
    
    int PP = XX.n_cols;
    int QQ = ZZ.n_cols;
    
    double soglia_p;
    
    double lnu, ln_acp;
    
    arma::mat loggamma_p(n, JJ);
    
    double log_full_soglia, log_full_soglia_p;
    
    log_full_soglia = 0.0;
    for(jj = 0; jj < JJ; jj ++){
      for(ii = 0; ii < n; ii++){
        log_full_soglia = log_full_soglia - lgamma(exp(loggamma(ii, jj)));
        log_full_soglia = log_full_soglia + exp(loggamma(ii, jj))*log(SS(ii, jj));
      }
    }
    
    log_full_soglia = log_full_soglia + ::Rf_dunif(soglia, lb, ub, 1);
    
    soglia_p = soglia + ::Rf_rnorm(0.0, .01);
    if(soglia_p >= ub){
      soglia_p = ub - (soglia_p - ub);
    }
    if(soglia_p <= lb){
      soglia_p = lb + (lb - soglia_p);
    }
    
    arma::mat bXZ(PP, QQ);
    arma::mat bZZ(QQ, QQ);
    
    for(jj = 0; jj < JJ; jj++){
      for(int qq2 = 0; qq2 < QQ; qq2++){
        for(int pp2 = 0; pp2 < PP; pp2++){
          int hh2 = qq2+((pp2)*QQ);
          bXZ(pp2, qq2)=bi_xz(jj, hh2);
        }
        
        for(int pp2 = 0; pp2 < QQ; pp2++){
          int hh2 = qq2 + ((pp2) * QQ);
          bZZ(pp2, qq2) = bi_zz(jj, hh2);
        }
      }
      
      loggamma_p.col(jj)=calculate_gamma_ri(XX, ZZ, Xl, Xnl, alpha, ran_int_exp, theta, theta2, bXZ,
                     bXXl.slice(jj), bXXnl.slice(jj),soglia_p, bZZ, soglia2, jj);
    }
    
    log_full_soglia_p = 0.0;
    for(jj = 0; jj < JJ; jj ++){
      for(ii = 0; ii <n; ii++){
        log_full_soglia_p = log_full_soglia_p - lgamma(exp(loggamma(ii, jj)));
        log_full_soglia_p = log_full_soglia_p + exp(loggamma(ii, jj))*log(SS(ii, jj));
      }
    }
    
    log_full_soglia_p = log_full_soglia_p + ::Rf_dunif(soglia_p, lb, ub, 1);
    
    
    ln_acp = log_full_soglia_p - log_full_soglia;
    lnu = log(arma::randu());
    
    if(lnu < ln_acp){
      soglia = soglia_p;
      acc_soglia_flag(pos) = 1;
      for(jj = 0; jj < JJ; jj++){
        loggamma.col(jj)= loggamma_p.col(jj);
      }
    }
  } else {
    
    int ii, jj;
    int JJ = SS.n_cols;
    int n = XX.n_rows;
    
    int PP = XX.n_cols;
    int QQ = ZZ.n_cols;
    
    double soglia_p;
    
    double lnu, ln_acp;
    
    arma::mat loggamma_p(n, JJ);
    
    double log_full_soglia, log_full_soglia_p;
    
    log_full_soglia = 0.0;
    for(jj = 0; jj < JJ; jj ++){
      for(ii = 0; ii < n; ii++){
        log_full_soglia = log_full_soglia - lgamma(exp(loggamma(ii, jj)));
        log_full_soglia = log_full_soglia + exp(loggamma(ii, jj))*log(SS(ii, jj));
      }
    }
    
    log_full_soglia = log_full_soglia + ::Rf_dunif(soglia2, lb, ub, 1);
    
    soglia_p = soglia2 + ::Rf_rnorm(0.0, .01);
    if(soglia_p >= ub){
      soglia_p = ub - (soglia_p - ub);
    }
    if(soglia_p <= lb){
      soglia_p = lb + (lb - soglia_p);
    }
    
    arma::mat bXZ(PP, QQ);
    arma::mat bZZ(QQ, QQ);
    
    for(jj = 0; jj < JJ; jj++){
      for(int qq2 = 0; qq2 < QQ; qq2++){
        for(int pp2 = 0; pp2 < PP; pp2++){
          int hh2 = qq2+((pp2)*QQ);
          bXZ(pp2, qq2)=bi_xz(jj, hh2);
        }
        
        for(int pp2 = 0; pp2 < QQ; pp2++){
          int hh2 = qq2 + ((pp2) * QQ);
          bZZ(pp2, qq2) = bi_zz(jj, hh2);
        }
      }
      
      loggamma_p.col(jj)=calculate_gamma_ri(XX, ZZ, Xl, Xnl, alpha, ran_int_exp, theta, theta2, bXZ,
                     bXXl.slice(jj), bXXnl.slice(jj), soglia, bZZ, soglia_p, jj);
    }
    
    log_full_soglia_p = 0.0;
    for(jj = 0; jj < JJ; jj ++){
      for(ii = 0; ii <n; ii++){
        log_full_soglia_p = log_full_soglia_p - lgamma(exp(loggamma(ii, jj)));
        log_full_soglia_p = log_full_soglia_p + exp(loggamma(ii, jj))*log(SS(ii, jj));
      }
    }
    
    log_full_soglia_p = log_full_soglia_p + ::Rf_dunif(soglia_p, lb, ub, 1);
    
    
    ln_acp = log_full_soglia_p - log_full_soglia;
    lnu = log(arma::randu());
    
    if(lnu < ln_acp){
      soglia = soglia_p;
      acc_soglia_flag(pos) = 1;
      for(jj = 0; jj < JJ; jj++){
        loggamma.col(jj)= loggamma_p.col(jj);
      }
    }
  }
  return soglia;
}

void update_bi_ri(arma::mat XX, arma::mat ZZ, arma::mat XXl, arma::mat XXT,
                  arma::mat SS, arma::mat &loggamma,
                  arma::mat lambda_hs, arma::vec tau_hs, arma::vec &alpha,
                  arma::vec ran_int_exp,
                  arma::vec theta, arma::vec theta2, arma::mat &bi_xz,
                  arma::mat &bi_zz, arma::mat &bi_p_mat, arma::mat bXXl,
                  arma::mat bXXnl, double soglia, double soglia2,
                  arma::mat &prop_per_bi, arma::mat &acc_bi_flag,
                  double mu_hs, int pp, int qq, int jj, int hered){
  
  int PP = XX.n_cols;
  int QQ = ZZ.n_cols;
  int JJ = SS.n_cols;
  int n = ZZ.n_rows;
  
  arma::vec loggamma_p(n);
  
  int ll, hh;
  
  ll = qq + jj * QQ;//indice per tau
  hh = qq + pp * QQ;//indice per bi e lambda
  
  //indici hereditariety
  int hh_theta = jj + pp * JJ;
  int ll_theta = jj + qq * JJ;
  
  arma::mat bXZ(PP, QQ);
  
  for(int qq2 = 0; qq2 < QQ; qq2++){
    for(int pp2 = 0; pp2 < PP; pp2++){
      int hh2 = qq2 + ((pp2) * QQ);
      bXZ(pp2, qq2) = bi_xz(jj, hh2);
    }
  }
  
  arma::mat bZZ(QQ, QQ);
  
  for(int qq2 = 0; qq2 < QQ; qq2++){
    for(int pp2 = 0; pp2 < QQ; pp2++){
      int hh2 = qq2 + ((pp2) * QQ);
      bZZ(pp2, qq2) = bi_zz(jj, hh2);
    }
  }
  
  //if hereditariety
  if(hered == 2){
    if((theta(hh_theta) == 0) | (theta2(ll_theta) == 0.0)){
      acc_bi_flag(jj, hh) = 0.0;
      bi_p_mat(pp, qq) = 0.0;
      
      loggamma.col(jj) = calculate_gamma_ri(XX, ZZ, XXl, XXT, alpha, ran_int_exp,
                   theta, theta2, bi_p_mat,
                   bXXl, bXXnl, soglia, bZZ, soglia2, jj);
      
      bi_xz(jj, hh) = bi_p_mat(pp, qq);
    } else {
      
      double lnu, ln_acp;
      
      double log_full_bi, log_full_bi_p;
      
      log_full_bi = 0.0;
      log_full_bi = log_full_bi - arma::sum(lgamma(exp(loggamma.col(jj)))) +
        arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)))
        + ::Rf_dnorm4(bi_xz(jj, hh), mu_hs, lambda_hs(jj, hh)*tau_hs(ll), 1);
      
      double propbipm = prop_per_bi(jj, hh);
      propbipm = adap_prop(propbipm, PP, JJ, QQ);
      
      bi_p_mat(pp, qq) += propbipm;
      
      loggamma.col(jj) = calculate_gamma_ri(XX, ZZ, XXl, XXT, alpha, ran_int_exp,
                   theta, theta2, bi_p_mat,
                   bXXl, bXXnl, soglia, bZZ, soglia2, jj);
      
      log_full_bi_p = 0.0;
      log_full_bi_p = log_full_bi_p - arma::sum(lgamma(exp(loggamma_p))) +
        arma::sum(exp(loggamma_p)%log(SS.col(jj)))
        + ::Rf_dnorm4(bi_p_mat(pp, qq), mu_hs, lambda_hs(jj, hh)*tau_hs(ll), 1);
      
      ln_acp = log_full_bi_p - log_full_bi;
      
      lnu = log(arma::randu());
      
      if(lnu < ln_acp){
        bi_xz(jj, hh) = bi_p_mat(pp, qq);
        acc_bi_flag(jj, hh) = 1;
        loggamma.col(jj)= loggamma_p;
      } else {
        bi_p_mat(pp, qq) -= propbipm;
      }//chiude ratio
    }
  }
  if(hered == 1){
    if((theta(hh_theta) == 0) & (theta2(ll_theta) == 0.0)){
      acc_bi_flag(jj, hh) = 0.0;
      bi_p_mat(pp, qq) = 0.0;
      
      loggamma.col(jj) = calculate_gamma_ri(XX, ZZ, XXl, XXT, alpha, ran_int_exp, theta, theta2, bi_p_mat,
                   bXXl, bXXnl, soglia, bZZ, soglia2, jj);
      
      bi_xz(jj, hh) = bi_p_mat(pp, qq);
    } else {
      
      double lnu, ln_acp;
      
      double log_full_bi, log_full_bi_p;
      
      log_full_bi = 0.0;
      log_full_bi = log_full_bi - arma::sum(lgamma(exp(loggamma.col(jj)))) +
        arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)))
        + ::Rf_dnorm4(bi_xz(jj, hh), mu_hs, lambda_hs(jj, hh)*tau_hs(ll), 1);
      
      double propbipm = prop_per_bi(jj, hh);
      propbipm = adap_prop(propbipm, PP, JJ, QQ);
      
      bi_p_mat(pp, qq) += propbipm;
      
      loggamma.col(jj) = calculate_gamma_ri(XX, ZZ, XXl, XXT, alpha, ran_int_exp, theta, theta2, bi_p_mat,
                   bXXl, bXXnl, soglia, bZZ, soglia2, jj);
      
      log_full_bi_p = 0.0;
      log_full_bi_p = log_full_bi_p - arma::sum(lgamma(exp(loggamma_p))) +
        arma::sum(exp(loggamma_p)%log(SS.col(jj)))
        + ::Rf_dnorm4(bi_p_mat(pp, qq), mu_hs, lambda_hs(jj, hh)*tau_hs(ll), 1);
      
      ln_acp = log_full_bi_p - log_full_bi;
      
      lnu = log(arma::randu());
      
      if(lnu < ln_acp){
        bi_xz(jj, hh) = bi_p_mat(pp, qq);
        acc_bi_flag(jj, hh) = 1;
        loggamma.col(jj)= loggamma_p;
      } else {
        bi_p_mat(pp, qq) -= propbipm;
      }//chiude ratio
    }
  }
  if(hered == 0){
    double lnu, ln_acp;
    
    double log_full_bi, log_full_bi_p;
    
    log_full_bi = 0.0;
    log_full_bi = log_full_bi - arma::sum(lgamma(exp(loggamma.col(jj)))) +
      arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)))
      + ::Rf_dnorm4(bi_xz(jj, hh), mu_hs, lambda_hs(jj, hh)*tau_hs(ll), 1);
    
    double propbipm = prop_per_bi(jj, hh);
    propbipm = adap_prop(propbipm, PP, JJ, QQ);
    
    bi_p_mat(pp, qq) += propbipm;
    
    loggamma.col(jj) = calculate_gamma_ri(XX, ZZ, XXl, XXT, alpha, ran_int_exp, theta, theta2, bi_p_mat,
                 bXXl, bXXnl, soglia, bZZ, soglia2, jj);
    
    log_full_bi_p = 0.0;
    log_full_bi_p = log_full_bi_p - arma::sum(lgamma(exp(loggamma_p))) +
      arma::sum(exp(loggamma_p)%log(SS.col(jj)))
      + ::Rf_dnorm4(bi_p_mat(pp, qq), mu_hs, lambda_hs(jj, hh)*tau_hs(ll), 1);
    
    ln_acp = log_full_bi_p - log_full_bi;
    
    lnu = log(arma::randu());
    
    if(lnu < ln_acp){
      bi_xz(jj, hh) = bi_p_mat(pp, qq);
      acc_bi_flag(jj, hh) = 1;
      loggamma.col(jj)= loggamma_p;
    } else {
      bi_p_mat(pp, qq) -= propbipm;
    }//chiude ratio
  }
}

void update_bi2_ri(arma::mat XX, arma::mat ZZ, arma::mat SS, arma::mat &loggamma,
                   arma::mat lambda_hs2, arma::vec tau_hs2, arma::vec &alpha,
                   arma::vec theta, arma::vec theta2, arma::mat &bi_xz,
                   arma::mat &bi_zz, arma::mat &bi_p_mat2, double soglia,
                   double soglia2, arma::mat &prop_per_bi2, arma::mat &acc_bi_flag2,
                   double mu_hs, int pp, int qq, int jj, int hered){
  
  int QQ = ZZ.n_cols;
  int JJ = SS.n_cols;
  int n = ZZ.n_rows;
  
  arma::vec loggamma_p(n);
  
  int ll, hh;
  ll = qq + jj * QQ;//indice per tau
  hh = qq + pp * QQ;//indice per bi e lambda QQ è il numero di colonne
  
  //indici hereditariety
  int hh_theta = jj + pp * JJ;
  int ll_theta = jj + qq * JJ;
  
  arma::mat bZZ(QQ, QQ);
  
  for(int qq2 = 0; qq2 < QQ; qq2++){
    for(int pp2 = 0; pp2 < QQ; pp2++){
      int hh2 = qq2 + ((pp2) * QQ);
      bZZ(pp2, qq2) = bi_zz(jj, hh2);
    }
  }
  
  //if hereditariety
  if(hered == 2){
    if((theta2(hh_theta) == 0) | (theta2(ll_theta) == 0.0)){
      acc_bi_flag2(jj, hh) = 0;
      bi_p_mat2(pp, qq) = 0.0;
      
      loggamma.col(jj) = loggamma.col(jj) - trans(threshold_mat(bZZ(pp, qq) *
        ZZ.col(qq).t(), soglia2)) % ZZ.col(qq) + trans(threshold_mat(bi_p_mat2(pp, qq) *
        ZZ.col(qq).t(), soglia2)) % ZZ.col(qq);
      
      bi_zz(jj, hh) = bi_p_mat2(pp, qq);
    } else {
      double lnu, ln_acp;
      
      double log_full_bi, log_full_bi_p;
      
      log_full_bi = 0.0;
      log_full_bi = log_full_bi - arma::sum(lgamma(exp(loggamma.col(jj)))) +
        arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)))
        + ::Rf_dnorm4(bi_zz(jj, hh), mu_hs, lambda_hs2(jj, hh)*tau_hs2(ll), 1);
      
      double propbipm = prop_per_bi2(jj, hh);
      propbipm = adap_prop(propbipm, QQ, JJ, QQ);
      bi_p_mat2(pp, qq) += propbipm;
      
      loggamma_p = loggamma.col(jj) - trans(threshold_mat(bZZ(pp, qq) *
        ZZ.col(qq).t(), soglia2)) % ZZ.col(qq) + trans(threshold_mat(bi_p_mat2(pp, qq) *
        ZZ.col(qq).t(), soglia2)) % ZZ.col(qq);
      
      log_full_bi_p = 0.0;
      log_full_bi_p = log_full_bi_p - arma::sum(lgamma(exp(loggamma_p))) +
        arma::sum(exp(loggamma_p)%log(SS.col(jj)))
        + ::Rf_dnorm4(bi_p_mat2(pp, qq), mu_hs, lambda_hs2(jj, hh)*tau_hs2(ll), 1);
      
      ln_acp = log_full_bi_p - log_full_bi;
      
      lnu = log(arma::randu());
      
      if(lnu < ln_acp){
        bi_zz(jj, hh) = bi_p_mat2(pp, qq);
        acc_bi_flag2(jj, hh) = 1;
        loggamma.col(jj)= loggamma_p;
      } else {
        bi_p_mat2(pp, qq) -= propbipm;
      }//chiude ratio
    }
  }
  if(hered == 1){
    if((theta2(hh_theta) == 0) & (theta2(ll_theta) == 0.0)){
      acc_bi_flag2(jj, hh) = 0;
      bi_p_mat2(pp, qq) = 0.0;
      
      loggamma.col(jj) = loggamma.col(jj) - trans(threshold_mat(bZZ(pp, qq) *
        ZZ.col(qq).t(), soglia2)) % ZZ.col(qq) + trans(threshold_mat(bi_p_mat2(pp, qq) *
        ZZ.col(qq).t(), soglia2)) % ZZ.col(qq);
      
      bi_zz(jj, hh) = bi_p_mat2(pp, qq);
    } else {
      double lnu, ln_acp;
      
      double log_full_bi, log_full_bi_p;
      
      log_full_bi = 0.0;
      log_full_bi = log_full_bi - arma::sum(lgamma(exp(loggamma.col(jj)))) +
        arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)))
        + ::Rf_dnorm4(bi_zz(jj, hh), mu_hs, lambda_hs2(jj, hh)*tau_hs2(ll), 1);
      
      double propbipm = prop_per_bi2(jj, hh);
      propbipm = adap_prop(propbipm, QQ, JJ, QQ);
      bi_p_mat2(pp, qq) += propbipm;
      
      loggamma_p = loggamma.col(jj) - trans(threshold_mat(bZZ(pp, qq) *
        ZZ.col(qq).t(), soglia2)) % ZZ.col(qq) + trans(threshold_mat(bi_p_mat2(pp, qq) *
        ZZ.col(qq).t(), soglia2)) % ZZ.col(qq);
      
      log_full_bi_p = 0.0;
      log_full_bi_p = log_full_bi_p - arma::sum(lgamma(exp(loggamma_p))) +
        arma::sum(exp(loggamma_p)%log(SS.col(jj)))
        + ::Rf_dnorm4(bi_p_mat2(pp, qq), mu_hs, lambda_hs2(jj, hh)*tau_hs2(ll), 1);
      
      ln_acp = log_full_bi_p - log_full_bi;
      
      lnu = log(arma::randu());
      
      if(lnu < ln_acp){
        bi_zz(jj, hh) = bi_p_mat2(pp, qq);
        acc_bi_flag2(jj, hh) = 1;
        loggamma.col(jj)= loggamma_p;
      } else {
        bi_p_mat2(pp, qq) -= propbipm;
      }//chiude ratio
    }
  }
  if(hered == 0){
    double lnu, ln_acp;
    
    double log_full_bi, log_full_bi_p;
    
    log_full_bi = 0.0;
    log_full_bi = log_full_bi - arma::sum(lgamma(exp(loggamma.col(jj)))) +
      arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)))
      + ::Rf_dnorm4(bi_zz(jj, hh), mu_hs, lambda_hs2(jj, hh)*tau_hs2(ll), 1);
    
    double propbipm = prop_per_bi2(jj, hh);
    propbipm = adap_prop(propbipm, QQ, JJ, QQ);
    bi_p_mat2(pp, qq) += propbipm;
    
    loggamma_p = loggamma.col(jj) - trans(threshold_mat(bZZ(pp, qq) *
      ZZ.col(qq).t(), soglia2)) % ZZ.col(qq) + trans(threshold_mat(bi_p_mat2(pp, qq) *
      ZZ.col(qq).t(), soglia2)) % ZZ.col(qq);
    
    log_full_bi_p = 0.0;
    log_full_bi_p = log_full_bi_p - arma::sum(lgamma(exp(loggamma_p))) +
      arma::sum(exp(loggamma_p)%log(SS.col(jj)))
      + ::Rf_dnorm4(bi_p_mat2(pp, qq), mu_hs, lambda_hs2(jj, hh)*tau_hs2(ll), 1);
    
    ln_acp = log_full_bi_p - log_full_bi;
    
    lnu = log(arma::randu());
    
    if(lnu < ln_acp){
      bi_zz(jj, hh) = bi_p_mat2(pp, qq);
      acc_bi_flag2(jj, hh) = 1;
      loggamma.col(jj)= loggamma_p;
    } else {
      bi_p_mat2(pp, qq) -= propbipm;
    }//chiude ratio
  }
}


arma::mat update_alpha_xx_ri(arma::mat XX, arma::mat ZZ, arma::mat XXl, arma::mat XXT,
                             arma::mat SS, arma::mat &loggamma, arma::vec &alpha, arma::vec ran_int_exp,
                             arma::vec theta,
                             arma::vec theta2, arma::mat &bi_xz, double soglia,
                             arma::mat &bi_zz, double soglia2, arma::vec dj, arma::mat alpha_lin,
                             arma::mat alpha_nl, arma::mat tausq_lin, arma::mat tausq_nl,
                             arma::mat gamma_lin, arma::mat gamma_nl, arma::mat ksi_lin, arma::mat ksi_nl,
                             bool linear, int jj, int hered, bool constraint){
  
  if(linear == true){
    int PP = XX.n_cols;
    int QQ = ZZ.n_cols;
    int KK = XXT.n_cols;
    int JJ = SS.n_cols;
    int n = ZZ.n_rows;
    
    arma::vec loggamma_p(n);
    
    int hh_theta, ll_theta;//indici x hereditariety
    
    arma::mat bXZ(PP, QQ);
    
    for(int qq2 = 0; qq2 < QQ; qq2++){
      for(int pp2 = 0; pp2 < PP; pp2++){
        int hh2 = qq2 + ((pp2) * QQ);
        bXZ(pp2, qq2) = bi_xz(jj, hh2);
      }
    }
    
    arma::mat bZZ(QQ, QQ);
    
    for(int qq2 = 0; qq2 < QQ; qq2++){
      for(int pp2 = 0; pp2 < QQ; pp2++){
        int hh2 = qq2 + ((pp2) * QQ);
        bZZ(pp2, qq2) = bi_zz(jj, hh2);
      }
    }
    
    double mhr_num, mhr_den, mhr_diff, mhr_u;
    
    arma::mat alpha_p(PP, PP);
    alpha_p.fill(0.0);
    
    for(int pp2 = 0; pp2 < PP; pp2++){
      for(int qq2 = pp2 + 1; qq2 < PP; qq2++){
        hh_theta = jj + pp2 * JJ;
        ll_theta = jj + qq2 * JJ;
        //cambia |/& a seconda dell'assunzione sull ereditarietà
        if(hered == 2){
          if((theta(hh_theta) == 0.0) | (theta(ll_theta) == 0.0)){
            alpha_lin(pp2, qq2) = 0.0;
            arma::mat beta_lin_j(PP, PP, arma::fill::zeros);
            arma::mat beta_nl_j(PP, KK, arma::fill::zeros);
            
            beta_lin_j = alpha_ksi_prod(alpha_lin, ksi_lin, dj, true);
            beta_nl_j = alpha_ksi_prod(alpha_nl, ksi_nl, dj, false);
            loggamma_p = calculate_gamma_ri(XX, ZZ, XXl, XXT, alpha, ran_int_exp, theta, theta2, bXZ,
                                            beta_lin_j, beta_nl_j, soglia, bZZ, soglia2, jj);
            
          } else{
            alpha_p(pp2, qq2) = alpha_lin(pp2, qq2) + ::Rf_rnorm(0, .01);
            
            mhr_den = 0.0;
            mhr_den = mhr_den - arma::sum(lgamma(exp(loggamma.col(jj)))) +
              arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)));
            
            double gammatau;
            gammatau = pow(tausq_lin(pp2, qq2), .5) * gamma_lin(pp2, qq2);
            mhr_den += ::Rf_dnorm4(alpha_lin(pp2, qq2), 0.0, gammatau, 1);
            
            arma::mat beta_lin_j_p(PP, PP, arma::fill::zeros);
            arma::mat beta_nl_j(PP, KK, arma::fill::zeros);
            
            beta_lin_j_p = alpha_ksi_prod(alpha_p, ksi_lin, dj, true);
            beta_nl_j = alpha_ksi_prod(alpha_nl, ksi_nl, dj, false);
            loggamma_p = calculate_gamma_ri(XX, ZZ, XXl, XXT, alpha, ran_int_exp, theta, theta2, bXZ,
                                            beta_lin_j_p, beta_nl_j, soglia, bZZ, soglia2, jj);
            
            mhr_num = 0.0;
            mhr_num = mhr_num - arma::sum(lgamma(exp(loggamma_p))) +
              arma::sum(exp(loggamma_p)%log(SS.col(jj)));
            
            gammatau = pow(tausq_lin(pp2, qq2), .5) * gamma_lin(pp2, qq2);
            mhr_num += ::Rf_dnorm4(alpha_p(pp2, qq2), 0.0, gammatau, 1);
            
            mhr_diff = mhr_num - mhr_den;
            
            mhr_u = log(arma::randu());
            
            if (mhr_u < mhr_diff) {
              alpha_lin(pp2, qq2) = alpha_p(pp2, qq2);
              loggamma.col(jj) = loggamma_p;
            }//chiude ratio
          }//chiude l else dell ereditarietà
        }
        if(hered == 1){
          if((theta(hh_theta) == 0.0) & (theta(ll_theta) == 0.0)){
            alpha_lin(pp2, qq2) = 0.0;
            arma::mat beta_lin_j(PP, PP, arma::fill::zeros);
            arma::mat beta_nl_j(PP, KK, arma::fill::zeros);
            
            beta_lin_j = alpha_ksi_prod(alpha_lin, ksi_lin, dj, true);
            beta_nl_j = alpha_ksi_prod(alpha_nl, ksi_nl, dj, false);
            loggamma_p = calculate_gamma_ri(XX, ZZ, XXl, XXT, alpha, ran_int_exp, theta, theta2, bXZ,
                                            beta_lin_j, beta_nl_j, soglia, bZZ, soglia2, jj);
            
          } else{
            alpha_p(pp2, qq2) = alpha_lin(pp2, qq2) + ::Rf_rnorm(0, .1);
            
            mhr_den = 0.0;
            mhr_den = mhr_den - arma::sum(lgamma(exp(loggamma.col(jj)))) +
              arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)));
            
            double gammatau;
            gammatau = pow(tausq_lin(pp2, qq2), .5) * gamma_lin(pp2, qq2);
            mhr_den += ::Rf_dnorm4(alpha_lin(pp2, qq2), 0.0, gammatau, 1);
            
            arma::mat beta_lin_j_p(PP, PP, arma::fill::zeros);
            arma::mat beta_nl_j(PP, KK, arma::fill::zeros);
            
            beta_lin_j_p = alpha_ksi_prod(alpha_p, ksi_lin, dj, true);
            beta_nl_j = alpha_ksi_prod(alpha_nl, ksi_nl, dj, false);
            loggamma_p = calculate_gamma_ri(XX, ZZ, XXl, XXT, alpha, ran_int_exp, theta, theta2, bXZ,
                                            beta_lin_j_p, beta_nl_j, soglia, bZZ, soglia2, jj);
            
            mhr_num = 0.0;
            mhr_num = mhr_num - arma::sum(lgamma(exp(loggamma_p))) +
              arma::sum(exp(loggamma_p)%log(SS.col(jj)));
            
            gammatau = pow(tausq_lin(pp2, qq2), .5) * gamma_lin(pp2, qq2);
            mhr_num += ::Rf_dnorm4(alpha_p(pp2, qq2), 0.0, gammatau, 1);
            
            mhr_diff = mhr_num - mhr_den;
            
            mhr_u = log(arma::randu());
            
            if (mhr_u < mhr_diff) {
              alpha_lin(pp2, qq2) = alpha_p(pp2, qq2);
              loggamma.col(jj) = loggamma_p;
            }//chiude ratio
          }//chiude l else dell ereditarietà
        }
        if(hered == 0){
          alpha_p(pp2, qq2) = alpha_lin(pp2, qq2) + ::Rf_rnorm(0, .1);
          
          mhr_den = 0.0;
          mhr_den = mhr_den - arma::sum(lgamma(exp(loggamma.col(jj)))) +
            arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)));
          
          double gammatau;
          gammatau = pow(tausq_lin(pp2, qq2), .5) * gamma_lin(pp2, qq2);
          mhr_den += ::Rf_dnorm4(alpha_lin(pp2, qq2), 0.0, gammatau, 1);
          
          arma::mat beta_lin_j_p(PP, PP, arma::fill::zeros);
          arma::mat beta_nl_j(PP, KK, arma::fill::zeros);
          
          beta_lin_j_p = alpha_ksi_prod(alpha_p, ksi_lin, dj, true);
          beta_nl_j = alpha_ksi_prod(alpha_nl, ksi_nl, dj, false);
          loggamma_p = calculate_gamma_ri(XX, ZZ, XXl, XXT, alpha, ran_int_exp, theta, theta2, bXZ,
                                          beta_lin_j_p, beta_nl_j, soglia, bZZ, soglia2, jj);
          
          mhr_num = 0.0;
          mhr_num = mhr_num - arma::sum(lgamma(exp(loggamma_p))) +
            arma::sum(exp(loggamma_p)%log(SS.col(jj)));
          
          gammatau = pow(tausq_lin(pp2, qq2), .5) * gamma_lin(pp2, qq2);
          mhr_num += ::Rf_dnorm4(alpha_p(pp2, qq2), 0.0, gammatau, 1);
          
          mhr_diff = mhr_num - mhr_den;
          
          mhr_u = log(arma::randu());
          
          if (mhr_u < mhr_diff) {
            alpha_lin(pp2, qq2) = alpha_p(pp2, qq2);
            loggamma.col(jj) = loggamma_p;
          }//chiude ratio
        }//chiude l else dell ereditarietà
      }
    }//
    return alpha_lin;
  } else {
    int PP = XX.n_cols;
    int KK = XXT.n_cols;
    int QQ = ZZ.n_cols;
    int JJ = SS.n_cols;
    int n = ZZ.n_rows;
    
    arma::vec loggamma_p(n);
    
    arma::mat bXZ(PP, QQ);
    
    for(int qq2 = 0; qq2 < QQ; qq2++){
      for(int pp2 = 0; pp2 < PP; pp2++){
        int hh2 = qq2 + ((pp2) * QQ);
        bXZ(pp2, qq2) = bi_xz(jj, hh2);
      }
    }
    
    arma::mat bZZ(QQ, QQ);
    
    for(int qq2 = 0; qq2 < QQ; qq2++){
      for(int pp2 = 0; pp2 < QQ; pp2++){
        int hh2 = qq2 + ((pp2) * QQ);
        bZZ(pp2, qq2) = bi_zz(jj, hh2);
      }
    }
    
    double mhr_num, mhr_den, mhr_diff, mhr_u;
    
    arma::mat alpha_p(PP, PP);
    alpha_p.fill(0.0);
    
    //indici x hereditariety
    int hh_theta;
    int ll_theta;
    
    for(int pp2 = 0; pp2 < PP; pp2++){
      for(int qq2 = pp2 + 1; qq2 < PP; qq2++){
        hh_theta = jj + pp2 * JJ;
        ll_theta = jj + qq2 * JJ;
        //unica condizione sul constraint
        if(constraint == true){
          if(gamma_lin(pp2, qq2) != 1.0){
            alpha_nl(pp2, qq2) = 0.0;
            arma::mat beta_lin_j(PP, PP, arma::fill::zeros);
            arma::mat beta_nl_j(PP, KK, arma::fill::zeros);
            
            beta_lin_j = alpha_ksi_prod(alpha_lin, ksi_lin, dj, true);
            beta_nl_j = alpha_ksi_prod(alpha_nl, ksi_nl, dj, false);
            loggamma_p = calculate_gamma_ri(XX, ZZ, XXl, XXT, alpha, ran_int_exp, theta, theta2, bXZ,
                                            beta_lin_j, beta_nl_j, soglia, bZZ, soglia2, jj);
          } else{
            alpha_p(pp2, qq2) = alpha_nl(pp2, qq2) + ::Rf_rnorm(0, .1);
            
            mhr_den = 0.0;
            mhr_den = mhr_den - arma::sum(lgamma(exp(loggamma.col(jj)))) +
              arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)));
            
            double gammatau;
            gammatau = pow(tausq_nl(pp2, qq2), .5) * gamma_nl(pp2, qq2);
            mhr_den += ::Rf_dnorm4(alpha_nl(pp2, qq2), 0.0, gammatau, 1);
            arma::mat beta_lin_j(PP, PP, arma::fill::zeros);
            arma::mat beta_nl_j_p(PP, KK, arma::fill::zeros);
            
            beta_lin_j = alpha_ksi_prod(alpha_lin, ksi_lin, dj, true);
            beta_nl_j_p = alpha_ksi_prod(alpha_p, ksi_nl, dj, false);
            loggamma_p = calculate_gamma_ri(XX, ZZ, XXl, XXT, alpha, ran_int_exp, theta, theta2, bXZ,
                                            beta_lin_j, beta_nl_j_p, soglia, bZZ, soglia2, jj);
            
            mhr_num = 0.0;
            mhr_num = mhr_num - arma::sum(lgamma(exp(loggamma_p))) +
              arma::sum(exp(loggamma_p)%log(SS.col(jj)));
            
            gammatau = pow(tausq_nl(pp2, qq2), .5) * gamma_nl(pp2, qq2);
            mhr_num += ::Rf_dnorm4(alpha_p(pp2, qq2), 0.0, gammatau, 1);
            
            
            mhr_diff = mhr_num - mhr_den;
            
            mhr_u = log(arma::randu());
            
            if (mhr_u < mhr_diff) {
              alpha_nl(pp2, qq2) = alpha_p(pp2, qq2);
              loggamma.col(jj) = loggamma_p;
            }//chiude ratio
          }//chiude else hereditariety
        } else {//chiude if(constraint == true)
          if(gamma_lin(pp2, qq2) != 1.0) {
            alpha_nl(pp2, qq2) = 0.0;
            arma::mat beta_lin_j(PP, PP, arma::fill::zeros);
            arma::mat beta_nl_j(PP, KK, arma::fill::zeros);
            
            beta_lin_j = alpha_ksi_prod(alpha_lin, ksi_lin, dj, true);
            beta_nl_j = alpha_ksi_prod(alpha_nl, ksi_nl, dj, false);
            loggamma_p = calculate_gamma_ri(XX, ZZ, XXl, XXT, alpha, ran_int_exp, theta, theta2, bXZ,
                                            beta_lin_j, beta_nl_j, soglia, bZZ, soglia2, jj);
          } else{
            alpha_p(pp2, qq2) = alpha_nl(pp2, qq2) + ::Rf_rnorm(0, .1);
            
            
            mhr_den = 0.0;
            mhr_den = mhr_den - arma::sum(lgamma(exp(loggamma.col(jj)))) +
              arma::sum(exp(loggamma.col(jj))%log(SS.col(jj)));
            
            double gammatau;
            gammatau = pow(tausq_nl(pp2, qq2), .5) * gamma_nl(pp2, qq2);
            mhr_den += ::Rf_dnorm4(alpha_nl(pp2, qq2), 0.0, gammatau, 1);
            
            arma::mat beta_lin_j(PP, PP, arma::fill::zeros);
            arma::mat beta_nl_j_p(PP, KK, arma::fill::zeros);
            
            beta_lin_j = alpha_ksi_prod(alpha_lin, ksi_lin, dj, true);
            beta_nl_j_p = alpha_ksi_prod(alpha_p, ksi_nl, dj, false);
            loggamma_p = calculate_gamma_ri(XX, ZZ, XXl, XXT, alpha, ran_int_exp, theta, theta2, bXZ,
                                            beta_lin_j, beta_nl_j_p, soglia, bZZ, soglia2, jj);
            
            mhr_num = 0.0;
            mhr_num = mhr_num - arma::sum(lgamma(exp(loggamma_p))) +
              arma::sum(exp(loggamma_p)%log(SS.col(jj)));
            
            gammatau = pow(tausq_nl(pp2, qq2), .5) * gamma_nl(pp2, qq2);
            mhr_num += ::Rf_dnorm4(alpha_p(pp2, qq2), 0.0, gammatau, 1);
            
            mhr_diff = mhr_num - mhr_den;
            
            mhr_u = log(arma::randu());
            
            if (mhr_u < mhr_diff) {
              alpha_nl(pp2, qq2) = alpha_p(pp2, qq2);
              loggamma.col(jj) = loggamma_p;
            }//chiude ratio
          }//chiude else gammalin =! 1.0
        }//chiude else constraint
      }//chiude qq
    }//chiude pp
    return alpha_nl;
  }//chiude else
}//chiude la funzione update_alpha_xx()

arma::mat update_ksi_ri(arma::mat XX, arma::mat ZZ, arma:: mat XXl, arma::mat XXT, arma::mat SS,
                        arma::mat &loggamma, arma::vec &alpha, arma::vec ran_int_exp, arma::vec theta, arma::vec theta2,
                        arma::mat &bi_xz, double soglia, arma::mat &bi_zz, double soglia2,
                        arma::vec dj, arma::mat alpha_lin,
                        arma::mat alpha_nl, arma::mat tausq, arma::mat tausq_nl, arma::mat gamma_lin,
                        arma::mat gamma_nl, arma::mat ksi_lin, arma::mat ksi_nl, arma::mat m_lin,
                        arma::mat m_nl, bool linear, int jj){
  
  if(linear == true){
    int PP = XX.n_cols;
    int KK = XXT.n_cols;
    int QQ = ZZ.n_cols;
    int n = ZZ.n_rows;
    
    arma::vec loggamma_p(n);
    
    arma::mat bXZ(PP, QQ);
    
    for(int qq2 = 0; qq2 < QQ; qq2++){
      for(int pp2 = 0; pp2 < PP; pp2++){
        int hh2 = qq2 + ((pp2) * QQ);
        bXZ(pp2, qq2) = bi_xz(jj, hh2);
      }
    }
    
    arma::mat bZZ(QQ, QQ);
    
    for(int qq2 = 0; qq2 < QQ; qq2++){
      for(int pp2 = 0; pp2 < QQ; pp2++){
        int hh2 = qq2 + ((pp2) * QQ);
        bZZ(pp2, qq2) = bi_zz(jj, hh2);
      }
    }
    
    double mhr_num, mhr_den, mhr_diff, mhr_u;
    
    arma::mat ksi_p(PP, PP, arma::fill::zeros);
    
    for(int pp2 = 0; pp2 < PP; pp2++){
      for(int qq2 = pp2 + 1; qq2 < PP; qq2++){
        ksi_p(pp2, qq2) = ksi_lin(pp2, qq2) + ::Rf_rnorm(0, .1);
        
        mhr_den = 0.0;
        mhr_den = mhr_den - arma::sum(lgamma(exp(loggamma.col(jj)))) +
          arma::sum(exp(loggamma.col(jj)) % log(SS.col(jj)));
        
        mhr_den += ::Rf_dnorm4(ksi_lin(pp2, qq2), m_lin(pp2, qq2), 1.0, 1);
        
        arma::mat beta_lin_j_p(PP, PP, arma::fill::zeros);
        arma::mat beta_nl_j(PP, KK, arma::fill::zeros);
        beta_lin_j_p = alpha_ksi_prod(alpha_lin, ksi_p, dj, true);
        beta_nl_j = alpha_ksi_prod(alpha_nl, ksi_nl, dj, false);
        
        loggamma_p = calculate_gamma_ri(XX, ZZ, XXl, XXT, alpha, ran_int_exp, theta, theta2, bXZ,
                                        beta_lin_j_p, beta_nl_j, soglia, bZZ, soglia2, jj);
        
        mhr_num = 0.0;
        mhr_num = mhr_num - arma::sum(lgamma(exp(loggamma_p))) +
          arma::sum(exp(loggamma_p)%log(SS.col(jj)));
        
        mhr_num += ::Rf_dnorm4(ksi_p(pp2, qq2), m_lin(pp2, qq2), 1.0, 1);
        
        mhr_diff = mhr_num - mhr_den;
        
        mhr_u = log(arma::randu());
        
        if (mhr_u < mhr_diff) {
          ksi_lin(pp2, qq2) = ksi_p(pp2, qq2);
          loggamma.col(jj) = loggamma_p;
        }//chiude ratio
      }
    }
    return ksi_lin;
  } else {//chiude (linear == true)
    int PP = XX.n_cols;
    int QQ = ZZ.n_cols;
    int KK = XXT.n_cols;
    int n = ZZ.n_rows;
    
    arma::vec loggamma_p(n);
    
    arma::mat bXZ(PP, QQ);
    
    for(int qq2 = 0; qq2 < QQ; qq2++){
      for(int pp2 = 0; pp2 < PP; pp2++){
        int hh2 = qq2 + ((pp2) * QQ);
        bXZ(pp2, qq2) = bi_xz(jj, hh2);
      }
    }
    
    arma::mat bZZ(QQ, QQ);
    
    for(int qq2 = 0; qq2 < QQ; qq2++){
      for(int pp2 = 0; pp2 < QQ; pp2++){
        int hh2 = qq2 + ((pp2) * QQ);
        bZZ(pp2, qq2) = bi_zz(jj, hh2);
      }
    }
    
    double mhr_num, mhr_den, mhr_diff, mhr_u;
    
    arma::mat ksi_p(PP, KK, arma::fill::zeros);
    
    arma::vec cdj(PP-1);
    cdj = cumsum(dj);
    double rs = 0.0;
    for(int pp2 = 0; pp2 < PP; pp2++){
      for(int qq2 = pp2 + 1; qq2 < PP; qq2++){
        rs = ::Rf_rnorm(0, .1);
        for(int ii = cdj(qq2-1); ii < cdj(qq2); ii++){
          ksi_p(pp2, ii) = ksi_nl(pp2, ii) + rs;
        }
      }
    }
    
    mhr_den = 0.0;
    mhr_den = mhr_den - arma::sum(lgamma(exp(loggamma.col(jj)))) +
      arma::sum(exp(loggamma.col(jj)) % log(SS.col(jj)));
    
    for(int pp2 = 0; pp2 < PP; pp2++){
      for(int qq2 = pp2 + 1; qq2 < PP; qq2++){
        for(int ii = cdj(qq2-1); ii < cdj(qq2); ii++){
          mhr_den += ::Rf_dnorm4(ksi_nl(pp2, ii), m_nl(pp2, ii), 1.0, 1);
        }
      }
    }
    
    arma::mat beta_lin_j(PP, PP, arma::fill::zeros);
    arma::mat beta_nl_j_p(PP, KK, arma::fill::zeros);
    beta_lin_j = alpha_ksi_prod(alpha_lin, ksi_lin, dj, true);
    beta_nl_j_p = alpha_ksi_prod(alpha_nl, ksi_p, dj, false);
    loggamma_p = calculate_gamma_ri(XX, ZZ, XXl, XXT, alpha, ran_int_exp, theta, theta2, bXZ,
                                    beta_lin_j, beta_nl_j_p, soglia, bZZ, soglia2, jj);
    
    mhr_num = 0.0;
    mhr_num = mhr_num - arma::sum(lgamma(exp(loggamma_p))) +
      arma::sum(exp(loggamma_p)%log(SS.col(jj)));
    
    for(int pp2 = 0; pp2 < PP; pp2++){
      for(int qq2 = pp2 + 1; qq2 < PP; qq2++){
        for(int ii = cdj(qq2-1); ii < cdj(qq2); ii++){
          mhr_num += ::Rf_dnorm4(ksi_p(pp2, ii), m_nl(pp2, ii), 1.0, 1);
        }
      }
    }
    
    mhr_diff = mhr_num - mhr_den;
    
    mhr_u = log(arma::randu());
    
    if (mhr_u < mhr_diff) {
      ksi_nl = ksi_p;
      loggamma.col(jj) = loggamma_p;
    }
    return ksi_nl;
  }
}//chiude funzione update_ksi()

void update_rand_int(arma::mat XX, arma::mat ZZ, arma::mat SS, arma::mat &loggamma, arma::vec &ran_int_sing,
                     arma::vec &ran_int_exp, arma::vec grplabel, int ngroup, arma::vec sig_int){
  
  int n = SS.n_rows;
  int JJ = SS.n_cols;
  
  arma::vec prop_ran_int_sing = ran_int_sing;
  arma::vec prop_ran_int_exp = ran_int_exp;
  arma::mat loggamma_p(n, JJ);
  
  double lnu, ln_acp;
  
  // Prepare the current and proposed full conditional values
  double log_full_int, log_full_int_p;
  
  // Calculate the full conditional for the current value
  log_full_int = 0.0;
  for(int jj = 0; jj < JJ; jj++){
    for(int ii = 0; ii < n; ii++){
      log_full_int = log_full_int - lgamma(exp(loggamma(ii, jj))) +
        (exp(loggamma(ii, jj)) * log(SS(ii, jj)));// +
      //::Rf_dnorm4(ran_int_exp(ii), 0.0, pow(sig_int(grplabel(ii)), .5), 1);
    }
  }
  
  for(int ii = 0; ii < n; ii++){
    log_full_int = log_full_int  +
      ::Rf_dnorm4(ran_int_exp(ii), 0.0, pow(sig_int(grplabel(ii)), .5), 1);
  }
  
  // Propose a new value for alpha[jj] using a random walk proposal centered on
  // the current value of alpha[jj]
  for(int gg = 0; gg < ngroup; gg++){
    prop_ran_int_sing(gg) = ran_int_sing(gg) + ::Rf_rnorm(0.0, .01);
    for(int ii = 0; ii < n; ii++){
      if(grplabel(ii) == (gg)){
        prop_ran_int_exp(ii) = prop_ran_int_sing(gg);
      }
    }
  }
  
  for(int jj = 0; jj < JJ; jj++){
    for(int ii = 0; ii < n; ii++){
      loggamma_p(ii, jj) = loggamma(ii, jj) - ran_int_exp(ii) + prop_ran_int_exp(ii);
    }
  }
  log_full_int_p = 0.0;
  for(int jj = 0; jj < JJ; jj++){
    for(int ii = 0; ii < n; ii++){
      log_full_int_p = log_full_int_p - lgamma(exp(loggamma_p(ii, jj))) +
        (exp(loggamma_p(ii,jj)) * log(SS(ii,jj)));// +
      //::Rf_dnorm4(prop_ran_int_exp(ii), 0, pow(sig_int(grplabel(ii)), .5), 1);
    }
  }
  for(int ii = 0; ii < n; ii++){
    log_full_int_p = log_full_int_p  +
      ::Rf_dnorm4(prop_ran_int_exp(ii), 0.0, pow(sig_int(grplabel(ii)), .5), 1);
  }
  ln_acp = log_full_int_p - log_full_int;
  lnu = log(arma::randu());
  if(lnu < ln_acp){
    //std::cout<<"accepted"<<std::endl;
    for(int gg = 0; gg < ngroup; gg++){
      ran_int_sing(gg) = prop_ran_int_sing(gg);
      for(int ii = 0; ii < n; ii++){
        if(grplabel(ii) == gg){
          ran_int_exp(ii) = ran_int_sing(gg);
        }
      }
      loggamma = loggamma_p;
    }
  }
} // Close function

//[[Rcpp::export]]
Rcpp::List sampler_randint(arma::mat YY, arma::mat XX, arma::mat ZZ, arma::mat XXl, arma::mat XXT,
                           arma::vec dj, int Niter, int burn, int thin, arma::vec hyper_theta_x,
                           arma::vec hyper_theta_z, arma::vec penmig_lin, arma::vec penmig_nl,
                           arma::vec tx, arma::vec tz, arma::vec prior_int, bool upsv,
                           int hereditariety, bool conlnl,
                           arma::vec grplabel, int ngroup, arma::vec theta_init,
                           arma::vec theta_init2){
  
  int Eff = (Niter-burn)/thin;
  int JJ = YY.n_cols;
  int PP = XX.n_cols;
  int KK = XXT.n_cols;
  int QQ = ZZ.n_cols;
  int n = XX.n_rows;
  
  int hh, ii, pp, jj, qq, ll, s;
  
  bool vincolo = conlnl;
  
  // adaptive proposal values for \thetaX & \thetaZ
  double last_mean;
  double last_var;
  double last_mean2;
  double last_var2;
  double last_meani;
  double last_vari;
  double last_meani2;
  double last_vari2;
  
  //MH acc ratios
  Rcpp::NumericVector acc_theta_flag0(PP*JJ);
  arma::vec acc_theta_flag(acc_theta_flag0.begin(), PP*JJ, false);
  acc_theta_flag.fill(0);
  
  Rcpp::NumericVector acc_theta_flag20(QQ*JJ);
  arma::vec acc_theta_flag2(acc_theta_flag20.begin(), QQ*JJ, false);
  acc_theta_flag2.fill(0);
  
  arma::vec acc_theta(PP*JJ);
  acc_theta.fill(0);
  
  arma::vec acc_theta2(QQ*JJ);
  acc_theta2.fill(0);
  
  Rcpp::NumericVector acc_alpha_flag0(JJ);
  arma::vec acc_alpha_flag(acc_alpha_flag0.begin(), JJ, false);
  acc_alpha_flag.fill(0);
  
  Rcpp::NumericVector acc_soglia_flag0(2);
  arma::vec acc_soglia_flag(acc_soglia_flag0.begin(), 2, false);
  acc_soglia_flag.fill(0);
  
  arma::vec acc_soglia(2);
  acc_soglia.fill(0);
  
  arma::vec acc_alpha(JJ);
  acc_alpha.fill(0);
  
  Rcpp::NumericMatrix acc_bi_flag0(JJ, PP*QQ);
  arma::mat acc_bi_flag(acc_bi_flag0.begin(), JJ, PP*QQ, false);
  acc_bi_flag.fill(0);
  
  Rcpp::NumericMatrix acc_bi_flag20(JJ, QQ*QQ);
  arma::mat acc_bi_flag2(acc_bi_flag20.begin(), JJ, QQ*QQ, false);
  acc_bi_flag2.fill(0);
  
  arma::mat acc_bi(JJ, PP*QQ);
  acc_bi.fill(0);
  
  arma::mat acc_bi2(JJ, QQ*QQ);
  acc_bi2.fill(0);
  
  //keep track of mean & var \tehtaX & \thetaZ
  arma::vec curr_mean(PP * JJ);
  curr_mean.fill(0.0);
  
  arma::vec curr_var(PP * JJ);
  curr_var.fill(0.5);
  
  arma::vec curr_mean2(QQ*JJ);
  curr_mean2.fill(0.0);
  
  arma::vec curr_var2(QQ*JJ);
  curr_var2.fill(0.5);
  
  arma::mat curr_meani(JJ, QQ*PP);
  curr_meani.fill(0.0);
  
  arma::mat curr_vari(JJ, QQ*PP);
  curr_vari.fill(.5);
  
  arma::mat curr_meani2(JJ, QQ * QQ);
  curr_meani2.fill(0.0);
  
  arma::mat curr_vari2(JJ, QQ * QQ);
  curr_vari2.fill(.5);
  
  //initialization for \thetaX & \thetaZ
  arma::vec theta(PP * JJ, arma::fill::zeros);
  
  arma::vec theta2(QQ * JJ, arma::fill::zeros);
  
  theta = theta_init;
  theta2 = theta_init2;
  
  
  //initialization for alpha
  Rcpp::NumericVector alpha0(JJ);
  arma::vec alpha(alpha0.begin(), JJ, false);
  alpha.fill(1.0);
  
  // variable inclusion indicator X & Z
  Rcpp::NumericVector inclusion_indicator0(PP*JJ);
  arma::vec inclusion_indicator(inclusion_indicator0.begin(), PP*JJ, false);
  for(jj = 0; jj < JJ; jj++){
    for(pp = 0; pp < PP; pp++){
      hh = jj + pp * JJ;
      if(theta(hh) != 0.0){
        inclusion_indicator(hh) = 1;
      } else {
        inclusion_indicator(hh) = 0;
      }
    }
  }
  
  Rcpp::NumericVector inclusion_indicator20(QQ*JJ);
  arma::vec inclusion_indicator2(inclusion_indicator20.begin(), QQ*JJ, false);
  for(jj = 0; jj < JJ; jj++){
    for(qq = 0; qq < QQ; qq++){
      hh = jj + qq * JJ;
      if(theta2(hh) != 0.0){
        inclusion_indicator2(hh) = 1;
      } else {
        inclusion_indicator2(hh) = 0;
      }
    }
  }
  
  //peNMIG
  arma::cube alpha_lin_xx(PP, PP, JJ, arma::fill::zeros);
  for(int jj = 0; jj < JJ; jj++){
    for(int pp = 0; pp < PP; pp++){
      for(int qq = pp + 1; qq < PP; qq++){
        alpha_lin_xx(pp, qq, jj) = 0.01;
      }
    }
  }
  
  arma::cube alpha_nl_xx(PP, PP, JJ, arma::fill::zeros);
  for(int jj = 0; jj < JJ; jj++){
    for(int pp = 0; pp < PP; pp++){
      for(int qq = pp + 1; qq < PP; qq++){
        alpha_nl_xx(pp, qq, jj) = 0.01;
      }
    }
  }
  
  arma::cube tausq_lin(PP, PP, JJ, arma::fill::zeros);
  for(int jj = 0; jj < JJ; jj++){
    for(int pp = 0; pp < PP; pp++){
      for(int qq = pp + 1; qq < PP; qq++){
        tausq_lin(pp, qq, jj) = 1.0;
      }
    }
  }
  
  arma::cube tausq_nl(PP, PP, JJ, arma::fill::zeros);
  for(int jj = 0; jj < JJ; jj++){
    for(int pp = 0; pp < PP; pp++){
      for(int qq = pp + 1; qq < PP; qq++){
        tausq_nl(pp, qq, jj) = 1.0;
      }
    }
  }
  
  arma::cube gamma_lin(PP, PP, JJ, arma::fill::zeros);
  for(int jj = 0; jj < JJ; jj++){
    for(int pp = 0; pp < PP; pp++){
      for(int qq = pp + 1; qq < PP; qq++){
        gamma_lin(pp, qq, jj) = 1.0;
      }
    }
  }
  
  arma::cube gamma_nl(PP, PP, JJ, arma::fill::zeros);
  for(int jj = 0; jj < JJ; jj++){
    for(int pp = 0; pp < PP; pp++){
      for(int qq = pp + 1; qq < PP; qq++){
        gamma_nl(pp, qq, jj) = 1.0;
      }
    }
  }
  
  arma::cube ksi_lin(PP, PP, JJ, arma::fill::zeros);
  for(int jj = 0; jj < JJ; jj++){
    for(int pp = 0; pp < PP; pp++){
      for(int qq = pp + 1; qq < PP; qq++){
        ksi_lin(pp, qq, jj) = 1.0;
      }
    }
  }
  
  arma::cube ksi_nl(PP, KK, JJ, arma::fill::zeros);
  arma::vec cdj(PP-1);
  cdj = cumsum(dj);
  for(int jj = 0; jj < JJ; jj++){
    for(int pp = 0; pp < PP; pp++){
      for(int qq = pp + 1; qq < PP; qq++){
        for(int ii = cdj(qq-1); ii < cdj(qq); ii++){
          ksi_nl(pp, ii, jj) = 1.0;
        }
      }
    }
  }
  
  arma::cube m_lin(PP, PP, JJ, arma::fill::zeros);
  for(int jj = 0; jj < JJ; jj++){
    for(int pp = 0; pp < PP; pp++){
      for(int qq = pp + 1; qq < PP; qq++){
        m_lin(pp, qq, jj) = 1.0;
      }
    }
  }
  
  arma::cube m_nl(PP, KK, JJ);
  for(int jj = 0; jj < JJ; jj++){
    for(int pp = 0; pp < PP; pp++){
      for(int qq = pp + 1; qq < PP; qq++){
        for(int ii = cdj(qq-1); ii < cdj(qq); ii++){
          m_nl(pp, qq, jj) = 1.0;
        }
      }
    }
  }
  
  arma::mat omega_lin(PP, JJ, arma::fill::zeros);
  for(int jj = 0; jj < JJ; jj++){
    for(int pp = 0; pp < PP; pp++){
      omega_lin(pp, jj) = 0.5;
    }
  }
  
  arma::mat omega_nl(PP, JJ, arma::fill::zeros);
  for(int jj = 0; jj < JJ; jj++){
    for(int pp = 0; pp < PP; pp++){
      omega_nl(pp, jj) = 0.5;
    }
  }
  
  Rcpp::List rescale_out;
  
  arma::cube bi_lin_xx(PP, PP, JJ);
  arma::cube bi_nl_xx(PP, KK, JJ);
  for(int jj = 0; jj < JJ; jj++){
    bi_lin_xx.slice(jj) = alpha_ksi_prod(alpha_lin_xx.slice(jj), ksi_lin.slice(jj), dj, true);
    bi_nl_xx.slice(jj) = alpha_ksi_prod(alpha_nl_xx.slice(jj), ksi_nl.slice(jj), dj, false);
  }
  
  double atau_lin = penmig_lin(0);
  double btau_lin = penmig_lin(1);
  
  double v0_lin = penmig_lin(2);
  
  double aomega = penmig_lin(3);
  double bomega = penmig_lin(4);
  
  double atau_nl = penmig_nl(0);
  double btau_nl = penmig_nl(1);
  
  double v0_nl = penmig_nl(2);
  
  double aomega_nl = penmig_nl(3);
  double bomega_nl = penmig_nl(4);
  
  //Y row sums
  arma::vec Ypiu(n);
  Ypiu.fill(0);
  
  for(int ii = 0; ii < n; ii ++){
    for(jj = 0; jj < JJ; jj ++){
      Ypiu(ii)=Ypiu(ii) + YY(ii, jj);
    }
  }
  
  //theta temp updates coeff matrix Ycategory-wise
  Rcpp::NumericVector theta_temp0(PP*JJ);
  arma::vec theta_temp(theta_temp0.begin(), PP*JJ, false);
  theta_temp.fill(0.0);
  
  Rcpp::NumericVector theta_temp20(QQ*JJ);
  arma::vec theta_temp2(theta_temp20.begin(), QQ*JJ, false);
  theta_temp2.fill(0.0);
  
  //matrici delle interazioni
  Rcpp::NumericMatrix bi0(PP*QQ*JJ);
  arma::mat bi_xz(bi0.begin(), JJ, PP*QQ, false);
  bi_xz.fill(0.0);
  arma::mat bXZ(PP, QQ);
  
  for(jj = 0 ; jj < JJ; jj ++){
    for(int qq2 = 0; qq2 < QQ; qq2 ++){
      for(int pp2 = 0; pp2 < PP; pp2 ++){
        int hh2 = qq2+((pp2)*QQ);
        bXZ(pp2, qq2)=bi_xz(jj, hh2);
      }
    }
  }
  
  //matrice delle interazioni temporanea
  Rcpp::NumericMatrix bi_p_mat0(PP, QQ);
  arma::mat bi_p_mat(bi_p_mat0.begin(), PP, QQ, false);
  bi_p_mat.fill(0.0);
  
  Rcpp::NumericMatrix bi20(JJ, QQ * QQ);
  arma::mat bi_zz(bi20.begin(), JJ, QQ * QQ, false);
  bi_zz.fill(0.0);
  arma::mat bZZ(QQ, QQ);
  
  for(jj = 0 ; jj < JJ; jj++){
    for(int qq2 = 0; qq2 < QQ; qq2++){
      for(int pp2 = 0; pp2 < QQ; pp2++){
        int hh2 = qq2+((pp2)*QQ);
        bZZ(pp2, qq2)=bi_zz(jj, hh2);
      }
    }
  }
  
  //matrice delle interazioni temporanea
  Rcpp::NumericMatrix bi_p_mat20(QQ, QQ);
  arma::mat bi_p_mat2(bi_p_mat20.begin(), QQ, QQ, false);
  bi_p_mat2.fill(0.0);
  
  //matrici per HS
  arma::mat lambda_hs(JJ, (PP*QQ));
  lambda_hs.fill(1.0);
  
  arma::vec tau_hs((QQ*JJ));
  tau_hs.fill(.1);
  
  arma::mat lambda_hs2(JJ, (QQ * QQ));
  lambda_hs2.fill(1.0);
  
  arma::vec tau_hs2((QQ*JJ));
  tau_hs2.fill(.1);
  
  double mu_hs = 0.0;
  
  //threshold parameters
  double soglia = 0.01;
  
  double lb = tx(0);
  double ub = tx(1);
  
  double soglia2 = 0.01;
  
  double lb2 = tz(0);
  double ub2 = tz(1);
  
  //initialization for random intercept
  Rcpp::NumericVector ran_int_sing0(ngroup);
  arma::vec ran_int_sing(ran_int_sing0.begin(), ngroup, false);
  ran_int_sing.fill(1.0);
  
  Rcpp::NumericVector ran_int_exp0(n);
  arma::vec ran_int_exp(ran_int_exp0.begin(), n, false);
  ran_int_exp.fill(1.0);
  
  arma::vec sig_int(ngroup);
  sig_int.fill(1.0);
  double asi = 1.0;
  double bsi = 1.0;
  
  //linear predictor  LOG SCALE
  Rcpp::NumericMatrix loggamma0(n, JJ);
  arma::mat loggamma(loggamma0.begin(), n, JJ, false);
  
  
  for(jj = 0 ; jj < JJ; jj ++){
    for(int qq2 = 0; qq2 < QQ; qq2 ++){
      for(int pp2 = 0; pp2 < PP; pp2 ++){
        int hh2 = qq2+((pp2)*QQ);
        bXZ(pp2, qq2)=bi_xz(jj, hh2);
      }
      
      for(int pp2 = 0; pp2 < QQ; pp2 ++){
        int hh2 = qq2+((pp2)*QQ);
        bZZ(pp2, qq2)=bi_zz(jj, hh2);
      }
    }
    loggamma.col(jj) = calculate_gamma_ri(XX, ZZ, XXl, XXT, alpha, ran_int_exp, theta, theta2, bXZ,
                 bi_lin_xx.slice(jj), bi_nl_xx.slice(jj), soglia, bZZ, soglia2, jj);
  }
  
  //latent variables
  arma::mat SS(n, JJ);
  SS.fill(0.0);
  arma::vec TT(n);
  TT.fill(0.0);
  for(ii = 0; ii < n; ii ++){
    for(jj = 0; jj < JJ; jj ++){
      SS(ii, jj) = (YY(ii, jj)/1.0);
      if(SS(ii, jj) < pow(10.0, -100.0)){
        SS(ii, jj) = pow(10.0, -100.0);
      }
      TT(ii)=TT(ii)+SS(ii, jj);
    }
  }
  
  arma::vec uu(n);
  for(ii = 0 ; ii < n ; ii ++){
    uu(ii) = ::Rf_rgamma(Ypiu(ii), 1.0/TT(ii));
  }
  
  //beta proposal variance
  Rcpp::NumericVector prop_per_theta0(PP*JJ);
  arma::vec prop_per_theta(prop_per_theta0.begin(), PP*JJ, false);
  prop_per_theta.fill(0.5);
  
  Rcpp::NumericVector prop_per_theta20(QQ*JJ);
  arma::vec prop_per_theta2(prop_per_theta20.begin(), QQ*JJ, false);
  prop_per_theta2.fill(0.5);
  
  //alpha proposal variance
  Rcpp::NumericVector prop_per_alpha0(JJ);
  arma::vec prop_per_alpha(prop_per_alpha0.begin(), JJ, false);
  prop_per_alpha.fill(0.5);
  
  //bi proposal variance
  Rcpp::NumericMatrix prop_per_bi0(JJ, PP*QQ);
  arma::mat prop_per_bi(prop_per_bi0.begin(), JJ, PP*QQ, false);
  prop_per_bi.fill(0.5);
  
  Rcpp::NumericMatrix prop_per_bi20(JJ, QQ * QQ);
  arma::mat prop_per_bi2(prop_per_bi20.begin(), JJ, QQ * QQ, false);
  prop_per_bi2.fill(0.5);
  
  //mean_slab
  arma::mat mu_theta(PP, JJ);
  mu_theta.fill(0.0);
  arma::mat mu_theta2(QQ, JJ);
  mu_theta2.fill(0.0);
  
  double tau_slab = hyper_theta_x(1);
  double tau_slab2 = hyper_theta_z(1);
  //slab_variance
  arma::mat sigma_theta(PP, JJ);
  sigma_theta.fill(tau_slab);
  arma::mat sigma_theta2(QQ, JJ);
  sigma_theta2.fill(tau_slab2);
  
  //Inverse Gamma Updating 1
  arma::mat sigma_theta_temp(PP, JJ);
  sigma_theta_temp = sigma_theta;
  
  double aig_slab = 1.0;
  double big_slab = 1.0;
  
  arma::vec idx_tmp(PP);
  arma::vec thetavec_tmp(PP);
  double thetavecTthetavec;
  int sum_idx_tmp;
  
  //Inverse Gamma Updating 2
  arma::mat sigma_theta_temp2(QQ, JJ);
  sigma_theta_temp2 = sigma_theta2;
  
  double aig_slab2 = 3.0;
  double big_slab2 = .5;
  
  arma::vec idx_tmp2(QQ);
  arma::vec thetavec_tmp2(QQ);
  double thetavecTthetavec2;
  int sum_idx_tmp2;
  
  // mean and standard deviation of the independent normal priors on alpha and beta
  arma:: vec mu_al(JJ);
  mu_al.fill(prior_int(0));
  arma::vec sig_al(JJ);
  sig_al.fill(prior_int(1));
  
  double aa_hp = hyper_theta_x(0);
  double bb_hp = 2 - aa_hp;
  
  double aa_hp2 = hyper_theta_z(0);
  double bb_hp2 = 2 - aa_hp2;
  
  arma::cube THETApost(PP, JJ, Eff);
  arma::cube THETApost2(QQ, JJ, Eff);
  arma::cube BIpost(JJ, (PP*QQ), Eff);
  arma::cube hBIpost(n, PP, JJ, arma::fill::zeros);
  arma::cube hBIpost_ppi(n, PP, JJ, arma::fill::zeros);
  arma::cube BIpost2(JJ, (QQ * QQ), Eff);
  arma::cube hBIpost2(n, QQ, JJ, arma::fill::zeros);
  arma::cube LOGLINPRED(n, JJ, Eff);
  arma::vec llvec(Eff);
  arma::vec SOGLIApost(Eff);
  arma::vec SOGLIApost2(Eff);
  arma::mat ALPHApost(Eff, JJ);
  
  arma::field<arma::cube> LINCOEF(Eff);
  LINCOEF.fill(arma::cube(PP, PP, JJ, arma::fill::zeros));
  arma::field<arma::cube> ALPHA_LIN(Eff);
  ALPHA_LIN.fill(arma::cube(PP, PP, JJ, arma::fill::zeros));
  arma::field<arma::cube> M_LIN(Eff);
  M_LIN.fill(arma::cube(PP, PP, JJ, arma::fill::zeros));
  arma::field<arma::cube> KSI_LIN(Eff);
  KSI_LIN.fill(arma::cube(PP, PP, JJ, arma::fill::zeros));
  arma::field<arma::cube> TAUSQ_LIN(Eff);
  TAUSQ_LIN.fill(arma::cube(PP, PP, JJ, arma::fill::zeros));
  arma::field<arma::cube> GAMMA_LIN(Eff);
  GAMMA_LIN.fill(arma::cube(PP, PP, JJ, arma::fill::zeros));
  arma::cube OMEGA_LIN(PP, JJ, Eff, arma::fill::zeros);
  
  arma::field<arma::cube> NONLINCOEF(Eff);
  NONLINCOEF.fill(arma::cube(PP, KK, JJ, arma::fill::zeros));
  arma::field<arma::cube> ALPHA_NL(Eff);
  ALPHA_NL.fill(arma::cube(PP, PP, JJ, arma::fill::zeros));
  arma::field<arma::cube> M_NL(Eff);
  M_NL.fill(arma::cube(PP, KK, JJ, arma::fill::zeros));
  arma::field<arma::cube> KSI_NL(Eff);
  KSI_NL.fill(arma::cube(PP, KK, JJ, arma::fill::zeros));
  arma::field<arma::cube> TAUSQ_NL(Eff);
  TAUSQ_NL.fill(arma::cube(PP, PP, JJ, arma::fill::zeros));
  arma::field<arma::cube> GAMMA_NL(Eff);
  GAMMA_NL.fill(arma::cube(PP, PP, JJ, arma::fill::zeros));
  arma::cube OMEGA_NL(PP, JJ, Eff, arma::fill::zeros);
  arma::mat RANDINTpost(Eff, n);
  arma::mat VARRANDINTpost(Eff, ngroup);
  
  Rcpp::List listainter;
  int curr_iter = 0;
  int zz = 0;
  
  for( s = 0; s < Niter+1; s++){
    
    ////////////////////////
    // AGGIORNO I THETA X //
    ////////////////////////
    
    //scelgo un taxa a caso per fare lo swap
    double u;
    u = arma::randu();
    u = u*(JJ/1.0);
    jj = Rcpp::as<int>(Rcpp::wrap(std::floor(u)));
    
    for(pp = 0; pp < PP; pp ++){
      int ci = jj + pp * JJ;
      theta_temp(ci) = theta(ci);
    }
    
    swap(XX, ZZ, SS, loggamma, alpha, theta_temp, theta_temp2, theta, theta2,
         bi_xz, soglia, bi_zz, soglia2, acc_theta_flag, acc_theta_flag2, inclusion_indicator,
         inclusion_indicator2, mu_theta, mu_theta2, sigma_theta, sigma_theta2,
         aa_hp, bb_hp, aa_hp2, bb_hp2, jj, true);
    
    // update theta with theta_temp
    for(pp = 0 ; pp < PP; pp ++){
      int ci = jj + pp * JJ;
      theta(ci) = theta_temp(ci);
    }
    //}
    
    //within model
    for(jj = 0; jj < JJ; jj ++){
      for(pp = 0; pp < PP; pp ++){
        int ci = jj + pp * JJ;
        theta_temp(ci) = theta(ci);
        
        if(s > (std::floor(burn/double(4)))){
          last_mean = curr_mean(ci);
          last_var = curr_var(ci);
          curr_mean(ci) = online_mean(s, last_mean, theta(ci));
          curr_var(ci) = online_var(s, last_mean, last_var, curr_mean(ci), theta(ci));
          
          //update proposal variance
          prop_per_theta(ci)=curr_var(ci);
        }
      }//chiude for pp
      
      //aggiorno i thetaX
      update_theta(XX, ZZ, SS, loggamma, alpha, theta_temp, theta_temp2, theta,
                   theta2, bi_xz, soglia, bi_zz, soglia2, inclusion_indicator,
                   inclusion_indicator2, prop_per_theta, prop_per_theta2,
                   mu_theta, mu_theta2, sigma_theta, sigma_theta2, aa_hp,
                   bb_hp, aa_hp2, bb_hp2, jj, true);
      
      for(pp = 0; pp < PP; pp ++){
        int ci = jj + pp * JJ;
        theta(ci) = theta_temp(ci);
        if((s > burn) && (s % thin == 0)){
          acc_theta(ci) = acc_theta(ci) + acc_theta_flag(ci);
          acc_theta_flag(ci) = 0;
        }
      }
    }//chiude for sui jj
    
    //qui salvo output
    if((s > burn) && (s % thin == 0)){
      zz = 0;
      for(pp = 0; pp < PP; pp++){
        for(jj = 0; jj < JJ; jj++){
          THETApost(pp, jj, curr_iter) = theta(zz);
          zz++;
        }
      }
    }
    
    if(upsv == true){
      //aggiorno inverse gamma
      for(jj = 0; jj < JJ; jj++){
        for(pp = 0; pp < PP; pp++){
          
          hh = jj + pp * JJ;
          idx_tmp(pp) = inclusion_indicator(hh);
          thetavec_tmp(pp) = theta(hh);
        }
        
        thetavecTthetavec = inner(thetavec_tmp, thetavec_tmp);
        sum_idx_tmp = arma::sum(idx_tmp);
        double rtauj;
        if(sum_idx_tmp == 0.0){
          rtauj = 1.0/::Rf_rgamma(aig_slab, big_slab);
        } else {
          rtauj = 1.0/::Rf_rgamma(0.5 * aig_slab + 0.5 * sum_idx_tmp, 0.5 * pow(big_slab, 2.0) + 0.5 * thetavecTthetavec);
        }
        sigma_theta_temp.col(jj).fill(rtauj);
      }
      sigma_theta = sigma_theta_temp;
    }
    
    ////////////////////////
    // AGGIORNO I THETA Z //
    ////////////////////////
    
    //scelgo un taxa a caso per fare lo swap
    u = arma::randu();
    u = u*(JJ/1.0);
    jj = Rcpp::as<int>(Rcpp::wrap(std::floor(u)));
    for(qq = 0; qq < QQ; qq ++){
      int ci = jj + qq * JJ;
      theta_temp2(ci)=theta2(ci);
    }
    
    swap(XX, ZZ, SS, loggamma, alpha, theta_temp, theta_temp2, theta, theta2,
         bi_xz, soglia, bi_zz, soglia2, acc_theta_flag, acc_theta_flag2, inclusion_indicator,
         inclusion_indicator2, mu_theta, mu_theta2, sigma_theta, sigma_theta2,
         aa_hp, bb_hp, aa_hp2, bb_hp2, jj, false);
    
    // update theta with theta_temp
    for(qq = 0 ; qq < QQ ; qq++){
      int ci = jj + qq * JJ;
      theta2(ci) = theta_temp2(ci);
    }
    
    //within model
    for(jj = 0; jj < JJ; jj ++){
      for(qq = 0; qq < QQ; qq ++){
        int ci = jj + qq * JJ;
        theta_temp2(ci)= theta2(ci);
        //aggiorno media e varianza adaptive
        if(s > (std::floor(burn/double(4)))){
          last_mean2 = curr_mean2(ci);
          last_var2 = curr_var2(ci);
          curr_mean2(ci) = online_mean(s, last_mean2, theta2(ci));
          curr_var2(ci) = online_var(s, last_mean2, last_var2, curr_mean2(ci),
                    theta2(ci));
          
          //update proposal variance
          prop_per_theta2(ci)=curr_var2(ci);
        }
      }//chiude for pp
      
      //aggiorno i thetaZ
      update_theta(XX, ZZ, SS, loggamma, alpha, theta_temp, theta_temp2, theta,
                   theta2, bi_xz, soglia, bi_zz, soglia2, inclusion_indicator, inclusion_indicator2,
                   prop_per_theta, prop_per_theta2, mu_theta, mu_theta2,
                   sigma_theta, sigma_theta2, aa_hp, bb_hp, aa_hp2, bb_hp2, jj,
                   false);
      
      for(qq = 0; qq < QQ; qq++){
        int ci = jj + qq * JJ;
        theta2(ci) = theta_temp2(ci);
        if((s>burn) && (s%thin==0)){
          acc_theta2(ci) = acc_theta2(ci) + acc_theta_flag2(ci);
          acc_theta_flag2(ci) = 0;
        }
      }
    }//chiude for sui jj
    
    //qui salvo output
    if((s > burn) && (s % thin == 0)){
      zz = 0;
      for(qq = 0; qq < QQ; qq++){
        for(jj = 0; jj < JJ; jj++){
          THETApost2(qq, jj, curr_iter) = theta2(zz);
          zz++;
        }
      }
    }
    
    if(upsv == true){
      //aggiorno inverse gamma
      for(jj = 0; jj < JJ; jj++){
        for(qq = 0; qq < QQ; qq++){
          
          hh = jj + qq * JJ;
          idx_tmp2(qq) = inclusion_indicator2(hh);
          thetavec_tmp2(qq) = theta2(hh);
        }
        
        thetavecTthetavec2 = inner(thetavec_tmp2, thetavec_tmp2);
        sum_idx_tmp2 = arma::sum(idx_tmp2);
        
        double rtau2j;
        if(sum_idx_tmp2 == 0.0){
          rtau2j = 1.0/::Rf_rgamma(aig_slab2, big_slab2);
        } else {
          rtau2j = 1.0/::Rf_rgamma(0.5 * aig_slab2 + 0.5 * sum_idx_tmp2, 0.5 * pow(big_slab2, 2.0) + 0.5 * thetavecTthetavec2);
        }
        sigma_theta_temp2.col(jj).fill(rtau2j);
      }
      sigma_theta2 = sigma_theta_temp2;
    }
    
    ///////////////////////////
    // AGGIORNO INTERAZIONI1 //
    ///////////////////////////
    
    u = arma::randu();
    u = u*(JJ/1.0);
    jj = Rcpp::as<int>(Rcpp::wrap(std::floor(u)));
    
    for(int qq2 = 0; qq2 < QQ; qq2++){
      for(int pp2 = 0; pp2 < PP; pp2++){
        int hh2 = qq2 + ((pp2) * QQ);
        bi_p_mat(pp2, qq2)= bi_xz(jj, hh2);// + ::Rf_rnorm(0.0, .1);
      }
    }
    
    // ADAPTIVE PROPOSAL
    if(s > 2000){
      for(pp = 0; pp < PP; pp ++){
        for(qq = 0; qq < QQ; qq++){
          
          int hh2 = qq + pp * QQ;
          last_meani = curr_meani(jj, hh2);
          last_vari = curr_vari(jj, hh2);
          curr_meani(jj, hh2) = online_mean(s, last_meani, bi_xz(jj, hh2));
          curr_vari(jj, hh2) = online_var(s, last_meani, last_vari, curr_meani(jj, hh2), bi_xz(jj, hh2));
          //update proposal variance
          prop_per_bi(jj, hh2)=curr_vari(jj, hh2);
          
          
        }//chiude ciclo su qq
      }
    }//chiude if
    
    for(pp = 0; pp < PP; pp ++){
      for(qq = 0; qq < QQ; qq++){
        update_bi_ri(XX, ZZ, XXl, XXT, SS, loggamma, lambda_hs, tau_hs, alpha, ran_int_exp, theta, theta2,
                     bi_xz, bi_zz, bi_p_mat, bi_lin_xx.slice(jj), bi_nl_xx.slice(jj), soglia, soglia2,
                     prop_per_bi, acc_bi_flag, mu_hs, pp, qq, jj, hereditariety);//hereditariety
        
        hh = qq + pp * QQ;
        if((s>burn) && (s%thin==0)){
          acc_bi(jj, hh)= acc_bi(jj, hh)+acc_bi_flag(jj, hh);
          acc_bi_flag(jj, hh) = 0;
        }
      }//chiude il ciclo sulla qq
    }
    
    if((s > burn)&& (s % thin == 0)){
      for(jj = 0; jj < JJ; jj ++){
        for(hh = 0; hh < (PP * QQ); hh ++){
          BIpost(jj, hh, curr_iter) = bi_xz(jj, hh);
        }
      }
    }
    
    ///////////////////////////
    // AGGIORNO HYP HS PRIOR //
    ///////////////////////////
    
    arma::vec bi_work(QQ);
    
    u = arma::randu();
    u = u*(double)JJ;
    jj = Rcpp::as<int>(Rcpp::wrap(std::floor(u)));
    
    u = arma::randu();
    u = u*(double)PP;
    pp = Rcpp::as<int>(Rcpp::wrap(std::floor(u)));
    
    for(qq = 0; qq < QQ; qq++){
      hh = qq + pp * QQ;
      bi_work(qq) = bi_xz(jj, hh);
    }
    
    u = arma::randu();
    u = u*(double)QQ;
    qq = Rcpp::as<int>(Rcpp::wrap(std::floor(u)));
    ll = jj + qq * JJ;
    
    lambda_hs.row(jj) = uplam(lambda_hs.row(jj), bi_work, pp, tau_hs(ll));
    tau_hs(ll) = uptau(lambda_hs.row(jj).t(), bi_work, pp, QQ, tau_hs(ll), true);
    
    //peNMIG parte lineare
    u = arma::randu();
    u = u*(double)JJ;
    jj = Rcpp::as<int>(Rcpp::wrap(std::floor(u)));
    
    alpha_lin_xx.slice(jj) = update_alpha_xx_ri(XX, ZZ, XXl, XXT, SS, loggamma, alpha, ran_int_exp,
                       theta, theta2, bi_xz, soglia, bi_zz, soglia2, dj, alpha_lin_xx.slice(jj), alpha_nl_xx.slice(jj),
                       tausq_lin.slice(jj), tausq_nl.slice(jj), gamma_lin.slice(jj),
                       gamma_nl.slice(jj), ksi_lin.slice(jj), ksi_nl.slice(jj), true, jj, hereditariety, vincolo);
    
    m_lin.slice(jj) = update_m(m_lin.slice(jj), ksi_lin.slice(jj));
    
    ksi_lin.slice(jj) = update_ksi_ri(XX, ZZ, XXl, XXT, SS, loggamma, alpha, ran_int_exp, theta, theta2,
                  bi_xz, soglia, bi_zz, soglia2, dj, alpha_lin_xx.slice(jj), alpha_nl_xx.slice(jj),
                  tausq_lin.slice(jj), tausq_nl.slice(jj), gamma_lin.slice(jj),
                  gamma_nl.slice(jj), ksi_lin.slice(jj), ksi_nl.slice(jj), m_lin.slice(jj),
                  m_nl.slice(jj), true, jj);
    
    rescale_out = rescale(alpha_lin_xx.slice(jj), ksi_lin.slice(jj), dj, true);
    alpha_lin_xx.slice(jj) = Rcpp::as<arma::mat>(rescale_out[0]);
    ksi_lin.slice(jj) = Rcpp::as<arma::mat>(rescale_out[1]);
    
    bi_lin_xx.slice(jj) = alpha_ksi_prod(alpha_lin_xx.slice(jj),
                    ksi_lin.slice(jj), dj, true);
    
    tausq_lin.slice(jj) = update_tau(tausq_lin.slice(jj), alpha_lin_xx.slice(jj),
                    atau_lin, btau_lin, gamma_lin.slice(jj), jj);
    
    gamma_lin.slice(jj) = update_gamma(alpha_lin_xx.slice(jj),
                    tausq_lin.slice(jj), gamma_lin.slice(jj),
                    omega_lin.col(jj), v0_lin, jj, JJ, theta, hereditariety);
    
    omega_lin.col(jj) = update_omega(gamma_lin.slice(jj), omega_lin.col(jj),
                  v0_lin, aomega, bomega);
    
    if((s>burn) && (s%thin==0)){
      ALPHA_LIN[curr_iter] = alpha_lin_xx;
      M_LIN[curr_iter] = m_lin;
      KSI_LIN[curr_iter] = ksi_lin;
      LINCOEF[curr_iter] = bi_lin_xx;
      TAUSQ_LIN[curr_iter] = tausq_lin;
      GAMMA_LIN[curr_iter] = gamma_lin;
      OMEGA_LIN.slice(curr_iter) = omega_lin;
    }
    
    //peNMIG non linear
    u = arma::randu();
    u = u*(double)JJ;
    jj = Rcpp::as<int>(Rcpp::wrap(std::floor(u)));
    alpha_nl_xx.slice(jj) = update_alpha_xx_ri(XX, ZZ, XXl, XXT, SS, loggamma, alpha, ran_int_exp,
                      theta, theta2, bi_xz, soglia, bi_zz, soglia2, dj, alpha_lin_xx.slice(jj), alpha_nl_xx.slice(jj),
                      tausq_lin.slice(jj), tausq_nl.slice(jj), gamma_lin.slice(jj),
                      gamma_nl.slice(jj), ksi_lin.slice(jj), ksi_nl.slice(jj), false, jj, hereditariety, vincolo);
    
    m_nl.slice(jj) = update_m(m_nl.slice(jj), ksi_nl.slice(jj));
    
    ksi_nl.slice(jj) = update_ksi_ri(XX, ZZ, XXl, XXT, SS, loggamma, alpha,ran_int_exp, theta, theta2,
                 bi_xz, soglia, bi_zz, soglia2, dj, alpha_lin_xx.slice(jj), alpha_nl_xx.slice(jj),
                 tausq_lin.slice(jj), tausq_nl.slice(jj), gamma_lin.slice(jj),
                 gamma_nl.slice(jj), ksi_lin.slice(jj), ksi_nl.slice(jj), m_lin.slice(jj),
                 m_nl.slice(jj), false, jj);
    
    rescale_out = rescale(alpha_nl_xx.slice(jj), ksi_nl.slice(jj), dj, false);
    alpha_nl_xx.slice(jj) = Rcpp::as<arma::mat>(rescale_out[0]);
    
    ksi_nl.slice(jj) = Rcpp::as<arma::mat>(rescale_out[1]);
    bi_nl_xx.slice(jj) = alpha_ksi_prod(alpha_nl_xx.slice(jj),
                   ksi_nl.slice(jj), dj, false);
    
    tausq_nl.slice(jj) = update_tau(tausq_nl.slice(jj), alpha_nl_xx.slice(jj),
                   atau_nl, btau_nl, gamma_nl.slice(jj), jj);
    
    gamma_nl.slice(jj) = update_gamma(alpha_nl_xx.slice(jj),
                   tausq_nl.slice(jj), gamma_nl.slice(jj), omega_nl.col(jj),
                   v0_nl, jj, JJ, theta, hereditariety);
    
    omega_nl.col(jj) = update_omega(gamma_nl.slice(jj), omega_nl.col(jj),
                 v0_nl, aomega_nl, bomega_nl);
    
    if((s>burn) && (s%thin==0)){
      ALPHA_NL[curr_iter] = alpha_nl_xx;
      M_NL[curr_iter] = m_nl;
      KSI_NL[curr_iter] = ksi_nl;
      NONLINCOEF[curr_iter] = bi_nl_xx;
      TAUSQ_NL[curr_iter] = tausq_nl;
      GAMMA_NL[curr_iter] = gamma_nl;
      OMEGA_NL.slice(curr_iter) = omega_nl;
    }
    
    ///////////////////////////
    //    AGGIORNO SOGLIA1   //
    ///////////////////////////
    
    soglia = update_soglia_ri(XX, ZZ, XXl, XXT, SS, loggamma, alpha, ran_int_exp, theta, theta2, bi_xz,
                              bi_lin_xx, bi_nl_xx, soglia, bi_zz, soglia2, acc_soglia_flag, lb, ub);
    
    if((s>burn) && (s%thin==0)){
      acc_soglia(0)=acc_soglia(0)+acc_soglia_flag(0);
      acc_soglia_flag(0)=0;
      SOGLIApost(curr_iter) = soglia;
    }
    
    ///////////////////////////
    //    Eff. Sub-Spec. 1   //
    ///////////////////////////
    
    if((s > burn)&& (s % thin == 0)){
      arma::mat bXZ(PP, QQ);
      arma::mat work(n, PP);
      
      for(jj = 0; jj < JJ; jj++){
        for(int qq2 = 0; qq2 < QQ; qq2++){
          for(int pp2 = 0; pp2 < PP; pp2++){
            int hh2 = qq2 + ((pp2) * QQ);
            bXZ(pp2, qq2) = bi_xz(jj, hh2);
          }
        }
        work = (bXZ * ZZ.t()).t() + (bi_lin_xx.slice(jj) * XXl.t()).t()
          + (bi_nl_xx.slice(jj) * XXT.t()).t();
        hBIpost.slice(jj) = hBIpost.slice(jj) + threshold_mat(work, soglia);
      }//nPJ
    }
    
    ///////////////////////////
    // AGGIORNO INTERAZIONI2 //
    ///////////////////////////
    
    u = arma::randu();
    u = u*(JJ/1.0);
    jj = Rcpp::as<int>(Rcpp::wrap(std::floor(u)));
    
    for(int qq2 = 0; qq2 < QQ; qq2++){
      for(int pp2 = 0; pp2 < QQ; pp2++){
        int hh2 = qq2 + ((pp2) * QQ);
        bi_p_mat2(pp2, qq2)= bi_zz(jj, hh2);
      }
    }
    
    // ADAPTIVE PROPOSAL
    if(s > 2000){
      for(pp = 0; pp < QQ; pp ++){
        for(qq = (pp + 1); qq < QQ; qq++){
          
          int hh2 = qq + pp * QQ;
          last_meani2 = curr_meani2(jj, hh2);
          last_vari2 = curr_vari2(jj, hh2);
          curr_meani2(jj, hh2) = online_mean(s, last_meani2, bi_zz(jj, hh2));
          curr_vari2(jj, hh2) = online_var(s, last_meani2, last_vari2, curr_meani2(jj, hh2), bi_zz(jj, hh2));
          //update proposal variance
          prop_per_bi2(jj, hh2)=curr_vari2(jj, hh2);
          
          
        }//chiude ciclo su qq
      }
    }//chiude if
    
    for(pp = 0; pp < QQ; pp ++){
      for(qq = (pp + 1); qq < QQ; qq++){
        update_bi2_ri(XX, ZZ, SS, loggamma, lambda_hs2, tau_hs2, alpha, theta, theta2,
                      bi_xz, bi_zz, bi_p_mat2, soglia, soglia2, prop_per_bi2,
                      acc_bi_flag2, mu_hs, pp, qq, jj, hereditariety);//hereditariety
      }
    }
    
    for(pp = 0; pp < QQ; pp ++){
      for(qq = 0; qq < QQ; qq++){
        hh = qq + pp * QQ;
        if((s>burn) && (s%thin==0)){
          acc_bi2(jj, hh)= acc_bi2(jj, hh)+acc_bi_flag2(jj, hh);
          acc_bi_flag2(jj, hh) = 0;
        }
      }//chiude il ciclo sulla qq
    }
    
    if((s > burn)&& (s % thin == 0)){
      for(jj = 0; jj < JJ; jj ++){
        for(hh = 0; hh < (QQ * QQ); hh ++){
          BIpost2(jj, hh, curr_iter) = bi_zz(jj, hh);
        }
      }
    }
    
    ///////////////////////////
    // AGGIORNO HYP HS PRIOR //
    ///////////////////////////
    
    arma::vec bi_work2(QQ);
    
    u = arma::randu();
    u = u*(double)JJ;
    jj = Rcpp::as<int>(Rcpp::wrap(std::floor(u)));
    
    u = arma::randu();
    u = u*(double)QQ;
    pp = Rcpp::as<int>(Rcpp::wrap(std::floor(u)));
    
    for(qq = 0; qq < QQ; qq++){
      hh = qq + pp * QQ;
      bi_work2(qq) = bi_zz(jj, hh);
    }
    
    u = arma::randu();
    u = u*(double)QQ;
    qq = Rcpp::as<int>(Rcpp::wrap(std::floor(u)));
    ll = jj + qq * JJ;
    
    lambda_hs2.row(jj) = uplam(lambda_hs2.row(jj), bi_work2, pp, tau_hs2(ll));
    tau_hs2(ll) = uptau(lambda_hs2.row(jj).t(), bi_work2, qq, QQ, tau_hs2(ll), true);
    
    ///////////////////////////
    //    AGGIORNO SOGLIA2   //
    ///////////////////////////
    
    soglia2 = update_soglia_ri(XX, ZZ, XXl, XXT, SS, loggamma, alpha, ran_int_exp, theta, theta2,
                               bi_xz, bi_lin_xx, bi_nl_xx, soglia, bi_zz, soglia2, acc_soglia_flag, lb2, ub2, 1);
    
    if((s>burn) && (s%thin==0)){
      acc_soglia(1)=acc_soglia(1)+acc_soglia_flag(1);
      acc_soglia_flag(1)=0;
      SOGLIApost2(curr_iter) = soglia2;
    }
    
    ///////////////////////////
    //    Eff. Sub-Spec. 2   //
    ///////////////////////////
    
    if((s > burn)&& (s % thin == 0)){
      arma::mat bZZ(QQ, QQ);
      arma::mat work(n, QQ);
      
      for(jj = 0; jj < JJ; jj++){
        for(int qq2 = 0; qq2 < QQ; qq2++){
          for(int pp2 = 0; pp2 < QQ; pp2++){
            int hh2 = qq2 + ((pp2) * QQ);
            bZZ(pp2, qq2) = bi_zz(jj, hh2);
          }
        }
        work = (bZZ * ZZ.t()).t();
        hBIpost2.slice(jj) = hBIpost2.slice(jj) + threshold_mat(work, soglia2);
      }
    }
    
    //////////////////////////
    // AGGIORNO INTERCETTA ///
    //////////////////////////
    
    for(jj = 0 ; jj < JJ ; jj++){
      update_alpha(XX, ZZ, SS, loggamma, alpha, theta, theta2, bi_xz, bi_zz, soglia, soglia2, prop_per_alpha, acc_alpha_flag,
                   mu_al, sig_al, jj);
      
      if((s>burn) && (s%thin==0)){
        acc_alpha(jj)=acc_alpha(jj)+acc_alpha_flag(jj);
        acc_alpha_flag(jj)=0;
      }
    }
    
    //qui salvo output
    if((s > burn) && (s % thin == 0)){
      ALPHApost.row(curr_iter) = alpha.t();
    }
    
    /////////////////////////////////
    // AGGIORNO RANDOM INTERCEPTS ///
    /////////////////////////////////
    
    update_rand_int(XX, ZZ, SS, loggamma, ran_int_sing, ran_int_exp, grplabel, ngroup, sig_int);
    
    for(int gg = 0; gg < ngroup; gg++){
      sig_int(gg) = 1.0/::Rf_rgamma(asi + 0.5, bsi + (pow(ran_int_sing(gg), 2.0)/2));
    }
    if((s > burn) && (s % thin == 0)){
      RANDINTpost.row(curr_iter) = ran_int_exp.t();
    }
    
    if((s > burn) && (s % thin == 0)){
      VARRANDINTpost.row(curr_iter) = sig_int.t();
    }
    
    llvec(curr_iter) = logliksimpleC(YY, exp(loggamma));
    
    //
    if((s > burn) && (s % thin == 0)){
      LOGLINPRED.slice(curr_iter) = loggamma;
      curr_iter++;
    }
    
    
    // update JJ and consequently TT
    for(ii = 0 ; ii < n ; ii++){
      TT(ii) = 0.0;
      for(jj = 0 ; jj < JJ ; jj++){
        SS(ii,jj) = ::Rf_rgamma(YY(ii, jj)+ exp(loggamma(ii,jj)), 1.0/(uu(ii) + 1.0));
        if(SS(ii,jj) < pow(10.0, -100.0)){
          SS(ii,jj) = pow(10.0, -100.0);
        }
        TT(ii) = TT(ii) + SS(ii, jj);
      }
    }
    
    // update latent variables uu
    for(ii = 0 ; ii < n ; ii++){
      uu(ii) = ::Rf_rgamma(Ypiu(ii), 1.0/TT(ii));
    }
  }//chiude le iterazioni della catena
  
  arma::vec tasso_alpha(JJ);
  arma::mat tasso_thetax(PP, JJ);
  arma::mat tasso_bi(JJ, PP * QQ);
  arma::mat tasso_bi2(JJ, QQ * QQ);
  arma::mat tasso_prop_thetax(PP, JJ);
  arma::mat tasso_thetaz(QQ, JJ);
  arma::mat tasso_prop_thetaz(QQ, JJ);
  
  for(jj = 0; jj < JJ; jj++){
    tasso_alpha(jj) = (double)acc_alpha(jj)/Eff;
    for(pp=0;pp<PP;pp++){
      tasso_thetax(pp, jj) = (double)acc_theta(pp+ jj*PP)/Eff;
      tasso_prop_thetax(pp, jj) = prop_per_theta(pp+jj*PP)/Eff;
    }
    
    for(qq = 0; qq < QQ; qq ++){
      tasso_thetaz(qq, jj) = (double)acc_theta2(qq+ jj*QQ)/Eff;
      tasso_prop_thetaz(qq, jj) = prop_per_theta2(qq+jj*QQ)/Eff;
    }
    
    for(hh = 0; hh < (PP*QQ); hh ++){
      tasso_bi(jj, hh) = (double)acc_bi(jj, hh)/Eff;
    }
    for(hh = 0; hh < (QQ * QQ); hh ++){
      tasso_bi2(jj, hh) = (double)acc_bi2(jj, hh)/Eff;
    }
  }
  
  Rcpp::List lista_acc_rates =
    Rcpp::List::create(Rcpp::Named("tassoalpha") = tasso_alpha,
                       Rcpp::Named("tassothetax") = tasso_thetax,
                       Rcpp::Named("tassothetaz") = tasso_thetaz,
                       Rcpp::Named("tassobi") = tasso_bi,
                       Rcpp::Named("tassobi2") = tasso_bi2,
                       Rcpp::Named("tasso_prop_thetax") = tasso_prop_thetax,
                       Rcpp::Named("tasso_prop_thetaz") = tasso_prop_thetaz);
  
  return Rcpp::List::create(Rcpp::Named("thetaxposterior") = THETApost,
                            Rcpp::Named("thetazposterior") = THETApost2,
                            Rcpp::Named("muposterior") = ALPHApost,
                            Rcpp::Named("thrxposterior") = SOGLIApost,
                            Rcpp::Named("thrzposterior") = SOGLIApost2,
                            Rcpp::Named("loglinpred") = LOGLINPRED,
                            Rcpp::Named("loglik") = llvec,
                            Rcpp::Named("alpha0posterior") = LINCOEF,
                            Rcpp::Named("alphastarposterior") = NONLINCOEF,
                            Rcpp::Named("eta0posterior") = ALPHA_LIN,
                            Rcpp::Named("etastarposterior") = ALPHA_NL,
                            //Rcpp::Named("KSI_LIN") = KSI_LIN,
                            //Rcpp::Named("KSI_NL") = KSI_NL,
                            //Rcpp::Named("M_LIN") = M_LIN,
                            //Rcpp::Named("M_NL") = M_NL,
                            //Rcpp::Named("TAUSQ_LIN") = TAUSQ_LIN,
                            //Rcpp::Named("TAUSQ_NL") = TAUSQ_NL,
                            Rcpp::Named("xi0posterior") = GAMMA_LIN,
                            Rcpp::Named("xistarposterior") = GAMMA_NL,
                            //Rcpp::Named("OMEGA_LIN") = OMEGA_LIN,
                            //Rcpp::Named("OMEGA_NL") = OMEGA_NL,
                            Rcpp::Named("bxzposterior") = BIpost,
                            Rcpp::Named("bzzposterior") = BIpost2,
                            Rcpp::Named("betaxz") = hBIpost,
                            //Rcpp::Named("hbipost_ppi") = hBIpost_ppi,
                            Rcpp::Named("betaz") = hBIpost2,
                            Rcpp::Named("rint") = RANDINTpost,
                            Rcpp::Named("varrint") = VARRANDINTpost);
}//chiude funzione SAMPLER

//[[Rcpp::export]]
Rcpp::List sampler(arma::mat YY, arma::mat XX, arma::mat ZZ, arma::mat XXl, arma::mat XXT,
                   arma::vec dj, int Niter, int burn, int thin, arma::vec hyper_theta_x,
                   arma::vec hyper_theta_z, arma::vec penmig_lin, arma::vec penmig_nl,
                   arma::vec tx, arma::vec tz, arma::vec prior_int,
                   bool upsv, int hereditariety, bool conlnl, arma::vec theta_init,
                   arma::vec theta_init2){
  
  int Eff = (Niter-burn)/thin;
  int JJ = YY.n_cols;
  int PP = XX.n_cols;
  int KK = XXT.n_cols;
  int QQ = ZZ.n_cols;
  int n = XX.n_rows;
  
  int hh, ii, pp, jj, qq, ll, s;
  
  bool vincolo = conlnl;
  
  // adaptive proposal values for \thetaX & \thetaZ
  double last_mean;
  double last_var;
  double last_mean2;
  double last_var2;
  double last_meani;
  double last_vari;
  double last_meani2;
  double last_vari2;
  
  //MH acc ratios
  Rcpp::NumericVector acc_theta_flag0(PP*JJ);
  arma::vec acc_theta_flag(acc_theta_flag0.begin(), PP*JJ, false);
  acc_theta_flag.fill(0);
  
  Rcpp::NumericVector acc_theta_flag20(QQ*JJ);
  arma::vec acc_theta_flag2(acc_theta_flag20.begin(), QQ*JJ, false);
  acc_theta_flag2.fill(0);
  
  arma::vec acc_theta(PP*JJ);
  acc_theta.fill(0);
  
  arma::vec acc_theta2(QQ*JJ);
  acc_theta2.fill(0);
  
  Rcpp::NumericVector acc_alpha_flag0(JJ);
  arma::vec acc_alpha_flag(acc_alpha_flag0.begin(), JJ, false);
  acc_alpha_flag.fill(0);
  
  Rcpp::NumericVector acc_soglia_flag0(2);
  arma::vec acc_soglia_flag(acc_soglia_flag0.begin(), 2, false);
  acc_soglia_flag.fill(0);
  
  arma::vec acc_soglia(2);
  acc_soglia.fill(0);
  
  arma::vec acc_alpha(JJ);
  acc_alpha.fill(0);
  
  Rcpp::NumericMatrix acc_bi_flag0(JJ, PP*QQ);
  arma::mat acc_bi_flag(acc_bi_flag0.begin(), JJ, PP*QQ, false);
  acc_bi_flag.fill(0);
  
  Rcpp::NumericMatrix acc_bi_flag20(JJ, QQ*QQ);
  arma::mat acc_bi_flag2(acc_bi_flag20.begin(), JJ, QQ*QQ, false);
  acc_bi_flag2.fill(0);
  
  arma::mat acc_bi(JJ, PP*QQ);
  acc_bi.fill(0);
  
  arma::mat acc_bi2(JJ, QQ*QQ);
  acc_bi2.fill(0);
  
  //keep track of mean & var \tehtaX & \thetaZ
  arma::vec curr_mean(PP * JJ);
  curr_mean.fill(0.0);
  
  arma::vec curr_var(PP * JJ);
  curr_var.fill(0.5);
  
  arma::vec curr_mean2(QQ*JJ);
  curr_mean2.fill(0.0);
  
  arma::vec curr_var2(QQ*JJ);
  curr_var2.fill(0.5);
  
  arma::mat curr_meani(JJ, QQ*PP);
  curr_meani.fill(0.0);
  
  arma::mat curr_vari(JJ, QQ*PP);
  curr_vari.fill(.5);
  
  arma::mat curr_meani2(JJ, QQ * QQ);
  curr_meani2.fill(0.0);
  
  arma::mat curr_vari2(JJ, QQ * QQ);
  curr_vari2.fill(.5);
  
  //initialization for \thetaX & \thetaZ
  arma::vec theta(PP * JJ, arma::fill::zeros);
  
  arma::vec theta2(QQ * JJ, arma::fill::zeros);
  
  theta = theta_init;
  theta2 = theta_init2;
  
  //initialization for alpha
  Rcpp::NumericVector alpha0(JJ);
  arma::vec alpha(alpha0.begin(), JJ, false);
  alpha.fill(1.0);
  
  // variable inclusion indicator X & Z
  Rcpp::NumericVector inclusion_indicator0(PP*JJ);
  arma::vec inclusion_indicator(inclusion_indicator0.begin(), PP*JJ, false);
  for(jj = 0; jj < JJ; jj++){
    for(pp = 0; pp < PP; pp++){
      hh = jj + pp * JJ;
      if(theta(hh) != 0.0){
        inclusion_indicator(hh) = 1;
      } else {
        inclusion_indicator(hh) = 0;
      }
    }
  }
  
  Rcpp::NumericVector inclusion_indicator20(QQ*JJ);
  arma::vec inclusion_indicator2(inclusion_indicator20.begin(), QQ*JJ, false);
  for(jj = 0; jj < JJ; jj++){
    for(qq = 0; qq < QQ; qq++){
      hh = jj + qq * JJ;
      if(theta2(hh) != 0.0){
        inclusion_indicator2(hh) = 1;
      } else {
        inclusion_indicator2(hh) = 0;
      }
    }
  }
  
  //peNMIG
  arma::cube alpha_lin_xx(PP, PP, JJ, arma::fill::zeros);
  for(int jj = 0; jj < JJ; jj++){
    for(int pp = 0; pp < PP; pp++){
      for(int qq = pp + 1; qq < PP; qq++){
        alpha_lin_xx(pp, qq, jj) = 0.01;
      }
    }
  }
  
  arma::cube alpha_nl_xx(PP, PP, JJ, arma::fill::zeros);
  for(int jj = 0; jj < JJ; jj++){
    for(int pp = 0; pp < PP; pp++){
      for(int qq = pp + 1; qq < PP; qq++){
        alpha_nl_xx(pp, qq, jj) = 0.01;
      }
    }
  }
  
  arma::cube tausq_lin(PP, PP, JJ, arma::fill::zeros);
  for(int jj = 0; jj < JJ; jj++){
    for(int pp = 0; pp < PP; pp++){
      for(int qq = pp + 1; qq < PP; qq++){
        tausq_lin(pp, qq, jj) = 1.0;
      }
    }
  }
  
  arma::cube tausq_nl(PP, PP, JJ, arma::fill::zeros);
  for(int jj = 0; jj < JJ; jj++){
    for(int pp = 0; pp < PP; pp++){
      for(int qq = pp + 1; qq < PP; qq++){
        tausq_nl(pp, qq, jj) = 1.0;
      }
    }
  }
  
  arma::cube gamma_lin(PP, PP, JJ, arma::fill::zeros);
  for(int jj = 0; jj < JJ; jj++){
    for(int pp = 0; pp < PP; pp++){
      for(int qq = pp + 1; qq < PP; qq++){
        gamma_lin(pp, qq, jj) = 1.0;
      }
    }
  }
  
  arma::cube gamma_nl(PP, PP, JJ, arma::fill::zeros);
  for(int jj = 0; jj < JJ; jj++){
    for(int pp = 0; pp < PP; pp++){
      for(int qq = pp + 1; qq < PP; qq++){
        gamma_nl(pp, qq, jj) = 1.0;
      }
    }
  }
  
  arma::cube ksi_lin(PP, PP, JJ, arma::fill::zeros);
  for(int jj = 0; jj < JJ; jj++){
    for(int pp = 0; pp < PP; pp++){
      for(int qq = pp + 1; qq < PP; qq++){
        ksi_lin(pp, qq, jj) = 1.0;
      }
    }
  }
  
  arma::cube ksi_nl(PP, KK, JJ, arma::fill::zeros);
  arma::vec cdj(PP-1);
  cdj = cumsum(dj);
  for(int jj = 0; jj < JJ; jj++){
    for(int pp = 0; pp < PP; pp++){
      for(int qq = pp + 1; qq < PP; qq++){
        for(int ii = cdj(qq-1); ii < cdj(qq); ii++){
          ksi_nl(pp, ii, jj) = 1.0;
        }
      }
    }
  }
  
  arma::cube m_lin(PP, PP, JJ, arma::fill::zeros);
  for(int jj = 0; jj < JJ; jj++){
    for(int pp = 0; pp < PP; pp++){
      for(int qq = pp + 1; qq < PP; qq++){
        m_lin(pp, qq, jj) = 1.0;
      }
    }
  }
  
  arma::cube m_nl(PP, KK, JJ);
  for(int jj = 0; jj < JJ; jj++){
    for(int pp = 0; pp < PP; pp++){
      for(int qq = pp + 1; qq < PP; qq++){
        for(int ii = cdj(qq-1); ii < cdj(qq); ii++){
          m_nl(pp, qq, jj) = 1.0;
        }
      }
    }
  }
  
  arma::mat omega_lin(PP, JJ, arma::fill::zeros);
  for(int jj = 0; jj < JJ; jj++){
    for(int pp = 0; pp < PP; pp++){
      omega_lin(pp, jj) = 0.5;
    }
  }
  
  arma::mat omega_nl(PP, JJ, arma::fill::zeros);
  for(int jj = 0; jj < JJ; jj++){
    for(int pp = 0; pp < PP; pp++){
      omega_nl(pp, jj) = 0.5;
    }
  }
  
  Rcpp::List rescale_out;
  
  arma::cube bi_lin_xx(PP, PP, JJ);
  arma::cube bi_nl_xx(PP, KK, JJ);
  for(int jj = 0; jj < JJ; jj++){
    bi_lin_xx.slice(jj) = alpha_ksi_prod(alpha_lin_xx.slice(jj), ksi_lin.slice(jj), dj, true);
    bi_nl_xx.slice(jj) = alpha_ksi_prod(alpha_nl_xx.slice(jj), ksi_nl.slice(jj), dj, false);
  }
  
  double atau_lin = penmig_lin(0);
  double btau_lin = penmig_lin(1);
  
  double v0_lin = penmig_lin(2);
  
  double aomega = penmig_lin(3);
  double bomega = penmig_lin(4);
  
  double atau_nl = penmig_nl(0);
  double btau_nl = penmig_nl(1);
  
  double v0_nl = penmig_nl(2);
  
  double aomega_nl = penmig_nl(3);
  double bomega_nl = penmig_nl(4);
  
  //Y row sums
  arma::vec Ypiu(n);
  Ypiu.fill(0);
  
  for(int ii = 0; ii < n; ii ++){
    for(jj = 0; jj < JJ; jj ++){
      Ypiu(ii)=Ypiu(ii) + YY(ii, jj);
    }
  }
  
  //theta temp updates coeff matrix Ycategory-wise
  Rcpp::NumericVector theta_temp0(PP*JJ);
  arma::vec theta_temp(theta_temp0.begin(), PP*JJ, false);
  theta_temp.fill(0.0);
  
  Rcpp::NumericVector theta_temp20(QQ*JJ);
  arma::vec theta_temp2(theta_temp20.begin(), QQ*JJ, false);
  theta_temp2.fill(0.0);
  
  //matrici delle interazioni
  Rcpp::NumericMatrix bi0(PP*QQ*JJ);
  arma::mat bi_xz(bi0.begin(), JJ, PP*QQ, false);
  bi_xz.fill(0.0);
  arma::mat bXZ(PP, QQ);
  
  for(jj = 0 ; jj < JJ; jj ++){
    for(int qq2 = 0; qq2 < QQ; qq2 ++){
      for(int pp2 = 0; pp2 < PP; pp2 ++){
        int hh2 = qq2+((pp2)*QQ);
        bXZ(pp2, qq2)=bi_xz(jj, hh2);
      }
    }
  }
  
  //matrice delle interazioni temporanea
  Rcpp::NumericMatrix bi_p_mat0(PP, QQ);
  arma::mat bi_p_mat(bi_p_mat0.begin(), PP, QQ, false);
  bi_p_mat.fill(0.0);
  
  Rcpp::NumericMatrix bi20(JJ, QQ * QQ);
  arma::mat bi_zz(bi20.begin(), JJ, QQ * QQ, false);
  bi_zz.fill(0.0);
  arma::mat bZZ(QQ, QQ);
  
  for(jj = 0 ; jj < JJ; jj++){
    for(int qq2 = 0; qq2 < QQ; qq2++){
      for(int pp2 = 0; pp2 < QQ; pp2++){
        int hh2 = qq2+((pp2)*QQ);
        bZZ(pp2, qq2)=bi_zz(jj, hh2);
      }
    }
  }
  
  //matrice delle interazioni temporanea
  Rcpp::NumericMatrix bi_p_mat20(QQ, QQ);
  arma::mat bi_p_mat2(bi_p_mat20.begin(), QQ, QQ, false);
  bi_p_mat2.fill(0.0);
  
  //matrici per HS
  arma::mat lambda_hs(JJ, (PP*QQ));
  lambda_hs.fill(1.0);
  
  arma::vec tau_hs((QQ*JJ));
  tau_hs.fill(.1);
  
  arma::mat lambda_hs2(JJ, (QQ * QQ));
  lambda_hs2.fill(1.0);
  
  arma::vec tau_hs2((QQ*JJ));
  tau_hs2.fill(.1);
  
  double mu_hs = 0.0;
  
  //threshold parameters
  double soglia = 0.01;
  
  double lb = tx(0);
  double ub = tx(1);
  
  double soglia2 = 0.01;
  
  double lb2 = tz(0);
  double ub2 = tz(1);
  
  //linear predictor  LOG SCALE
  Rcpp::NumericMatrix loggamma0(n, JJ);
  arma::mat loggamma(loggamma0.begin(), n, JJ, false);
  
  
  for(jj = 0 ; jj < JJ; jj ++){
    for(int qq2 = 0; qq2 < QQ; qq2 ++){
      for(int pp2 = 0; pp2 < PP; pp2 ++){
        int hh2 = qq2+((pp2)*QQ);
        bXZ(pp2, qq2)=bi_xz(jj, hh2);
      }
      
      for(int pp2 = 0; pp2 < QQ; pp2 ++){
        int hh2 = qq2+((pp2)*QQ);
        bZZ(pp2, qq2)=bi_zz(jj, hh2);
      }
    }
    loggamma.col(jj) = calculate_gamma(XX, ZZ, XXl, XXT, alpha, theta, theta2, bXZ,
                 bi_lin_xx.slice(jj), bi_nl_xx.slice(jj), soglia, bZZ, soglia2, jj);
  }
  
  //latent variables
  arma::mat SS(n, JJ);
  SS.fill(0.0);
  arma::vec TT(n);
  TT.fill(0.0);
  for(ii = 0; ii < n; ii ++){
    for(jj = 0; jj < JJ; jj ++){
      SS(ii, jj) = (YY(ii, jj)/1.0);
      if(SS(ii, jj) < pow(10.0, -100.0)){
        SS(ii, jj) = pow(10.0, -100.0);
      }
      TT(ii)=TT(ii)+SS(ii, jj);
    }
  }
  
  arma::vec uu(n);
  for(ii = 0 ; ii < n ; ii ++){
    uu(ii) = ::Rf_rgamma(Ypiu(ii), 1.0/TT(ii));
  }
  
  //beta proposal variance
  Rcpp::NumericVector prop_per_theta0(PP*JJ);
  arma::vec prop_per_theta(prop_per_theta0.begin(), PP*JJ, false);
  prop_per_theta.fill(0.5);
  
  Rcpp::NumericVector prop_per_theta20(QQ*JJ);
  arma::vec prop_per_theta2(prop_per_theta20.begin(), QQ*JJ, false);
  prop_per_theta2.fill(0.5);
  
  //alpha proposal variance
  Rcpp::NumericVector prop_per_alpha0(JJ);
  arma::vec prop_per_alpha(prop_per_alpha0.begin(), JJ, false);
  prop_per_alpha.fill(0.5);
  
  //bi proposal variance
  Rcpp::NumericMatrix prop_per_bi0(JJ, PP*QQ);
  arma::mat prop_per_bi(prop_per_bi0.begin(), JJ, PP*QQ, false);
  prop_per_bi.fill(0.5);
  
  Rcpp::NumericMatrix prop_per_bi20(JJ, QQ * QQ);
  arma::mat prop_per_bi2(prop_per_bi20.begin(), JJ, QQ * QQ, false);
  prop_per_bi2.fill(0.5);
  
  //mean_slab
  arma::mat mu_theta(PP, JJ);
  mu_theta.fill(0.0);
  arma::mat mu_theta2(QQ, JJ);
  mu_theta2.fill(0.0);
  
  double tau_slab = hyper_theta_x(1);
  double tau_slab2 = hyper_theta_z(1);
  //slab_variance
  arma::mat sigma_theta(PP, JJ);
  sigma_theta.fill(tau_slab);
  arma::mat sigma_theta2(QQ, JJ);
  sigma_theta2.fill(tau_slab2);
  
  //Inverse Gamma Updating 1
  arma::mat sigma_theta_temp(PP, JJ);
  sigma_theta_temp = sigma_theta;
  
  double aig_slab = 3.0;
  double big_slab = .5;
  
  arma::vec idx_tmp(PP);
  arma::vec thetavec_tmp(PP);
  double thetavecTthetavec;
  int sum_idx_tmp;
  
  //Inverse Gamma Updating 2
  arma::mat sigma_theta_temp2(QQ, JJ);
  sigma_theta_temp2 = sigma_theta2;
  
  double aig_slab2 = 3.0;
  double big_slab2 = .5;
  
  arma::vec idx_tmp2(QQ);
  arma::vec thetavec_tmp2(QQ);
  double thetavecTthetavec2;
  int sum_idx_tmp2;
  
  // mean and standard deviation of the independent normal priors on alpha and beta
  arma:: vec mu_al(JJ);
  mu_al.fill(prior_int(0));
  arma::vec sig_al(JJ);
  sig_al.fill(prior_int(1));
  
  double aa_hp = hyper_theta_x(0);
  double bb_hp = 2 - aa_hp;
  
  double aa_hp2 = hyper_theta_z(0);
  double bb_hp2 = 2 - aa_hp2;
  
  arma::cube THETApost(PP, JJ, Eff);
  arma::cube THETApost2(QQ, JJ, Eff);
  arma::cube BIpost(JJ, (PP*QQ), Eff);
  arma::cube hBIpost(n, PP, JJ, arma::fill::zeros);
  arma::cube hBIpost_ppi(n, PP, JJ, arma::fill::zeros);
  arma::cube BIpost2(JJ, (QQ * QQ), Eff);
  arma::cube hBIpost2(n, QQ, JJ, arma::fill::zeros);
  arma::cube LOGLINPRED(n, JJ, Eff);
  arma::vec llvec(Eff);
  arma::vec SOGLIApost(Eff);
  arma::vec SOGLIApost2(Eff);
  arma::mat ALPHApost(Eff, JJ);
  
  arma::field<arma::cube> LINCOEF(Eff);
  LINCOEF.fill(arma::cube(PP, PP, JJ, arma::fill::zeros));
  arma::field<arma::cube> ALPHA_LIN(Eff);
  ALPHA_LIN.fill(arma::cube(PP, PP, JJ, arma::fill::zeros));
  arma::field<arma::cube> M_LIN(Eff);
  M_LIN.fill(arma::cube(PP, PP, JJ, arma::fill::zeros));
  arma::field<arma::cube> KSI_LIN(Eff);
  KSI_LIN.fill(arma::cube(PP, PP, JJ, arma::fill::zeros));
  arma::field<arma::cube> TAUSQ_LIN(Eff);
  TAUSQ_LIN.fill(arma::cube(PP, PP, JJ, arma::fill::zeros));
  arma::field<arma::cube> GAMMA_LIN(Eff);
  GAMMA_LIN.fill(arma::cube(PP, PP, JJ, arma::fill::zeros));
  arma::cube OMEGA_LIN(PP, JJ, Eff, arma::fill::zeros);
  
  arma::field<arma::cube> NONLINCOEF(Eff);
  NONLINCOEF.fill(arma::cube(PP, KK, JJ, arma::fill::zeros));
  arma::field<arma::cube> ALPHA_NL(Eff);
  ALPHA_NL.fill(arma::cube(PP, PP, JJ, arma::fill::zeros));
  arma::field<arma::cube> M_NL(Eff);
  M_NL.fill(arma::cube(PP, KK, JJ, arma::fill::zeros));
  arma::field<arma::cube> KSI_NL(Eff);
  KSI_NL.fill(arma::cube(PP, KK, JJ, arma::fill::zeros));
  arma::field<arma::cube> TAUSQ_NL(Eff);
  TAUSQ_NL.fill(arma::cube(PP, PP, JJ, arma::fill::zeros));
  arma::field<arma::cube> GAMMA_NL(Eff);
  GAMMA_NL.fill(arma::cube(PP, PP, JJ, arma::fill::zeros));
  arma::cube OMEGA_NL(PP, JJ, Eff, arma::fill::zeros);
  
  Rcpp::List listainter;
  int curr_iter = 0;
  int zz = 0;
  
  for( s = 0; s < Niter+1; s++){
    
    ////////////////////////
    // AGGIORNO I THETA X //
    ////////////////////////
    
    //scelgo un taxa a caso per fare lo swap
    double u;
    u = arma::randu();
    u = u*(JJ/1.0);
    jj = Rcpp::as<int>(Rcpp::wrap(std::floor(u)));
    
    for(pp = 0; pp < PP; pp ++){
      int ci = jj + pp * JJ;
      theta_temp(ci) = theta(ci);
    }
    
    swap(XX, ZZ, SS, loggamma, alpha, theta_temp, theta_temp2, theta, theta2,
         bi_xz, soglia, bi_zz, soglia2, acc_theta_flag, acc_theta_flag2, inclusion_indicator,
         inclusion_indicator2, mu_theta, mu_theta2, sigma_theta, sigma_theta2,
         aa_hp, bb_hp, aa_hp2, bb_hp2, jj, true);
    
    // update theta with theta_temp
    for(pp = 0 ; pp < PP; pp ++){
      int ci = jj + pp * JJ;
      theta(ci) = theta_temp(ci);
    }
    //}
    
    //within model
    for(jj = 0; jj < JJ; jj ++){
      for(pp = 0; pp < PP; pp ++){
        int ci = jj + pp * JJ;
        theta_temp(ci) = theta(ci);
        
        if(s > (std::floor(burn/double(4)))){
          last_mean = curr_mean(ci);
          last_var = curr_var(ci);
          curr_mean(ci) = online_mean(s, last_mean, theta(ci));
          curr_var(ci) = online_var(s, last_mean, last_var, curr_mean(ci), theta(ci));
          
          //update proposal variance
          prop_per_theta(ci)=curr_var(ci);
        }
      }//chiude for pp
      
      //aggiorno i thetaX
      update_theta(XX, ZZ, SS, loggamma, alpha, theta_temp, theta_temp2, theta,
                   theta2, bi_xz, soglia, bi_zz, soglia2, inclusion_indicator,
                   inclusion_indicator2, prop_per_theta, prop_per_theta2,
                   mu_theta, mu_theta2, sigma_theta, sigma_theta2, aa_hp,
                   bb_hp, aa_hp2, bb_hp2, jj, true);
      
      for(pp = 0; pp < PP; pp ++){
        int ci = jj + pp * JJ;
        theta(ci) = theta_temp(ci);
        if((s > burn) && (s % thin == 0)){
          acc_theta(ci) = acc_theta(ci) + acc_theta_flag(ci);
          acc_theta_flag(ci) = 0;
        }
      }
    }//chiude for sui jj
    
    //qui salvo output
    if((s > burn) && (s % thin == 0)){
      zz = 0;
      for(pp = 0; pp < PP; pp++){
        for(jj = 0; jj < JJ; jj++){
          THETApost(pp, jj, curr_iter) = theta(zz);
          zz++;
        }
      }
    }
    
    if(upsv == true){
      //aggiorno inverse gamma
      for(jj = 0; jj < JJ; jj++){
        for(pp = 0; pp < PP; pp++){
          
          hh = jj + pp * JJ;
          idx_tmp(pp) = inclusion_indicator(hh);
          thetavec_tmp(pp) = theta(hh);
        }
        
        thetavecTthetavec = inner(thetavec_tmp, thetavec_tmp);
        sum_idx_tmp = arma::sum(idx_tmp);
        double rtauj;
        if(sum_idx_tmp == 0.0){
          rtauj = 1.0/::Rf_rgamma(aig_slab, big_slab);
        } else {
          rtauj = 1.0/::Rf_rgamma(0.5 * aig_slab + 0.5 * sum_idx_tmp, 0.5 * pow(big_slab, 2.0) + 0.5 * thetavecTthetavec);
        }
        sigma_theta_temp.col(jj).fill(rtauj);
      }
      sigma_theta = sigma_theta_temp;
    }
    
    ////////////////////////
    // AGGIORNO I THETA Z //
    ////////////////////////
    
    //scelgo un taxa a caso per fare lo swap
    u = arma::randu();
    u = u*(JJ/1.0);
    jj = Rcpp::as<int>(Rcpp::wrap(std::floor(u)));
    for(qq = 0; qq < QQ; qq ++){
      int ci = jj + qq * JJ;
      theta_temp2(ci)=theta2(ci);
    }
    
    swap(XX, ZZ, SS, loggamma, alpha, theta_temp, theta_temp2, theta, theta2,
         bi_xz, soglia, bi_zz, soglia2, acc_theta_flag, acc_theta_flag2, inclusion_indicator,
         inclusion_indicator2, mu_theta, mu_theta2, sigma_theta, sigma_theta2,
         aa_hp, bb_hp, aa_hp2, bb_hp2, jj, false);
    
    // update theta with theta_temp
    for(qq = 0 ; qq < QQ ; qq++){
      int ci = jj + qq * JJ;
      theta2(ci) = theta_temp2(ci);
    }
    
    //within model
    for(jj = 0; jj < JJ; jj ++){
      for(qq = 0; qq < QQ; qq ++){
        int ci = jj + qq * JJ;
        theta_temp2(ci)= theta2(ci);
        //aggiorno media e varianza adaptive
        if(s > (std::floor(burn/double(4)))){
          last_mean2 = curr_mean2(ci);
          last_var2 = curr_var2(ci);
          curr_mean2(ci) = online_mean(s, last_mean2, theta2(ci));
          curr_var2(ci) = online_var(s, last_mean2, last_var2, curr_mean2(ci),
                    theta2(ci));
          
          //update proposal variance
          prop_per_theta2(ci)=curr_var2(ci);
        }
      }//chiude for pp
      
      //aggiorno i thetaZ
      update_theta(XX, ZZ, SS, loggamma, alpha, theta_temp, theta_temp2, theta,
                   theta2, bi_xz, soglia, bi_zz, soglia2, inclusion_indicator, inclusion_indicator2,
                   prop_per_theta, prop_per_theta2, mu_theta, mu_theta2,
                   sigma_theta, sigma_theta2, aa_hp, bb_hp, aa_hp2, bb_hp2, jj,
                   false);
      
      for(qq = 0; qq < QQ; qq++){
        int ci = jj + qq * JJ;
        theta2(ci) = theta_temp2(ci);
        if((s>burn) && (s%thin==0)){
          acc_theta2(ci) = acc_theta2(ci) + acc_theta_flag2(ci);
          acc_theta_flag2(ci) = 0;
        }
      }
    }//chiude for sui jj
    
    //qui salvo output
    if((s > burn) && (s % thin == 0)){
      zz = 0;
      for(qq = 0; qq < QQ; qq++){
        for(jj = 0; jj < JJ; jj++){
          THETApost2(qq, jj, curr_iter) = theta2(zz);
          zz++;
        }
      }
    }
    
    if(upsv == true){
      //aggiorno inverse gamma
      for(jj = 0; jj < JJ; jj++){
        for(qq = 0; qq < QQ; qq++){
          
          hh = jj + qq * JJ;
          idx_tmp2(qq) = inclusion_indicator2(hh);
          thetavec_tmp2(qq) = theta2(hh);
        }
        
        thetavecTthetavec2 = inner(thetavec_tmp2, thetavec_tmp2);
        sum_idx_tmp2 = arma::sum(idx_tmp2);
        
        double rtau2j;
        if(sum_idx_tmp2 == 0.0){
          rtau2j = 1.0/::Rf_rgamma(aig_slab2, big_slab2);
        } else {
          rtau2j = 1.0/::Rf_rgamma(0.5 * aig_slab2 + 0.5 * sum_idx_tmp2, 0.5 * pow(big_slab2, 2.0) + 0.5 * thetavecTthetavec2);
        }
        sigma_theta_temp2.col(jj).fill(rtau2j);
      }
      sigma_theta2 = sigma_theta_temp2;
    }
    
    ///////////////////////////
    // AGGIORNO INTERAZIONI1 //
    ///////////////////////////
    
    u = arma::randu();
    u = u*(JJ/1.0);
    jj = Rcpp::as<int>(Rcpp::wrap(std::floor(u)));
    
    for(int qq2 = 0; qq2 < QQ; qq2++){
      for(int pp2 = 0; pp2 < PP; pp2++){
        int hh2 = qq2 + ((pp2) * QQ);
        bi_p_mat(pp2, qq2)= bi_xz(jj, hh2);// + ::Rf_rnorm(0.0, .1);
      }
    }
    
    // ADAPTIVE PROPOSAL
    if(s > 2000){
      for(pp = 0; pp < PP; pp ++){
        for(qq = 0; qq < QQ; qq++){
          
          int hh2 = qq + pp * QQ;
          last_meani = curr_meani(jj, hh2);
          last_vari = curr_vari(jj, hh2);
          curr_meani(jj, hh2) = online_mean(s, last_meani, bi_xz(jj, hh2));
          curr_vari(jj, hh2) = online_var(s, last_meani, last_vari, curr_meani(jj, hh2), bi_xz(jj, hh2));
          //update proposal variance
          prop_per_bi(jj, hh2)=curr_vari(jj, hh2);
          
          
        }//chiude ciclo su qq
      }
    }//chiude if
    
    for(pp = 0; pp < PP; pp ++){
      for(qq = 0; qq < QQ; qq++){
        update_bi(XX, ZZ, XXl, XXT, SS, loggamma, lambda_hs, tau_hs, alpha, theta, theta2,
                  bi_xz, bi_zz, bi_p_mat, bi_lin_xx.slice(jj), bi_nl_xx.slice(jj), soglia, soglia2,
                  prop_per_bi, acc_bi_flag, mu_hs, pp, qq, jj, hereditariety);//hereditariety
        
        hh = qq + pp * QQ;
        if((s>burn) && (s%thin==0)){
          acc_bi(jj, hh)= acc_bi(jj, hh)+acc_bi_flag(jj, hh);
          acc_bi_flag(jj, hh) = 0;
        }
      }//chiude il ciclo sulla qq
    }
    
    if((s > burn)&& (s % thin == 0)){
      for(jj = 0; jj < JJ; jj ++){
        for(hh = 0; hh < (PP * QQ); hh ++){
          BIpost(jj, hh, curr_iter) = bi_xz(jj, hh);
        }
      }
    }
    
    ///////////////////////////
    // AGGIORNO HYP HS PRIOR //
    ///////////////////////////
    
    arma::vec bi_work(QQ);
    
    u = arma::randu();
    u = u*(double)JJ;
    jj = Rcpp::as<int>(Rcpp::wrap(std::floor(u)));
    
    u = arma::randu();
    u = u*(double)PP;
    pp = Rcpp::as<int>(Rcpp::wrap(std::floor(u)));
    
    for(qq = 0; qq < QQ; qq++){
      hh = qq + pp * QQ;
      bi_work(qq) = bi_xz(jj, hh);
    }
    
    u = arma::randu();
    u = u*(double)QQ;
    qq = Rcpp::as<int>(Rcpp::wrap(std::floor(u)));
    ll = jj + qq * JJ;
    
    lambda_hs.row(jj) = uplam(lambda_hs.row(jj), bi_work, pp, tau_hs(ll));
    tau_hs(ll) = uptau(lambda_hs.row(jj).t(), bi_work, pp, QQ, tau_hs(ll), true);
    
    //peNMIG parte lineare
    u = arma::randu();
    u = u*(double)JJ;
    jj = Rcpp::as<int>(Rcpp::wrap(std::floor(u)));
    
    alpha_lin_xx.slice(jj) = update_alpha_xx(XX, ZZ, XXl, XXT, SS, loggamma, alpha,
                       theta, theta2, bi_xz, soglia, bi_zz, soglia2, dj, alpha_lin_xx.slice(jj), alpha_nl_xx.slice(jj),
                       tausq_lin.slice(jj), tausq_nl.slice(jj), gamma_lin.slice(jj),
                       gamma_nl.slice(jj), ksi_lin.slice(jj), ksi_nl.slice(jj), true, jj, hereditariety, vincolo);
    m_lin.slice(jj) = update_m(m_lin.slice(jj), ksi_lin.slice(jj));
    
    ksi_lin.slice(jj) = update_ksi(XX, ZZ, XXl, XXT, SS, loggamma, alpha, theta, theta2,
                  bi_xz, soglia, bi_zz, soglia2, dj, alpha_lin_xx.slice(jj), alpha_nl_xx.slice(jj),
                  tausq_lin.slice(jj), tausq_nl.slice(jj), gamma_lin.slice(jj),
                  gamma_nl.slice(jj), ksi_lin.slice(jj), ksi_nl.slice(jj), m_lin.slice(jj),
                  m_nl.slice(jj), true, jj);
    
    rescale_out = rescale(alpha_lin_xx.slice(jj), ksi_lin.slice(jj), dj, true);
    alpha_lin_xx.slice(jj) = Rcpp::as<arma::mat>(rescale_out[0]);
    
    ksi_lin.slice(jj) = Rcpp::as<arma::mat>(rescale_out[1]);
    
    bi_lin_xx.slice(jj) = alpha_ksi_prod(alpha_lin_xx.slice(jj),
                    ksi_lin.slice(jj), dj, true);
    
    tausq_lin.slice(jj) = update_tau(tausq_lin.slice(jj), alpha_lin_xx.slice(jj),
                    atau_lin, btau_lin, gamma_lin.slice(jj), jj);
    
    gamma_lin.slice(jj) = update_gamma(alpha_lin_xx.slice(jj),
                    tausq_lin.slice(jj), gamma_lin.slice(jj),
                    omega_lin.col(jj), v0_lin, jj, JJ, theta, hereditariety);
    
    omega_lin.col(jj) = update_omega(gamma_lin.slice(jj), omega_lin.col(jj),
                  v0_lin, aomega, bomega);
    
    if((s>burn) && (s%thin==0)){
      ALPHA_LIN[curr_iter] = alpha_lin_xx;
      M_LIN[curr_iter] = m_lin;
      KSI_LIN[curr_iter] = ksi_lin;
      LINCOEF[curr_iter] = bi_lin_xx;
      TAUSQ_LIN[curr_iter] = tausq_lin;
      GAMMA_LIN[curr_iter] = gamma_lin;
      OMEGA_LIN.slice(curr_iter) = omega_lin;
    }
    
    //peNMIG non linear
    u = arma::randu();
    u = u*(double)JJ;
    jj = Rcpp::as<int>(Rcpp::wrap(std::floor(u)));
    alpha_nl_xx.slice(jj) = update_alpha_xx(XX, ZZ, XXl, XXT, SS, loggamma, alpha,
                      theta, theta2, bi_xz, soglia, bi_zz, soglia2, dj, alpha_lin_xx.slice(jj), alpha_nl_xx.slice(jj),
                      tausq_lin.slice(jj), tausq_nl.slice(jj), gamma_lin.slice(jj),
                      gamma_nl.slice(jj), ksi_lin.slice(jj), ksi_nl.slice(jj), false, jj, hereditariety, vincolo);
    
    m_nl.slice(jj) = update_m(m_nl.slice(jj), ksi_nl.slice(jj));
    
    ksi_nl.slice(jj) = update_ksi(XX, ZZ, XXl, XXT, SS, loggamma, alpha, theta, theta2,
                 bi_xz, soglia, bi_zz, soglia2, dj, alpha_lin_xx.slice(jj), alpha_nl_xx.slice(jj),
                 tausq_lin.slice(jj), tausq_nl.slice(jj), gamma_lin.slice(jj),
                 gamma_nl.slice(jj), ksi_lin.slice(jj), ksi_nl.slice(jj), m_lin.slice(jj),
                 m_nl.slice(jj), false, jj);
    
    rescale_out = rescale(alpha_nl_xx.slice(jj), ksi_nl.slice(jj), dj, false);
    alpha_nl_xx.slice(jj) = Rcpp::as<arma::mat>(rescale_out[0]);
    ksi_nl.slice(jj) = Rcpp::as<arma::mat>(rescale_out[1]);
    bi_nl_xx.slice(jj) = alpha_ksi_prod(alpha_nl_xx.slice(jj),
                   ksi_nl.slice(jj), dj, false);
    
    tausq_nl.slice(jj) = update_tau(tausq_nl.slice(jj), alpha_nl_xx.slice(jj),
                   atau_nl, btau_nl, gamma_nl.slice(jj), jj);
    
    gamma_nl.slice(jj) = update_gamma(alpha_nl_xx.slice(jj),
                   tausq_nl.slice(jj), gamma_nl.slice(jj), omega_nl.col(jj),
                   v0_nl, jj, JJ, theta, hereditariety);
    
    omega_nl.col(jj) = update_omega(gamma_nl.slice(jj), omega_nl.col(jj),
                 v0_nl, aomega_nl, bomega_nl);
    
    if((s>burn) && (s%thin==0)){
      ALPHA_NL[curr_iter] = alpha_nl_xx;
      M_NL[curr_iter] = m_nl;
      KSI_NL[curr_iter] = ksi_nl;
      NONLINCOEF[curr_iter] = bi_nl_xx;
      TAUSQ_NL[curr_iter] = tausq_nl;
      GAMMA_NL[curr_iter] = gamma_nl;
      OMEGA_NL.slice(curr_iter) = omega_nl;
    }
    
    ///////////////////////////
    //    AGGIORNO SOGLIA1   //
    ///////////////////////////
    
    soglia = update_soglia(XX, ZZ, XXl, XXT, SS, loggamma, alpha, theta, theta2, bi_xz,
                           bi_lin_xx, bi_nl_xx, soglia, bi_zz, soglia2, acc_soglia_flag, lb, ub);
    
    if((s>burn) && (s%thin==0)){
      acc_soglia(0)=acc_soglia(0)+acc_soglia_flag(0);
      acc_soglia_flag(0)=0;
      SOGLIApost(curr_iter) = soglia;
    }
    
    ///////////////////////////
    //    Eff. Sub-Spec. 1   //
    ///////////////////////////
    
    if((s > burn)&& (s % thin == 0)){
      arma::mat bXZ(PP, QQ);
      arma::mat work(n, PP);
      
      for(jj = 0; jj < JJ; jj++){
        for(int qq2 = 0; qq2 < QQ; qq2++){
          for(int pp2 = 0; pp2 < PP; pp2++){
            int hh2 = qq2 + ((pp2) * QQ);
            bXZ(pp2, qq2) = bi_xz(jj, hh2);
          }
        }
        work = (bXZ * ZZ.t()).t() + (bi_lin_xx.slice(jj) * XXl.t()).t()
          + (bi_nl_xx.slice(jj) * XXT.t()).t();
        hBIpost.slice(jj) = hBIpost.slice(jj) + threshold_mat(work, soglia);
      }//nPJ
    }
    
    ///////////////////////////
    // AGGIORNO INTERAZIONI2 //
    ///////////////////////////
    
    u = arma::randu();
    u = u*(JJ/1.0);
    jj = Rcpp::as<int>(Rcpp::wrap(std::floor(u)));
    
    for(int qq2 = 0; qq2 < QQ; qq2++){
      for(int pp2 = 0; pp2 < QQ; pp2++){
        int hh2 = qq2 + ((pp2) * QQ);
        bi_p_mat2(pp2, qq2)= bi_zz(jj, hh2);
      }
    }
    
    // ADAPTIVE PROPOSAL
    if(s > 2000){
      for(pp = 0; pp < QQ; pp ++){
        for(qq = (pp + 1); qq < QQ; qq++){
          
          int hh2 = qq + pp * QQ;
          last_meani2 = curr_meani2(jj, hh2);
          last_vari2 = curr_vari2(jj, hh2);
          curr_meani2(jj, hh2) = online_mean(s, last_meani2, bi_zz(jj, hh2));
          curr_vari2(jj, hh2) = online_var(s, last_meani2, last_vari2, curr_meani2(jj, hh2), bi_zz(jj, hh2));
          //update proposal variance
          prop_per_bi2(jj, hh2)=curr_vari2(jj, hh2);
          
          
        }//chiude ciclo su qq
      }
    }//chiude if
    
    for(pp = 0; pp < QQ; pp ++){
      for(qq = (pp + 1); qq < QQ; qq++){
        update_bi2(XX, ZZ, SS, loggamma, lambda_hs2, tau_hs2, alpha, theta, theta2,
                   bi_xz, bi_zz, bi_p_mat2, soglia, soglia2, prop_per_bi2,
                   acc_bi_flag2, mu_hs, pp, qq, jj, hereditariety);//hereditariety
      }
    }
    
    for(pp = 0; pp < QQ; pp ++){
      for(qq = 0; qq < QQ; qq++){
        hh = qq + pp * QQ;
        if((s>burn) && (s%thin==0)){
          acc_bi2(jj, hh)= acc_bi2(jj, hh)+acc_bi_flag2(jj, hh);
          acc_bi_flag2(jj, hh) = 0;
        }
      }//chiude il ciclo sulla qq
    }
    
    if((s > burn)&& (s % thin == 0)){
      for(jj = 0; jj < JJ; jj ++){
        for(hh = 0; hh < (QQ * QQ); hh ++){
          BIpost2(jj, hh, curr_iter) = bi_zz(jj, hh);
        }
      }
    }
    
    ///////////////////////////
    // AGGIORNO HYP HS PRIOR //
    ///////////////////////////
    
    arma::vec bi_work2(QQ);
    
    u = arma::randu();
    u = u*(double)JJ;
    jj = Rcpp::as<int>(Rcpp::wrap(std::floor(u)));
    
    u = arma::randu();
    u = u*(double)QQ;
    pp = Rcpp::as<int>(Rcpp::wrap(std::floor(u)));
    
    for(qq = 0; qq < QQ; qq++){
      hh = qq + pp * QQ;
      bi_work2(qq) = bi_zz(jj, hh);
    }
    
    u = arma::randu();
    u = u*(double)QQ;
    qq = Rcpp::as<int>(Rcpp::wrap(std::floor(u)));
    ll = jj + qq * JJ;
    
    lambda_hs2.row(jj) = uplam(lambda_hs2.row(jj), bi_work2, pp, tau_hs2(ll));
    tau_hs2(ll) = uptau(lambda_hs2.row(jj).t(), bi_work2, qq, QQ, tau_hs2(ll), true);
    
    ///////////////////////////
    //    AGGIORNO SOGLIA2   //
    ///////////////////////////
    
    soglia2 = update_soglia(XX, ZZ, XXl, XXT, SS, loggamma, alpha, theta, theta2,
                            bi_xz, bi_lin_xx, bi_nl_xx, soglia, bi_zz, soglia2, acc_soglia_flag, lb2, ub2, 1);
    
    if((s>burn) && (s%thin==0)){
      acc_soglia(1)=acc_soglia(1)+acc_soglia_flag(1);
      acc_soglia_flag(1)=0;
      SOGLIApost2(curr_iter) = soglia2;
    }
    
    ///////////////////////////
    //    Eff. Sub-Spec. 2   //
    ///////////////////////////
    
    if((s > burn)&& (s % thin == 0)){
      arma::mat bZZ(QQ, QQ);
      arma::mat work(n, QQ);
      
      for(jj = 0; jj < JJ; jj++){
        for(int qq2 = 0; qq2 < QQ; qq2++){
          for(int pp2 = 0; pp2 < QQ; pp2++){
            int hh2 = qq2 + ((pp2) * QQ);
            bZZ(pp2, qq2) = bi_zz(jj, hh2);
          }
        }
        work = (bZZ * ZZ.t()).t();
        hBIpost2.slice(jj) = hBIpost2.slice(jj) + threshold_mat(work, soglia);
      }
    }
    
    //////////////////////////
    // AGGIORNO INTERCETTE ///
    //////////////////////////
    
    for(jj = 0 ; jj < JJ ; jj++){
      update_alpha(XX, ZZ, SS, loggamma, alpha, theta, theta2, bi_xz, bi_zz, soglia, soglia2, prop_per_alpha, acc_alpha_flag,
                   mu_al, sig_al, jj);
      
      if((s>burn) && (s%thin==0)){
        acc_alpha(jj)=acc_alpha(jj)+acc_alpha_flag(jj);
        acc_alpha_flag(jj)=0;
      }
    }
    
    //qui salvo output
    if((s > burn) && (s % thin == 0)){
      ALPHApost.row(curr_iter) = alpha.t();
    }
    
    llvec(curr_iter) = logliksimpleC(YY, exp(loggamma));
    
    if((s > burn) && (s % thin == 0)){
      LOGLINPRED.slice(curr_iter) = loggamma;
      curr_iter++;
    }
    
    // update JJ and consequently TT
    for(ii = 0 ; ii < n ; ii++){
      TT(ii) = 0.0;
      for(jj = 0 ; jj < JJ ; jj++){
        SS(ii,jj) = ::Rf_rgamma(YY(ii, jj)+ exp(loggamma(ii,jj)), 1.0/(uu(ii) + 1.0));
        if(SS(ii,jj) < pow(10.0, -100.0)){
          SS(ii,jj) = pow(10.0, -100.0);
        }
        TT(ii) = TT(ii) + SS(ii, jj);
      }
    }
    
    // update latent variables uu
    for(ii = 0 ; ii < n ; ii++){
      uu(ii) = ::Rf_rgamma(Ypiu(ii), 1.0/TT(ii));
    }
  }//chiude le iterazioni della catena
  
  arma::vec tasso_alpha(JJ);
  arma::mat tasso_thetax(PP, JJ);
  arma::mat tasso_bi(JJ, PP * QQ);
  arma::mat tasso_bi2(JJ, QQ * QQ);
  arma::mat tasso_prop_thetax(PP, JJ);
  arma::mat tasso_thetaz(QQ, JJ);
  arma::mat tasso_prop_thetaz(QQ, JJ);
  
  for(jj = 0; jj < JJ; jj++){
    tasso_alpha(jj) = (double)acc_alpha(jj)/Eff;
    for(pp=0;pp<PP;pp++){
      tasso_thetax(pp, jj) = (double)acc_theta(pp+ jj*PP)/Eff;
      tasso_prop_thetax(pp, jj) = prop_per_theta(pp+jj*PP)/Eff;
    }
    
    for(qq = 0; qq < QQ; qq ++){
      tasso_thetaz(qq, jj) = (double)acc_theta2(qq+ jj*QQ)/Eff;
      tasso_prop_thetaz(qq, jj) = prop_per_theta2(qq+jj*QQ)/Eff;
    }
    
    for(hh = 0; hh < (PP*QQ); hh ++){
      tasso_bi(jj, hh) = (double)acc_bi(jj, hh)/Eff;
    }
    for(hh = 0; hh < (QQ * QQ); hh ++){
      tasso_bi2(jj, hh) = (double)acc_bi2(jj, hh)/Eff;
    }
  }
  
  Rcpp::List lista_acc_rates =
    Rcpp::List::create(Rcpp::Named("tassoalpha") = tasso_alpha,
                       Rcpp::Named("tassothetax") = tasso_thetax,
                       Rcpp::Named("tassothetaz") = tasso_thetaz,
                       Rcpp::Named("tassobi") = tasso_bi,
                       Rcpp::Named("tassobi2") = tasso_bi2,
                       Rcpp::Named("tasso_prop_thetax") = tasso_prop_thetax,
                       Rcpp::Named("tasso_prop_thetaz") = tasso_prop_thetaz);
  
  return Rcpp::List::create(Rcpp::Named("thetaxposterior") = THETApost,
                            Rcpp::Named("thetazposterior") = THETApost2,
                            Rcpp::Named("muposterior") = ALPHApost,
                            Rcpp::Named("thrxposterior") = SOGLIApost,
                            Rcpp::Named("thrzposterior") = SOGLIApost2,
                            Rcpp::Named("loglinpred") = LOGLINPRED,
                            Rcpp::Named("loglik") = llvec,
                            Rcpp::Named("alpha0posterior") = LINCOEF,
                            Rcpp::Named("alphastarposterior") = NONLINCOEF,
                            Rcpp::Named("eta0posterior") = ALPHA_LIN,
                            Rcpp::Named("etastarposterior") = ALPHA_NL,
                            //Rcpp::Named("KSI_LIN") = KSI_LIN,
                            //Rcpp::Named("KSI_NL") = KSI_NL,
                            //Rcpp::Named("M_LIN") = M_LIN,
                            //Rcpp::Named("M_NL") = M_NL,
                            //Rcpp::Named("TAUSQ_LIN") = TAUSQ_LIN,
                            //Rcpp::Named("TAUSQ_NL") = TAUSQ_NL,
                            Rcpp::Named("xi0posterior") = GAMMA_LIN,
                            Rcpp::Named("xistarposterior") = GAMMA_NL,
                            //Rcpp::Named("OMEGA_LIN") = OMEGA_LIN,
                            //Rcpp::Named("OMEGA_NL") = OMEGA_NL,
                            Rcpp::Named("bxzposterior") = BIpost,
                            Rcpp::Named("bzzposterior") = BIpost2,
                            Rcpp::Named("betaxz") = hBIpost,
                            //Rcpp::Named("hbipost_ppi") = hBIpost_ppi,
                            Rcpp::Named("betaz") = hBIpost2);
  
}//chiude funzione SAMPLER
