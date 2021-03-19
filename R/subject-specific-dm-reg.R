#' Subject Specific Dirichlet Multinomial Regression from Microbiota analysis
#'
#' The implementation has been done in \code{C++} through the use of \code{Rcpp} and \code{RcppArmadillo}.
#' @author
#' Matteo Pedone
#'
#' @docType package
#' @name subject-specific-dm-reg
#'
#' @useDynLib subject-specific-dm-reg, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @importFrom stats cor.test
#' @importFrom stats cor
#' @importFrom stats p.adjust
#' @importFrom spikeSlabGAM sm
#' @importFrom dplyr between 
#' @importFrom plotrix plotCI
#' @useDynLib subject-specific-dm-reg, .registration = TRUE
NULL
