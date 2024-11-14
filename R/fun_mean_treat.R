#' Estimate treatment effect on the treatment arm using different covariates adjustment methods.
#'
#' @param X The n by p covariates matrix
#' @param Y Vector of n dimensional observed response.
#' @param A Vector of n dimensional treatment assignment.
#' @param H The n by n hat matrix corresponding to X.
#'
#' @return
#' \item{tau_unadj}{The point estimation of unadjusted estimator tau_unadj.}
#' \item{tau_adj}{The point estimation of Lin’s (2013) adjusted estimator with linear working model.}
#' \item{tau_db}{The point estimation of debiased estimator tau_db by Lu et al.(2023).}
#' \item{tau_adj2}{The point estimation of HOIF-motivated adjusted estimator.}
#' \item{tau_adj2_dag}{The point estimation of another form of debiased estimator tau_db.}
#' \item{tau_adj3}{The point estimation of bias-free adjusted estimator based on tau_adj2.}
#' @export
#'
#' @references
#' Winston Lin. Agnostic notes on regression adjustments to experimental data: Reexamining Freedman’s critique. The Annals of Statistics, 7(1):295–318, 2013.
#'
#' Xin Lu, Fan Yang, and Yuhao Wang. Debiased regression adjustment in completely randomized experiments with moderately high-dimensional covariates. arXiv preprint arXiv:2309.02073, 2023.
#'
#' @examples
#' set.seed(100)
#' n <- 500
#' p <- n * 0.3
#' beta <- runif(p, -1 / sqrt(p), 1 / sqrt(p))
#'
#' X <- mvtnorm::rmvt(n, sigma = diag(1, p), df = 2)
#' Y1 <- as.numeric(X %*% beta)
#' Y0 <- rep(0, n)
#'
#' pi1 <- 0.45
#' n1 <- n * pi1
#' ind <- sample(n, size = n1)
#' A <- rep(0, n)
#' A[ind] <- 1
#' Y <- Y1 * A + Y0 * (1 - A)
#'
#' Xc_svd <- svd(X)
#' H <- Xc_svd$u %*% t(Xc_svd$u)
#'
#' result_vec <- esti_mean_treat(X, Y, A, H)
#' print(paste0('True mean treat:', round(mean(Y1), digits = 3), '.'))
#' print('Absolute bias:')
#' print(abs(result_vec - mean(Y1)))
#'
#'
esti_mean_treat <- function(X,Y,A,H=NULL){
  n <- nrow(X); p <- ncol(X)
  pi1 <- mean(A)
  n1 <- sum(A)
  n0 <- n-n1

  if(is.null(H)){
    Xc <- scale(X,scale=FALSE)
    Xc_svd <- svd(Xc)
    H <- Xc_svd$u %*% t(Xc_svd$u)
  }

  tau_unadj <- mean(A * Y) / mean(A)

  YcA <- as.numeric(A * (Y - tau_unadj))
  tau_adj <- tau_unadj - mean((A / pi1 - 1) * H %*% YcA *  n / n1 )

  tau_db <- tau_adj + (1 - pi1) / pi1 * mean(A / pi1 * diag(H) * (Y - tau_unadj))

  ps_resid <- A / pi1 - 1
  or_resid <- A * Y / pi1
  or_resid_c <- YcA / pi1
  IF22 <- (t(ps_resid) %*% H %*% or_resid - sum(ps_resid * diag(H) * or_resid)) / n
  tau_adj2 <- tau_unadj - IF22

  IF22_dag <- - (t(ps_resid) %*% H %*% or_resid_c - sum(ps_resid * diag(H) * or_resid_c)) / n
  tau_adj2_dag <- tau_unadj + IF22_dag

  tau_adj3 <- tau_adj2 + (1-pi1)/pi1/n/(n-1)*sum(diag(H)*or_resid)

  point_est <- c(tau_unadj,tau_adj,tau_db,tau_adj2,tau_adj2_dag,tau_adj3)
  names(point_est) <- c('unadj','adj','db','adj2','adj2_dag','adj3')
  return(point_est)

}


if(F){

  set.seed(100)
  n <- 500
  p <- n*0.3
  beta <- runif(p, - 1 / sqrt(p), 1 / sqrt(p))

  X <- mvtnorm::rmvt(n, sigma = diag(1, p), df = 2)
  Y1 <- as.numeric(X %*% beta)
  Y0 <- rep(0,n)

  pi1 <- 0.45
  n1 <- n*pi1
  ind <- sample(n, size = n1)
  A <- rep(0, n)
  A[ind] <- 1
  Y <- Y1*A + Y0*(1-A)

  source('./fun_estimate_mean_treat.R')
  Xc_svd <- svd(Xc)
  H <- Xc_svd$u %*% t(Xc_svd$u)

  result_vec <- esti_mean_treat(X,Y,A,H)
  print(paste0('True mean treat:', round(mean(Y1),digits=3),'.'))
  print('Absolute bias:')
  print(abs(result_vec - mean(Y1)))



}
