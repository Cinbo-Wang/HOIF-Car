#' Estimate the oracle bias, the oracle variance of the unadjusted estimator, the Lin’s (2013) adjusted estimator and the debiased estimator tau_db by Lu et al.(2023).
#'
#' @param X The n by p covariates matrix.
#' @param Y1 Vector of n dimensional potential response Y(1).
#' @param pi1 The proportion of subjects in the treatment group.
#'
#' @return A list of variance_unadj, variance_adj and variance_db.
#' \item{variance_unadj}{The oracle variance of the unadjusted estimator.}
#' \item{variance_adj}{The oracle variance of Lin’s (2013) adjusted estimator tau_adj with linear working model.}
#' \item{variance_db}{The oracle variance of the debiased estimator tau_db by Lu et al.(2023).}
#' @export
#'
#' @references
#' Winston Lin. Agnostic notes on regression adjustments to experimental data: Reexamining Freedman’s critique. The Annals of Statistics, 7(1):295–318, 2013.
#'
#' Xin Lu, Fan Yang, and Yuhao Wang. Debiased regression adjustment in completely randomized experiments with moderately high-dimensional covariates. arXiv preprint arXiv:2309.02073, 2023.
#'
#' @examples
get_oracle_var_adj_db <- function(X,Y1,pi1=0.5){
  n <- nrow(X); p <- ncol(X)
  n1 <- n*pi1; n0 <- n-n1

  Xc <- scale(X,scale=FALSE)
  Xc_svd <- svd(Xc)
  H <- Xc_svd$u %*% t(Xc_svd$u)

  get_S2 <- function(x){
    n <- length(x)
    s2 <- sum((x-mean(x))**2)/(n-1)
    return(s2)
  }
  SY1 <- get_S2(Y1)
  Stau <- get_S2(Y1)
  sigma2_unadj <- SY1/pi1 - Stau


  e1 <- Y1 - mean(Y1) - H %*% (Y1-mean(Y1))
  Se1 <- get_S2(e1)
  Stau_e <- get_S2(e1)
  var_adj <- Se1/pi1 - Stau_e


  s1 <- diag(H)*(Y1-mean(Y1)) - mean(diag(H)*(Y1-mean(Y1)))
  Se1_s1 <- get_S2(e1+s1)
  Stau_e_s <- get_S2(e1+s1)
  sigma2_hd_l <- Se1_s1 / pi1 - Stau_e_s

  Q <- H^2
  diag(Q) <- diag(H) - diag(H)^2
  get_mat_S2 <- function(Q,x){
    tmp_mat <- Q * ((x-mean(x)) %*% t(x-mean(x)))
    return(sum(tmp_mat)/(nrow(Q)-1))
  }

  sigma2_hd_q <- get_mat_S2(Q=Q,x=Y1/pi1**2) * (pi1*(1-pi1))^2
  sigma2_hd <- sigma2_hd_l + sigma2_hd_q


  return(list(
    variance_unadj = sigma2_unadj/n,
    variance_adj = var_adj/n,
    variance_db = sigma2_hd/n
  ))

}




#' Estimate the oracle bias, the exact variance and approximated variance of the debiased estimator motivated by HOIF.
#'
#' @param X The n by p covariates matrix.
#' @param Y1 Vector of n dimensional potential response Y(1).
#' @param pi1 The proportion of subjects in the treatment group.
#'
#' @return A list of oracle bias and variance.
#' \item{bias_adj2_dag}{The oracle bias of the estimator.}
#' \item{variance_exact_adj2_dag}{The oracle exact variance of the estimator.}
#' \item{variance_approx_adj2_dag}{The oracle approximated variance of the estimator which omits the term of order o(1/n).}
#' @export
#'
#' @references
#' Zhao, S., Wang, X., Liu, L., & Zhang, X. (2024). Covariate adjustment in randomized experiments motivated by higher-order influence functions. arXiv preprint. https://arxiv.org/abs/2411.08491
#'
#' @examples
get_oracle_bias_var_adj2_dag <- function(X,Y1,pi1=0.5){
  n <- nrow(X); p <- ncol(X)
  n1 <- n*pi1; n0 <- n-n1

  Xc <- scale(X,scale=FALSE)
  Xc_svd <- svd(Xc)
  H <- Xc_svd$u %*% t(Xc_svd$u)
  h <- diag(H)
  all_1s <- rep(1, n)
  tau <- mean(Y1)
  tau2 <- mean(Y1**2)

  Sigmay <- t(Xc) %*% diag(Y1) %*% Xc

  bias_tau_adj2_dag <- 2 * (1 - pi1) / pi1 * (n1 - 1) / (n1 - pi1) * (p / n * tau + (sum(H %*% Y1) - sum(diag(H) * Y1)) / n) / (n - 2)

  var_tau_unadj <- (1 - pi1) / pi1 / n * var(Y1)


  Sigma_hat <- t(Xc) %*% Xc
  pi0 <- 1-pi1

  tmp_term <- 0
  for(i in 1:n){
    for(j in 1:n){
      if(i!=j){
        tmp_term <- tmp_term + H[i,j]**2*Y1[i]*Y1[j]
      }
    }
  }
  term12_right <- pi0/pi1^3*(pi1-pi0/(n-1))*(
    (n^2*pi1^2*(2-pi1)/(n-2)/(n-3) + 2*n*pi1*(4*pi1-5)/(n-2)/(n-3) + 12*pi0/(n-2)/(n-3))*tau^2 +
      (n*pi1*(pi1-1)*(pi1-2)/(n-2)/(n-3) + (-5*pi1^2+8*pi1-4)/(n-2)/(n-3))*tau2
  )*sum(diag(H)*(1-diag(H))) +
    2*pi0/pi1^3*(pi1-pi0/(n-1))*(
      n^2*pi1^2*(pi1-2)/(n-2)/(n-3) + n*pi1*(-pi1^2-5*pi1+8)/(n-2)/(n-3) + (5*pi1^2+4*pi1-8)/(n-2)/(n-3)
    )*tau*sum(diag(H)*(1-diag(H))*Y1) +
    (pi0-pi0/pi1/n + 2*pi0/pi1^3*(pi1-pi0/(n-1))*(
      -n^2*pi1*(pi1^2-4*pi1+2)/n/(n-2)/(n-3) + (4*pi1^2-4*pi1+2)/n/(n-2)/(n-3) +
        n*(pi1**3-4*pi1**2-2*pi1+2)/n/(n-2)/(n-3)
    ))*sum(diag(H)*(1-diag(H))*Y1**2) +
    pi0/pi1^3*(pi1-pi0/(n-1))*(n^3*pi1^2*pi0/n/(n-2)/(n-3) + n^2*pi1*(3*pi1^2+3*pi1-4)/n/(n-2)/(n-3)+
                                 n*(-2*pi1**3-8*pi1**2+4)/n/(n-2)/(n-3)+(-2*pi1^2+8*pi1-4)/n/(n-2)/(n-3))*tmp_term

  tmp_term1 <- tmp_term2 <- tmp_term3 <- tmp_term4 <- 0
  for(i in 1:n){
    for(j in 1:n){
      if(i!=j){
        tmp_term1 <- tmp_term1 + H[i,i]*H[i,j]*Y1[j]
        tmp_term2 <- tmp_term2 + H[i,i]*H[i,j]*Y1[j]**2
        tmp_term3 <- tmp_term3 + H[i,i]*H[i,j]*Y1[i]*Y1[j]
        tmp_term4 <- tmp_term4 + H[i,j]*(1-2*H[j,j])*Y1[i]*Y1[j]
      }
    }
  }
  term34_right <- pi0/pi1^3*(pi1-pi0/(n-1))*(pi1-2*pi0/(n-2))*(
    (pi1-3*pi0/(n-3))*((n+8)*pi1/(n-4)-8/(n-4))*tau^2 + (n*pi1*pi0/(n-3)/(n-4)+(-8*pi1^2+13*pi1-6)/(n-3)/(n-4))*tau2
  )*sum(diag(H)*(2*diag(H)-1)) +
    pi0/pi1^3*(pi1-pi0/(n-1))*(pi1-2*pi0/(n-2))*(
      n*pi1*(-7*pi1+5)/(n-3)/(n-4) + (4*pi1^2+14*pi1-12)/(n-3)/(n-4)
    )*tau*sum(diag(H)*(2*diag(H)-1)*Y1) +
    pi0/pi1^3*(pi1-pi0/(n-1))*(pi1-2*pi0/(n-2))*(
      -2*n^2*pi1^2/(n-3)/(n-4) + n*pi1*(-7*pi1+15)/(n-3)/(n-4) + (12*pi1^2+8*pi1-24)/(n-3)/(n-4)
    )*tau*(sum(diag(H)*(diag(H)-1)*Y1)-tmp_term1) +
    pi0/pi1^3*(pi1-pi0/(n-1))*(n^2*pi1*(7*pi1^2-9*pi1+3)/n/(n-2)/(n-3)/(n-4) +
                                 n*(-4*pi1**3-9*pi1**2+13*pi1-4)/n/(n-2)/(n-3)/(n-4)+
                                 (-4*pi1**2+16*pi1-8)/n/(n-2)/(n-3)/(n-4))*sum(diag(H)*(2*diag(H)-1)*Y1**2) +
    pi0/pi1^3*(pi1-pi0/(n-1))*(-2*n^3*pi1**2*pi0/n/(n-2)/(n-3)/(n-4) + n**2*pi1*(7*pi1^2-14*pi1+9)/n/(n-2)/(n-3)/(n-4) +
                                 n*(-12*pi1^3-pi1^2+15*pi1-8)/n/(n-2)/(n-3)/(n-4) - 4*pi0*(1+3*pi0)/n/(n-2)/(n-3)/(n-4)
    )*(sum(diag(H)*(diag(H)-1)*Y1**2)-tmp_term2) +
    pi0/pi1^3*(pi1-pi0/(n-1))*(
      n^3*pi1^2*(7*pi1-5)/n/(n-2)/(n-3)/(n-4) + n^2*pi1*(-18*pi1**2-10*pi1+16)/n/(n-2)/(n-3)/(n-4) +
        n*(8*pi1^3+26*pi1^2+2*pi1-16)/n/(n-2)/(n-3)/(n-4)+(8*pi1^2-32*pi1+16)/n/(n-2)/(n-3)/(n-4)
    )*(-tmp_term3 - sum(diag(Sigmay %*% solve(Sigma_hat) %*% Sigmay %*% solve(Sigma_hat))) + sum(diag(H)**2*Y1**2)) +
    pi0/pi1^3*(pi1-pi0/(n-1))*(
      pi1^3 + n^3*pi1^2*(7*pi1-5)/n/(n-2)/(n-3)/(n-4) + n^2*pi1*(-30*pi1^2+8*pi1+10)/n/(n-2)/(n-3)/(n-4) +
        n*(32*pi1**3-8*pi1-8)/n/(n-2)/(n-3)/(n-4) + (8*pi1**2-12*pi1+8)/n/(n-2)/(n-3)/(n-4)
    )*tmp_term4


  tmp_term1 <- tmp_term2 <- tmp_term3 <- tmp_term4 <- 0
  for(i in 1:n){
    for(j in 1:n){
      if(i!=j){
        tmp_term1 <- tmp_term1 + H[i,i]*H[i,j]*Y1[j]
        tmp_term2 <- tmp_term2 + H[i,i]*H[i,j]*Y1[j]**2
        tmp_term3 <- tmp_term3 + H[i,i]*H[i,j]*Y1[i]*Y1[j]
        tmp_term4 <- tmp_term4 + H[i,j]*(1-2*H[j,j])*Y1[i]*Y1[j]
      }
    }
  }
  term56_right <- pi0/pi1^3*(pi1-pi0/(n-1)) *(
    (pi1-2*pi0/(n-2))*(n*pi1*(4*pi1-5)/(n-3)/(n-4)+24*pi0/(n-3)/(n-4))*tau^2 +
      (-4*n*pi1*pi0**2/(n-2)/(n-3)/(n-4)+(20*pi1^2-30*pi1+12)/(n-2)/(n-3)/(n-4))*tau2
  )*sum(diag(H)*(2*diag(H)-1)) +
    pi0/pi1^3*(pi1-pi0/(n-1))*(pi1-2*pi0/(n-2))*(5*n*pi1*pi0/(n-3)/(n-4)+(4*pi1^2+24*pi1-24)/(n-3)/(n-4))*tau*(
      sum(diag(H)*(diag(H)-1)*Y1) - tmp_term1
    ) +
    pi0/pi1^3*(pi1-pi0/(n-1))*(n^2*pi1^2*(-3*pi1+5)/(n-2)/(n-3)/(n-4) +
                                 n*pi1*(4*pi1^2+14*pi1-26)/(n-2)/(n-3)/(n-4) +
                                 (-32*pi1^2+12*pi1+24)/(n-2)/(n-3)/(n-4))*tau*sum(diag(H)*(2*diag(H)-1)*Y1) +
    pi0/pi1^3*(pi1-pi0/(n-1))*(n^2*pi1*(5*pi1^2-8*pi1+3)/n/(n-2)/(n-3)/(n-4) +
                                 n*(-4*pi1^3-19*pi1^2+29*pi1-8)/n/(n-2)/(n-3)/(n-4) +
                                 (-4*pi1^2+20*pi1-16)/n/(n-2)/(n-3)/(n-4))*(
                                   sum(diag(H)*(diag(H)-1)*Y1**2) - tmp_term2
                                 ) +
    (pi0/pi1*(1-n*pi1)/n/(n-1) + pi0/pi1^3*(pi1-pi0/(n-1))*(
      n^2*pi1*(3*pi1^2 -11 *pi1+5)/n/(n-2)/(n-3)/(n-4) + n*(-4*pi1^3+25*pi1^2-pi1-4)/n/(n-2)/(n-3)/(n-4)+
        (-28*pi1^2+16*pi1-8)/n/(n-2)/(n-3)/(n-4)
    ))*sum(diag(H)*(2*diag(H)-1)*Y1**2) +
    pi0/pi1^3*(pi1-pi0/(n-1))*(
      -3*n^3*pi1^2*pi0/n/(n-2)/(n-3)/(n-4) + n^2*pi1*(-10*pi1^2-10*pi1+16)/n/(n-2)/(n-3)/(n-4) +
        n*(8*pi1^3+34*pi1^2-10*pi1-16)/n/(n-2)/(n-3)/(n-4) + (8*pi1^2-32*pi1+16)/n/(n-2)/(n-3)/(n-4)
    )*(-tmp_term3 - sum(diag(Sigmay %*% solve(Sigma_hat) %*% Sigmay %*% solve(Sigma_hat))) + sum(diag(H)**2*Y1**2)) +
    pi0/pi1^3*(pi1-pi0/(n-1))*(
      -n^3*pi1^2*pi0/n/(n-2)/(n-3)/(n-4) + n^2*pi1*(-2*pi1^2-4*pi1+6)/n/(n-2)/(n-3)/(n-4)+
        n*(6*pi1^2-8)/n/(n-2)/(n-3)/(n-4)+(-4*pi1+8)/n/(n-2)/(n-3)/(n-4)
    )*tmp_term4


  tmp_term1 <- tmp_term2 <- tmp_term3 <- tmp_term4 <- 0
  for(i in 1:n){
    for(j in 1:n){
      if(i!=j){
        tmp_term1 <- tmp_term1 + H[i,j]*H[j,j]*Y1[i]
        tmp_term2 <- tmp_term2 + H[i,j]*H[j,j]*Y1[i]**2
        tmp_term3 <- tmp_term3 + H[i,j]*(4*H[j,j]-p)*Y1[i]*Y1[j]

        tmp_term4 <- tmp_term4 + H[i,j]*(4*H[j,j]-1)*Y1[i]*Y1[j]
      }
    }
  }
  term7_right <- pi0/pi1^3*(pi1-pi0/(n-1))*(pi1-2*pi0/(n-2))*
    ((-n**2*pi1**2+n*pi1*(-20*pi1+23)-60*pi0)*tau**2/(n-3)/(n-4)/(n-5) +
       (-n*pi1*pi0+20*pi1**2-30*pi1+12)*tau2/(n-3)/(n-4)/(n-5))*(
         p*(2+p) - 6*sum(diag(H)**2)
       ) +
    pi0/pi1^3*(pi1-pi0/(n-1))*(pi1-2*pi0/(n-2))*((n**2*pi1**2+n*pi1*(19*pi1-22) -
                                                    20*pi1^2-30*pi1+48))/(n-3)/(n-4)/(n-5)*tau*2*(
                                                      sum(diag(H)*(2+p-4*diag(H))*Y1)+2*tmp_term1
                                                    ) +
    pi0/pi1^3*(pi1-pi0/(n-1))/n/(n-2)*(
      n^3*pi1^2*pi0/(n-3)/(n-4)/(n-5) + n^2*(-19*pi1^3+26*pi1^2-10*pi1)/(n-3)/(n-4)/(n-5) +
        n*(20*pi1^3+21*pi1^2-42*pi1+12)/(n-3)/(n-4)/(n-5) +
        (20*pi1^2-60*pi1+36)/(n-3)/(n-4)/(n-5)
    )*2*(
      sum(diag(H)*(2+p-4*diag(H))*Y1**2) + 2*tmp_term2
    ) + pi0/pi1^3*(pi1-pi0/(n-1))*(pi1-2*pi0/(n-2))*(
      n^2*(3*pi1-4*pi1**2)/n/(n-3)/(n-4)/(n-5) +
        (n*(pi1+10*pi1^2-6)-10*pi1+6)/n/(n-3)/(n-4)/(n-5)
    )*2*(
      -2*sum(diag(H)**2*Y1**2) + 2*sum(diag(Sigmay %*% solve(Sigma_hat) %*% Sigmay %*% solve(Sigma_hat))) +
        tmp_term3
    ) + pi0/pi1^3*(pi1-pi0/(n-1))*(
      -n^4*pi1^3+n^3*pi1^2*(-9*pi1+16)+n^2*pi1*(38*pi1**2-2*pi1-48) +
        n*(-40*pi1^3-22*pi1^2+16*pi1+48) + (-40*pi1^2+80*pi1-48)
    )/n/(n-2)/(n-3)/(n-4)/(n-5)*(
      -2*sum(diag(H)**2*Y1**2) + sum(diag(Sigmay %*% solve(Sigma_hat)))**2 +
        sum(diag(Sigmay %*% solve(Sigma_hat) %*% Sigmay %*% solve(Sigma_hat))) +
        tmp_term4
    )

  Expectation_square <- term12_right + term34_right + term56_right + term7_right

  tr_half <- sum(diag(Sigmay %*% solve(Sigma_hat)))

  Square_expectation <- (pi0/pi1*(pi1-pi0/(n-1))/(n-2))^2*(
    4*tau^2*(p**2) - 8*tau*p*sum(h*Y1)+ 4*tr_half**2
  )

  var_IF_22_dag_exact <- (Expectation_square - Square_expectation) / n**2/pi1**2



  U12 <- (t(Y1)%*% H %*% (Y1) - sum(h*Y1**2))/n^2/pi1^2
  U11 <- -tau/n/pi1^2*sum(h*Y1) + sum(h*Y1**2)/n^2/pi1^2 - U12
  U13 <- - sum(h*Y1**2)/n^2/pi1^2

  cov11 <- - n * pi1 * (1 - pi1) / (n - 1) / (n - 2) + 2*(1 - pi1)^2 / (n - 1) / (n - 2)
  cov12 <- pi1 * (1 - pi1) - (1 - pi1) * (1 - 2 * pi1) / (n - 1)
  cov13 <- - (1 - pi1)^2 / (n - 1)

  cov_unadj_IF_22_dag_term1 <- U11*cov11 + U12*cov12 + U13*cov13


  U21 <- pi0/pi1^2/n*U11
  U22 <- pi0/pi1**2/n*(U12+U13)
  cov21 <- -2*pi1*(1-pi1)*(n*pi1-1)/(n-1)/(n-2)
  cov22 <- pi1*pi0*(n*pi1-1) / (n-1)
  cov_unadj_IF_22_dag_term2 <- U21*cov21 + U22*cov22


  U31 <- U11/n/pi1
  U32 <- U12/n/pi1
  U33 <- U13/n/pi1
  cov31 <- - n * pi1 * (1 - pi1) / (n - 1) / (n - 2) + 2*(1 - pi1)^2 / (n - 1) / (n - 2)
  cov32 <- pi1 * (1 - pi1) - (1 - pi1) * (1 - 2 * pi1) / (n - 1)
  cov33 <- - (1 - pi1)^2 / (n - 1)
  cov_unadj_IF_22_dag_term3 <- U31*cov31 + U32*cov32 + U33*cov33


  U41 <- -tau**2*p/n/pi1**3 + tau2*p/pi1^3/n^2 + 4*tau/pi1^3/n^2*sum(h*Y1) -
    4/pi1^3/n^3*sum(h*Y1**2) + 2/pi1^3/n^3*(t(Y1)%*%H%*%Y1 - sum(h*Y1**2))
  U43 <- U42 <- U11/n/pi1
  U44 <- -tau2*p/pi1^3/n^2 + 2/pi1^3/n^3*sum(h*Y1**2)
  cov41 <- -(1-pi1)*(n*pi1-1)*(n*pi1-6*pi0)/(n-1)/(n-2)/(n-3)
  cov42 <- (1-2*pi1)*(n*pi1-1)*(n*pi1-2)/(n-1)/(n-2) + pi1^2*(n*pi1-1)/(n-1)
  cov43 <- -2*pi0^2*(n*pi1-1)/(n-1)/(n-2)
  cov44 <- -2*pi0^2*(n*pi1-1)/(n-1)/(n-2)
  cov_unadj_IF_22_dag_term4 <- U41*cov41 + U42*cov42 + U43*cov43 + U44*cov44

  cov_unadj_IF_22_dag <- cov_unadj_IF_22_dag_term1 - (cov_unadj_IF_22_dag_term2 +
                                                        cov_unadj_IF_22_dag_term3 + cov_unadj_IF_22_dag_term4)


  cov_unadj_IF_22_dag_approx <- tau^2*pi0*p/pi1/n^2+pi0/pi1/n^2*(t(Y1) %*% H %*% Y1 - sum(h*Y1**2))

  var_exact_adj2_dag <- var_tau_unadj + var_IF_22_dag_exact - 2*cov_unadj_IF_22_dag


  tr <- sum(diag(Sigmay %*% solve(Sigma_hat) %*% Sigmay %*% solve(Sigma_hat)))
  term_approx <- (2*pi0/pi1**2*p/n**2 - pi0/pi1*p**2/n**3)*tau**2 +
    (-2*pi0/pi1**2/n**2 + 3*pi0/pi1/n**2)*tau**2*sum(h**2) +
    (-4*pi0/pi1**2/n**2 + 2*pi0*p / pi1/n**3)*tau*sum(h*Y1) +
    4*(pi0/pi1)**2/n**2*tau*sum(h**2*Y1) +
    2*pi0/pi1/n**2*tau*(sum(h*H%*%Y1) - sum(h**2*Y1)) +
    pi0/pi1**2/n**2*sum(h*Y1**2) -
    (pi0/pi1**2/n**2 + (pi0**2/pi1**2)/n**2)*sum(h**2*Y1**2) -
    pi0/pi1/n**3*tr_half**2 +
    (pi0/pi1)^2/n**2*tr -
    pi0/pi1/n**2*(t(Y1)%*%H%*%Y1-sum(h*Y1**2)) -
    2*pi0/pi1/n**2*sum(t(h*Y1) %*% H %*% Y1 - sum(h*Y1*h*Y1))

  var_approx_adj2_dag <- var_tau_unadj + term_approx

  return(list(
    bias_adj2_dag = as.numeric(bias_tau_adj2_dag),
    variance_exact_adj2_dag = var_exact_adj2_dag,
    variance_approx_adj2_dag = var_approx_adj2_dag,
    variance_unadj = var_tau_unadj
  ))


}



#' Estimate the bias, the exact variance and approximated variance of the HOIF-motivated adjusted estimators.
#'
#' @param X The n by p covariates matrix.
#' @param Y1 Vector of n dimensional potential response Y(1).
#' @param pi1 The proportion of subjects in the treatment group.
#'
#' @return A list of bias_adj2, variance_exact_adj2, variance_approx_adj2 and variance_exact_adj3.
#' \item{bias_adj2}{The bias of the HOIF-motivated estimator tau_adj_2.}
#' \item{variance_exact_adj2}{The exact vairance of the HOIF-motivated estimator tau_adj_2.}
#' \item{variance_approx_adj2}{The approximated vairance of the HOIF-motivated estimator tau_adj_2.}
#' \item{variance_exact_adj3}{The exact vairance of the HOIF-motivated estimator tau_adj_3.}
#' @export
#'
#' @references
#' Zhao, S., Wang, X., Liu, L., & Zhang, X. (2024). Covariate adjustment in randomized experiments motivated by higher-order influence functions. arXiv preprint. https://arxiv.org/abs/2411.08491
#'
#' @examples
#' # Linear setting
#' set.seed(100)
#' n <- 500
#' p <- 50
#' beta <- rt(p,3)
#'
#' X <- mvtnorm::rmvt(n, sigma = diag(1, p), df = 2)
#' Y1 <- as.numeric(X %*% beta)
#' pi1 <- 0.45
#'
#' result_adj_db <- get_oracle_var_adj_db(X = X,Y1=Y1,pi1=pi1) # from LYW paper
#' result_adj2_dag <- get_oracle_bias_var_adj2_dag(X = X,Y1=Y1,pi1=pi1)
#' result_adj2_3 <- get_oracle_bias_var_adj_2_3(X=X, Y1=Y1,pi1=pi1)
#' unlist(result_adj_db)
#' unlist(result_adj2_dag)
#' unlist(result_adj2_3)
#'
#'
#'
#' # Nonlinear setting
#' n <- 500;
#' alpha <- 0.2;
#' set.seed(1000)
#' p <- ceiling(n*alpha)
#' Sigma_true <- matrix(0,nrow=p,ncol=p)
#' for(i in 1:p){
#'   for(j in 1:p){
#'     Sigma_true[i,j] <- 0.3**(abs(i-j))
#'   }
#' }
#'
#' X <- mvtnorm::rmvt(n, sigma = Sigma_true, df = 3)
#' beta <- rt(p,3)
#' or_baseline <- sign(X %*% beta) * abs(X %*% beta)^(1/2) + sin(X %*% beta)
#' epsilon1 <- epsilon0 <- rt(n,3)
#' Y1 <- 1 + as.numeric(or_baseline) + epsilon1
#'
#'
#' pi1 <- 0.45
#' result_adj_db <- get_oracle_var_adj_db(X = X,Y1=Y1,pi1=pi1)
#' result_adj2_dag <- get_oracle_bias_var_adj2_dag(X = X,Y1=Y1,pi1=pi1)
#' result_adj2_3 <- get_oracle_bias_var_adj_2_3(X=X, Y1=Y1,pi1=pi1)
#' unlist(result_adj_db)
#' unlist(result_adj2_dag)
#' unlist(result_adj2_3)
#'
get_oracle_bias_var_adj_2_3 <- function(X,Y1,pi1=0.5){

  n <- nrow(X); p <- ncol(X)
  pi0 <- 1-pi1
  n1 <- n*pi1; n0 <- n - n1

  Xc <- scale(X,scale=FALSE)
  Xc_svd <- svd(Xc)

  H <- Xc_svd$u %*% t(Xc_svd$u)
  h <- diag(H)
  all_1s <- rep(1, n)
  Sigma_hat <- t(Xc) %*% Xc

  if(ncol(X)==1){
    Omega <- Xc_svd$d^(-2) * Xc_svd$v  %*% t(Xc_svd$v)
  }else{
    Omega <- Xc_svd$v %*% diag(Xc_svd$d^(-2)) %*% t(Xc_svd$v)
  }

  Sigmay <- t(Xc) %*% diag(Y1) %*% Xc #


  tau <- mean(Y1)
  bias_tau_adj2 <- (1 - pi1) / pi1 * (t(all_1s) %*% H %*% Y1 - sum(all_1s * diag(H) * Y1)) / (n * (n - 1))

  var_tau_unadj <- (1 - pi1) / pi1 / n * var(Y1)
  V1 <- sum(diag(H) * (1 - diag(H)) * Y1^2) / n^2
  V2 <- sum(diag(Sigmay %*% Omega %*% Sigmay %*% Omega)) / n^2 - sum(diag(H)^2 * Y1^2) / n^2
  V3 <- (t(Y1) %*% H %*% diag(1 - 2 * diag(H)) %*% Y1 - sum(diag(H) * (1 - 2 * diag(H)) * Y1^2)) / n^2
  V4 <- V5 <- - t(Y1) %*% H %*% diag(diag(H)) %*% Y1 / n^2 + sum(diag(H)^2 * Y1^2) / n^2 - sum(diag(Sigmay %*% Omega %*% Sigmay %*% Omega)) / n^2 + sum(diag(H)^2 * Y1^2) / n^2
  V6 <- sum(diag(H) * (2 * diag(H) - 1) * Y1^2) / n^2
  V7 <- t(Y1) %*% (diag(4 * diag(H) - 1)) %*% H %*% Y1 / n^2 - sum(Y1^2 * diag(H) * (4 * diag(H) - 1)) / n^2 + sum(diag(Sigmay %*% Omega))^2 / n^2 + sum(diag(Sigmay %*% Omega %*% Sigmay %*% Omega)) / n^2 - 2 * sum(diag(H)^2 * Y1^2) / n^2

  cov1 <- (1 - pi1) - (1 - pi1) * (1 - 2 * pi1) / pi1 / (n - 1) - (1 - pi1)^2 / (n - 1)^2
  cov2 <- (1 - pi1)^2 - (1 - pi1)^3 / pi1 / (n - 1) - (1 - pi1)^2 / (n - 1)^2
  cov3 <- pi1 * (1 - pi1) - (1 - pi1)^2 / (n - 1) - 2 * (1 - pi1) * (1 - 2 * pi1) / (n - 2) + (2 / pi1 - 5 + 1 / (n - 1)) * (1 - pi1)^2 / (n - 1) / (n - 2)
  cov4 <- cov5 <- - 2 * (1 - pi1)^2 / (n - 2) + (2 / pi1 - 3 + 1 / (n - 1)) * (1 - pi1)^2 / (n - 1) / (n - 2)
  cov6 <- - (1 - pi1) * n / (n - 1) / (n - 2) - (1 - pi1)^2 / (n - 1)^2 +  2 * (1 - pi1)^2 / pi1 / (n - 1) / (n - 2)
  cov7 <- - pi1 * (1 - pi1) / (n - 2) + 3 * (2 - 3 * pi1) * (1 - pi1) / (n - 2) / (n - 3) - 3 * (2 - 3 * pi1) * (1 - pi1)^2 / pi1 / (n - 1) / (n - 2) / (n - 3) + (1 - pi1)^2 / (n - 1)^2 / (n - 2)

  cov8 <- - n * pi1 * (1 - pi1) / (n - 1) / (n - 2) + 2*(1 - pi1)^2 / (n - 1) / (n - 2)
  cov9 <- pi1 * (1 - pi1) - (1 - pi1) * (1 - 2 * pi1) / (n - 1)
  cov10 <- - (1 - pi1)^2 / (n - 1)


  cov1_app <- 1 - pi1
  cov2_app <- (1 - pi1)^2
  cov3_app <- pi1 * (1 - pi1)
  cov4_app <- cov5_app <- - 2 * (1 - pi1)^2 / n
  cov6_app <- - (1 - pi1) / n
  cov7_app <- - pi1 * (1 - pi1) / n

  cov8_app <- - pi1 * (1 - pi1) / n
  cov9_app <- pi1 * (1 - pi1)
  cov10_app <- - (1 - pi1)^2 / n

  var_IF22_unadj <- (V1 * cov1 + V2 * cov2 + V3 * cov3 + V4 * cov4 + V5 * cov5 + V6 * cov6 + V7 * cov7) / pi1^2
  var_IF22_unadj_app <- (V1 * cov1_app + V2 * cov2_app + V3 * cov3_app + V4 * cov4_app + V5 * cov5_app + V6 * cov6_app + V7 * cov7_app) / pi1^2
  var_IF22_unadj_app_form <- (1 - pi1) / pi1 * sum(diag(H) * (1 - (2 - pi1) * diag(H)) * Y1^2 / pi1) / n^2 +
    ((1 - pi1) / pi1)^2 * sum(diag(Sigmay %*% Omega %*% Sigmay %*% Omega)) / n^2 -
    (1 - pi1) / pi1 * sum(diag(Sigmay %*% Omega))^2 / n^3 +
    (1 - pi1) / pi1 * (t(Y1) %*% H %*% diag(1 - 2 * diag(H)) %*% Y1 - sum(Y1^2 * diag(H) * (1 - 2 * diag(H)))) / n^2

  U1 <- - tau * sum(diag(H) * Y1) / n + 2 * sum(diag(H) * Y1^2) / n^2 - t(Y1) %*% H %*% Y1 / n^2
  U2 <- t(Y1) %*% H %*% Y1 / n^2 - sum(diag(H) * Y1^2) / n^2
  U3 <- - sum(diag(H) * Y1^2) / n^2

  cov_tau_IF22_unadj <- (U1 * cov8 + U2 * cov9 + U3 * cov10) / pi1^2
  cov_tau_IF22_unadj_app <- (U1 * cov8_app + U2 * cov9_app + U3 * cov10_app) / pi1^2
  cov_tau_IF22_unadj_app_form <- tau * (1 - pi1) / pi1 * sum(diag(Sigmay %*% Omega)) / n^2 + (1 - pi1) / pi1 * U2

  var_tau_adj2 <- var_tau_unadj + var_IF22_unadj - 2 * cov_tau_IF22_unadj
  var_tau_adj2_app <- var_tau_unadj + var_IF22_unadj_app - 2 * cov_tau_IF22_unadj_app
  var_tau_adj2_app_form <- var_tau_unadj - (1 - pi1) / pi1 * (t(Y1) %*% H %*% diag((1 + 2 * diag(H))) %*% Y1 - sum(Y1^2 * diag(H) * (1 + 2 * diag(H)))) / n^2 +
    (1 - pi1) / pi1 * sum(diag(H) * (1 - (2 - pi1) * diag(H)) * Y1^2 / pi1) / n^2 +
    (1 - pi1)^2 / pi1^2 * sum(diag(Sigmay %*% Omega %*% Sigmay %*% Omega)) / n^2 -
    (1 - pi1) / pi1 * sum(diag(Sigmay %*% Omega))^2 / n^3 -
    2 * tau * pi0/pi1*sum(diag(Sigmay %*% Omega)) / n^2

  var_tau_adj2_app_form_2 <- (1 - pi1) / pi1 / n^2 * sum(((1 + h) * (Y1) - H %*% (Y1) - mean((1 + h) * (Y1)))^2) +
    (1 - pi1)^2 / pi1^2 / n * (t(Y1) %*% H^2 %*% (Y1) / n + sum(h * (1 - 2 * h) * (Y1)^2) / n)

  tr <- sum(diag(Sigmay %*% solve(Sigma_hat)))
  var_delta <- (pi0/pi1)^3/n/(n-1)^3*(
    sum(h**2*Y1**2) - tr**2/n
  )
  cov_tau_delta <- (pi0/pi1)^2/n/(n-1)**2 * sum(h*Y1*(Y1-tau))

  cov_IF22_delta <- (pi0/pi1)^2/n^2/(n-1)^2*(pi1*n^2-n+pi0)/(n-2)/pi1*(t(Y1*h)%*%H%*%Y1 - sum(Y1*h*h*Y1)) +
    (pi0/pi1)^2/n/(n-1)^2/(n-2)*(1-pi0/pi1/n)*tr**2 +
    (pi0/pi1)^2/n/(n-1)^2*(pi0/pi1/n - 1/(n-2)*(1-pi0/pi1/n))*sum(h^2*Y1^2)
  var_tau_adj3 <- var_tau_adj2 + var_delta + cov_tau_delta - cov_IF22_delta

  return(list(
    bias_adj2 = as.numeric(bias_tau_adj2),
    variance_exact_adj2 = as.numeric(var_tau_adj2),
    variance_approx_adj2 = as.numeric(var_tau_adj2_app_form),
    variance_exact_adj3 = as.numeric(var_tau_adj3),
    variance_unadj = var_tau_unadj
  ))
}



if(F){
  # Linear setting
  set.seed(100)
  n <- 500
  p <- 50
  beta <- rt(p,3)

  X <- mvtnorm::rmvt(n, sigma = diag(1, p), df = 2)
  Y1 <- as.numeric(X %*% beta)
  pi1 <- 0.45

  result_adj_db <- get_oracle_var_adj_db(X = X,Y1=Y1,pi1=pi1)
  result_adj2_dag <- get_oracle_bias_var_adj2_dag(X = X,Y1=Y1,pi1=pi1)
  result_adj2_3 <- get_oracle_bias_var_adj_2_3(X=X, Y1=Y1,pi1=pi1)
  unlist(result_adj_db)
  unlist(result_adj2_dag)
  unlist(result_adj2_3)



  # Nonlinear setting
  n <- 500;
  alpha <- 0.2;
  set.seed(1000)
  p <- ceiling(n*alpha)
  Sigma_true <- matrix(0,nrow=p,ncol=p)
  for(i in 1:p){
    for(j in 1:p){
      Sigma_true[i,j] <- 0.3**(abs(i-j))
    }
  }

  X <- mvtnorm::rmvt(n, sigma = Sigma_true, df = 3)
  beta <- rt(p,3)
  or_baseline <- sign(X %*% beta) * abs(X %*% beta)^(1/2) + sin(X %*% beta)
  epsilon1 <- epsilon0 <- rt(n,3)
  Y1 <- 1 + as.numeric(or_baseline) + epsilon1


  pi1 <- 0.45
  result_adj_db <- get_oracle_var_adj_db(X = X,Y1=Y1,pi1=pi1)
  result_adj2_dag <- get_oracle_bias_var_adj2_dag(X = X,Y1=Y1,pi1=pi1)
  result_adj2_3 <- get_oracle_bias_var_adj_2_3(X=X, Y1=Y1,pi1=pi1)
  unlist(result_adj_db)
  unlist(result_adj2_dag)
  unlist(result_adj2_3)

}
