#' Estimate treatment effect and the corresponding variance estimation on the treatment arm using different covariates adjustment methods.
#'
#'
#' @param X The n by p covariates matrix.
#' @param Y Vector of n dimensional observed response.
#' @param A Vector of n dimensional treatment assignment.
#' @param H The n by n hat projection matrix corresponding to X.
#'
#' @return the list of point estimation and variance estimation for all kinds of estimators.
#' \item{tau_unadj}{The point estimation of the unadjusted estimator tau_unadj.}
#' \item{tau_hat_lin}{The point estimation of Lin’s (2013) adjusted estimator.}
#' \item{tau_hat_lin_db}{The point estimation of Lei’s (2020) debiased adjusted estimator using population leverage scores.}
#' \item{tau_db}{The point estimation of the debiased estimator by Lu et al.(2023).}
#' \item{tau_adj2c}{The point estimation of another form of the debiased estimator inspired by HOIF (Zhao et al.(2024)).}
#' \item{tau_adj2}{The point estimation of the adjusted estimator motivatied by HOIF (Zhao et al.(2024)).}
#' \item{tau_adj3}{The point estimation of the bias-free adjusted estimator based on tau_adj2.}
#' \item{var_hat_unadj}{The oracle variance of the unadjusted estimator tau_unadj.}
#' \item{var_hat_lin}{The HC3 type variance estimation of Lin’s (2013) adjusted estimator tau_hat_lin.}
#' \item{var_hat_lin_db}{The HC3 type variance estimation of Lei’s (2020) adjusted estimator tau_hat_lin_db.}
#' \item{var_hat_db}{The variance estimation of the debiased estimator tau_db.}
#' \item{var_hat_adj2c_wang}{The variance estimation of the adjusted estimator tau_adj2c where the estimation method follows Lu et al.(2023).}
#' \item{var_hat_adj2c_v2}{The variance estimation of the adjusted estimator tau_adj2c using the conservative method introduced in Zhao et al.(2024).}
#' \item{var_hat_adj2_wang}{The variance estimation of the adjusted estimator tau_adj2 where the estimation method follows Lu et al.(2023).}
#' \item{var_hat_adj2_v2}{The variance estimation of the adjusted estimator tau_adj2 using using the conservative method introduced in Zhao et al.(2024).}
#' @export
#'
#' @references
#' Winston Lin. Agnostic notes on regression adjustments to experimental data: Reexamining Freedman’s critique. The Annals of Statistics, 7(1):295–318, 2013.
#'
#' Lihua Lei, Peng Ding. Regression adjustment in completely randomized experiments with a diverging number of covariates. Biometrika, 815–828, 2020.
#'
#' Sihui Zhao, Xinbo Wang, Lin Liu, & Xin Zhang. Covariate adjustment in randomized experiments motivated by higher-order influence functions. arXiv preprint. https://arxiv.org/abs/2411.08491.
#'
#' Xin Lu, Fan Yang, and Yuhao Wang. Debiased regression adjustment in completely randomized experiments with moderately high-dimensional covariates. arXiv preprint arXiv:2309.02073, 2023.
#'
#' @examples
#' set.seed(100)
#' n <- 500
#' p <- n * 0.3
#' beta <- runif(p, -1 / sqrt(p), 1 / sqrt(p))
#'
#' X <- mvtnorm::rmvt(n, sigma = diag(1, p), df = 3)
#' Y1 <- as.numeric(X %*% beta)
#' Y0 <- rep(0, n)
#'
#' pi1 <- 2/3
#' n1 <- ceiling(n * pi1)
#' ind <- sample(n, size = n1)
#' A <- rep(0, n)
#' A[ind] <- 1
#' Y <- Y1 * A + Y0 * (1 - A)
#'
#' Xc_svd <- svd(X)
#' H <- Xc_svd$u %*% t(Xc_svd$u)
#'
#' result_ls <- esti_mean_treat(X, Y, A, H)
#' point_est <- result_ls$point_est
#' var_est <- result_ls$var_est
#' print(paste0('True mean treat:', round(mean(Y1), digits = 3), '.'))
#' print('Absolute bias:')
#' print(abs(point_est - mean(Y1)))
#' print('Estimate variance:')
#' print(var_est)
#'
esti_mean_treat <- function(X,Y,A,H=NULL){

  n <- nrow(X); p <- ncol(X)
  pi1 <- mean(A)
  pi0 <- 1-pi1

  n1 <- sum(A)
  n0 <- n-n1

  Xc <- scale(X,scale=FALSE)
  if(is.null(H)){
    Xc_svd <- svd(Xc)
    H <- Xc_svd$u %*% t(Xc_svd$u)
  }

  tau_unadj <- mean(A * Y) / mean(A)


  if(n1 > p){
    fit1 <- lm(Y~Xc,subset = (A==1))
    tau_hat_lin <- fit1$coefficients[1]
    e1 <- resid(fit1)
    lscores <- hat(Xc)
    lscores1 <- lscores[A==1]
    Delta1 <- sum(lscores1 * e1) / n1
    tau_hat_lin_db <- tau_hat_lin + n0 / n1 * Delta1

    point_vec_lin <- c(tau_hat_lin,tau_hat_lin_db)
  }else{
    point_vec_lin <- c(NA,NA)
  }
  names(point_vec_lin) <- c('lin','lin_db')


  YcA <- as.numeric(A * (Y - tau_unadj))
  tau_adj <- tau_unadj - mean((A / pi1 - 1) * H %*% YcA * n / n1)


  tau_db <- tau_adj + (1 - pi1) / pi1 * mean(A / pi1 * diag(H) * (Y - tau_unadj))


  ps_resid <- A / pi1 - 1
  or_resid <- A * Y / pi1
  or_resid_c <- YcA / pi1
  IF22 <- (t(ps_resid) %*% H %*% or_resid - sum(ps_resid * diag(H) * or_resid)) / n
  tau_adj2 <- tau_unadj - IF22


  IF22_c <- - (t(ps_resid) %*% H %*% or_resid_c - sum(ps_resid * diag(H) * or_resid_c)) / n
  tau_adj2c <- tau_unadj + IF22_c


  tau_adj3 <- tau_adj2 + (1-pi1)/pi1/n/(n-1)*sum(diag(H)*or_resid)

  point_est <- c(tau_unadj,tau_db,tau_adj2c,tau_adj2,tau_adj3,point_vec_lin)
  names(point_est) <- c('unadj','db','adj2c','adj2','adj3','lin','lin_db')



  h <- diag(H)

  var_hat_unadj <- pi0/pi1/n*sum((Y[A==1] - tau_unadj)^2)/(n1-1)


  if(n1 > p){
    mod1 <- lm(Y[A==1] ~ Xc[A==1, ])
    e1_hat <- as.numeric(resid(mod1))
    lscores1 <- hat(Xc[A==1, ])
    s1_HC3 <- sqrt(mean(e1_hat^2 / (1 - pmin(lscores1, 0.99))^2))
    var_hat_lin <- (1/n1 - 1/n)*s1_HC3^2

  }else{
    var_hat_lin <- NA
  }


  Q <- H^2
  diag(Q) <- diag(H) - diag(H)^2
  ones <- rep(1,n)
  M <- diag(n) - 1/n * ones %*% t(ones) - H + (diag(n) - 1/n * ones %*% t(ones)) %*% diag(h)
  B <- t(M) %*% M
  get_vec_s2 <- function(q,x,z,z0=1){
    q_sub <- q[z==z0]
    x_sub <- x[z==z0]
    s2 <- mean(q_sub*(x_sub-mean(x_sub))**2)
    return(s2)
  }
  I1_hat <- pi1*pi0*(get_vec_s2(q=pi1*pi0*diag(Q),x=Y/pi1^2,z=A,z0=1) +
                       get_vec_s2(q=diag(B)/pi1^2,x=Y,z=A,z0=1))

  get_mat_s2 <- function(Q,x,z,z0=1){
    Q_sub <- Q[z==z0,z==z0]
    x_sub <- x[z==z0] - mean(x[z==z0])

    s2 <- 1/pi1/n1*(t(x_sub) %*% Q_sub %*% x_sub - sum(x_sub * diag(Q_sub) * x_sub))
    return(s2)
  }
  I2_hat <- pi1*pi0*(get_mat_s2(Q = pi1*pi0*Q,x=Y/pi1**2,z=A,z0=1) +
                       get_mat_s2(Q=1/pi1^2*B,x=Y,z=A,z0=1))
  I3_hat <- 0
  I4_hat <- 0
  var_hat_db <- (I1_hat + I2_hat + I3_hat + I4_hat)/n


  Y1_hat <- A*Y/pi1
  tr_half_square_hat <- sum(h^2*A*Y^2/pi1) + sum((h*Y1_hat) %*% t(h*Y1_hat)) - sum((h*Y1_hat)^2)
  tr_hat <- sum(h^2*A*Y^2/pi1) + t(Y1_hat) %*% H^2 %*% Y1_hat - sum(Y1_hat*h^2*Y1_hat)
  Y1_sub <- Y[A==1]; H1_sub <- H[A==1,A==1]
  P <- diag(n) - 1/n*rep(1,n) %*% t(rep(1,n))
  M <- P - H + P%*%diag(diag(H))
  B <- t(M) %*% M
  var_hat_adj2c_wang <- (pi0/pi1)/n/n*(
    sum(diag(B)*A*(Y - tau_unadj)^2/pi1) + t(A*(Y-tau_unadj)/pi1) %*% B %*% (A*(Y-tau_unadj)/pi1) - sum((A*(Y-tau_unadj)/pi1) * diag(B) * (A*(Y-tau_unadj)/pi1))
  ) + (pi0/pi1)^2/n*(
    1/n1*sum(A*h*(1-h)*(Y-tau_unadj)^2) + 1/n1/pi1*(t(Y1_sub-tau_unadj) %*% H1_sub^2 %*% (Y1_sub - tau_unadj) - sum(diag(H1_sub)^2 * (Y1_sub - tau_unadj)^2))
  )


  resi <- Y - tau_unadj - (H - diag(diag(H))) %*% (A*(Y - tau_unadj)/pi1) - mean((1+diag(H))*A*(Y-tau_unadj)/pi1)
  var_hat_adj2c_v2 <- (pi0/pi1)/n/n*sum(A*resi^2/pi1) +
    (pi0/pi1)^2/n*(
      1/n1*sum(A*h*(1-h)*(Y-tau_unadj)^2) + 1/n1/pi1*(t(Y1_sub-tau_unadj) %*% H1_sub^2 %*% (Y1_sub - tau_unadj) - sum(diag(H1_sub)^2 * (Y1_sub - tau_unadj)^2))
    )


  P <- diag(n) - 1/n*rep(1,n) %*% t(rep(1,n))
  M <- P - H + P%*%diag(diag(H))
  B <- t(M) %*% M
  Y1_sub <- Y[A==1]; H1_sub <- H[A==1,A==1]
  var_hat_adj2_wang <- pi0/pi1/n/n*(
    sum(diag(B)*A*Y^2/pi1) + t(A*Y/pi1) %*% B %*% (A*Y/pi1) - sum((A*Y/pi1) * diag(B) * (A*Y/pi1))
  ) + (pi0/pi1)^2/n*(
    1/n1*sum(A*h*(1-h)*Y^2) + 1/n1/pi1*(t(Y1_sub) %*% H1_sub^2 %*% Y1_sub - sum(diag(H1_sub)^2 * Y1_sub^2))
  )


  resi <- Y - (H - diag(diag(H))) %*% (A*(Y)/pi1) - mean((1+diag(H))*A*(Y)/pi1)
  var_hat_adj2_v2 <- (pi0/pi1)/n/n*sum(A*resi^2/pi1) +
    (pi0/pi1)^2/n*(
      1/n1*sum(A*h*(1-h)*(Y)^2) + 1/n1/pi1*(t(Y1_sub) %*% H1_sub^2 %*% (Y1_sub) - sum(diag(H1_sub)^2 * (Y1_sub)^2))
    )


  var_hat_vec <- c(var_hat_unadj,var_hat_db,var_hat_adj2c_wang,var_hat_adj2c_v2,
                   var_hat_adj2_wang,var_hat_adj2_v2,
                   var_hat_adj2_wang,var_hat_adj2_v2,
                   var_hat_lin,var_hat_lin)
  names(var_hat_vec) <- c('unadj','db','adj2c','adj2c_v2',
                          'adj2','adj2_v2',
                          'adj3','adj3_v2',
                          'lin','lin_db')


  return(list(
    point_est = point_est,
    var_est = var_hat_vec
  ))

}

